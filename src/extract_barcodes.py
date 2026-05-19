from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import reverse_complement
import numpy as np
import gzip
import regex
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import nullcontext
import os
import logging

logger = logging.getLogger(__name__)
# per-process worker state and helpers (used by ProcessPoolExecutor initializer)
_worker_processor = None


class BarcodeBatchProcessor:
    """
    Worker class that processes batches of FastQ records.
    """

    def __init__(self, polyt_min_length: int, barcode_len: int, umi_len: int, search_range: int,
                 min_qual: int, whitelist: set, exact_anchor: bool, write_upstream_region: bool = False):
        """
        Initialize the processor with parameters for barcode extraction.

        Args:
            polyt_min_length (int):       The minimum size required for a polyT(A) to be identified in the read
            barcode_len (int):            The length of the cell barcode sequences
            umi_len (int):                The length of the UMI used
            search_range (int):           The range around the UMI length to look for exact matches
            min_qual (int):               The minimum average quality in a barcode for it to be considered
            whitelist (set):              The whitelist of cell barcodes
            exact_anchor (bool):          Whether or not to enforce anchors must be one per read and at the exact position as expected
            write_upstream_region (bool): Whether to emit trimmed upstream sequence lines for output
        """
        self.barcode_len = barcode_len
        self.umi_len = umi_len
        self.search_range = search_range
        self.min_qual = min_qual
        self.whitelist = whitelist
        self.exact_anchor = exact_anchor
        self.write_upstream_region = write_upstream_region
        self.polyt_pattern = f'TTTT(?:T{{{polyt_min_length - 4},}}){{e<=1}}'

    def process(self, records_batch):
        """
        Process a batch of FastQ records to extract barcodes.

        Args:
            records_batch (list): A list of tuples (title, sequence) for each read in the batch.

        Returns:
            list: Lines for putative barcodes to write to file
            list: Lines for anchor barcodes to write to file
            list: Lines for trimmed upstream region to write to file
            Counter: Frequency of each putative barcode in the batch
            Counter: Frequency of each anchor barcode in the batch
            int: Number of reads with polyA found in the batch
            int: Total number of reads processed in the batch
        """
        local_putative_lines = []
        local_anchor_lines = []
        local_upstream_region_lines = []
        local_putative_freq = Counter()
        local_anchor_freq = Counter()
        local_polya = 0
        local_no_polya = 0
        local_processed = 0
        local_anchor = 0
        local_quality_filtered_reads = 0
        local_quality_filtered_barcodes = 0
        local_only_putative_reads = 0

        for record in records_batch:
            upstream_tail, upstream_qual, search_to_polyt_region = self._find_upstream(record[1], record[2])
            if upstream_tail is not None:
                local_polya += 1
                if self.write_upstream_region:
                    local_upstream_region_lines.append(f'{record[0]}\t{search_to_polyt_region}\n')
                putative_barcodes, putative_umis, putative_position, anchor, anchor_umi, quality_filtered = self._extract_putatives(
                    upstream_tail,
                    upstream_qual)

                for barcode, umi, pos in zip(putative_barcodes, putative_umis, putative_position):
                    local_putative_lines.append(f'{record[0]}\t{barcode}\t{umi}\t{pos}\n')
                    local_putative_freq[barcode] += 1
                if anchor is not None:
                    local_anchor_lines.append(f'{record[0]}\t{anchor}\t{anchor_umi}\n')
                    local_anchor_freq[anchor] += 1
                    local_anchor += 1
                if quality_filtered > 0:
                    local_quality_filtered_reads += 1
                    local_quality_filtered_barcodes += quality_filtered
                if anchor is None and putative_barcodes:
                    local_only_putative_reads += 1

            else:
                local_no_polya += 1
            local_processed += 1

        return {
                'local_putative_lines': local_putative_lines,
                'local_anchor_lines': local_anchor_lines,
                'local_upstream_region_lines': local_upstream_region_lines,
                'local_putative_freq':local_putative_freq,
                'local_anchor_freq': local_anchor_freq,
                'local_polya': local_polya,
                'local_no_polya': local_no_polya,
                'local_anchor': local_anchor,
                'local_quality_filtered_reads': local_quality_filtered_reads,
                'local_quality_filtered_barcodes': local_quality_filtered_barcodes,
                'local_only_putative_reads': local_only_putative_reads,
                'local_processed': local_processed
                }

    @staticmethod
    def _avg_qual(qual_str: str):
        """
        Determine the average quality of a Phred encoded quality string

        Args:
            qual_str (str): Phred encoded quality string

        Returns:
            int: average quality of the string
        """
        total = 0
        for qchar in qual_str:
            total += ord(qchar) - 33
        return total / len(qual_str)

    def _extract_putatives(self, upstream_tail: str, upstream_qual: str):
        """Collect all putative barcode/UMI pairs and mark which are anchors.

        Args:
            upstream_tail (str): The sequence segment upstream of the polyT to search.
            upstream_qual (str): Quality string for upstream_tail.

        Returns:
            tuple: (putative_barcodes (list), putative_umis (list), is_anchors (list of bool))
        """
        putative_barcodes = []
        putative_umis = []
        putative_position = []
        anchor = None
        anchor_umi = None
        quality_filtered = 0
        upstream_tail_len = len(upstream_tail)

        expected_pos = upstream_tail_len - self.umi_len - self.barcode_len
        lower_search = expected_pos - self.search_range
        upper_search = expected_pos + self.search_range + 1

        for i in range(lower_search, upper_search):
            checker = upstream_tail[i:i + self.barcode_len]
            checker_qual = upstream_qual[i:i + self.barcode_len]
            umi = upstream_tail[i + self.barcode_len:i + self.barcode_len + self.umi_len]
            if len(checker) != self.barcode_len:
                quality_filtered += 1
                continue
            if self._avg_qual(checker_qual) <= self.min_qual:
                quality_filtered += 1
                continue
            if not self.exact_anchor:
                if checker in self.whitelist:
                    anchor = checker
                    anchor_umi = umi
            else:
                if i == expected_pos and checker in self.whitelist:
                    anchor = checker
                    anchor_umi = umi
            putative_barcodes.append(checker)
            putative_position.append(i - expected_pos)
            # add T's to the UMI when too short from running the window upstream and including the identified polyT in the UMI
            while len(umi) != self.umi_len:
                umi = umi + 'T'
            putative_umis.append(umi)

        return putative_barcodes, putative_umis, putative_position, anchor, anchor_umi, quality_filtered

    def _find_upstream(self, sequence: str, qual: str):
        """Return the subsequence upstream of the polyT using the processor's pattern and buffer.

        Args:
            sequence (str): Full read sequence.
            qual (str):     Quality string for the read

        Returns:
            str|None: Portion of the sequence upstream of the polyT (including buffer) or None.
            str|None: Quality of the sequence, or None.
            str|None: Sequence used for optional TSV output, or None.
        """
        reverse_seq = reverse_complement(sequence)
        reverse_qual = qual[::-1]

        forward = regex.findall(self.polyt_pattern, sequence[:200])
        reverse = regex.findall(self.polyt_pattern, reverse_seq[:200])

        if forward:
            polyt_start = sequence.find(forward[0])
            # include some of the polyT in case barcode/UMI ends in T
            upstream_tail = sequence[:polyt_start]
            upstream_qual = qual[:polyt_start]
            search_to_polyt_region = upstream_tail[-(self.barcode_len + self.umi_len + self.search_range):]

        elif reverse:
            polyt_start = reverse_seq.find(reverse[0])
            # include some of the polyT in case barcode/UMI ends in T
            upstream_tail = reverse_seq[:polyt_start]
            upstream_qual = reverse_qual[:polyt_start]
            search_to_polyt_region = upstream_tail[-(self.barcode_len + self.umi_len + self.search_range):]

        else:
            upstream_tail = None
            upstream_qual = None
            search_to_polyt_region = None

        return upstream_tail, upstream_qual, search_to_polyt_region


def _worker_init(polyt_min_length, barcode_len, umi_len, search_range, min_qual, whitelist, exact_anchor,
                 write_upstream_region):
    """
    Initializer executed in each worker process to set up a local processor.
    Using an initializer avoids pickling large shared objects per-task

    Args:
        polyt_min_length (int): The minimum size required for a polyT(A) to be identified in the read
        barcode_len (int):      The length of the cell barcode sequences
        umi_len (int):          The length of the UMI used
        search_range (int):     The range around the UMI length to look for exact matches
        min_qual (int):         The minimum average quality in a barcode for it to be considered
        whitelist (set):        The whitelist of cell barcodes
        exact_anchor (bool):    Whether or not to enforce anchors must be one per read and at the exact position as expected
        write_upstream_region (bool): Whether to emit trimmed upstream sequence lines for output
    """
    global _worker_processor
    _worker_processor = BarcodeBatchProcessor(polyt_min_length, barcode_len, umi_len, search_range,
                                              min_qual, whitelist, exact_anchor,
                                              write_upstream_region=write_upstream_region)


def _worker_process_batch(records_batch):
    """
    Module-level entry point used by ProcessPoolExecutor workers.

    Args:
        records_batch (list): A list of tuples (title, sequence, qual) for each read in the batch.

    Returns:
        list: Lines for putative barcodes to write to file
        list: Lines for anchor barcodes to write to file
        Counter: Frequency of each putative barcode in the batch
        Counter: Frequency of each anchor barcode in the batch
        int: Number of reads with polyA found in the batch
        int: Total number of reads processed in the batch
    """
    if _worker_processor is None:
        raise RuntimeError(
            'Worker not initialized. This function must be run inside a worker process with _worker_init called.')
    return _worker_processor.process(records_batch)


def _safe_future_result(future):
    """Return future.result() but ensure the exception is logged with traceback.

    Worker exceptions are logged here then re-raised so callers still fail-fast.
    KeyboardInterrupt/SystemExit propagate normally.
    """
    try:
        return future.result()
    except Exception:
        logger.exception('Worker task raised an exception (see traceback) — re-raising')
        raise


def extract_barcodes(fastq: str, polyt_min_length: int, barcode_len: int, umi_len: int,
                     search_range: int, min_qual: int, min_anchor_percentile: int, exact_anchor: bool,
                     whitelist: set = None, processes: int = None, batch_size: int = 10000,
                     upstream_region_tsv: str = None):
    """
    Given a fastq file, extract all the barcodes.

    Args:
        fastq (str):                 The path to the fastq file
        polyt_min_length (int):      The minimum size required for a polyT(A) to be identified in the read
        barcode_len (int):           The length of the cell barcode sequenceS
        umi_len (int):               The length of the UMI used
        search_range (int):          The range around the UMI length to look for exact matches
        min_qual (int):              The minimum average quality in a barcode for it to be considered
        min_anchor_percentile (int): The minimum number of exact matches an anchor requires to be kept
        exact_anchor (bool):         Whether or not to enforce anchors must be one per read and at the exact position as expected
        whitelist (set):             The whitelist of cell barcodes that are expected as a set for quick accession
        processes (int):             Number of worker processes (defaults to cpu_count() - 1)
        batch_size (int):            Number of reads per worker task (default 1000)
        upstream_region_tsv (str):   Output path for optional TSV of search-window-to-polyT sequences

    Returns:
        defaultdict: anchor barcode: frequency
        defaultdict: putative barcode: frequency
    """
    logger.section('BARCODE EXTRACTATION')
    logger.info('Beginning extraction')
    putative_anchors_freq = defaultdict(int)
    putative_barcodes_freq = defaultdict(int)
    polya_found_reads = 0
    no_polya_found_reads = 0
    reads_with_putative_anchor_count = 0
    processed = 0
    reads_with_quality_filtered = 0
    total_quality_filtered_barcodes = 0
    only_putative_barcodes = 0

    max_workers = processes if processes and processes > 0 else max(1, (os.cpu_count() or 2) - 1)
    logger.info('Using %d worker processes (batch_size=%d)', max_workers, batch_size)

    # stream through fastq file and submit batches to worker processes (decompression/streaming stays in main process)
    reader = gzip.open(fastq, 'rt') if fastq.endswith('gz') else open(fastq, 'rt')

    # track reads with multiple anchors
    reads_with_anchors = set()
    reads_with_multiple_anchors = 0

    # start a process pool — use an initializer so each worker creates its own BarcodeBatchProcessor
    executor_cm = ProcessPoolExecutor(max_workers=max_workers, initializer=_worker_init,
                                      initargs=(polyt_min_length, barcode_len, umi_len, search_range,
                                                min_qual, whitelist, exact_anchor, upstream_region_tsv is not None))

    optional_upstream_cm = open(upstream_region_tsv, 'w') if upstream_region_tsv is not None else nullcontext(None)

    with (open('.putative_barcodes', 'w') as putative_set_file, open('.anchors', 'w') as anchor_set_file,
          optional_upstream_cm as upstream_region_file, executor_cm as executor):
        if upstream_region_file is not None:
            upstream_region_file.write('read_id\tsearch_start_to_polyt\n')

        pending = []
        batch = []
        with reader as file:
            for record in FastqGeneralIterator(file):  # returns tuple: title, sequence, qual
                batch.append(record)
                if len(batch) >= batch_size:
                    pending.append(executor.submit(_worker_process_batch, batch))
                    batch = []

                # keep pending futures bounded and flush completed ones
                if len(pending) >= max_workers * 3:
                    done_future = next(as_completed(pending))
                    pending.remove(done_future)
                    res = _safe_future_result(done_future)
                    result_dict = res

                    for line in result_dict['local_putative_lines']:
                        putative_set_file.write(line)
                    for line in result_dict['local_anchor_lines']:
                        anchor_set_file.write(line)
                        read_id = line.split('\t')[0]
                        if read_id in reads_with_anchors:
                            reads_with_multiple_anchors += 1
                        reads_with_anchors.add(read_id)
                    if upstream_region_file is not None:
                        for line in result_dict['local_upstream_region_lines']:
                            upstream_region_file.write(line)
                    for k, v in result_dict['local_putative_freq'].items():
                        putative_barcodes_freq[k] += v
                    for k, v in result_dict['local_anchor_freq'].items():
                        putative_anchors_freq[k] += v

                    polya_found_reads += result_dict['local_polya']
                    no_polya_found_reads += result_dict['local_no_polya']
                    reads_with_putative_anchor_count += result_dict['local_anchor']
                    processed += result_dict['local_processed']
                    reads_with_quality_filtered += result_dict['local_quality_filtered_reads']
                    total_quality_filtered_barcodes += result_dict['local_quality_filtered_barcodes']
                    only_putative_barcodes += result_dict['local_only_putative_reads']
                    if processed % 500000 == 0:
                        logger.info('Extracted barcodes from %d reads', processed)

            # submit any remaining records
            if batch:
                pending.append(executor.submit(_worker_process_batch, batch))

            # collect remaining futures
            for future in as_completed(pending):
                res = _safe_future_result(future)

                result_dict = res

                for line in result_dict['local_putative_lines']:
                    putative_set_file.write(line)
                for line in result_dict['local_anchor_lines']:
                    anchor_set_file.write(line)
                    read_id = line.split('\t')[0]
                    if read_id in reads_with_anchors:
                        reads_with_multiple_anchors += 1
                    reads_with_anchors.add(read_id)
                if upstream_region_file is not None:
                    for line in result_dict['local_upstream_region_lines']:
                        upstream_region_file.write(line)
                for k, v in result_dict['local_putative_freq'].items():
                    putative_barcodes_freq[k] += v
                for k, v in result_dict['local_anchor_freq'].items():
                    putative_anchors_freq[k] += v

                polya_found_reads += result_dict['local_polya']
                no_polya_found_reads += result_dict['local_no_polya']
                reads_with_putative_anchor_count += result_dict['local_anchor']
                processed += result_dict['local_processed']
                reads_with_quality_filtered += result_dict['local_quality_filtered_reads']
                total_quality_filtered_barcodes += result_dict['local_quality_filtered_barcodes']
                only_putative_barcodes += result_dict['local_only_putative_reads']
                if processed % 500000 == 0:
                    logger.info('Extracted barcodes from %d reads', processed)

    # filter putative anchors to remove spurious matches
    logger.info('Number of putative anchor barcodes before filter: %d', len(putative_anchors_freq))
    anchor_freq_sizes = set(putative_anchors_freq.values())
    cut = np.floor(np.percentile(list(anchor_freq_sizes), min_anchor_percentile))
    logger.info(
        'The minimum number of exact anchors matches required is %d calculated as the %d percentile of the anchor'
            ' frequency set.',
        cut, min_anchor_percentile)
    putative_anchors_freq = {key: value for key, value in putative_anchors_freq.items() if value > cut}

    # log information
    logger.info('Total reads processed: %d', processed)
    logger.info('polyA+ reads found: %d', polya_found_reads)
    logger.info('polyA- reads found: %d', no_polya_found_reads)
    logger.info('Reads that contained quality filtered barcodes: %d', reads_with_quality_filtered)
    logger.info('Total quality filtered barcodes: %d', total_quality_filtered_barcodes)
    logger.info('polyA+ reads with anchors: %d', reads_with_putative_anchor_count)
    logger.info('polyA+ reads without anchors: %d', only_putative_barcodes)
    logger.info('Maximum recoverable reads: %d', reads_with_putative_anchor_count + only_putative_barcodes)
    logger.info('Number of anchor barcodes after filter: %d', len(putative_anchors_freq))
    logger.info('Number of unique putative barcodes: %d', len(putative_barcodes_freq))

    putative_anchors_freq = dict(sorted(putative_anchors_freq.items()))
    putative_barcodes_freq = dict(sorted(putative_barcodes_freq.items()))

    return putative_anchors_freq, putative_barcodes_freq

