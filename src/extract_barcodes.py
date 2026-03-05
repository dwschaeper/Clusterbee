from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import reverse_complement
import numpy as np
import gzip
import regex
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import logging

logger = logging.getLogger(__name__)
# per-process worker state and helpers (used by ProcessPoolExecutor initializer)
_worker_processor = None


class BarcodeBatchProcessor:
    """
    Worker class that processes batches of FastQ records.
    """

    def __init__(self, polyt_min_length: int, barcode_len: int, umi_len: int, polyt_buffer: int, search_range: int, min_qual: int):
        """
        Initialize the processor with parameters for barcode extraction.

        Args:
            polyt_min_length (int): The minimum size required for a polyT(A) to be identified in the read
            barcode_len (int):      The length of the cell barcode sequences
            umi_len (int):          The length of the UMI used
            polyt_buffer (int):     How many bases of the polyT to include when searching for UMI/barcode in case the sequence ends in T
            search_range (int):     The range around the UMI length to look for exact matches
            min_qual (int):         The minimum average quality in a barcode for it to be considered
        """
        self.barcode_len = barcode_len
        self.umi_len = umi_len
        self.polyt_buffer = polyt_buffer
        self.search_range = search_range
        self.min_qual = min_qual
        self.polyt_pattern = f'TTTT(?:T{{{polyt_min_length - 4},}}){{e<=1}}'

    def process(self, records_batch):
        """
        Process a batch of FastQ records to extract barcodes.

        Args:
            records_batch (list): A list of tuples (title, sequence) for each read in the batch.

        Returns:
            list: Lines for putative barcodes to write to file
            Counter: Frequency of each putative barcode in the batch
            int: Number of reads with polyA found in the batch
            int: Total number of reads processed in the batch
        """
        local_putative_lines = []
        local_putative_freq = Counter()
        local_polya = 0
        local_processed = 0

        for record in records_batch:
            upstream_tail, upstream_qual = self._find_upstream(record[1], record[2])
            if upstream_tail is not None:
                local_polya += 1
                putative_barcodes, putative_umis = self._extract_putatives(upstream_tail, upstream_qual)

                for barcode, umi in zip(putative_barcodes, putative_umis):
                    local_putative_lines.append(f'{record[0]}\t{barcode}\t{umi}\n')
                    local_putative_freq[barcode] += 1
            local_processed += 1

        return local_putative_lines, local_putative_freq, local_polya, local_processed
    
    def _avg_qual(self, qual_str: str):
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
        """Collect all putative barcode/UMI pairs using search and quality filters.

        Args:
            upstream_tail (str): The sequence segment upstream of the polyT to search.

        Returns:
            tuple: (putative_barcodes (list), putative_umis (list))
        """
        putative_barcodes = []
        putative_umis = []
        lower_search = len(upstream_tail) - self.umi_len - (2 * self.search_range) - (2 * self.polyt_buffer) - self.barcode_len
        upper_search = len(upstream_tail) - self.umi_len - self.barcode_len - 1  # adjust for 0-based indexing

        for i in range(lower_search, upper_search):
            checker = upstream_tail[i:i + self.barcode_len]
            checker_qual = upstream_qual[i:i + self.barcode_len]
            umi = upstream_tail[i + self.barcode_len:i + self.barcode_len + self.umi_len]
            if len(checker) != self.barcode_len:
                continue
            if self._avg_qual(checker_qual) <= self.min_qual:
                continue
            putative_barcodes.append(checker)
            putative_umis.append(umi)

        return putative_barcodes, putative_umis

    def _find_upstream(self, sequence: str, qual: str):
        """Return the subsequence upstream of the polyT using the processor's pattern and buffer.

        Args:
            sequence (str): Full read sequence.
            qual (str):     Quality string for the read

        Returns:
            str|None: Portion of the sequence upstream of the polyT (including buffer) or None.
            str|None: Quality of the sequence, or None
        """
        reverse_seq = reverse_complement(sequence)
        reverse_qual = qual[::-1]

        forward = regex.findall(self.polyt_pattern, sequence[:200])
        reverse = regex.findall(self.polyt_pattern, reverse_seq[:200])

        if forward:
            polyt_start = sequence.find(forward[0])
            upstream_tail = sequence[:polyt_start + self.search_range + self.polyt_buffer]  # include some of the polyT in case barcode/UMI ends in T
            upstream_qual = qual[:polyt_start + self.search_range + self.polyt_buffer]

        elif reverse:
            polyt_start = reverse_seq.find(reverse[0])
            upstream_tail = reverse_seq[:polyt_start + self.search_range + self.polyt_buffer]  # include some of the polyT in case barcode/UMI ends in T
            upstream_qual = reverse_qual[:polyt_start + self.search_range + self.polyt_buffer]

        else:
            upstream_tail = None
            upstream_qual = None

        return upstream_tail, upstream_qual


def _worker_init(polyt_min_length, barcode_len, umi_len, polyt_buffer, search_range, min_qual):
    """
    Initializer executed in each worker process to set up a local processor.
    Using an initializer avoids pickling large shared objects per-task

    Args:
        polyt_min_length (int): The minimum size required for a polyT(A) to be identified in the read
        barcode_len (int):      The length of the cell barcode sequences
        umi_len (int):          The length of the UMI used
        polyt_buffer (int):     How many bases of the polyT to include when searching for UMI/barcode in case the sequence ends in T
        search_range (int):     The range around the UMI length to look for exact matches
        min_qual (int):         The minimum average quality in a barcode for it to be considered
    """
    global _worker_processor
    _worker_processor = BarcodeBatchProcessor(polyt_min_length, barcode_len, umi_len, polyt_buffer, search_range, min_qual)


def _worker_process_batch(records_batch):
    """
    Module-level entry point used by ProcessPoolExecutor workers.

    Args:
        records_batch (list): A list of tuples (title, sequence, qual) for each read in the batch.

    Returns:
        list: Lines for putative barcodes to write to file
        Counter: Frequency of each putative barcode in the batch
        int: Number of reads with polyA found in the batch
        int: Total number of reads processed in the batch
    """
    if _worker_processor is None:
        raise RuntimeError('Worker not initialized. This function must be run inside a worker process with _worker_init called.')
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

def extract_barcodes(fastq: str, polyt_min_length: int, barcode_len: int, umi_len: int, polyt_buffer: int, search_range: int, min_qual: int, min_anchor_percentile: int, whitelist: set=None, processes: int=None, batch_size: int=10000):
    """
    Given a fastq file, extract all the barcodes.

    Args:
        fastq (str):                 The path to the fastq file
        polyt_min_length (int):      The minimum size required for a polyT(A) to be identified in the read
        barcode_len (int):           The length of the cell barcode sequenceS
        umi_len (int):               The length of the UMI used
        polyt_buffer (int):          How many bases of the polyT to include when searching for UMI/barcode in case the sequence ends in T
        search_range (int):          The range around the UMI length to look for exact matches
        min_qual (int):              The minimum average quality in a barcode for it to be considered
        min_anchor_percentile (int): The minimum number of exact matches an anchor requires to be kept
        whitelist (set):             The whitelist of cell barcodes that are expected as a set for quick accession
        processes (int):             Number of worker processes (defaults to cpu_count() - 1)
        batch_size (int):            Number of reads per worker task (default 1000)

    Returns:
        defaultdict: anchor barcode: frequency
        defaultdict: putative barcode: frequency
    """
    logger.info('Beginning extraction')
    putative_anchors_freq = defaultdict(int)
    putative_barcodes_freq = defaultdict(int)
    polya_found_reads = 0
    reads_with_putative_anchor_count = 0
    processed = 0

    max_workers = processes if processes and processes > 0 else max(1, (os.cpu_count() or 2) - 1)
    logger.info('Using %d worker processes (batch_size=%d)', max_workers, batch_size)

    # stream through fastq file and submit batches to worker processes (decompression/streaming stays in main process)
    reader = gzip.open(fastq, 'rt') if fastq.endswith('gz') else open(fastq, 'rt')

    # track reads with multiple anchors
    reads_with_anchors = set()
    reads_with_multiple_anchors = 0

    # start a process pool — use an initializer so each worker creates its own BarcodeBatchProcessor
    executor_cm = ProcessPoolExecutor(max_workers=max_workers, initializer=_worker_init,
                                      initargs=(polyt_min_length, barcode_len, umi_len, polyt_buffer, search_range, min_qual))

    with open('.putative_barcodes', 'w') as putative_set_file, open('.anchors', 'w') as anchor_set_file, executor_cm as executor:
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
                    putative_lines, putative_cnt, local_polya, local_processed = res

                    local_anchor_reads = set()
                    for line in putative_lines:
                        putative_set_file.write(line)
                        read_id, barcode, _ = line.rstrip('\n').split('\t')
                        if barcode in whitelist:
                            anchor_set_file.write(line)
                            local_anchor_reads.add(read_id)
                            if read_id in reads_with_anchors:
                                reads_with_multiple_anchors += 1
                            reads_with_anchors.add(read_id)
                    for k, v in putative_cnt.items():
                        putative_barcodes_freq[k] += v
                        if k in whitelist:
                            putative_anchors_freq[k] += v

                    polya_found_reads += local_polya
                    reads_with_putative_anchor_count += len(local_anchor_reads)
                    processed += local_processed
                    if processed % 500000 == 0:
                        logger.info('Extracted barcodes from %d reads', processed)

            # submit any remaining records
            if batch:
                pending.append(executor.submit(_worker_process_batch, batch))

            # collect remaining futures
            for future in as_completed(pending):
                res = _safe_future_result(future)

                putative_lines, putative_cnt, local_polya, local_processed = res

                local_anchor_reads = set()
                for line in putative_lines:
                    putative_set_file.write(line)
                    read_id, barcode, _ = line.rstrip('\n').split('\t')
                    if barcode in whitelist:
                        anchor_set_file.write(line)
                        local_anchor_reads.add(read_id)
                        if read_id in reads_with_anchors:
                            reads_with_multiple_anchors += 1
                        reads_with_anchors.add(read_id)
                for k, v in putative_cnt.items():
                    putative_barcodes_freq[k] += v
                    if k in whitelist:
                        putative_anchors_freq[k] += v

                polya_found_reads += local_polya
                reads_with_putative_anchor_count += len(local_anchor_reads)
                processed += local_processed
                if processed % 500000 == 0:
                    logger.info('Extracted barcodes from %d reads', processed)

    # filter putative anchors to remove spurious matches
    anchor_freq_sizes = set(putative_anchors_freq.values())
    cut = np.floor(np.percentile(list(anchor_freq_sizes), min_anchor_percentile))
    logger.info('The minimum number of exact anchors matches required is %d calculated as the %d percentile of the anchor frequncy set.', cut, min_anchor_percentile)
    putative_anchors_freq = {key: value for key, value in putative_anchors_freq.items() if value > cut}

    logger.info('polyA read+ reads found: %d', polya_found_reads)
    logger.info('Reads with putative anchors: %d', reads_with_putative_anchor_count)
    logger.info('Reads with multiple anchors: %d', reads_with_multiple_anchors)
    logger.info('Reads with polyA but no anchor: %d', polya_found_reads - reads_with_putative_anchor_count)
    logger.info('Number of putative anchor barcodes: %d', len(putative_anchors_freq))
    logger.info('Number of putative barcodes: %d', len(putative_barcodes_freq))

    return putative_anchors_freq, putative_barcodes_freq
