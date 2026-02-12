from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import reverse_complement
import gzip
import regex
from collections import defaultdict, Counter
import numpy as np
import logging

logger = logging.getLogger(__name__)


def whitelist_match(whitelist: set, wiggle: int, umi_len: int, barcode_len: int, polyt_include: int, upstream_tail: str):
    """
    Use the provided whitelist to identify exact barcode matches in the 5' end of the reads prior to the polyT(A)
    By default this scans 7 windows.

    Args:
        whitelist (set):            The whitelist set of cell barcodes
        wiggle (int):               The range around the UMI length to look for exact matches
        umi_len (int):              The length of the UMI if used
        barcode_len (int):          The length of the cell barcode
        polyt_include (int):        How much of the polyT to include in the search
        upstream_tail (str):        The bases extracted from the read to search

    Returns:
        str:  The barcode anchor if valid, otherwise None
        str:  The UMI associated with the anchor
        list: The putative barcodes
        list: The putative UMIs
    """
    anchor = None
    anchor_umi = None
    putative_barcodes = []
    putative_umis = []
    lower_search = len(upstream_tail) - umi_len - wiggle - barcode_len
    upper_search = len(upstream_tail) - umi_len + wiggle + polyt_include - barcode_len + 1 # adjust for 0-based indexing

    for i in range(lower_search, upper_search):
        checker = upstream_tail[i:i + barcode_len]
        umi = upstream_tail[i + barcode_len:i + barcode_len + umi_len]
        if len(checker) != barcode_len:
            continue
        if checker in whitelist:
            anchor = checker
            anchor_umi = umi
        else:
            putative_barcodes.append(checker)
            putative_umis.append(umi)

    return anchor, anchor_umi, putative_barcodes, putative_umis


def find_upstream(sequence: str, polyt_include: int, polyt_min_threshold: int):
    """
    Extract the cell barcode from a single long read scRNA-seq sequencing read.

    Args:
        sequence (str): The full sequencing read, ensuring that it corresponds to the ID given during initialization
    """
    polyt_pattern = f'TTTT(?:T{{{polyt_min_threshold - 4},}}){{e<=1}}'
    reverse_seq = reverse_complement(sequence)

    forward = regex.findall(polyt_pattern, sequence[:200])
    reverse = regex.findall(polyt_pattern, reverse_seq[:200])

    if forward:
        polyt_start = sequence.find(forward[0])
        upstream_tail = sequence[:polyt_start + polyt_include]  # include some of the polyT in case barcode/UMI ends in T

    elif reverse:
        polyt_start = reverse_seq.find(reverse[0])
        upstream_tail = reverse_seq[:polyt_start + polyt_include]  # include some of the polyT in case barcode/UMI ends in T

    else:
        upstream_tail = None

    return upstream_tail


def extract_barcodes(fastq: str, polyt_min_threshold: int, barcode_len: int, umi_len: int, polyt_include: int, wiggle: int, whitelist: set=None):
    """
    Given a fastq file, extract all the barcodes.

    Args:
        fastq (str):               The path to the fastq file
        polyt_min_threshold (int): The minimum size required for a polyT(A) to be identified in the read
        barcode_len (int):         The length of the cell barcode sequenceS
        umi_len (int):             The length of the UMI used
        polyt_include (int):       How many bases of the polyT to include when searchig for UMI/barcode in case the sequence ends in T
        wiggle (int):              The range around the UMI length to look for exact matches
        whitelist (set):           The whitelist of cell barcodes that are expected as a set for quick accession

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
    
    # stream through fastq file
    reader = gzip.open(fastq, 'rt') if fastq.endswith('gz') else open(fastq, 'rt')
    
    with open('.putative_barcodes', 'w') as putative_set_file:
        with open('.anchors', 'w') as anchor_set_file:
            with reader as file:
                for record in FastqGeneralIterator(file):  # record is tuple of title, sequence, quality
                    upstream_tail = find_upstream(record[1], polyt_include, polyt_min_threshold)
            
                    # perform exact matching
                    if upstream_tail:
                        polya_found_reads += 1
                        anchor, anchor_umi, putative_barcodes, putative_umis = whitelist_match(whitelist, wiggle, umi_len, barcode_len, polyt_include, upstream_tail)
                        # append results to collection
                        if anchor:
                            anchor_set_file.write(f'{record[0]}\t{anchor}\t{anchor_umi}\n')
                            putative_anchors_freq[anchor] += 1
                            reads_with_putative_anchor_count += 1

                        for barcode, umi in zip(putative_barcodes, putative_umis):
                            putative_set_file.write(f'{record[0]}\t{barcode}\t{umi}\n')
                            putative_barcodes_freq[barcode] += 1
                    processed += 1
                    if processed % 500000 == 0:
                        logger.info('Extracted barcodes from %d reads', processed)
    
    logger.info('polyA read+ reads found: %d', polya_found_reads)
    logger.info('Reads with putative anchors: %d', reads_with_putative_anchor_count)
    logger.info('Reads with polyA but no anchor: %d', polya_found_reads - reads_with_putative_anchor_count)
    logger.info('Number of putative anchor barcodes: %d', len(putative_anchors_freq))
    logger.info('Number of putative barcodes: %d', len(putative_barcodes_freq))

    return putative_anchors_freq, putative_barcodes_freq
