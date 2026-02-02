from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import reverse_complement
import gzip
import regex
from collections import defaultdict, Counter
from typing import TextIO
import matplotlib.pyplot as plt


def whitelist_match(readID: str, whitelist: set, wiggle: int, umi_len: int, barcode_len: int, polyt_include: int, no_umi: bool, upstream_tail: str):
    """
    Use the provided whitelist to identify exact barcode matches in the 5' end of the reads prior to the polyT(A)
    By default this scans 7 windows.

    Args:
        readID (str):               The ID of the read
        whitelist (set):            The whitelist set of cell barcodes
        wiggle (int):               The range around the UMI length to look for exact matches
        umi_len (int):              The length of the UMI if used
        barcode_len (int):          The length of the cell barcode
        polyt_include (int):        How much of the polyT to include in the search
        no_umi (bool):              Whether or not a UMI was used
        upstream_tail (str):        The bases extracted from the read to search

    Returns:
        str:  The barcode anchor if valid, otherwise None
        list: The putative barcodes if valid, otherwise None
    """
    anchor = None
    umi = None
    putative_barcodes = []
    putative_umis = []
    lower_search = len(upstream_tail) - umi_len - wiggle - polyt_include - barcode_len
    upper_search = len(upstream_tail) - umi_len + wiggle + polyt_include - barcode_len + 1 # adjust for 0-based indexing

    if no_umi:
        for i in range(lower_search, upper_search):
            checker = upstream_tail[i:i + barcode_len]
            if checker in whitelist:
                anchor = checker
                break
            else:
                putative_barcodes.append(checker)

        #if putative_barcodes:        
        #    for barcode in putative_barodes:
        #        if len(barcode) == barcode_len:
        #            putative_set.write(f'{readID}\t{barcode}\n')

    else:
        for i in range(lower_search, upper_search):
            checker = upstream_tail[i:i + barcode_len]
            umi = upstream_tail[i + barcode_len:i + barcode_len + umi_len]
            if checker in whitelist:
                anchor = checker
                # anchor_umi_holder.write(f'{readID}\t{anchor}\t{umi}\n')
            else:
                putative_barcodes.append(checker)
                putative_umis.append(umi)
        
        #if putative_barcodes:
        #    for barcode, umi in zip(putative_barcodes, putative_umis):
        #        if len(barcode) == barcode_len and len(umi) == umi_len:
        #            putative_set.write(f'{readID}\t{barcode}\t{umi}\n')

    return anchor, putative_barcodes


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


def find_supported_anchors(all_anchors: list[tuple[str, str]]):
    """
    Finds realistic anchor barcodes from the set of exact match barcodes.

    Args:
        all_anchors (list): A list of all the barcodes
    """
    barcodes = defaultdict(int)
    pass_barcodes = set()
    final_barcodes = []
    final_names = set()

    for name, barcode in all_anchors:
        barcodes[barcode] += 1
    
    # plot barcode usage frequency
    potential_anchor_freq = Counter(barcodes.values())
    plt.figure(figsize=(8,5))
    plt.scatter(potential_anchor_freq.keys(), potential_anchor_freq.values())
    plt.xlabel('Barcodes supporting anchor')
    plt.ylabel('Frequency')
    plt.savefig('anchor_elbow_plot.png')
    

    for barcode, freq in barcodes.items():
        if freq >= 35:  # simple cutoff for development purposes
            pass_barcodes.add(barcode)

    for name, barcode in all_anchors:
        if barcode in pass_barcodes:
            final_barcodes.append(barcode)
            final_names.add(name)

    return final_names, final_barcodes


def extract_barcodes(fastq: str, polyt_min_threshold: int, barcode_len: int, no_umi: bool, umi_len: int, polyt_include: int, wiggle: int, whitelist: set=None):
    """
    Given a fastq file, extract all the barcodes.

    Args:
        fastq (str):               The path to the fastq file
        polyt_min_threshold (int): The minimum size required for a polyT(A) to be identified in the read
        barcode_len (int):         The length of the cell barcode sequence
        no_umi (bool):             Whether or not UMIs are used
        umi_len (int):             The length of the UMI used
        polyt_include (int):       How many bases of the polyT to include when searchig for UMI/barcode in case the sequence ends in T
        wiggle (int):              The range around the UMI length to look for exact matches
        whitelist (set):           The whitelist of cell barcodes that are expected as a set for quick accession

    Returns:
        defaultdict: anchor barcode: frequency
        defaultdict: putative barcode: frequency
    """
    anchor_barcodes_freq = defaultdict(int)
    putative_barcodes_freq = defaultdict(int)
    putative_anchors = []
    all_putative_barcodes = []
    reads_with_anchors = 0
    reads_without_anchors = 0

    reader = gzip.open(fastq, 'rt') if fastq.endswith('gz') else open(fastq, 'rt')
    #with open('.putative_barcodes', 'w') as putative_set_file:
    #    with open('.anchor_umis', 'w') as umi_set_file:
    with reader as file:
        for record in FastqGeneralIterator(file):  # record is tuple of title, sequence, quality
            upstream_tail = find_upstream(record[1], polyt_include, polyt_min_threshold)
            
            # perform exact matching if a whitelist is provided
            if whitelist and upstream_tail:
                anchor, putatives = whitelist_match(record[0], whitelist, wiggle, umi_len, barcode_len, polyt_include, no_umi, upstream_tail)
                if anchor:
                    putative_anchors.append((record[0], anchor))
                all_putative_barcodes.append((record[0], putatives))

    anchor_names, anchor_barcodes = find_supported_anchors(putative_anchors)
    
    # write out the anchors and putative barcodes
    with open('.putative_barcodes', 'w') as putative_set_file:
        with open('.anchors', 'w') as anchor_set_file:
            for name, barcode in zip(anchor_names, anchor_barcodes):
                anchor_barcodes_freq[barcode] += 1
                reads_with_anchors += 1
                anchor_set_file.write(f'{name}\t{barcode}')
            for name, barcodes in all_putative_barcodes:
                if name not in anchor_names:
                    reads_without_anchors += 1
                    for barcode in barcodes:
                        putative_barcodes_freq[barcode] += 1
                        putative_set_file.write(f'{name}\t{barcode}')
                
            
                #if anchor:
                #    reads_with_anchors += 1
                #    anchor_barcodes[anchor] += 1
                #if putatives:
                #    reads_without_anchors += 1
                #    for putative in putatives:
                #        putative_barcodes[putative] += 1
    print(f'Reads with anchors: {reads_with_anchors}')
    print(f'Reads without anchors: {reads_without_anchors}')
    return anchor_barcodes_freq, putative_barcodes_freq
