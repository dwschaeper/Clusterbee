from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import reverse_complement
import gzip
import regex
from collections import defaultdict, Counter
from typing import TextIO
import matplotlib.pyplot as plt


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


#def find_supported_anchors(all_anchors: list[tuple[str, str]], used_umi: bool):
#    """
#    Finds realistic anchor barcodes from the set of exact match barcodes.
#
#    Args:
#        all_anchors (list): A list of all the barcodes
#        used_umi (bool):    Whether or not a UMI was used
#
#    Returns:
#        set:  IDs for the reads with anchors
#        list: The supported anchor barcodes
#        list: The UMIs for the supported anchor barcodes, only returned if used_umi=True
#    """
#    barcodes = defaultdict(int)
#    pass_barcodes = set()
#    final_barcodes = []
#    final_names = set()
#    final_umis = []

#    if used_umi:
#        for name, barcode, umi in all_anchors:
#            barcodes[barcode] += 1
#    else:
#        for name, barcode in all_anchors:
#            barcodes[barcode] += 1
#    
#    # plot barcode usage frequency
#    potential_anchor_freq = Counter(barcodes.values())
#    plt.figure(figsize=(8,5))
#    plt.scatter(potential_anchor_freq.keys(), potential_anchor_freq.values())
#    plt.xlabel('Barcodes supporting anchor')
#    plt.ylabel('Frequency')
#    plt.savefig('anchor_elbow_plot.png')
#   
#
#    for barcode, freq in barcodes.items():
#        if freq >= 40:  # simple cutoff for development purposes
#            pass_barcodes.add(barcode)
#
#    if used_umi:
#        for name, barcode, umi in all_anchors:
#            if barcode in pass_barcodes:
#                final_barcodes.append(barcode)
#                final_names.add(name)
#                final_umis.append(umi)
#
#        return final_names, final_barcodes, final_umis
#
#    else:
#        for name, barcode in all_anchors:
#            if barcode in pass_barcodes:
#                final_barcodes.append(barcode)
#                final_names.add(name)
#
#        return final_names, final_barcodes


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
    putative_anchors_freq = defaultdict(int)
    putative_barcodes_freq = defaultdict(int)
    polya_found_reads = 0
    reads_with_putative_anchor_count = 0
    
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
                            anchor_set_file.write(f'{record[0]}\t{barcode}\t{umi}\n')
                            putative_anchors_freq[anchor] += 1
                            reads_with_putative_anchor_count += 1

                        for barcode, umi in zip(putative_barcodes, putative_umis):
                            putative_set_file.write(f'{record[0]}\t{barcode}\t{umi}\n')
                            putative_barcodes_freq[barcode] += 1
    
    # THIS SHOULD BE DONE AT THE GRAPH LEVEL
    # identify which anchors are well supported and keep as final
    # anchor_names, anchor_barcodes, anchor_umis = find_supported_anchors(putative_anchors, True)
    
    # write out the final anchors and putative barcodes
    #with open('.putative_barcodes', 'w') as putative_set_file:
    #    with open('.anchors', 'w') as anchor_set_file:
    #        for name, barcode, umi in zip(anchor_names, anchor_barcodes, anchor_umis):
    #            anchor_barcodes_freq[barcode] += 1
    #            reads_with_anchors += 1
    #            anchor_set_file.write(f'{name}\t{barcode}\t{umi}\n')
    #        for name, barcodes, umis in all_putative_barcodes:
    #            if name not in anchor_names:
    #                reads_without_anchors += 1
    #                for barcode, umi in zip(barcodes, umis):
    #                    putative_barcodes_freq[barcode] += 1
    #                    putative_set_file.write(f'{name}\t{barcode}\t{umi}\n')

    print(f'polyA read+ reads found: {polya_found_reads}')            
    print(f'Reads with putative anchors: {reads_with_putative_anchor_count}')
    return putative_anchors_freq, putative_barcodes_freq
