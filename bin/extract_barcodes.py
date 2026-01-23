from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import reverse_complement
import gzip
import regex

import matplotlib.pyplot as plt


translate = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


class BarcodeExtractor:
    def __init__(self, readID: str):
        """
        Initialize the BarcodeExtractor class. This class functions to extract the barcode from a sequence and
        maintain easily accessed attributes.

        Args:
            name (str): The ID of the sequence read
        """
        self.readID = readID
        self.putative_barcodes = []
        self.barcode = ''
        self.umi = ''


    def whitelist_match(self, whitelist: set, wiggle: int, umi_len: int, barcode_len: int, polyt_include: int, no_umi: bool, upstream_tail: str):
        """
        Use the provided whitelist to identify exact barcode matches in the 5' end of the reads prior to the polyT(A)
        By default this scans 7 windows.

        Args:
            wiggle (int):    The range around the UMI length to look for exact matches
            whitelist (set): The whitelist set of cell barcodes
        """
        lower_search = len(upstream_tail) - umi_len - wiggle - polyt_include - barcode_len
        upper_search = len(upstream_tail) - umi_len + wiggle - polyt_include - barcode_len + 1 # adjust for 0-based indexing
        if no_umi:
            for i in range(lower_search, upper_search):
                checker = upstream_tail[i:i + barcode_len]
                if checker in whitelist:
                    self.barcode = checker
                    distance_from_polyt = len(upstream_tail) - (i + barcode_len)
                    self.putative_barcode = []
                    break
                else:
                    self.putative_barcodes.append(checker)
        else:
            for i in range(lower_search, upper_search):
                checker = upstream_tail[i:i + barcode_len]
                if checker in whitelist:
                    self.barcode = checker
                    distance_from_polyt = len(upstream_tail) - (i + barcode_len)
                    self.putative_barcodes = []
                    self.umi = upstream_tail[i + barcode_len:i + barcode_len + umi_len]
                    break
                else:
                    self.putative_barcodes.append(checker)


    def find_upstream(self, sequence: str, polyt_include: int, polyt_min_threshold: int):
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
            self.barcode = None
            upstream_tail = None

        return upstream_tail


def extract_barcodes(fastq: str, polyt_min_threshold: int, barcode_len: int, no_umi: bool, umi_len: int, polyt_include: int, wiggle: int=None, whitelist: set=None):
    """
    Given a fastq file, extract all the barcodes.

    Args:
        fastq (str):               The path to the fastq file
        polyt_min_threshold (int): The minimum size required for a polyT(A) to be identified in the read
        barcode_len (int):         The length of the cell barcode sequence
        no_umi (bool):             Whether or not UMIs are used
        umi_len (int):             The length of the UMI used
        wiggle (int):              The range around the UMI length to look for exact matches
        polyt_include (int):       How many bases of the polyT to include when searchig for UMI/barcode in case the sequence ends in T
        whitelist (set):           The whitelist of cell barcodes that are expected as a set for quick accession

    Returns:
        list: A list of BarcodeExtractor object. Each entry in the list corresponds to a singl read
    """
    polyt_count = 0
    polya_count = 0
    missing_count = 0
    extractor_objects = []
    barcode_dict = {}
    barcode_start = {}
    distance_from_polyt = {}
    #barcode_usage_per_location = {}

    reader = gzip.open(fastq, 'rt') if fastq.endswith('gz') else open(fastq, 'rt')

    with reader as file:
        for record in FastqGeneralIterator(file):  # record is tuple of title, sequence, quality
            extractor = BarcodeExtractor(record[0])
            upstream_tail = extractor.find_upstream(record[1], polyt_include, polyt_min_threshold)
            
            # perform exact matching if a whitelist is provided
            if whitelist and upstream_tail:
                extractor.whitelist_match(whitelist, wiggle, umi_len, barcode_len, polyt_include, no_umi, upstream_tail)
                extractor_objects.append(extractor)
   
    return extractor_objects
