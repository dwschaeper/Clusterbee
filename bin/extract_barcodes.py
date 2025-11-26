from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import reverse_complement
import gzip
import regex


class BarcodeExtractor:
    def __init__(self, name: str, polyt_min_threshold: int, barcode_len: int, umi_len: int):
        """
        Initialize the BarcodeExtractor class. This class functions to extract the barcode from a sequence and
        maintain easily accessed attributes.

        Args:
            name (str): The ID of the sequence read
            polyt_min_threshold (int): The minimum size required for a polyT(A) to be identified in the read
            barcode_len (int): The length of the cell barcode
        """
        self.name = name
        self.polyt_min_threshold = polyt_min_threshold
        self.barcode_len = barcode_len
        self.umi_len = umi_len
        self.polyt_start = 0
        self.polyt_end = 0
        self.strand = ''
        self.upstream_tail = ''
        self.barcode = ''
        self.barcode_start = -1
        self.distance_from_polyt = -100

    def whitelist_match(self, whitelist: set, wiggle: int):
        """
        Use the provided whitelist to identify exact barcode matches in the 5' end of the reads prior to the polyT(A)
        By default this scans 7 windows.

        Args:
            wiggle (int): The range around the UMI length to look for exact matches
            whitelist (set): The whitelist set of cell barcodes
        """

        lower_search = len(self.upstream_tail) - self.umi_len - wiggle - self.barcode_len  
        upper_search = len(self.upstream_tail) - self.umi_len + wiggle + 2 - self.barcode_len + 1  # adjust for 0-based and the case barcode/UMI ends in T and is captured in polyT

        for i in range(lower_search, upper_search):
            checker = self.upstream_tail[i:i + self.barcode_len]
            if checker in whitelist:
                self.barcode = checker
                self.barcode_start = i
                self.distance_from_polyt = len(self.upstream_tail) - (i + self.barcode_len)
                break
        #for i in range():

    def find_upstream(self, sequence):
        """
        Extract the cell barcode from a single long read scRNA-seq sequencing read.

        Args:
            sequence (str): The full sequencing read, ensuring that it corresponds to the ID given during initialization
        """
        polyt_pattern = f'T(?:T{{{self.polyt_min_threshold - 1},}}){{e<=1}}'
        reverse_seq = reverse_complement(sequence)

        forward = regex.findall(polyt_pattern, sequence[:200])
        reverse = regex.findall(polyt_pattern, reverse_seq[:200])

        if forward:
            self.strand = '+'
            self.polyt_start = sequence.find(forward[0])
            self.polyt_end = self.polyt_start + len(forward[0])
            self.upstream_tail = sequence[:self.polyt_start]  # include some of the polyT in case barcode ends in T

        elif reverse:
            self.strand = '-'
            self.polyt_start = reverse_seq.find(reverse[0])
            self.polyt_end = self.polyt_start + len(reverse[0])
            self.upstream_tail = reverse_seq[:self.polyt_start]  # include some of the polyT in case barcode ends in T

        else:
            self.strand, self.polyt_start, self.polyt_end, self.barcode = None, None, None, None


def extract_barcodes(fastq: str, polyt_min_threshold: int, barcode_len: int, umi_len: int, wiggle: int=None, whitelist: set=None):
    """
    Given a fastq file, extract all the barcodes.

    Args:
        fastq (str): The path to the fastq file
        polyt_min_threshold (int): The minimum size required for a polyT(A) to be identified in the read
        barcode_len (int): The length of the cell barcode sequence
        umi_len (int): The length of the UMI used
        wiggle (int): The range around the UMI length to look for exact matches
        whitelist (set): The whitelist of cell barcodes that are expected as a set for quick accession

    Returns:
        dictionary: The read names and barcodes. Unique barcodes are the keys, the reads
                    that have that barcode are held in a list as the value
    """
    polyt_count = 0
    polya_count = 0
    missing_count = 0
    barcode_dict = {}
    barcode_start = {}
    distance_from_polyt = {}
    #barcode_usage_per_location = {}

    reader = gzip.open(fastq, 'rt') if fastq.endswith('gz') else open(fastq, 'rt')

    with reader as file:
        for record in FastqGeneralIterator(file):  # tuple of title, sequence, quality
            extractor = BarcodeExtractor(record[0], polyt_min_threshold, barcode_len, umi_len)
            extractor.find_upstream(record[1])
            
            if whitelist and extractor.strand:
                extractor.whitelist_match(whitelist, wiggle)
            if extractor.strand == '+':
                polyt_count += 1
                barcode_dict[extractor.barcode] = barcode_dict.get(extractor.barcode, 0) + 1
                barcode_start[extractor.barcode_start] = barcode_start.get(extractor.barcode_start, 0) + 1
                distance_from_polyt[extractor.distance_from_polyt] = distance_from_polyt.get(extractor.distance_from_polyt, 0) + 1
                #barcode_usage_per_location[extractor.barcode_start] = barcode_usage_per_location.get(extractor.barcode, set()).add(extractor.barcode)
            elif extractor.strand == '-':
                polya_count += 1
                barcode_dict[extractor.barcode] = barcode_dict.get(extractor.barcode, 0) + 1
                barcode_start[extractor.barcode_start] = barcode_start.get(extractor.barcode_start, 0) + 1
                distance_from_polyt[extractor.distance_from_polyt] = distance_from_polyt.get(extractor.distance_from_polyt, 0) + 1
            else:
                missing_count += 1
                barcode_dict[extractor.barcode] = barcode_dict.get(extractor.barcode, 0) + 1

    print(f'poly(T)+ reads: {polyt_count}')
    print(f'poly(A)+ reads: {polya_count}')
    print(f'No polyT(A): {missing_count}')
    print(f"Empty barcodes: {barcode_dict['']}")
    sorted_start = dict(sorted(barcode_start.items()))
    for key, value in sorted_start.items():
        print(f'{key}: total barcodes: {value}')
    
    sorted_dist = dict(sorted(distance_from_polyt.items()))
    for key, value in sorted_dist.items():
        print(f'{key}: distance from polyt: {value}')
    #print(sum(sorted_start.values()))
    #print(barcode_dict)
    #del barcode_dict['']
    #sorted_vals = sorted(barcode_dict.values())
    #print(f'Top 20 barcode counts: {sorted_vals[-20:]}')
