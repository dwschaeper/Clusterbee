from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import reverse_complement
import gzip
import regex

import matplotlib.pyplot as plt


translate = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


class BarcodeExtractor:
    def __init__(self, readID: str, polyt_min_threshold: int, barcode_len: int, umi_len: int, polyt_include: int, no_umi: bool):
        """
        Initialize the BarcodeExtractor class. This class functions to extract the barcode from a sequence and
        maintain easily accessed attributes.

        Args:
            name (str):                The ID of the sequence read
            polyt_min_threshold (int): The minimum size required for a polyT(A) to be identified in the read
            barcode_len (int):         The length of the cell barcode
            umi_len (int):             The length of the UMI sequence
            polyt_include (int):       The number of nucleotides labeled within the polyT to include in the search section
            no_umi (bool):             True if there are no UMIs used
        """
        self.readID = readID
        self.polyt_min_threshold = polyt_min_threshold
        self.barcode_len = barcode_len
        self.umi_len = umi_len
        self.polyt_include = polyt_include
        self.no_umi = no_umi
        self.polyt_start = 0
        self.polyt_end = 0
        self.strand = ''
        self.upstream_tail = ''
        self.putative_barcodes = []
        self.barcode = ''
        self.barcod_encoding = ''
        self.umi = ''
        self.barcode_start = -1
        self.distance_from_polyt = -100


    def whitelist_match(self, whitelist: set, wiggle: int):
        """
        Use the provided whitelist to identify exact barcode matches in the 5' end of the reads prior to the polyT(A)
        By default this scans 7 windows.

        Args:
            wiggle (int):    The range around the UMI length to look for exact matches
            whitelist (set): The whitelist set of cell barcodes
        """
        lower_search = len(self.upstream_tail) - self.umi_len - wiggle - self.polyt_include - self.barcode_len
        upper_search = len(self.upstream_tail) - self.umi_len + wiggle - self.polyt_include - self.barcode_len + 1 # adjust for 0-based indexing
        if self.no_umi:
            for i in range(lower_search, upper_search):
                checker = self.upstream_tail[i:i + self.barcode_len]
                if checker in whitelist:
                    self.barcode = checker
                    self.barcode_start = i
                    self.distance_from_polyt = len(self.upstream_tail) - (i + self.barcode_len)
                    self.putative_barcode = []
                else:
                    self.putative_barcodes.append(checker)
        else:
            for i in range(lower_search, upper_search):
                checker = self.upstream_tail[i:i + self.barcode_len]
                if checker in whitelist:
                    self.barcode = checker
                    self.barcode_start = i
                    self.distance_from_polyt = len(self.upstream_tail) - (i + self.barcode_len)
                    self.putative_barcodes = []
                    self.umi = self.upstream_tail[i + self.barcode_len:i + self.barcode_len + self.umi_len]
                    break
                else:
                    self.putative_barcodesappend(checker)


    def find_upstream(self, sequence):
        """
        Extract the cell barcode from a single long read scRNA-seq sequencing read.

        Args:
            sequence (str): The full sequencing read, ensuring that it corresponds to the ID given during initialization
        """
        polyt_pattern = f'TTTT(?:T{{{self.polyt_min_threshold - 4},}}){{e<=1}}'
        reverse_seq = reverse_complement(sequence)

        forward = regex.findall(polyt_pattern, sequence[:200])
        reverse = regex.findall(polyt_pattern, reverse_seq[:200])

        if forward:
            self.strand = '+'
            self.polyt_start = sequence.find(forward[0])
            self.polyt_end = self.polyt_start + len(forward[0])
            self.upstream_tail = sequence[:self.polyt_start + self.polyt_include]  # include some of the polyT in case barcode/UMI ends in T

        elif reverse:
            self.strand = '-'
            self.polyt_start = reverse_seq.find(reverse[0])
            self.polyt_end = self.polyt_start + len(reverse[0])
            self.upstream_tail = reverse_seq[:self.polyt_start + self.polyt_include]  # include some of the polyT in case barcode/UMI ends in T

        else:
            self.strand, self.polyt_start, self.polyt_end, self.barcode = None, None, None, None


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
            extractor = BarcodeExtractor(record[0], polyt_min_threshold, barcode_len, umi_len, polyt_include, no_umi)
            extractor.find_upstream(record[1])
            
            # perform exact matching if a whitelist is provided
            if whitelist and extractor.strand:
                extractor.whitelist_match(whitelist, wiggle)

            # don't hold on to reads that don't have a polyT(A)
            if extractor.polyt_start == None:  
                continue

            extractor_objects.append(extractor)

            # quick counters for development sanity checks
            if extractor.strand == '+':
                polyt_count += 1
                #barcode_dict[extractor.barcode] = barcode_dict.get(extractor.barcode, 0) + 1
                #barcode_start[extractor.barcode_start] = barcode_start.get(extractor.barcode_start, 0) + 1
                #distance_from_polyt[extractor.distance_from_polyt] = distance_from_polyt.get(extractor.distance_from_polyt, 0) + 1
                #barcode_usage_per_location[extractor.barcode_start] = barcode_usage_per_location.get(extractor.barcode, set()).add(extractor.barcode)
            elif extractor.strand == '-':
                polya_count += 1
                #barcode_dict[extractor.barcode] = barcode_dict.get(extractor.barcode, 0) + 1
                #barcode_start[extractor.barcode_start] = barcode_start.get(extractor.barcode_start, 0) + 1
                #distance_from_polyt[extractor.distance_from_polyt] = distance_from_polyt.get(extractor.distance_from_polyt, 0) + 1
            else:
                missing_count += 1
                #barcode_dict[extractor.barcode] = barcode_dict.get(extractor.barcode, 0) + 1

    #print(f'poly(T)+ reads: {polyt_count}')
    #print(f'poly(A)+ reads: {polya_count}')
    #print(f'No polyT(A): {missing_count}')
    #print(f"Empty barcodes: {barcode_dict['']}")
    
    #print(barcode_dict)
    #del barcode_dict[None]
    #del barcode_dict['']
    #sorted_vals = sorted(list(barcode_dict.values()))
    #print(f'Number of barcodes: {len(sorted_vals)}')
    #print(f'Number of reads assigned to barcodes: {sum(sorted_vals)}')
    #print(f'Top 20 barcode counts: {sorted_vals[-20:]}')
    #print(f'Average barcode count: {sum(sorted_vals)/len(sorted_vals)}')
    #print(f'Number of barcodes with less than 5 reads: {sum([i for i in sorted_vals if i < 5 ])}')
    #print(f'Number of barcodes with less than 4 reads: {sum([i for i in sorted_vals if i < 4 ])}')
    #print(f'Number of barcodes with less than 3 reads: {sum([i for i in sorted_vals if i < 3 ])}')
    #print(f'Number of barcodes with less than 2 reads: {sum([i for i in sorted_vals if i < 2 ])}')

    #sorted_dict = {}
    #for key in sorted(barcode_dict, key=barcode_dict.get):
    #    sorted_dict[key] = barcode_dict[key]
    
    #print('Barcode: Reads supporting')
    #for key, value in sorted_dict.items():
    #    print(key, ': ', value)

	
    # plots
    #fig, axes = plt.subplots(1,3, figsize=(12,4))
    
    # histogram
    #axes[0].hist(sorted_vals, bins=30)
    #axes[0].set_title('Histogram of Barcode Counts')

    # density
    #axes[1].hist(sorted_vals, bins=30, density=True)
    #axes[1].set_title('Density of Barcode Counts')

    # violin
    #axes[2].violinplot(sorted_vals, showmeans=True)
    #axes[2].set_title('Violin of Barcode Counts')

    #plt.tight_layout()
    #plt.savefig('exact_match_distribution_plots.pdf')
    #plt.close()
    

    return extractor_objects
