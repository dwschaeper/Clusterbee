from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import reverse_complement
import gzip
import regex


class BarcodeExtractor:
    def __init__(self, name: str, polyt_min_threshold: int, barcode_len: int):
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
        self.polyt_start = 0
        self.polyt_end = 0
        self.strand = ''
        self.barcode = ''

    def find_barcode(self, sequence):
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
            self.barcode = sequence[self.polyt_start - self.barcode_len: self.polyt_start]

        elif reverse:
            self.strand = '-'
            self.polyt_start = reverse_seq.find(reverse[0])
            self.polyt_end = self.polyt_start + len(reverse[0])
            self.barcode = reverse_seq[self.polyt_start - self.barcode_len: self.polyt_start]

        else:
            self.strand, self.polyt_start, self.polyt_end, self.barcode = None, None, None, None


def extract_barcodes(fastq: str, polyt_min_threshold: int, barcode_len: int):
    """
    Given a fastq file, extract all the barcodes.

    Args:
        fastq (str): The path to the fastq file
        polyt_min_threshold (int): The minimum size required for a polyT(A) to be identified in the read
        barcode_len (int): The length of the cell barcode sequence

    Returns:
        dictionary: The read names and barcodes. Unique barcodes are the keys, the reads
                    that have that barcode are held in a list as the value
    """
    polyt_count = 0
    polya_count = 0
    missing_count = 0
    barcode_dict = {}

    reader = gzip.open(fastq, 'rt') if fastq.endswith('gz') else open(fastq, 'rt')

    with reader as file:
        for record in FastqGeneralIterator(file):  # tuple of title, sequence, quality
            extractor = BarcodeExtractor(record[0], polyt_min_threshold, barcode_len)
            extractor.find_barcode(record[1])
            if extractor.strand == '+':
                polyt_count += 1
                barcode_dict[extractor.barcode] = barcode_dict.get(extractor.barcode, 0) + 1
            elif extractor.strand == '-':
                polya_count += 1
                barcode_dict[extractor.barcode] = barcode_dict.get(extractor.barcode, 0) + 1
            else:
                missing_count += 1

    print(polyt_count)
    print(polya_count)
    print(missing_count)
    print(sorted(barcode_dict.values())[:-20])
