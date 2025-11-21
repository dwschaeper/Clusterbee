from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip


def find_barcode(seq: str):
    """
    Extract the cell barcode from a single long read scRNA-seq sequencing read.

    Args:
        seq (str): The full sequencing read
    
    Returns:
        str | None: The barcode for the read, or None if there is not one identified
    """
    pass


def extract_barcodes(fastq: str):
    """
    Given a fastq file, extract all of the barcodes.

    Args:
        fastq (str): The path to the fastq file

    Returns:
        dictionary: The read names and barcodes. Unique barcodes are the keys, the reads
                    that have that barcode are held in a list as the value
    """
    if fastq.endswith('.gz'):
        with gzip.open(fastq, 'rt') as file:
            for record in FastqGeneralIterator(file):  # tuple of title, sequence, quality
                barcode = find_barcode(record[1])
    else:
        with open(fastq, 'rt') as file:
            for record in FastqGeneralIterator(file):  # tuple of title, sequence, quality
                barcode = find_barcode(record[1])
