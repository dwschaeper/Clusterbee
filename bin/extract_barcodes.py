from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import reverse_complement
import gzip
import regex


def find_barcode(seq: str, polyT_min_threshold: int, barcode_len: int):
    """
    Extract the cell barcode from a single long read scRNA-seq sequencing read.

    Args:
        seq (str): The full sequencing read
        polyT_min_threshold (int): The minimum size required for a polyT(A) to be identified in the read
    
    Returns:
        tuple : The strand for the read and barcode, or empty tuple if there is not one identified
    """
    polyT_pattern = f'T(?:T{{{polyT_min_threshold - 1},}}){{e<=1}}'
    reverse_seq = reverse_complement(seq)
    
    forward = regex.findall(polyT_pattern, seq[:200])
    reverse = regex.findall(polyT_pattern, reverse_seq[:200])
    
    if forward:
        polyT_start = seq.find(forward[0])
        polyT_end = polyT_start + len(forward[0])
        barcode = seq[polyT_start - barcode_len: polyT_start]
        return 'polyT', barcode
    elif reverse:
        polyT_start = reverse_seq.find(reverse[0])
        polyT_end = polyT_start + len(reverse[0])
        barcode = reverse_seq[polyT_start - barcode_len: polyT_start]
        return 'polyA', barcode
    else:
        return ()
    

def extract_barcodes(fastq: str, polyT_min_threshold: int, barcode_len: int):
    """
    Given a fastq file, extract all of the barcodes.

    Args:
        fastq (str): The path to the fastq file
        polyT_min_threshold (int): The minimum size required for a polyT(A) to be identified in the read

    Returns:
        dictionary: The read names and barcodes. Unique barcodes are the keys, the reads
                    that have that barcode are held in a list as the value
    """
    polyT_count = 0
    polyA_count = 0
    missing_count = 0

    if fastq.endswith('.gz'):
        with gzip.open(fastq, 'rt') as file:
            for record in FastqGeneralIterator(file):  # tuple of title, sequence, quality
                barcode = find_barcode(record[1], polyT_min_threshold, barcode_len)
                if barcode == ():
                    missing_count += 1
                elif barcode[0] == 'polyT':
                    polyT_count += 1
                    print(barcode[1])
                elif barcode[0] == 'polyA':
                    polyA_count += 1
                    print(barcode[1])
    else:
        with open(fastq, 'rt') as file:
            for record in FastqGeneralIterator(file):  # tuple of title, sequence, quality
                barcode = find_barcode(record[1], polyT_min_threshold, barcode_len)
                if barcode == ():
                    missing_count += 1
                elif barcode[0] == 'polyT':
                    polyT_count += 1
                    print(barcode[1])
                elif barcode[0] == 'polyA':
                    polyA_count += 1
                    print(barcode[1])
    print(polyT_count)
    print(polyA_count)
    print(missing_count)
