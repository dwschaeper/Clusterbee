import argparse

from bin.extract_barcodes import extract_barcodes

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, type=str, help='The input fastq file.')
    parser.add_argument('--polyT_min_threshold', default=10, type=int, 
                        help='The minimum length required to identify a polyT(A) in the read.')
    parser.add_argument('--barcode_len', default=16, type=int,
                        help='The length of the barcode used to tag the cells.')

    args = parser.parse_args()

    return args

if __name__ == '__main__':
    arguments = parse_args()
    
    extract_barcodes(arguments.input, arguments.polyT_min_threshold, arguments.barcode_len)
