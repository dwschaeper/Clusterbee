import argparse

from bin.extract_barcodes import extract_barcodes

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='The input fastq file.')

    args = parser.parse_args()

    return args

if __name__ == '__main__':
    arguments = parse_args()
    
    extract_barcodes(arguments.input)
