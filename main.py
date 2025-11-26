import argparse

from bin.extract_barcodes import extract_barcodes

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, type=str, help='The input fastq file.')
    parser.add_argument('--polyT_min_threshold', default=10, type=int, 
                        help='The minimum length required to identify a polyT(A) in the read.')
    parser.add_argument('--barcode_len', default=16, type=int,
                        help='The length of the barcode used to tag the cells.')
    parser.add_argument('--umi_len', default=12, type=int,
                        help='The length of the UMI sequencie used.')
    parser.add_argument('--cell_barcode_whitelist', default=None, type=str,
                        help='The whitelist of cell barcodes where each line in the whitelist is a barcode.')
    parser.add_argument('--wiggle', default=2, type=int,
                        help='The range around the UMI length to look for exact barcode matches if whitelist is provided.')

    args = parser.parse_args()

    return args

if __name__ == '__main__':
    arguments = parse_args()

    # read in the whitelist if present
    if arguments.cell_barcode_whitelist:
        with open(arguments.cell_barcode_whitelist, 'r') as file:
            lines = file.readlines()
        lines = [line.strip() for line in lines]
        whitelist = set(lines)
        extract_barcodes(fastq=arguments.input,
                         polyt_min_threshold=arguments.polyT_min_threshold,
                         barcode_len=arguments.barcode_len,
                         umi_len=arguments.umi_len,
                         wiggle=arguments.wiggle,
                         whitelist=whitelist
                         )
    
    else:
        extract_barcodes(fastq=arguments.input, 
                         polyt_min_threshold=arguments.polyT_min_threshold, 
                         barcode_len=arguments.barcode_len, 
                         umi_len=arguments.umi_len
                         )
