import argparse
import logging

from bin.logging_config import setup_logging
from bin.extract_barcodes import extract_barcodes
from bin.edge_construction import construct_graph_edges
from bin.make_graph import make_graph


setup_logging(level=logging.INFO, logfile="GBBD.log")
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, type=str, help='The input fastq file.')
    parser.add_argument('--polyT_min_threshold', default=10, type=int, 
                        help='The minimum length required to identify a polyT(A) in the read.')
    parser.add_argument('--barcode_len', default=16, type=int,
                        help='The length of the barcode used to tag the cells.')
    parser.add_argument('--no_umi', action='store_true',
                        help='Use this if there are no UMIs used.')
    parser.add_argument('--umi_len', default=12, type=int,
                        help='The length of the UMI sequencie used.')
    parser.add_argument('--polyt_include', default=2, type=int,
                        help='How many bases within the found polyT to include when searching for barcodes and UMIs to account for the cases where these sequences end in T and merge with the polyT.')
    parser.add_argument('--cell_barcode_whitelist', default=None, type=str, required=True,
                        help='The whitelist of cell barcodes where each line in the whitelist is a barcode.')
    parser.add_argument('--wiggle', default=2, type=int,
                        help='The range around the UMI length to look for exact barcode matches if whitelist is provided.')
    parser.add_argument('--max_dist', default=2, type=int,
                        help='The maximum hamming distance between a putative barcode and an anchor for it to be considered.')

    args = parser.parse_args()

    return args

if __name__ == '__main__':
    logger.info("Starting program")

    arguments = parse_args()
    
    # read in the whitelist
    with open(arguments.cell_barcode_whitelist, 'r') as file:
        lines = file.readlines()
    lines = [line.strip() for line in lines]
    whitelist = set(lines)
    putative_anchors, putative_barcodes = extract_barcodes(fastq=arguments.input,
                     polyt_min_threshold=arguments.polyT_min_threshold,
                     barcode_len=arguments.barcode_len,
                     umi_len=arguments.umi_len,
                     polyt_include=arguments.polyt_include,
                     wiggle=arguments.wiggle,
                     whitelist=whitelist
                     )

    construct_graph_edges(putative_anchors, putative_barcodes, arguments.max_dist)

    make_graph(putative_anchors, putative_barcodes)

