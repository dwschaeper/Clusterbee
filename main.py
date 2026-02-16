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
    parser.add_argument('--polyT_min_length', default=10, type=int, 
                        help='The minimum length required to identify a polyT(A) in the read.')
    parser.add_argument('--barcode_len', default=16, type=int,
                        help='The length of the barcode used to tag the cells.')
    parser.add_argument('--no_umi', action='store_true',
                        help='Use this if there are no UMIs used.')
    parser.add_argument('--umi_len', default=12, type=int,
                        help='The length of the UMI sequencie used.')
    parser.add_argument('--polyt_buffer', default=2, type=int,
                        help='How many bases within the found polyT to include when searching for barcodes and UMIs to account for the cases where these sequences end in T and merge with the polyT.')
    parser.add_argument('--cell_barcode_whitelist', default=None, type=str, required=True,
                        help='The whitelist of cell barcodes where each line in the whitelist is a barcode.')
    parser.add_argument('--search_range', default=2, type=int,
                        help='The range around the UMI length to look for exact barcode matches if whitelist is provided.')
    parser.add_argument('--processes', default=None, type=int,
                        help='Number of worker processes to use for barcode extraction (default: cpu_count()-1).')
    parser.add_argument('--batch_size', default=10000, type=int,
                        help='Number of reads per worker batch (default: 1000).')
    parser.add_argument('--max_dist', default=2, type=int,
                        help='The maximum hamming distance between a putative barcode and an anchor for it to be considered.')

    # Controls how strict the knee cutoff is. Values < 1.0 are less restrictive (keep smaller communities).
    parser.add_argument('--knee_scale', default=1.0, type=float,
                        help='Multiplier applied to the knee cutoff (use <1.0 to be less restrictive).')
    parser.add_argument('--diagnostics', action='store_true',
                        help='Write community diagnostics (community_summary.tsv and extra plots).')

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
                     polyt_min_length=arguments.polyT_min_length,
                     barcode_len=arguments.barcode_len,
                     umi_len=arguments.umi_len,
                     polyt_buffer=arguments.polyt_buffer,
                     search_range=arguments.search_range,
                     whitelist=whitelist,
                     processes=arguments.processes,
                     batch_size=arguments.batch_size
                     )

    construct_graph_edges(putative_anchors, putative_barcodes, arguments.max_dist, arguments.barcode_len, processes=arguments.processes, batch_size=arguments.batch_size)

    make_graph(putative_anchors, putative_barcodes, knee_scale=arguments.knee_scale, diagnostics=arguments.diagnostics)
