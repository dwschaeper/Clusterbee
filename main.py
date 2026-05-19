import argparse
import logging

from src.logging_config import setup_logging
from src.extract_barcodes import extract_barcodes
from src.edge_construction_vectorized_levenshtein_fallback import construct_graph_edges
from src.make_graph import make_graph
from src.reads_to_cells import assignment


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
    parser.add_argument('--min_qual', default=10, type=int,
                        help='The minimum average quality across a barcode for it to be considered.')
    parser.add_argument('--min_anchor_percentile', default=1,
                        help='The percentile used to filter spurious anchors found during initial the search. By default, the anchors with exact matches in the bottom 1% are filtered as spurious.')
    parser.add_argument('--exact_anchor', default=False,
                        help='Whether or not to enforce anchors must be one per read and at the exact position as expected')
    parser.add_argument('--upstream_region_tsv', default=None, type=str,
                        help='Optional output TSV path with sequence from the first barcode-search base through the polyT match.')
    parser.add_argument('--keep_best_match_only', action='store_true',
                        help='If set, keep only matches at the minimum hamming distance (including ties). By default, keep all matches within max_dist.')
    parser.add_argument('--kmer_size', default=4, type=int,
                        help='K-mer size used to prefilter anchors before Myers fallback (default: 4).')
    parser.add_argument('--min_kmer_matches', default=3, type=int,
                        help='Minimum shared k-mers required to compare a putative/anchor pair with Myers (default: 3).')
    parser.add_argument('--max_myers_candidates', default=7, type=int,
                        help='Maximum number of candidate anchors to send into Myers after k-mer filtering (default: 7).')
    parser.add_argument('--resolution', default=0.05, type=float,
                        help='Resolution for leiden clustering.')
    parser.add_argument('--iterations', default=-1, type=int,
                        help='Iterations in leiden clustering.')

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
                     search_range=arguments.search_range,
                     min_qual=arguments.min_qual,
                     min_anchor_percentile=arguments.min_anchor_percentile,
                     exact_anchor = arguments.exact_anchor,
                     whitelist=whitelist,
                     processes=arguments.processes,
                     batch_size=arguments.batch_size,
                     upstream_region_tsv=arguments.upstream_region_tsv
                     )

    with open('.anchor_freqs', 'w') as out:
        for key, value in putative_anchors.items():
            out.write(f'{key},{value}\n')
    
    construct_graph_edges(
        putative_anchors,
        putative_barcodes,
        arguments.max_dist,
        arguments.barcode_len,
        processes=arguments.processes,
        batch_size=arguments.batch_size,
        keep_best_match_only=arguments.keep_best_match_only,
        kmer_size=5,
        min_kmer_matches=6,
        max_myers_candidates=arguments.max_myers_candidates
    )
    
    graph_data = make_graph(putative_anchors, putative_barcodes, arguments.resolution, arguments.iterations)

    assignment(graph_data, putative_barcodes)
