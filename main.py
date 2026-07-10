import argparse
import logging
import os

from src.logging_config import setup_logging
from src.extract_barcodes import extract_barcodes
from src.edge_construction_vectorized_levenshtein_fallback import construct_graph_edges
from src.make_graph import make_graph
from src.assignment import assignment

from config import set_global_prefix, set_global_outdir


setup_logging(level=logging.INFO, logfile="honeycomb.log")
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, type=str, 
                        help='The input fastq file.')
    parser.add_argument('--processes', default=None, type=int,
                        help='Number of worker processes to use for barcode extraction (default: cpu_count()-1).')
    parser.add_argument('--batch_size', default=5000, type=int,
                        help='Number of reads per worker batch (default: 1000).')
    parser.add_argument('--outdir', default='results', type=str,
                        help='The output directory used for all generated files.')
    parser.add_argument('--prefix', default='honeycomb', type=str,
                        help='The filename prefix used for all generated files.')

    # barcode extraction parameters
    parser.add_argument('--polyT_min_length', default=4, type=int, 
                        help='The minimum length required to identify a polyT(A) in the read.')
    parser.add_argument('--barcode_len', default=16, type=int,
                        help='The length of the barcode used to tag the cells.')
    parser.add_argument('--no_umi', action='store_true',
                        help='Use this if there are no UMIs used.')
    parser.add_argument('--umi_len', default=12, type=int,
                        help='The length of the UMI sequence used.')
    parser.add_argument('--cell_barcode_whitelist', default=None, type=str, required=True,
                        help='The whitelist of cell barcodes where each line in the whitelist is a barcode.')
    parser.add_argument('--search_range', default=5, type=int,
                        help='The range around the UMI length to look for exact barcode matches if whitelist is provided.')
    parser.add_argument('--min_qual', default=10, type=int,
                        help='The minimum average quality across a barcode for it to be considered.')
    parser.add_argument('--upstream_region_tsv', default=None, type=str,
                        help='Optional output TSV path with sequence from the first barcode-search base through the polyT match.')
    parser.add_argument('--probability_threshold', default=0.85, type=float,
                        help='The probabilty required from the NB mixture to be considered a true anchor.')
    parser.add_argument('--max_iter', default=200, type=int,
                        help='The maximum number of iterations for the EM step.')
    parser.add_argument('--adapter', default='CTACACGACGCTCTTCCGATCT', type=str,
                        help='The adapter used during sequencing that will be looked for to anchor the barcode search.')
    parser.add_argument('--adapter_tail_length', default=10, type=int,
                        help='How many bases of the adapter to look for.')
    parser.add_argument('--expected_cells', default=None, type=int,
                        help='The number cells that were expected to be sequenced.')
    
    # edge construction parameters
    parser.add_argument('--max_hamming_dist', default=2, type=int,
                        help='The maximum hamming distance between a putative barcode and an anchor for it to be considered.')
    parser.add_argument('--max_levenshtein_dist', default=3, type=int,
                        help='The max levenshtein distance between a putative barcode and an anchor for it io be considered.')
    parser.add_argument('--keep_best_match_only', action='store_true',
                        help='If set, keep only matches at the minimum hamming distance (including ties). By default, keep all matches within max_dist.')
    parser.add_argument('--kmer_size', default=4, type=int,
                        help='K-mer size used to prefilter anchors before Myers fallback (default: 4).')
    parser.add_argument('--min_kmer_matches', default=3, type=int,
                        help='Minimum shared k-mers required to compare a putative/anchor pair with Myers (default: 4).')
    parser.add_argument('--max_myers_candidates', default=10, type=int,
                        help='Maximum number of candidate anchors to send into Myers after k-mer filtering (default: 7).')
    
    # graph construction parameters
    parser.add_argument('--resolution', default=0.05, type=float,
                        help='Resolution for leiden clustering.')
    parser.add_argument('--iterations', default=-1, type=int,
                        help='Iterations in leiden clustering.')

    args = parser.parse_args()

    return args

if __name__ == '__main__':
    logger.info("Starting program")

    arguments = parse_args()

    # set prefix and create directory if needed for outputs
    set_global_outdir(arguments.outdir)
    set_global_prefix(arguments.prefix)
    os.makedirs(arguments.outdir, exist_ok=True)


    # read in the whitelist
    with open(arguments.cell_barcode_whitelist, 'r') as file:
        lines = file.readlines()
    lines = [line.strip() for line in lines]
    whitelist = set(lines)
    putative_anchors, putative_barcodes = extract_barcodes(fastq=arguments.input,
                     polyt_min_length=arguments.polyT_min_length,
                     adapter=arguments.adapter[-arguments.adapter_tail_length:],
                     barcode_len=arguments.barcode_len,
                     umi_len=arguments.umi_len,
                     search_range=arguments.search_range,
                     min_qual=arguments.min_qual,
                     expected_cells = arguments.expected_cells,
                     whitelist=whitelist,
                     processes=arguments.processes,
                     batch_size=arguments.batch_size,
                     upstream_region_tsv=arguments.upstream_region_tsv,
                     probability_threshold=arguments.probability_threshold,
                     max_iter=arguments.max_iter
                     )

    construct_graph_edges(
        putative_anchors,
        putative_barcodes,
        arguments.max_hamming_dist,
        arguments.max_levenshtein_dist,
        arguments.barcode_len,
        processes=arguments.processes,
        batch_size=arguments.batch_size,
        keep_best_match_only=arguments.keep_best_match_only,
        kmer_size=arguments.kmer_size,
        min_kmer_matches=arguments.min_kmer_matches,
        max_myers_candidates=arguments.max_myers_candidates
    )
    
    graph_data = make_graph(putative_anchors, putative_barcodes, arguments.resolution, arguments.iterations)

    assignment(graph_data, set(putative_anchors.keys()), arguments.processes, arguments.batch_size * 50)
