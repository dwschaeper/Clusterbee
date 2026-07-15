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
    parser = argparse.ArgumentParser(description='Honeycomb Barcode Deconvolution', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group('Required Parameters')
    required.add_argument('-i', '--input', required=True, type=str, 
                        help='The input fastq file.')

    required.add_argument('--cell_barcode_whitelist', default=None, type=str, required=True,
                        help='The whitelist of cell barcodes where each line in the whitelist is a barcode.')


    global_params = parser.add_argument_group('Global Parameters')

    global_params.add_argument('--processes', default=None, type=int,
                        help='Number of worker processes to use for barcode extraction. If default, uses the number of CPUs - 1.')

    global_params.add_argument('--batch_size', default=5000, type=int,
                        help='Number of reads per worker batch.')

    global_params.add_argument('--outdir', default='results', type=str,
                        help='The output directory used for all generated files.')

    global_params.add_argument('--prefix', default='honeycomb', type=str,
                        help='The filename prefix used for all generated files.')


    # barcode extraction parameters
    extraction = parser.add_argument_group('Barcode Extraction Parameters')

    extraction.add_argument('--polyT_min_length', default=4, type=int, 
                        help='The minimum length required to identify a polyT(A) in the read.')

    extraction.add_argument('--barcode_len', default=16, type=int,
                        help='The length of the barcode used to tag the cells.')

    extraction.add_argument('--umi_len', default=12, type=int,
                        help='The length of the UMI sequence used.')

    extraction.add_argument('--search_range', default=5, type=int,
                        help='The range of bases after the adapter to look for a barcode.')

    extraction.add_argument('--min_qual', default=10, type=int,
                        help='The minimum average quality across a barcode for it to be considered a high quality match.')

    extraction.add_argument('--probability_threshold', default=0.85, type=float,
                        help='The probabilty required from the NB mixture to be considered a true anchor.')

    extraction.add_argument('--max_iter', default=200, type=int,
                        help='The maximum number of iterations for the EM step.')

    extraction.add_argument('--adapter', default='CTACACGACGCTCTTCCGATCT', type=str,
                        help='The adapter used during sequencing that will be looked for to anchor the barcode search. The default is for 10X v3')

    extraction.add_argument('--adapter_tail_length', default=10, type=int,
                        help='How many bases of the adapter to look for.')

    extraction.add_argument('--expected_cells', default=None, type=int,
                        help='The number cells that were expected to be sequenced.')

    extraction.add_argument('--upstream_region_tsv', default=None, type=str,
                        help='Optional output TSV path with sequence from the first barcode-search base through the polyT match.')
    
    # edge construction parameters
    edge = parser.add_argument_group('Graph Edge Construction Parameters')

    edge.add_argument('--max_hamming_dist', default=2, type=int,
                        help='The maximum hamming distance between a putative barcode and an anchor for it to be considered.')

    edge.add_argument('--max_levenshtein_dist', default=3, type=int,
                        help='The max levenshtein distance between a putative barcode and an anchor for it io be considered.')

    edge.add_argument('--keep_best_match_only', action='store_true',
                        help='If set, keep only matches at the minimum hamming distance (including ties). By default, keep all matches within max_dist.')

    edge.add_argument('--kmer_size', default=4, type=int,
                        help='K-mer size used to prefilter anchors before fallback to calculate levenshtein distance.')

    edge.add_argument('--min_kmer_matches', default=3, type=int,
                        help='Minimum shared k-mers between a putative barcode and anchor barcode required to perform levenshtein distance fallback calculation between a putative/barcode pair.')

    edge.add_argument('--max_levenshtein_candidates', default=10, type=int,
                        help='Maximum number of candidate anchors to send to levenshtein distance fallback calculation after k-mer filtering.')
    
    # graph construction parameters
    graph = parser.add_argument_group('Graph Construction Parameters')
    graph.add_argument('--resolution', default=0.05, type=float,
                        help='Resolution for leiden clustering.')

    graph.add_argument('--iterations', default=-1, type=int,
                        help='Iterations in leiden clustering. -1 means to iterate until convergence.')

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
        max_levenshtein_candidates=arguments.max_levenshtein_candidates
    )
    
    graph_data = make_graph(putative_anchors, putative_barcodes, arguments.resolution, arguments.iterations)

    assignment(graph_data, set(putative_anchors.keys()), arguments.processes, arguments.batch_size * 50)
