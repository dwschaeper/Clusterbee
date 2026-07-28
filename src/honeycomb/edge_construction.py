"""
Construct barcode graph edges with vectorized Hamming and Levenshtein distance fallback.

This module assigns putative barcodes to anchors, then emits weighted edge lists
for downstream community detection.
"""


import logging
import os
from collections import defaultdict
from concurrent.futures import Future, ProcessPoolExecutor, as_completed
from itertools import combinations
from typing import Any

import config
import edlib
import matplotlib.pyplot as plt
import numpy as np

logger = logging.getLogger(__name__)
TRANSLATE = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
_worker_processor = None
# Byte-popcount lookup for NumPy versions without np.bit_count.
# This is over uint8 values (0..255), independent of barcode alphabet size.
_BYTE_POPCOUNT_TABLE = np.array([(i).bit_count() for i in range(256)], dtype=np.uint8)


class _VectorizedBatchProcessor:
    """
    Process putative barcodes in worker processes.

    Args:
        max_hamming_dist (int):                   Maximum accepted Hamming distance.
        max_levenshtein_dist (int):               Maximum accepted Levenshtein distance.
        anchor_barcodes (tuple[str, ...]):        Anchor barcode sequences.
        anchor_encoded (np.ndarray):              2-bit encoded anchors.
        anchor_kmer_index (dict[str, list[int]]): K-mer to anchor-index lookup table.
        keep_best_match_only (bool):              Keep only minimum-distance ties.
        kmer_size (int):                          K-mer size for fallback prefilter.
        min_kmer_matches (int):                   Minimum shared k-mers to evaluate a candidate.
        max_levenshtein_candidates (int):               Maximum fallback candidates per putative.
    """

    def __init__(self, max_hamming_dist: int, max_levenshtein_dist: int, anchor_barcodes: tuple[str, ...], 
                anchor_encoded: np.ndarray, anchor_kmer_index: dict[str, list[int]], keep_best_match_only: bool,
                kmer_size: int, min_kmer_matches: int, max_levenshtein_candidates: int) -> None:
        self.max_hamming_dist = max_hamming_dist
        self.max_levenshtein_dist = max_levenshtein_dist
        self.anchor_barcodes = anchor_barcodes
        self.anchor_encoded = anchor_encoded
        self.anchor_kmer_index = anchor_kmer_index
        self.keep_best_match_only = keep_best_match_only
        self.kmer_size = kmer_size
        self.min_kmer_matches = min_kmer_matches
        self.max_levenshtein_candidates = max_levenshtein_candidates


    @staticmethod
    def _bit_count_uint32(values: np.ndarray) -> np.ndarray:
        """
        Vectorized bit count for uint32 arrays with compatibility fallback.

        Args:
            values (np.ndarray): Comparison of barcodes array

        Returns:
            np.ndarray: Dist
        """
        values = values.astype(np.uint32, copy=False)
        if hasattr(np, 'bit_count'):
            return np.bit_count(values)
        byte_view = values.view(np.uint8).reshape(values.shape + (4,))
        return _BYTE_POPCOUNT_TABLE[byte_view].sum(axis=-1, dtype=np.uint16)
    

    def _hamming_all(self, query: int, barcodes: np.ndarray) -> np.ndarray:
        """
        Compute hamming distance from one encoded barcode to all encoded anchors.

        Args:
            query (int):           An encoded barcode
            barcodes (np.ndarray): All encoded barcodes

        Returns:
            np.ndarray: Dist
        """
        x = np.bitwise_xor(barcodes, np.uint32(query))
        x = np.bitwise_and(np.bitwise_or(x, np.right_shift(x, np.uint32(1))), np.uint32(0x55555555))
        return self._bit_count_uint32(x)
    

    def _get_kmers(self, putative_barcode: str) -> list[str]:
        """
        Collect the kmers for a barcode.

        Args:
            putative_barcode (str): The barcode to get the kmers from

        Returns:
            list[str]: The kmers
        """
        kmers = []
        for start in range(len(putative_barcode) - self.kmer_size + 1):
            kmers.append(putative_barcode[start: start + self.kmer_size])
        return kmers
    

    def _kmer_filter_levenshtein(self, putative_barcode: str) -> list[tuple[int, int]]:
        """
        Filter anchor set by minimum kmer matches to the putative barcode.

        Args:
            putative_barcode (str): The barcode of interest that the anchors should match
       
        Returns:
            list[tuple[int, int]]: The filtered set of anchors to use with kmer matches and anchor idx in the tuple
        """
        anchor_match_counts = defaultdict(int)
        for kmer in self._get_kmers(putative_barcode):
            for anchor_idx in self.anchor_kmer_index.get(kmer, ()):  # sparse list lookup
                anchor_match_counts[anchor_idx] += 1

        filtered = [(count, anchor_idx) for anchor_idx, count in anchor_match_counts.items() if count >= self.min_kmer_matches]
        filtered = sorted(filtered, key=lambda x: (-x[0], x[1]))

        return filtered[:self.max_levenshtein_candidates]


    def _edlib_alignment(self, query: str, target: str) -> int:
        """
        Calculate edit distance between a barcode and anchor.

        Args:
            query (str):  The barcode of interest
            target (str): An anchor

        Returns:
            int: The distance
        """
        aln = edlib.align(target, query, mode='NW', task='distance', k=self.max_levenshtein_dist)
        edit_dist = aln['editDistance']

        # check if no alignment
        if edit_dist == -1:
            return self.max_levenshtein_dist + 1

        return edit_dist

    
    def _process_levenshtein(self, putative_barcode: str) -> list[tuple[str, int]]:
        """
        Find the anchors that match the putative with edit distance.

        Args:
            putative_barcode (str): The barcode of interest

        Returns:
            list[tuple[str, int]]: List of (anchor, edit_distance)
        """
        candidates = self._kmer_filter_levenshtein(putative_barcode)

        results = []
        for _, anchor_idx in candidates:
            anchor_barcode = self.anchor_barcodes[anchor_idx]
            dist = self._edlib_alignment(putative_barcode, anchor_barcode)
            if dist <= self.max_levenshtein_dist:
                results.append((anchor_barcode, dist))
        
        return results


    def process(self, records_batch: list[tuple[str, int]]) -> list[tuple[str, int, int, list[str], str]]:
        """
        Process a batch of putative barcodes using Hamming first, then Levenshtein fallback.
        
        Args:
            records_batch (list[tuple[str, int]]): Batch of records with putative barcode and its frequency.

        Returns:
            list[tuple[str, int, int, list[str], str]]: Comparison result with barcode, barcode frequencies, encodings, matches, and method used.
        """
        local_out = []
        for putative_barcode, put_freq in records_batch:
            encoded = _encode_barcode(putative_barcode)
            distances = self._hamming_all(encoded, self.anchor_encoded)
            valid_idx = np.flatnonzero(distances <= self.max_hamming_dist)

            if valid_idx.size > 0:
                used = 'hamming'
                matches = [(self.anchor_barcodes[idx], int(distances[idx])) for idx in valid_idx]
            else:
                levenshtein_matches = self._process_levenshtein(putative_barcode)
                if levenshtein_matches:
                    used = 'levenshtein'
                    matches = levenshtein_matches
                else:
                    used = 'unmatched'
                    matches = []

            if matches and self.keep_best_match_only:
                    min_dist = min(dist for _anchor, dist in matches)
                    matches = [match for match in matches if match[1] == min_dist]

            if matches:
                local_out.append((putative_barcode, put_freq, encoded, matches, used))
            else:
                local_out.append((putative_barcode, put_freq, None, [], used))

        return local_out


def _worker_init(max_hamming_dist: int, max_levenshtein_dist: int, anchor_barcodes: tuple[str, ...], 
                anchor_encoded: np.ndarray, anchor_kmer_index: dict[str, list[int]], keep_best_match_only: bool,
                kmer_size: int, min_kmer_matches: int, max_levenshtein_candidates: int) -> None:
    """
    Initialize per-process _VectorizedBatchProcessor.

    Args:
        max_hamming_dist (int):           Maximum hamming distance allowed for matches.
        max_levenshtein_dist (int):       Maximum edit distance allowed for matches.
        anchor_barcodes (tuple[str]):     Anchor barcode strings.
        anchor_encoded (list[np.ndarray): Anchor barcode encodings (uint32).
        anchor_kmer_index (dict):         Index related kmer to anchor index.
        keep_best_match_only (bool):      Whether to keep only lowest-distance ties.
        kmer_size (int):                  K-mer size used for Levenshtein fallback candidate prefiltering.
        min_kmer_matches (int):           Minimum shared k-mers required to evaluate an anchor via Levenshtein distance.
        max_levenshtein_candidates (int): Maximum number of candidate anchors to send into Levenshtein distance calculation per putative.

    Returns:
        None
    """
    global _worker_processor
    _worker_processor = _VectorizedBatchProcessor(
        max_hamming_dist,
        max_levenshtein_dist,
        anchor_barcodes,
        anchor_encoded,
        anchor_kmer_index,
        keep_best_match_only,
        kmer_size,
        min_kmer_matches,
        max_levenshtein_candidates
        )


def _worker_process_batch(records_batch: list[tuple[str, int]]) -> list[tuple[str, int, int | None, list[tuple[str, int]], str]]:
    """
    Worker entry point used by ProcessPoolExecutor.

    Args:
        records_batch (list[tuple[str, int]]): Batch of (putative_barcode, freq) tuples.

    Returns:
        list[tuple[str, int, int, list[str], str]]: Comparison result with barcode, barcode frequencies, encodings, matches, and method used.
    """
    if _worker_processor is None:
        raise RuntimeError(
            'Worker not initialized. This function must be run inside a worker process with _worker_init called.')
    else:
        results = _worker_processor.process(records_batch)
        return results


def _safe_future_result(future: Future[Any]) -> Any:
    """
    Resolve a Future while logging exceptions with traceback.

    Args:
        future (concurrent.futures.Future): Future to resolve.

    Returns:
        Any: The value returned by future.result(); exceptions are logged and re-raised.
    """
    try:
        return future.result()
    except Exception:
        logger.exception('Worker task raised an exception (see traceback) — re-raising')
        raise


def _hamming(e1: int, e2: int) -> int:
    """
    Calculate the hamming distance between encoded barcodes

    Args:
        e1 (int): encoding of barcode 1
        e2 (int): encoding of barcode 2

    Returns:
        int: hamming distance between barcodes
    """
    x = e1 ^ e2
    x = (x | (x >> 1)) & 0x55555555
    return x.bit_count()


def _encode_barcode(barcode: str) -> int:
    """Encode a barcode string into a compact 2-bit-per-base integer.

    Args:
        barcode (str): The barcode to be encoded

    Returns:
        int: The encoded barcode
    """
    x = 0
    for base in barcode:
        x = (x << 2) | TRANSLATE[base]
    return x


def _calculate_anchor_neighborhoods(putative_barcodes: list[tuple[str, int, int]], max_hamming_dist: int) -> list[tuple[str, str, float]]:
    """
    Takes the putative barcodes associated with a single anchor and calculate the
    hamming distance between them all.

    Args:
        putative_barcodes (list): List of tuples, (barcode, encoded) for each putative for an anchor
        max_hamming_dist (int):          Maximum allowed hamming distance to create an edge

    Returns:
        list: graph edges for each barcode pair
    """
    edges = []
    for (barcode1, encoded1, freq1), (barcode2, encoded2, freq2) in combinations(putative_barcodes, 2):
        dist = _hamming(encoded1, encoded2)
        if dist > max_hamming_dist:
            continue
        # weight = (math.sqrt((freq1 + freq2)/2) * (2 ** -dist))
        #weight = (np.exp(-1*dist))
        #weight = (3 ** (-1.3 * dist)) * (1 + 0.05 * np.log1p(min(freq1, freq2)))
        weight = (2 ** -dist)
        edges.append((barcode1, barcode2, weight))

    return edges


def _canonical_edge(barcode1: str, barcode2: str) -> tuple[str, str]:
    """Return a stable ordering for undirected barcode-barcode edges.

    Args:
        barcode1 (str): The first barcode in the edge
        barcode2 (str): The second barcode in the edge

    Returns:
        tuple[str, str]: Tuple of (barcode1, barcode2)
    """
    if barcode1 <= barcode2:
        return barcode1, barcode2
    return barcode2, barcode1


def _build_anchor_kmer_index(anchor_barcodes: tuple[str, ...], kmer_size: int) -> dict[str, list[int]]:
    """
    Build an index of anchor k-mers to anchor positions.

    Args:
        anchor_barcodes (tuple[str, ...]): Anchor barcode strings.
        kmer_size (int): K-mer size.

    Returns:
        dict[str, list[int]]: Mapping from k-mer to list of anchor indices.
    """
    kmer_to_anchor = defaultdict(list)
    for anchor_idx, anchor in enumerate(anchor_barcodes):
        for start in range(len(anchor) - kmer_size + 1):
            kmer_to_anchor[anchor[start:start + kmer_size]].append(anchor_idx)
    return kmer_to_anchor


def _combine_batch_results(results_batch: list[tuple[str, int, int | None, list[tuple[str, int]], str]],
                          combined_result: dict[str, Any], anchors: dict[str, int],
                          keep_best_match_only: bool = False) -> dict[str, Any]:
    """Process a worker result batch and update mutable edge-collection state.

    Args:
        results_batch (list):        Results from _VectorizedBatchProcessor.process() for the batch
        combined_results (dict):     Dictionary containing all of the running totals
        anchors (dict[str, int]):    Dictionary of anchor freqs
        keep_best_match_only (bool): Whether matches were pre-filtered to lowest-distance ties.

    Returns:
        dict: Running total of all results
    """
    for putative_barcode, put_freq, encoded, matches, used in results_batch:
        combined_result['assigned_method'][used] += 1

        if not matches:
            combined_result['putatives_without_matches'] += 1
        else:
            combined_result['putatives_with_matches'] += 1
            combined_result['putative_barcodes_match_size'][len(matches)] += 1

            if keep_best_match_only:
                max_freq_match = (0, '', 0)
                for anchor_barcode, mismatches in matches:
                    if mismatches == 0:
                        continue
                    if anchors[anchor_barcode] > max_freq_match[0]:
                        max_freq_match = (anchors[anchor_barcode], anchor_barcode, mismatches)

                # Final tie-break among lowest-distance matches: keep highest-frequency anchor.
                if max_freq_match[0] != 0:
                    _anchor_freq, top_anchor, mismatches = max_freq_match
                    weight = (2 ** -mismatches) + (np.sqrt(anchors[top_anchor] + put_freq))
                    barcode1, barcode2 = _canonical_edge(top_anchor, putative_barcode)
                    combined_result['graph_edges'].append((barcode1, barcode2, weight))
                    combined_result['total_edges'] += 1
                    combined_result['anchor_neighborhood'][top_anchor].add((putative_barcode, encoded, put_freq))
            else:
                for anchor_barcode, mismatches in matches:
                    if mismatches == 0:
                        continue
                    #weight = math.sqrt((anchors[anchor_barcode] * put_freq)) * (2 ** -mismatches)
                    #weight = (np.exp(-1*mismatches)) + (0.01 * np.log(anchors[anchor_barcode]))
                    weight = (2 ** -mismatches) + (np.sqrt(anchors[anchor_barcode] + put_freq))
                    barcode1, barcode2 = _canonical_edge(anchor_barcode, putative_barcode)
                    combined_result['graph_edges'].append((barcode1, barcode2, weight))
                    combined_result['total_edges'] += 1
                    combined_result['anchor_neighborhood'][anchor_barcode].add((putative_barcode, encoded, put_freq))

        combined_result['processed'] += 1
        if combined_result['processed'] % 500000 == 0:
            logger.info('Queried %d putative barcodes', combined_result['processed'])

    return combined_result


def construct_graph_edges(anchors: dict[str, int], putatives: dict[str, int], max_hamming_dist: int, max_levenshtein_dist: int, barcode_len: int,
                          processes: int | None = None, batch_size: int = 10000, keep_best_match_only: bool = False,
                          kmer_size: int = 4, min_kmer_matches: int = 3, max_levenshtein_candidates: int = 7) -> None:
    """
    Build graph edges using vectorized Hamming matching with Levenshtein distance fallback for Hamming misses.

    Args:
        anchors (dict):                   Anchor barcode -> frequency mapping.
        putatives (dict):                 Putative barcode -> frequency mapping.
        max_hamming_dist (int):           Largest acceptable distance for Hamming distance calculation.
        max_levenshtein_dist (int):       Largest acceptable distance for Levenshtein distance calculation.
        barcode_len (int):                Length of the barcode.
        processes (int):                  Number of worker processes (defaults to cpu_count()-1).
        batch_size (int):                 Number of putative barcodes to send per worker task.
        keep_best_match_only (bool):      Keep only lowest-distance ties when True; otherwise keep all valid matches.
        kmer_size (int):                  K-mer size for Levenshtein distance fallback candidate filtering.
        min_kmer_matches (int):           Minimum shared k-mers required for Levenshtein distance candidate anchors.
        max_levenshtein_candidates (int): Maximum number of candidate anchors to consider with Levenshtein distance per putative.

    Returns:
        None
    """
    logger.section('EDGE CONSTRUCTION')
    logger.info('Begin construction of edges')
    
    # initialize file handles
    graph_edge_file = os.path.join(config.GLOBAL_OUTDIR, '.graph_edges')
    anchor_size_fig = os.path.join(config.GLOBAL_OUTDIR, f'{config.GLOBAL_PREFIX}_anchor_sizes.png')

    logger.info('Encoding anchor barcodes for vectorized hamming search')
    anchor_barcodes = tuple(anchors.keys())
    anchor_encoded = np.array([_encode_barcode(barcode) for barcode in anchor_barcodes], dtype=np.uint32)
    # create the anchor kmer index kmer -> anchor
    anchor_kmer_index = _build_anchor_kmer_index(anchor_barcodes, kmer_size)

    logger.info('Searching putative barcodes against encoded anchors')
    combined_result = {'total_edges': 0,
              'graph_edges': [],
              'putatives_with_matches': 0,
              'putatives_without_matches': 0,
              'putative_barcodes_match_size': defaultdict(int),
              'processed': 0,
              'anchor_neighborhood': defaultdict(set),
              'assigned_method': defaultdict(int)}

    anchor_size_freq = defaultdict(int)
    neighborhood_edge_count = 0

    max_workers = processes if processes and processes > 0 else max(1, (os.cpu_count() or 2) - 1)

    executor_cm = ProcessPoolExecutor(max_workers=max_workers, initializer=_worker_init,
                                      initargs=(
                                          max_hamming_dist,
                                          max_levenshtein_dist,
                                          anchor_barcodes,
                                          anchor_encoded,
                                          anchor_kmer_index,
                                          keep_best_match_only,
                                          kmer_size,
                                          min_kmer_matches,
                                          max_levenshtein_candidates
                                      ))

    with open(graph_edge_file, "w") as graph_edges_f, executor_cm as executor:
        pending = []
        batch = []

        for putative_barcode, put_freq in putatives.items():
            batch.append((putative_barcode, put_freq))
            if len(batch) >= batch_size:
                pending.append(executor.submit(_worker_process_batch, batch))
                batch = []

            # keep pending futures bounded and flush completed ones
            if len(pending) >= max_workers * 2:
                done_future = next(as_completed(pending))
                pending.remove(done_future)
                results = _safe_future_result(done_future)
                _combine_batch_results(results, combined_result, anchors, keep_best_match_only)

        # submit remaining records
        if batch:
            pending.append(executor.submit(_worker_process_batch, batch))
        for future in as_completed(pending):
            results = _safe_future_result(future)
            _combine_batch_results(results, combined_result, anchors, keep_best_match_only)

        for putatives_set in combined_result['anchor_neighborhood'].values():
            # Convert set to list for calculate_anchor_neighborhoods
            neighborhood_edges = _calculate_anchor_neighborhoods(list(putatives_set), max_hamming_dist)
            for neighborhood_edge in neighborhood_edges:
                barcode1, barcode2 = _canonical_edge(neighborhood_edge[0], neighborhood_edge[1])
                combined_result['graph_edges'].append((barcode1, barcode2, neighborhood_edge[2]))
                neighborhood_edge_count += 1

        # sort graph edges for reproducibility
        combined_result['graph_edges'] = sorted(combined_result['graph_edges'])

        # Write sorted edges
        for barcode1, barcode2, weight in combined_result['graph_edges']:
            graph_edges_f.write(f'{barcode1}\t{barcode2}\t{weight}\n')

    logger.info('Number of anchor-putative edges: %d', combined_result['total_edges'])
    logger.info('Number of putative barcodes with associated edges: %d', combined_result['putatives_with_matches'])
    logger.info('Number of putative barcodes without associated edges: %d', combined_result['putatives_without_matches'])
    logger.info('Assignment distribution: %s', str(combined_result['assigned_method']))
    logger.info('Distribution of putatives that match multiple anchors %s', str(combined_result['putative_barcodes_match_size']))

    anchor_sizes = [len(x) for x in combined_result['anchor_neighborhood'].values()]
    for val in anchor_sizes:
        anchor_size_freq[val] += 1
    plt.hist(anchor_sizes)
    plt.title('Putative Barcodes per Anchor Barcode')
    plt.xlabel('Number of Putative Barcodes')
    plt.ylabel('Frequency')
    plt.savefig(anchor_size_fig)

