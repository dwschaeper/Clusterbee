from collections import defaultdict
from itertools import combinations
import matplotlib.pyplot as plt
import numpy as np
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
import os

logger = logging.getLogger(__name__)
TRANSLATE = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
_worker_processor = None
# Byte-popcount lookup for NumPy versions without np.bit_count.
# This is over uint8 values (0..255), independent of barcode alphabet size.
_BYTE_POPCOUNT_TABLE = np.array([bin(i).count('1') for i in range(256)], dtype=np.uint8)


class VectorizedBatchProcessor:
    """Worker class that processes batches of putative barcodes against encoded anchors."""
    def __init__(self, max_dist, anchor_barcodes, anchor_encoded, anchor_kmer_index, anchor_myer_states, keep_best_match_only,
                 kmer_size, min_kmer_matches, max_myers_candidates):
        self.max_dist = max_dist
        self.anchor_barcodes = anchor_barcodes
        self.anchor_encoded = anchor_encoded
        self.anchor_kmer_index = anchor_kmer_index
        self.anchor_myer_states = anchor_myer_states
        self.keep_best_match_only = keep_best_match_only
        self.kmer_size = kmer_size
        self.min_kmer_matches = min_kmer_matches
        self.max_myers_candidates = max_myers_candidates


    @staticmethod
    def bit_count_uint32(values: np.ndarray) -> np.ndarray:
        """Vectorized bit count for uint32 arrays with compatibility fallback."""
        values = values.astype(np.uint32, copy=False)
        if hasattr(np, 'bit_count'):
            return np.bit_count(values)
        byte_view = values.view(np.uint8).reshape(values.shape + (4,))
        return _BYTE_POPCOUNT_TABLE[byte_view].sum(axis=-1, dtype=np.uint16)
    

    def hamming_all(self, query: int, barcodes: np.ndarray) -> np.ndarray:
        """Compute hamming distance from one encoded barcode to all encoded anchors."""
        x = np.bitwise_xor(barcodes, np.uint32(query))
        x = np.bitwise_and(np.bitwise_or(x, np.right_shift(x, np.uint32(1))), np.uint32(0x55555555))
        return self.bit_count_uint32(x)
    

    def get_kmers(self, putative_barcode):
        kmers = []
        for start in range(0, len(putative_barcode) - self.kmer_size + 1):
            kmers.append(putative_barcode[start: start + self.kmer_size])
        return kmers
    

    def kmer_filter_myers(self, putative_barcode):
        anchor_match_counts = defaultdict(int)
        for kmer in self.get_kmers(putative_barcode):
            for anchor_idx in self.anchor_kmer_index.get(kmer, ()):  # sparse list lookup
                anchor_match_counts[anchor_idx] += 1

        filtered = [(count, anchor_idx) for anchor_idx, count in anchor_match_counts.items() if count >= self.min_kmer_matches]
        filtered = sorted(filtered, key=lambda x: (-x[0], x[1]))

        return filtered[:self.max_myers_candidates]
    
    
    def myers_levenshtein(self, text: str, max_dist: int, peq: dict[str, int], last_bit: int, bit_mask: int, pattern_len: int) -> int | None:
        """Return Levenshtein distance if it is within cutoff, else None.

        The implementation uses the Myers bit-parallel algorithm and prunes as soon as
        the best possible final score can no longer reach the current cutoff.
        """
        text_len = len(text)

        pv = bit_mask
        mv = 0
        score = pattern_len

        for idx, base in enumerate(text):
            eq = peq.get(base, 0)
            xv = eq | mv
            xh = (((eq & pv) + pv) ^ pv) | eq
            ph = mv | ~(xh | pv)
            mh = pv & xh

            if ph & last_bit:
                score += 1
            elif mh & last_bit:
                score -= 1

            remaining = text_len - idx - 1
            if score - remaining > max_dist:
                return None

            ph = ((ph << 1) | 1) & bit_mask
            mh = (mh << 1) & bit_mask
            pv = (mh | ~(xv | ph)) & bit_mask
            mv = (ph & xv) & bit_mask

        return score if score <= max_dist else None

    
    def process_myers(self, putative_barcode):
        candidates = self.kmer_filter_myers(putative_barcode)

        results = []
        for _, anchor_idx in candidates:
            anchor_barcode, peq, last_bit, bit_mask, pattern_len = self.anchor_myer_states[anchor_idx]
            dist = self.myers_levenshtein(putative_barcode, self.max_dist, peq, last_bit, bit_mask, pattern_len)
            if dist:
                results.append((anchor_barcode, dist))
        
        return results


    def process(self, records_batch):
        """Process a batch of putative barcodes using Hamming first, then Myers fallback."""
        local_out = []
        for putative_barcode, put_freq in records_batch:
            encoded = encode_barcode(putative_barcode)
            distances = self.hamming_all(encoded, self.anchor_encoded)
            valid_idx = np.flatnonzero(distances <= self.max_dist)

            if valid_idx.size > 0:
                used = 'hamming'
                matches = [(self.anchor_barcodes[idx], int(distances[idx])) for idx in valid_idx]
            else:
                myers_matches = self.process_myers(putative_barcode)
                if myers_matches:
                    used = 'myers'
                    matches = myers_matches
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


def _worker_init(max_dist, anchor_barcodes, anchor_encoded, anchor_kmer_index, anchor_myer_states, keep_best_match_only,
                 kmer_size, min_kmer_matches, max_myers_candidates):
    """
    Initialize per-process VectorizedBatchProcessor.

    Args:
        max_dist (int):                    Maximum distance allowed for matches.
        anchor_barcodes (tuple[str, ...]): Anchor barcode strings.
        anchor_encoded (np.ndarray):       Anchor barcode encodings (uint32).
        anchor_kmer_index
        anchor_myer_states
        keep_best_match_only (bool): Whether to keep only lowest-distance ties.
        kmer_size (int):              K-mer size used for Myers fallback candidate prefiltering.
        min_kmer_matches (int):       Minimum shared k-mers required to evaluate an anchor via Myers.
        max_myers_candidates (int):    Maximum number of candidate anchors to send into Myers per putative.

    Returns:
        None
    """
    global _worker_processor
    _worker_processor = VectorizedBatchProcessor(
        max_dist,
        anchor_barcodes,
        anchor_encoded,
        anchor_kmer_index,
        anchor_myer_states,
        keep_best_match_only,
        kmer_size,
        min_kmer_matches,
        max_myers_candidates
        )


def _worker_process_batch(records_batch):
    """
    Worker entry point used by ProcessPoolExecutor.

    Args:
        records_batch (list): Batch of (putative_barcode, freq) tuples.

    Returns:
        list: Results from VectorizedBatchProcessor.process() for the batch.
    """
    if _worker_processor is None:
        raise RuntimeError(
            'Worker not initialized. This function must be run inside a worker process with _worker_init called.')
    else:
        results = _worker_processor.process(records_batch)
        return results


def _safe_future_result(future):
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


def build_peq(pattern: str) -> dict[str, int]:
    """Build Myers position masks for a fixed pattern barcode."""
    peq = defaultdict(int)
    for idx, base in enumerate(pattern):
        peq[base] |= 1 << idx
    return peq


def precompute_myers_pattern(pattern: str) -> tuple[dict[str, int], int, int, int]:
    """Precompute Myers pattern state to reuse across many text queries."""
    pattern_len = len(pattern)
    peq = build_peq(pattern)
    last_bit = 1 << (pattern_len - 1) if pattern_len > 0 else 0
    bit_mask = (1 << pattern_len) - 1 if pattern_len > 0 else 0
    return peq, last_bit, bit_mask, pattern_len


def hamming(e1: int, e2: int) -> int:
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


def encode_barcode(barcode: str) -> int:
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


def calculate_anchor_neighborhoods(putative_barcodes, max_dist):
    """
    Takes the putative barcodes associated with a single anchor and calculate the
    hamming distance between them all.

    Args:
        putative_barcodes (list): List of tuples, (barcode, encoded) for each putative for an anchor
        max_dist (int):          Maximum allowed hamming distance to create an edge

    Returns:
        list: graph edges for each barcode pair
    """
    edges = []
    for (barcode1, encoded1, freq1), (barcode2, encoded2, freq2) in combinations(putative_barcodes, 2):
        dist = hamming(encoded1, encoded2)
        if dist > max_dist:
            continue
        # weight = (math.sqrt((freq1 + freq2)/2) * (2 ** -dist))
        #weight = (np.exp(-1*dist))
        #weight = (3 ** (-1.3 * dist)) * (1 + 0.05 * np.log1p(min(freq1, freq2)))
        weight = (2 ** -dist)
        edges.append((barcode1, barcode2, weight))

    return edges


def canonical_edge(barcode1: str, barcode2: str) -> tuple[str, str]:
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


def build_anchor_kmer_index(anchor_barcodes, kmer_size):
    kmer_to_anchor = defaultdict(list)
    for anchor_idx, anchor in enumerate(anchor_barcodes):
        for start in range(0, len(anchor) - kmer_size + 1):
            kmer_to_anchor[anchor[start:start + kmer_size]].append(anchor_idx)
    return kmer_to_anchor


def combine_batch_results(results_batch, combined_result, anchors: dict[str, int], keep_best_match_only=False):
    """Process a worker result batch and update mutable edge-collection state.

    Args:
        results_batch (list):        Results from VectorizedBatchProcessor.process() for the batch
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
                    barcode1, barcode2 = canonical_edge(top_anchor, putative_barcode)
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
                    barcode1, barcode2 = canonical_edge(anchor_barcode, putative_barcode)
                    combined_result['graph_edges'].append((barcode1, barcode2, weight))
                    combined_result['total_edges'] += 1
                    combined_result['anchor_neighborhood'][anchor_barcode].add((putative_barcode, encoded, put_freq))

        combined_result['processed'] += 1
        if combined_result['processed'] % 500000 == 0:
            logger.info('Queried %d putative barcodes', combined_result['processed'])

    return combined_result


def construct_graph_edges(anchors: dict[str, int], putatives: dict[str, int], max_dist: int, barcode_len: int,
                          processes: int = None, batch_size: int = 10000, keep_best_match_only: bool = False,
                          kmer_size: int = 4, min_kmer_matches: int = 3, max_myers_candidates: int = 7):
    """
    Build graph edges using vectorized Hamming matching with Myers fallback for Hamming misses.

    Args:
        anchors (dict):         Anchor barcode -> frequency mapping.
        putatives (dict):       Putative barcode -> frequency mapping.
        max_dist (int):         Largest acceptable distance.
        barcode_len (int):      Length of the barcode.
        processes (int):        Number of worker processes (defaults to cpu_count()-1).
        batch_size (int):       Number of putative barcodes to send per worker task.
        keep_best_match_only (bool): Keep only lowest-distance ties when True; otherwise keep all valid matches.
        kmer_size (int):        K-mer size for Myers fallback candidate filtering.
        min_kmer_matches (int): Minimum shared k-mers required for Myers candidate anchors.
        max_myers_candidates (int): Maximum number of candidate anchors to consider with Myers per putative.

    Returns:
        None
    """
    logger.section('EDGE CONSTRUCTION')
    logger.info('Begin construction of edges')

    if barcode_len > 16:
        raise ValueError('uint32 2-bit barcode encoding supports barcode_len <= 16')
    if kmer_size <= 0:
        raise ValueError('kmer_size must be > 0')
    if min_kmer_matches <= 0:
        raise ValueError('min_kmer_matches must be > 0')
    if max_myers_candidates <= 0:
        raise ValueError('max_myers_candidates must be > 0')
    if kmer_size > barcode_len:
        raise ValueError('kmer_size must be <= barcode_len')

    logger.info('Encoding anchor barcodes for vectorized hamming search')
    anchor_barcodes = tuple(anchors.keys())
    anchor_encoded = np.array([encode_barcode(barcode) for barcode in anchor_barcodes], dtype=np.uint32)
    # create the anchor kmer index kmer -> anchor
    anchor_kmer_index = build_anchor_kmer_index(anchor_barcodes, kmer_size)
    anchor_myer_states = [(barcode, *precompute_myers_pattern(barcode)) for barcode in anchor_barcodes]

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

    logger.info('Using %d worker processes (batch_size=%d)', max_workers, batch_size)
    logger.info('Myers fallback k-mer prefilter: k=%d, min_shared_kmers=%d, max_candidates=%d',
                kmer_size, min_kmer_matches, max_myers_candidates)
    executor_cm = ProcessPoolExecutor(max_workers=max_workers, initializer=_worker_init,
                                      initargs=(
                                          max_dist,
                                          anchor_barcodes,
                                          anchor_encoded,
                                          anchor_kmer_index,
                                          anchor_myer_states,
                                          keep_best_match_only,
                                          kmer_size,
                                          min_kmer_matches,
                                          max_myers_candidates
                                      ))

    with open(".graph_edges", "w") as graph_edges_f, executor_cm as executor:
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
                combine_batch_results(results, combined_result, anchors, keep_best_match_only)

        # submit remaining records
        if batch:
            pending.append(executor.submit(_worker_process_batch, batch))
        for future in as_completed(pending):
            results = _safe_future_result(future)
            combine_batch_results(results, combined_result, anchors, keep_best_match_only)

        for anchor, putatives_set in combined_result['anchor_neighborhood'].items():
            # Convert set to list for calculate_anchor_neighborhoods
            neighborhood_edges = calculate_anchor_neighborhoods(list(putatives_set), max_dist)
            for neighborhood_edge in neighborhood_edges:
                barcode1, barcode2 = canonical_edge(neighborhood_edge[0], neighborhood_edge[1])
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
    logger.info('Neighborhood edge count: %d', neighborhood_edge_count)
    logger.info('Distribution of putatives that match multiple anchors %s', str(combined_result['putative_barcodes_match_size']))

    anchor_sizes = [len(x) for x in combined_result['anchor_neighborhood'].values()]
    for val in anchor_sizes:
        anchor_size_freq[val] += 1
    plt.hist(anchor_sizes)
    plt.savefig('anchor_sizes.png')
    logger.info('The distribution of the number of putative barcodes linked to an anchor barcode is saved to "anchor_sizes.png"')

