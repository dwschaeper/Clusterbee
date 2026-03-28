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


class TrieNode:
    """
    Trie node object
    """
    __slots__ = ("children", "anchors")

    def __init__(self):
        self.children = {}
        self.anchors = []


def insert(root: TrieNode, barcode: str):
    """
    Insert an anchor barcode into the Trie

    Args:
       root (TrieNode): TrieNode to insert into.
       barcode (str):   The barcode to insert.
    """
    node = root
    for base in barcode:
        if base not in node.children:
            node.children[base] = TrieNode()
        node = node.children[base]
    node.anchors.append(barcode)


class TrieBatchProcessor:
    """Worker class that processes batches of putative barcodes against an anchor trie."""

    def __init__(self, max_hamming_dist, trie, barcode_len, keep_all):
        """
        Initialize the TrieBatchProcessor with the anchor trie and search parameters.

        Args:
            max_hamming_dist (int): Maximum Hamming distance for matches.
            trie (TrieNode):        Pre-built anchor trie to search against.
            barcode_len (int):      Barcode length.
        """
        self.max_hamming_dist = max_hamming_dist
        self.trie = trie
        self.barcode_len = barcode_len
        self.keep_all = keep_all

    def process(self, records_batch):
        """
        Process a batch of putative barcodes against the local anchor trie.

        Args:
            records_batch (list): List of (putative_barcode, freq) tuples to query.

        Returns:
            list: Each element is (putative_barcode, freq, encoded_barcode, matches) where matches is
                  a list of (anchor_barcode, mismatches).
        """
        local_out = []
        for putative_barcode, put_freq in records_batch:
            matches = []
            # use the instance method (DFS owns max_hamming_dist)
            self._depth_first_search(self.trie, putative_barcode, 0, 0, matches)

            if not self.keep_all:
                if matches:
                    encoded = self._encode(putative_barcode)
                else:
                    encoded = None
                local_out.append((putative_barcode, put_freq, encoded, matches))

            else:
                hits = []
                if matches:
                    # take best match, or all those tied for it
                    min_dist = 10
                    for match in matches:
                        if match[1] < min_dist:
                            hits = [match]
                            min_dist = match[1]
                        elif match[1] == min_dist:
                            hits.append(match)
                    encoded = self._encode(putative_barcode)
                else:
                    encoded = None
                local_out.append((putative_barcode, put_freq, encoded, hits))

        return local_out

    def _depth_first_search(self, node: TrieNode, putative_barcode: str, pos: int, mismatches: int, results: list):
        """
        Depth-first search of the trie for a single putative barcode.

        Args:
            node (TrieNode): Current trie node being visited.
            putative_barcode (str): Barcode being queried.
            pos (int): Current index within the putative barcode.
            mismatches (int): Number of mismatches accumulated so far.
            results (list): Mutable list where (anchor, mismatches) matches are appended.

        Returns:
            None
        """
        if mismatches > self.max_hamming_dist:
            return

        if pos == self.barcode_len:
            for anchor in node.anchors:
                results.append((anchor, mismatches))
            return

        barcode_base = putative_barcode[pos]
        for trie_base, child in node.children.items():
            self._depth_first_search(child, putative_barcode, pos + 1,
                                     mismatches + (trie_base != barcode_base), results)

    def _encode(self, barcode: str) -> int:
        """
        Encode the barcode into a 32-bit integer

        Args:
            barcode (str): The barcode to be encoded

        Returns:
            int: The integer encoding
        """
        x = 0
        for base in barcode:
            x = (x << 2) | TRANSLATE[base]
        return x


def _worker_init(max_hamming_dist, trie, barcode_len, keep_all):
    """
    Initialize per-process TrieBatchProcessor.

    Args:
        max_hamming_dist (int): Maximum Hamming distance allowed for matches.
        trie (TrieNode):        Pre-built anchor trie to copy into the worker.
        barcode_len (int):      Expected barcode length (used by DFS termination).
        keep_all (bool):        Whether to keep all matches, or just the top match

    Returns:
        None
    """
    global _worker_processor
    _worker_processor = TrieBatchProcessor(max_hamming_dist, trie, barcode_len, keep_all)


def _worker_process_batch(records_batch):
    """
    Worker entry point used by ProcessPoolExecutor.

    Args:
        records_batch (list): Batch of (putative_barcode, freq) tuples.

    Returns:
        list: Results from TrieBatchProcessor.process() for the batch.
    """
    if _worker_processor is None:
        raise RuntimeError(
            'Worker not initialized. This function must be run inside a worker process with _worker_init called.')
    else:
        return _worker_processor.process(records_batch)


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


def calculate_anchor_neighborhoods(putative_barcodes):
    """
    Takes the putative barcodes associated with a single anchor and calculate the
    hamming distance between them all.

    Args:
        putative_barcodes (list): List of tuples, (barcode, encoded) for each putative for an anchor

    Returns:
        list: graph edges for each barcode pair
    """
    edges = []
    for (barcode1, encoded1, freq1), (barcode2, encoded2, freq2) in combinations(putative_barcodes, 2):
        dist = hamming(encoded1, encoded2)
        # weight = (math.sqrt((freq1 + freq2)/2) * (2 ** -dist))
        #weight = (np.exp(-1*dist))
        weight = (3 ** -dist)
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


def combine_batch_results(results_batch, anchors: dict[str, int], total_edges, graph_edges, putatives_with_matches,
                          putatives_without_matches, putative_barcodes_match_size, processed, anchor_neighborhood, keep_all=True):
    """Process a worker result batch and update mutable edge-collection state.

    Args:
        results_batch (list):                       Results from TrieBatchProcessor.process() for the batch
        anchors (dict[str, int]):                   Dictionary of anchor freqs
        total_edges (int):                          Running total number of edges
        graph_edges (list):                         List of tuples, (barcode1, barcode2, weight) for each edge in the graph
        putatives_with_matches (int):               Number of putative barcodes with matches
        putatives_without_matches (int):            Number of putative barcodes without matches
        putative_barcodes_match_size (defaultdict): Number of anchors a putative barcode matches
        anchor_neighborhood (defaultdict):          Dictionary of anchor neighborhoods
        processed (int):                            Running total of number of reads processed
        keep_all (bool):                            Whether to keep all matches, or just the top match

    Returns:
        int:         Running total number of edges
        list:        List of tuples, (barcode1, barcode2, weight) for each edge in the graph
        int:         Number of putative barcodes with matches
        int:         Number of putative barcodes without matches
        defaultdict: Number of anchors a putative barcode matches
        defaultdict: Dictionary of anchor neighborhoods
        int:         Running total of number of reads processed
    """
    if not keep_all:
        for putative_barcode, put_freq, encoded, matches in results_batch:
            if not matches:
                putatives_without_matches += 1
            else:
                putatives_with_matches += 1
                putative_barcodes_match_size[len(matches)] += 1

                for anchor_barcode, mismatches in matches:
                    if mismatches == 0:
                        continue
                    #weight = math.sqrt((anchors[anchor_barcode] * put_freq)) * (2 ** -mismatches)
                    #weight = (np.exp(-1*mismatches)) + (0.01 * np.log(anchors[anchor_barcode]))
                    weight = (3 ** -mismatches) + (np.log(anchors[anchor_barcode]))
                    barcode1, barcode2 = canonical_edge(anchor_barcode, putative_barcode)
                    graph_edges.append((barcode1, barcode2, weight))
                    total_edges += 1
                    anchor_neighborhood[anchor_barcode].add((putative_barcode, encoded, put_freq))

            processed += 1
            if processed % 500000 == 0:
                logger.info('Queried %d putative barcodes', processed)

        return total_edges, graph_edges, putatives_with_matches, putatives_without_matches, putative_barcodes_match_size, processed, anchor_neighborhood

    else:
        for putative_barcode, put_freq, encoded, matches in results_batch:
            if not matches:
                putatives_without_matches += 1
            else:
                putatives_with_matches += 1
                putative_barcodes_match_size[len(matches)] += 1
                max_freq_match = (0, '', 0.0)
                for anchor_barcode, mismatches in matches:
                    if mismatches == 0:
                        continue

                    if anchors[anchor_barcode] > max_freq_match[0]:
                        #weight = math.sqrt((anchors[anchor_barcode] * put_freq)) * (2 ** -mismatches)
                        #weight = (np.exp(-1*mismatches)) + (0.01 * np.log(anchors[anchor_barcode]))
                        weight = (3 ** -mismatches) + (np.log(anchors[anchor_barcode]))
                        max_freq_match = (anchors[anchor_barcode], anchor_barcode, weight)

                # keep the highest-frequency anchor among tied best-distance matches
                if max_freq_match[0] != 0:
                    _anchor_freq, top_anchor, weight = max_freq_match
                    barcode1, barcode2 = canonical_edge(top_anchor, putative_barcode)
                    graph_edges.append((barcode1, barcode2, weight))
                    total_edges += 1
                    anchor_neighborhood[top_anchor].add((putative_barcode, encoded, put_freq))

            processed += 1
            if processed % 500000 == 0:
                logger.info('Queried %d putative barcodes', processed)

        return total_edges, graph_edges, putatives_with_matches, putatives_without_matches, putative_barcodes_match_size, processed, anchor_neighborhood


def construct_graph_edges(anchors: dict[str, int], putatives: dict[str, int], max_hamming_dist: int, barcode_len: int,
                          processes: int = None, batch_size: int = 10000, keep_all: bool = True):
    """
    Build graph edges using trie-based search bounded by hamming distance.

    Args:
        anchors (dict):         Anchor barcode -> frequency mapping.
        putatives (dict):       Putative barcode -> frequency mapping.
        max_hamming_dist (int): Largest acceptable hamming distance.
        barcode_len (int):      Length of the barcode.
        processes (int):        Number of worker processes (defaults to cpu_count()-1).
        batch_size (int):       Number of putative barcodes to send per worker task.
        keep_all (bool):        Whether to keep all matches, or just the top match

    Returns:
        None
    """
    logger.section('EDGE CONSTRUCTION')
    logger.info('Begin construction of edges')

    logger.info('Creating anchor trie')
    anchor_barcodes = set(anchors.keys())
    trie = TrieNode()
    for barcode in anchor_barcodes:
        insert(trie, barcode)

    logger.info('Searching the putative barcodes in the anchor trie')
    total_edges = 0
    graph_edges = []
    putatives_without_matches = 0
    putatives_with_matches = 0
    putative_barcodes_match_size = defaultdict(int)
    processed = 0
    anchor_neighborhood = defaultdict(set)

    anchor_size_freq = defaultdict(int)
    neighborhood_edge_count = 0

    max_workers = processes if processes and processes > 0 else max(1, (os.cpu_count() or 2) - 1)

    logger.info('Using %d worker processes (batch_size=%d)', max_workers, batch_size)
    executor_cm = ProcessPoolExecutor(max_workers=max_workers, initializer=_worker_init,
                                      initargs=(max_hamming_dist, trie, barcode_len, keep_all))

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
                results = _safe_future_result(done_future)  # (putative_barcode, put_freq, encoded, matches)
                total_edges, graph_edges, putatives_with_matches, putatives_without_matches, putative_barcodes_match_size, processed, anchor_neighborhood = combine_batch_results(
                    results, anchors, total_edges, graph_edges, putatives_with_matches, putatives_without_matches,
                    putative_barcodes_match_size, processed, anchor_neighborhood, keep_all)

        # submit remaining records
        if batch:
            pending.append(executor.submit(_worker_process_batch, batch))
        for future in as_completed(pending):
            results = _safe_future_result(future)
            total_edges, graph_edges, putatives_with_matches, putatives_without_matches, putative_barcodes_match_size, processed, anchor_neighborhood = combine_batch_results(
                results, anchors, total_edges, graph_edges, putatives_with_matches, putatives_without_matches,
                putative_barcodes_match_size, processed, anchor_neighborhood, keep_all)

        for anchor, putatives_set in anchor_neighborhood.items():
            # Convert set to list for calculate_anchor_neighborhoods
            neighborhood_edges = calculate_anchor_neighborhoods(list(putatives_set))
            for neighborhood_edge in neighborhood_edges:
                barcode1, barcode2 = canonical_edge(neighborhood_edge[0], neighborhood_edge[1])
                graph_edges.append((barcode1, barcode2, neighborhood_edge[2]))
                neighborhood_edge_count += 1

        # sort graph edges for reproducibility
        graph_edges = sorted(graph_edges)

        # Write sorted edges
        for barcode1, barcode2, weight in graph_edges:
            graph_edges_f.write(f'{barcode1}\t{barcode2}\t{weight}\n')

    logger.info('Number of anchor-putative edges: %d', total_edges)
    logger.info('Number of putative barcodes with associated edges: %d', putatives_with_matches)
    logger.info('Number of putative barcodes without associated edges: %d', putatives_without_matches)
    logger.info('Neighborhood edge count: %d', neighborhood_edge_count)
    logger.info('Distribution of putatives that match multiple anchors %s', str(putative_barcodes_match_size))

    anchor_sizes = [len(x) for x in anchor_neighborhood.values()]
    for val in anchor_sizes:
        anchor_size_freq[val] += 1
    plt.hist(anchor_sizes)
    plt.savefig('anchor_sizes.png')
    logger.info(
        'The distrubition of the number of putative barcodes linked to an anchor barcode is saved to "anchor_sizes.png"')

