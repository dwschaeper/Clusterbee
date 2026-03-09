from collections import defaultdict
from itertools import combinations
import matplotlib.pyplot as plt
import math
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

    def __init__(self, max_hamming_dist, trie, barcode_len):
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
            hits = []
            if matches:
                # take best match, or all those tied for it
                min_dist = 10
                for match in matches:
                    #if match[1] == 0:
                    #    continue  # skip perfect matches (putative barcode is an anchor matching itself)
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

    def _encode(self,barcode: str) -> int:
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


def _worker_init(max_hamming_dist, trie, barcode_len):
    """
    Initialize per-process TrieBatchProcessor.

    Args:
        max_hamming_dist (int): Maximum Hamming distance allowed for matches.
        trie (TrieNode): Pre-built anchor trie to copy into the worker.
        barcode_len (int): Expected barcode length (used by DFS termination).

    Returns:
        None
    """
    global _worker_processor
    _worker_processor = TrieBatchProcessor(max_hamming_dist, trie, barcode_len)


def _worker_process_batch(records_batch):
    """
    Worker entry point used by ProcessPoolExecutor.

    Args:
        records_batch (list): Batch of (putative_barcode, freq) tuples.

    Returns:
        list: Results from TrieBatchProcessor.process() for the batch.
    """
    if _worker_processor is None:
        raise RuntimeError('Worker not initialized. This function must be run inside a worker process with _worker_init called.')
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
    """Encode a barcode string into a compact 2-bit-per-base integer."""
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
        weight = (math.sqrt((freq1 + freq2)/2) * (2 ** -dist))
        #weight = math.sqrt((freq1 * freq2)) * (2 ** -dist)
        edges.append((barcode1, barcode2, weight))

    return edges


def build_read_linkage_weights(used_barcodes: set[str], barcode_to_anchor: dict[str, str], anchors: dict[str, int], putatives: dict[str, int]) -> dict[str, list[tuple[str, int, int]]]:
    """
    Build pairwise linkage data for barcodes that co-occur in the same read.

    Returns barcodes ready to add to anchor neighborhoods, organized by anchor.
    
    Args:
        used_barcodes (set[str]): Barcodes already represented in the graph.
        barcode_to_anchor (dict[str, str]): Mapping of barcode -> anchor.
        anchors (dict[str, int]): Anchor barcode -> frequency mapping.
        putatives (dict[str, int]): Putative barcode -> frequency mapping.

    Returns:
        dict[str, list[tuple[str, int, int]]]: anchor -> [(barcode, encoded, freq), ...] tuples ready to add to neighborhood.
    """
    read_to_barcodes = defaultdict(set)
    pair_counts = defaultdict(int)

    with open('.putative_barcodes', 'r') as infile:
        for line in infile:
            row = line.rstrip('\n').split('\t')
            read_id, barcode = row[0], row[1]
            if barcode in used_barcodes:
                read_to_barcodes[read_id].add(barcode)

    for barcodes in read_to_barcodes.values():
        for barcode1, barcode2 in combinations(sorted(barcodes), 2):
            pair_counts[(barcode1, barcode2)] += 1

    # Organize results by anchor with actual frequencies and encoded values.
    # Only cross-anchor pairs are returned because same-anchor links are already
    # represented within each anchor neighborhood.
    results_by_anchor = defaultdict(list)
    encoded_cache = {}
    for (barcode1, barcode2), _co_occurrence_count in pair_counts.items():
        anchor1 = barcode_to_anchor.get(barcode1)
        anchor2 = barcode_to_anchor.get(barcode2)
        
        if not anchor1 or not anchor2:
            continue  # Skip if barcode anchor mapping not found

        if anchor1 == anchor2:
            continue
        
        # Get actual frequencies
        freq1 = anchors.get(barcode1, putatives.get(barcode1, 1))
        freq2 = anchors.get(barcode2, putatives.get(barcode2, 1))
        
        # Encode each barcode once.
        if barcode1 not in encoded_cache:
            encoded_cache[barcode1] = encode_barcode(barcode1)
        if barcode2 not in encoded_cache:
            encoded_cache[barcode2] = encode_barcode(barcode2)
        encoded1 = encoded_cache[barcode1]
        encoded2 = encoded_cache[barcode2]
        
        # Return tuples organized by anchor where they need to be added (barcode, encoded, freq)
        results_by_anchor[anchor1].append((barcode2, encoded2, freq2))
        results_by_anchor[anchor2].append((barcode1, encoded1, freq1))

    return results_by_anchor


def canonical_edge(barcode1: str, barcode2: str) -> tuple[str, str]:
    """Return a stable ordering for undirected barcode-barcode edges."""
    if barcode1 <= barcode2:
        return barcode1, barcode2
    return barcode2, barcode1


def consume_results_batch(results_batch, anchors: dict[str, int], multi_match, non_anchor_barcodes: set[str], edge_weights,
    used_barcodes: set[str], anchor_neighborhood, barcode_to_anchor: dict[str, str], edges: int, without_matches: int, with_matches: int, processed: int,) -> tuple[int, int, int, int, dict, dict]:
    """Process a worker result batch and update mutable edge-collection state."""
    for putative_barcode, put_freq, encoded, matches in results_batch:
        if not matches:
            without_matches += 1
        else:
            with_matches += 1
            multi_match[len(matches)] += 1
            max_freq_match = (0, '', '', 0.0)
            for anchor_barcode, mismatches in matches:
                if anchors[anchor_barcode] > max_freq_match[0]:
                    weight = math.sqrt(anchors[anchor_barcode]) * (2 ** -mismatches)
                    max_freq_match = (anchors[anchor_barcode], anchor_barcode, putative_barcode, weight)

            # keep the highest-frequency anchor among tied best-distance matches
            if max_freq_match[0] != 0:
                _anchor_freq, anchor, putative, weight = max_freq_match
                non_anchor_barcodes.add(putative)
                edge = canonical_edge(anchor, putative)
                edge_weights[edge] += weight
                used_barcodes.add(anchor)
                used_barcodes.add(putative)
                barcode_to_anchor[putative] = anchor
                barcode_to_anchor[anchor] = anchor
                edges += 1
                anchor_neighborhood[anchor].add((putative, encoded, put_freq))

        processed += 1
        if processed % 500000 == 0:
            logger.info('Queried %d putative barcodes', processed)

    return edges, without_matches, with_matches, processed, anchor_neighborhood, barcode_to_anchor


def construct_graph_edges(anchors: dict[str, int], putatives: dict[str, int], max_hamming_dist: int, barcode_len: int, processes: int = None, batch_size: int = 10000):
    """
    Build graph edges using trie-based search bounded by hamming distance.

    Args:
        anchors (dict):         Anchor barcode -> frequency mapping.
        putatives (dict):       Putative barcode -> frequency mapping.
        max_hamming_dist (int): Largest acceptable hamming distance.
        barcode_len (int):      Length of the barcode.
        processes (int):        Number of worker processes (defaults to cpu_count()-1).
        batch_size (int):       Number of putative barcodes to send per worker task.

    Returns:
        None
    """
    logger.info('Begin construction of edges')

    logger.info('Creating anchor trie')
    anchor_barcodes = set(anchors.keys())
    trie = TrieNode()
    for barcode in anchor_barcodes:
        insert(trie, barcode)

    logger.info('Searching the putative barcodes in the anchor trie')
    edges = 0
    without_matches = 0
    with_matches = 0
    neighborhood_edge_count = 0
    read_linkage_edge_count = 0
    anchor_neighborhood = defaultdict(set)
    anchor_size_freq = defaultdict(int)
    non_anchor_barcodes = set()
    non_anchor_reads = set()
    anchor_reads = set()
    used_barcodes = set()
    processed = 0
    multi_match = defaultdict(int)

    max_workers = processes if processes and processes > 0 else max(1, (os.cpu_count() or 2) - 1)

    logger.info('Using %d worker processes (batch_size=%d)', max_workers, batch_size)
    executor_cm = ProcessPoolExecutor(max_workers=max_workers, initializer=_worker_init,
                                      initargs=(max_hamming_dist, trie, barcode_len))

    with open(".graph_edges", "w") as graph_edges, executor_cm as executor:
        pending = []
        batch = []
        edge_weights = defaultdict(float)
        barcode_to_anchor = {}

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
                edges, without_matches, with_matches, processed, anchor_neighborhood, barcode_to_anchor = consume_results_batch(
                    results,
                    anchors,
                    multi_match,
                    non_anchor_barcodes,
                    edge_weights,
                    used_barcodes,
                    anchor_neighborhood,
                    barcode_to_anchor,
                    edges,
                    without_matches,
                    with_matches,
                    processed,
                )

        # submit remaining records
        if batch:
            pending.append(executor.submit(_worker_process_batch, batch))
        for future in as_completed(pending):
            results = _safe_future_result(future)
            edges, without_matches, with_matches, processed, anchor_neighborhood, barcode_to_anchor = consume_results_batch(
                results,
                anchors,
                multi_match,
                non_anchor_barcodes,
                edge_weights,
                used_barcodes,
                anchor_neighborhood,
                barcode_to_anchor,
                edges,
                without_matches,
                with_matches,
                processed,
            )

        # Process read-linkage edges to strengthen cross-anchor connectivity.
        linkage_results = build_read_linkage_weights(used_barcodes, barcode_to_anchor, anchors, putatives)
        cross_anchor_links = 0
        
        for anchor, cross_anchor_barcodes in linkage_results.items():
            for barcode, encoded, freq in cross_anchor_barcodes:
                # Add this barcode to the current anchor neighborhood only once.
                before = len(anchor_neighborhood[anchor])
                anchor_neighborhood[anchor].add((barcode, encoded, freq))
                if len(anchor_neighborhood[anchor]) != before:
                    cross_anchor_links += 1

        read_linkage_edge_count = cross_anchor_links
        logger.info('Read-linkage cross-anchor additions: %d', cross_anchor_links)

        for anchor, putatives_set in anchor_neighborhood.items():
            # Convert set to list for calculate_anchor_neighborhoods
            neighborhood_edges = calculate_anchor_neighborhoods(list(putatives_set))
            for neighborhood_edge in neighborhood_edges:
                edge = canonical_edge(neighborhood_edge[0], neighborhood_edge[1])
                edge_weights[edge] += neighborhood_edge[2]
                neighborhood_edge_count += 1


        # Write each undirected edge once with its cumulative weight.
        for (barcode1, barcode2), weight in edge_weights.items():
            graph_edges.write(f'{barcode1}\t{barcode2}\t{weight}\n')

    logger.info('Number of anchor-putative edges: %d', edges)
    logger.info('Number of barcodes without associated edges: %d', without_matches)
    logger.info('Number of barcodes with associated edges: %d', with_matches)
    logger.info('Total number of unique barcodes with edges: %d', len(used_barcodes))
    logger.info('Neighborhood edge count: %d', neighborhood_edge_count)
    logger.info('Read-linkage edge count: %d', read_linkage_edge_count)
    logger.info('Distribution of putatives that match multiple anchors %s', str(multi_match))

    # find out how many reads are covered by the non-anchor nodes
    with open('.putative_barcodes', 'r') as infile:
        for line in infile:
            row = line.strip().split('\t')
            if row[1] in non_anchor_barcodes:
                non_anchor_reads.add(row[0])
    logger.info('Number of reads associated with non-anchor barcodes: %d', len(non_anchor_reads))

    # find number of reads associated with anchor barcodes
    with open('.anchors', 'r') as infile:
        for line in infile:
            row = line.strip().split('\t')
            if row[1] in anchor_barcodes:
                anchor_reads.add(row[0])
    logger.info('Number of reads associated with anchor barcodes: %d', len(anchor_reads))
                
    anchor_sizes = [len(x) for x in anchor_neighborhood.values()]
    for val in anchor_sizes:
        anchor_size_freq[val] += 1
    plt.hist(anchor_sizes)
    plt.savefig('anchor_sizes.png')
    logger.info('The distrubition of the number of putative barcodes linked to an anchor barcode is saved to "anchor_sizes.png"')


