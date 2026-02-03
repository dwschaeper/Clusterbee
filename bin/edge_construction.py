from collections import defaultdict
from itertools import combinations
import matplotlib.pyplot as plt
import math

TRANSLATE = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

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


def depth_first_search(node: TrieNode, putative_barcode: str, pos: int, mismatches: int, max_hamming_dist: int, results: list):
    """
    Search the tree of anchor barcodes with a putative barcode ensuring that there is a hamming distance less than max_mismatches.

    Args:
        node (TrieNode):        Tree (subtree on recursive calls) to search.
        putative_barcode (str): Barcode to find an anchor for in the tree.
        pos (int):              The current position in the barcode to match.
        mismatches (int):       The current number of mismatches.
        max_mismatches (int):   The largest acceptable hamming distance.
        results (list):         The list of acceptable anchor matches.
    """
    if mismatches > max_hamming_dist:
        return

    if pos == len(putative_barcode):
        for anchor in node.anchors:
            results.append((anchor, mismatches))
        return

    barcode_base = putative_barcode[pos]

    for trie_base, child in node.children.items():
        depth_first_search(child, putative_barcode, pos + 1, mismatches + (trie_base != barcode_base), max_hamming_dist, results)


def encode(barcode: str) -> int:
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
        dist = 1 / hamming(encoded1, encoded2)
        weight = 0.3 * dist * math.sqrt(math.log(freq1) * math.log(freq2))
        edges.append((barcode1, barcode2, weight))

    return edges


def construct_graph_edges(anchors: dict[str, int], putatives: dict[str, int], max_hamming_dist: int):
    """
    Build graph edges using trie-based search bounded by hamming distance.

    Args:
        anchors (dict):         Anchor barcode: frequency
        putatives (dict):       Putative barcode: frequency
        max_hamming_dist (int): Largest acceptable hamming distance.
    """
    trie = TrieNode()
    for barcode in anchors.keys():
        insert(trie, barcode)

    edges = 0
    without_matches = 0
    with_matches = 0
    neighborhood_edge_count = 0
    anchor_candidates = defaultdict(list)
    anchor_size_freq = defaultdict(int)
    non_anchor_barcodes = set()
    non_anchor_reads = set()
    
    with open(".graph_edges", "w") as graph_edges:
        for putative_barcode, put_freq in putatives.items():
            matches = []
            depth_first_search(trie, putative_barcode, 0, 0, max_hamming_dist, matches)
            if not matches:
                without_matches += 1
            else:
                with_matches += 1
            for anchor_barcode, mismatches in matches:
                edges += 1
                anchor_candidates[anchor_barcode].append((putative_barcode, encode(putative_barcode), put_freq))
                try:
                    weight = (1 / mismatches) * math.sqrt(math.log(anchors[anchor_barcode]) * math.log(put_freq))
                except ZeroDivisionError:  # putative anchors are also stored as putative barcodes
                    weight = 2
                non_anchor_barcodes.add(putative_barcode)
                graph_edges.write(f'{anchor_barcode}\t{putative_barcode}\t{weight}\n')

        for anchor, putatives in anchor_candidates.items():
            neighborhood_edges = calculate_anchor_neighborhoods(putatives)
            for neighborhood_edge in neighborhood_edges:
                graph_edges.write(f'{neighborhood_edge[0]}\t{neighborhood_edge[1]}\t{neighborhood_edge[2]}\n')
                neighborhood_edge_count += 1

    # development metrics
    print(f'Number of anchor edges: {edges}')
    print(f'Number of barcodes without associated edges: {without_matches}')
    print(f'Number of barcodes with associated edges: {with_matches}')
    print(f'Neighborhood edge count: {neighborhood_edge_count}')

    # find out how many reads are covered by the non-anchor nodes
    with open('.putative_barcodes', 'r') as infile:
        for line in infile:
            row = line.strip().split('\t')
            if row[1] in non_anchor_barcodes:
                non_anchor_reads.add(row[0])
    print(f'Number of reads associated with non-anchor barcodes: {len(non_anchor_reads)}')
                
    anchor_sizes = [len(x) for x in anchor_candidates.values()]
    for val in anchor_sizes:
        anchor_size_freq[val] += 1
    print(anchor_size_freq)
    plt.hist(anchor_sizes)
    plt.savefig('anchor_sizes.png')


