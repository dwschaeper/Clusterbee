from collections import defaultdict, Counter

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
    with open(".graph_edges", "w") as graph_edges:
        for putative_barcode, put_freq in putatives.items():
            matches = []
            depth_first_search(trie, putative_barcode, 0, 0, max_hamming_dist, matches)

            for anchor_barcode, mismatches in matches:
                edges += 1
                graph_edges.write(f'{anchor_barcode}\t{putative_barcode}\t{1 / mismatches}\n')

    print(f'Number of edges: {edges}')
    return edges

