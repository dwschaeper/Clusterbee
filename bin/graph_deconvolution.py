from collections import defaultdict
from itertools import combinations

TRANSLATE = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


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


def build_anchor_index(anchors: set, barcode_len: int, k: int):
    """
    Create kmer index over the anchor barcodes

    Args:
        anchors (set):     The set of anchor barcodes
        barcode_len (int): Length of the barcodes
        k (int):           Size of the kmer

    Returns:
        defaultdict(list): Index of anchors kmer: (barcode, encoded barcode)
    """
    index = defaultdict(list)

    for barcode in anchors:
        encoded = encode(barcode)
        for i in range(barcode_len - k + 1):
            kmer = barcode[i:i + k]
            index[kmer].append((barcode,encoded))
    
    return index


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


def find_edges(barcode: str, barcode_encoding: int, index, barcode_len: int, k: int):
    """
    Identify close anchors as candidates, then find distance between them.

    Args:
        barcode: ambiguous barcode
        barcoding_encoding: integer encoding of the barcode
        index: index of the anchor barcodes
        barcode_len (int): Length of the barcodes
        k (int):           Size of the kmer

    Returns:
        list: edges for the graph
    """
    candidates = set()
    edges = []

    for i in range(barcode_len - k + 1):
        kmer = barcode[i:i + k]
        candidates.update(index.get(kmer, ()))

    for anchor in candidates:
        inverse_hamming = 1 / hamming(anchor[1], barcode_encoding)
        edges.append((anchor[0], barcode, inverse_hamming))
        
    return edges


def graph_deconvolute(extractors: list, barcode_len: int, k: int):
    """
    Takes a list of barcode extractor objects and builds the graph used to assign cell barcodes to reads.

    Args:
        extractors (list): List of extractor objects from each read that contains a polyT(A)
        barcode_len (int): Length of the barcodes
        k (int):           Size of the kmer
    """
    anchor_set = set()
    putative_set = set()
    
    for extractor in extractors:
        if extractor.barcode:
            anchor_set.add(extractor.barcode)
        else:
            putative_set.update(extractor.putative_barcode)

    anchor_index = build_anchor_index(anchor_set, barcode_len, k)
    
    edge_list = []
    for barcode in putative_set:
        edges = find_edges(barcode, encode(barcode), anchor_index, barcode_len, k)
        edge_list.extend(edges)
