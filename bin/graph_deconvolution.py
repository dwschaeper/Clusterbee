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


def build_anchor_index(anchors: dict, barcode_len: int, k: int):
    """
    Create kmer index over the anchor barcodes

    Args:
        anchors (dict):    The dictionary of anchors
        barcode_len (int): Length of the barcodes
        k (int):           Size of the kmer

    Returns:
        defaultdict(list): Index of anchors kmer: (barcode, encoded barcode)
    """
    index = defaultdict(list)

    for name, barcode in anchors.items():
        encoded = encode(barcode)
        for i in range(barcode_len - k + 1):
            kmer = barcode[i:i + k]
            index[kmer].append((name, encoded))
    
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


def find_edges(barcode: str, barcode_encoding: int, index, barcode_len: int, k: int, readID: str):
    """
    Identify close anchors as candidates, then find distance between them.

    Args:
        barcode: ambiguous barcode
        barcoding_encoding: integer encoding of the barcode
        index: index of the anchor barcodes
        barcode_len (int): Length of the barcodes
        k (int):           Size of the kmer
        readID (str):      ReadID corresponding to the barcode

    Returns:
        list: edges for the graph
    """
    candidates = set()
    edges = []

    for i in range(barcode_len - k + 1):
        kmer = barcode[i:i + k]
        candidates.update(index.get(kmer, ()))

    for anchor in candidates:
        dist = hamming(anchor[1], barcode_encoding)
        if dist < 3:
            edges.append((anchor[0], readID, 1 / dist))
        
    return edges


def graph_deconvolute(anchors: dict, barcode_len: int, k: int):
    """
    Takes a list of identified barcode anchors and uses them with a graph to assign ambiguous cell barcodes

    Args:
        anchors (dict): Dict of anchors
        barcode_len (int): Length of the barcodes
        k (int):           Size of the kmer used in indexing
    """
    anchor_index = build_anchor_index(anchors, barcode_len, k)

    edge_list = []
    with open('.putative_barcodes', 'r') as putatives:
        for line in putatives:
            data = line.strip().split('\t')
            readID = data[0]
            putative_barcode = data[1]
            edges = find_edges(putative_barcode, encode(putative_barcode), anchor_index, barcode_len, k, readID)
            edge_list.extend(edges)

    print(edge_list)

