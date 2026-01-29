from collections import defaultdict
from itertools import combinations
import math

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
        defaultdict(list): Index of anchors kmer: anchor_id
        dict:              anchor_id: encoded
    """
    index = defaultdict(list)
    anchor_encodings = {}

    for name, barcode in anchors.items():
        encoded = encode(barcode)
        anchor_encodings[name] = encoded

        for i in range(barcode_len - k + 1):
            kmer = barcode[i:i + k]
            index[kmer].append(name)
    
    return index, anchor_encodings


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


def find_edges(barcode: str, barcode_encoding: int, index, barcode_len: int, k: int, readID: str, min_hits: int, max_mismatches: int, anchor_encodings: list):
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
    hit_counts = defaultdict(int)
    edges = []
    index_local = index
    anchor_encodings_local = anchor_encodings

    for i in range(barcode_len - k + 1):
        kmer = barcode[i:i + k]
        for anchor_id in index_local.get(kmer, ()):
            hit_counts[anchor_id] += 1

    for anchor_id, hits in hit_counts.items():
        if hits < min_hits:
            continue

        min_possible_mismatches = math.ceil(barcode_len - hits * k) / k
        if min_possible_mismatches > max_mismatches:
            continue

        dist = hamming(anchor_encodings_local[anchor_id], barcode_encoding)
        if dist <= max_mismatches:
            edges.append((anchor_id, readID, dist))
        
    return edges


def graph_deconvolute(anchors: dict, barcode_len: int, k: int, min_hits: int, max_mismatches: int):
    """
    Takes a list of identified barcode anchors and uses them with a graph to assign ambiguous cell barcodes

    Args:
        anchors (dict): Dict of anchors
        barcode_len (int): Length of the barcodes
        k (int):           Size of the kmer used in indexing
    """
    anchor_index, anchor_encodings = build_anchor_index(anchors, barcode_len, k)

    with open('.graph_edges', 'w') as graph_edges:
        with open('.putative_barcodes', 'r') as putatives:
            for line in putatives:
                data = line.strip().split('\t')
                readID = data[0]
                putative_barcode = data[1]
                encoded = encode(putative_barcode)
                edges = find_edges(putative_barcode, encoded, anchor_index, barcode_len, k, readID, min_hits, max_mismatches, anchor_encodings)
                for edge in edges:
                    graph_edges.write(f'{edge[0]}\t{edge[1]}\t{edge[2]}\n')


