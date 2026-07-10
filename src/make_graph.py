"""
Build and annotate the barcode graph prior to read assignment.

This module consumes edge lists emitted during barcode matching and performs
Leiden community detection with anchor-presence filtering.
"""

from collections import Counter, defaultdict
import logging
import os

import igraph as ig
import leidenalg as la
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
import ruptures as rpt
from scipy.special import logsumexp
from scipy.stats import nbinom

import config

logger = logging.getLogger(__name__)


def _get_anchor_communities(graph: ig.Graph) -> tuple[set[int], set[int]]:
    """
    Split communities into those with anchors and those without.

    Args:
        graph (ig.Graph): Graph with ``community`` and ``type`` vertex
            attributes.

    Returns:
        tuple[set[int], set[int]]: ``(anchored_communities, anchorless_communities)``.
    """
    anchor_communities = set()
    anchorless_communities = set(graph.vs['community'])

    for community, node_type in zip(graph.vs['community'], graph.vs['type']):
        if node_type == 'anchor':
            anchor_communities.add(community)
    anchorless_communities -= anchor_communities

    return anchor_communities, anchorless_communities


def _write_info(graph: ig.Graph, header: str) -> None:
    """Log per-community summary statistics.

    Args:
        graph (ig.Graph): Graph with ``community``, ``type``, and ``freq``
            attributes.
        header (str): Header line emitted before statistics.
    """
    anchor_freq = dict.fromkeys(set(graph.vs['community']), 0)
    reads_freq = defaultdict(int)
    node_count = defaultdict(int)
    communities = set()

    for community, node_type, freq in zip(graph.vs['community'], graph.vs['type'], graph.vs['freq']):
        if node_type == 'anchor':
            anchor_freq[community] += 1
        else:
            node_count[community] += 1
        reads_freq[community] += freq
        communities.add(community)

    no_anchor_communities = 0
    one_anchor_communities = 0
    multi_anchor_communities = 0
    for community, anchor_count in anchor_freq.items():
        if anchor_count == 0:
            no_anchor_communities += 1
        elif anchor_count == 1:
            one_anchor_communities += 1
        else:
            multi_anchor_communities += 1

    logger.info(header + '\n'
               '  Number of communities              : %d\n'
               '  Number of zero anchor communities  : %d\n'
               '  Number of one anchor communities   : %d\n'
               '  Number of multi anchor communities : %d\n'
               '  Number of anchors per community    :\n\t%s\n'
               '  Barcodes per community distribution:\n\t%s',
               len(communities), no_anchor_communities, one_anchor_communities, 
               multi_anchor_communities, str(Counter(anchor_freq.values())), str(Counter(node_count.values())))   
        

def _add_node_attributes(graph: ig.Graph, put_anchor_freq: dict[str, int], put_barcode_freq: dict[str, int]) -> None:
    """
    Populate frequency/type/degree node attributes used downstream.

    Args:
        graph (ig.Graph): Graph containing barcode names in
            ``graph.vs['name']``.
        put_anchor_freq (dict[str, int]): Anchor barcode frequencies.
        put_barcode_freq (dict[str, int]): Putative barcode frequencies.
    """
    put_anchors = set(put_anchor_freq.keys())
    freqs = []

    # set barcode usage frequency
    for name in graph.vs['name']:
        if name in put_anchors:
            freqs.append(put_anchor_freq[name])
        else:
            freqs.append(put_barcode_freq[name])
    graph.vs['freq'] = freqs
            
    # set barcode type, anchor/putative
    graph.vs['type'] = ['anchor' if name in put_anchors else 'putative' for name in graph.vs['name']]

    # find log(frequency)
    graph.vs['log_freq'] = [np.log10(x) for x in graph.vs['freq']]

    # find degrees
    graph.vs['degree'] = graph.degree()


def _filter_anchorless(graph: ig.Graph, anchorless: set[int]) -> None:
    """Delete vertices belonging to communities that lack anchors.

    Args:
        graph (ig.Graph): Graph to mutate in place.
        anchorless (set[int]): Community identifiers to remove.
    """
    delete_nodes = []
    for name, community in zip(graph.vs['name'], graph.vs['community']):
        if community in anchorless:
            delete_nodes.append(name)
    logger.info('Number of nodes in anchorless communities that will be deleted: %d', len(delete_nodes))
    graph.delete_vertices(delete_nodes)


def make_graph(put_anchor_freq: dict[str, int], put_barcode_freq: dict[str, int], resolution: float, iterations: int) -> dict[str, list]:
    """Build graph, run Louvain with targeted Leiden refinement, apply knee cutoff, and optionally write diagnostics.

    Args:
        put_anchor_freq (dict):  anchor barcode -> read count
        put_barcode_freq (dict): putative barcode -> read count
        resolution (float):      resolution for leiden calculation
        iterations (int):        number of iterations used in leiden clustering
    """
    # initialize values
    seed = 357
    logger.section('GRAPH CONSTRUCTION')
    logger.info('making the graph')
    edges = []
    weights = []
    nodes = set()
    anchor_nodes = set(put_anchor_freq.keys())
    anchor_anchor_edges = 0

    # initialize file handles
    graph_edge_file = os.path.join(config.GLOBAL_OUTDIR, '.graph_edges')
    community_size_csv = os.path.join(config.GLOBAL_OUTDIR, f'{config.GLOBAL_PREFIX}_community_sizes.csv')

    with open(graph_edge_file) as graph_file:
        for line in graph_file:
            node1, node2, weight = line.strip().split('\t')
            edges.append((node1, node2))
            weights.append(float(weight))
            nodes.add(node1)
            nodes.add(node2)
            if node1 in anchor_nodes and node2 in anchor_nodes:
                anchor_anchor_edges += 1
    
    graph = ig.Graph.TupleList(edges, directed=False)
    graph.es['weight'] = weights
    logger.info('Number of anchor-anchor edges: %d', anchor_anchor_edges)

    _add_node_attributes(graph, put_anchor_freq, put_barcode_freq)

    logger.info('resolution parameter: ' + str(resolution))
    logger.info('iterations: ' + str(iterations))
    communities = la.find_partition(graph, la.CPMVertexPartition, seed=seed, weights='weight', n_iterations=iterations, resolution_parameter=resolution)
    graph.vs['community'] = communities.membership
    _anchored_communities, anchorless_communities = _get_anchor_communities(graph)
    if len(anchorless_communities) > 0:
        _filter_anchorless(graph, anchorless_communities)
    _write_info(graph, '\nCOMMUNITY STATS\n')

    # compute per-community total reads
    community_freq = defaultdict(int)
    for community, freq in zip(graph.vs['community'], graph.vs['freq']):
        community_freq[community] += freq

    with open(community_size_csv, 'w') as out:
        for name, community, node_type in zip(graph.vs['name'], graph.vs['community'], graph.vs['type']):
            if node_type == 'anchor':
                out.write(f'{name},{community},{community_freq[community]}\n')

    return {
        'name': list(graph.vs['name']),
        'community': list(graph.vs['community']),
        'type': list(graph.vs['type']),
        'freq': list(graph.vs['freq']),
    }
