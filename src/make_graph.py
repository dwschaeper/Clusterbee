from collections import defaultdict, Counter
import igraph as ig
import leidenalg as la
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
import ruptures as rpt

logger = logging.getLogger(__name__)
ig.set_random_seed = 357
random.seed = 357
np.random.seed = 357
seed = 357


def plot_knee(put_community_freq: dict[str, int]):
    """
    Finds realistic communities from the set of communities.
    Calculates knee on non-log (linear) scale.

    Args:
        put_community_freq (dict): A dict of all the communities: freq

    Returns:
        int: cutoff to use for community freq
    """
    freqs = np.array(sorted(put_community_freq.values(), reverse=True), dtype=float)
    ranks = np.arange(1, len(freqs) + 1, dtype=float)

    if len(freqs) == 0:
        return 0
    if len(freqs) < 3:
        return int(freqs[0])

    # Calculate knee using linear space (non-log)
    start_points = np.array([ranks[0], freqs[0]], dtype=float)
    end_points = np.array([ranks[-1], freqs[-1]], dtype=float)
    line_vector = end_points - start_points
    norm = np.linalg.norm(line_vector)
    
    if norm == 0.0:
        knee_index = 0
    else:
        line_unit = line_vector / norm
        points = np.column_stack((ranks, freqs))
        point_vectors = points - start_points
        projection_lengths = np.dot(point_vectors, line_unit)
        projected_points = np.outer(projection_lengths, line_unit) + start_points
        distances = np.linalg.norm(points - projected_points, axis=1)
        knee_index = int(np.argmax(distances))
    
    cutoff_freq = int(freqs[knee_index])

    # plot
    plt.figure(figsize=(6,5))
    plt.plot(ranks, freqs)
    plt.scatter(ranks[knee_index], cutoff_freq, s=80)
    plt.axvline(ranks[knee_index], linestyle='--')
    plt.axhline(cutoff_freq, linestyle='--')
    plt.xlabel('Barcode Rank')
    plt.ylabel('Read Count')
    plt.title('Anchor Barcode Knee Plot')
    plt.tight_layout()
    plt.savefig('anchor_barcode_knee_plot.png')
    plt.close()

    return cutoff_freq


def PELT_changepoints(put_community_freq):
    sizes = np.array(sorted(put_community_freq.values(), reverse=True), dtype=float)
    signal = sizes.reshape((-1, 1))

    # TRY MIXTURE MODELS..GULP

    model = rpt.Pelt(model='rbf').fit(signal)
    penalty = 15
    changepoints = model.predict(pen=penalty)
    logger.info('PELT breakpoints p=15: %s', str(changepoints))

    return changepoints[4]   # TODO: make method to grab a breakpoint


def filter_anchorless(graph, anchorless):
    delete_nodes = []
    for name, community in zip(graph.vs['name'], graph.vs['community']):
        if community in anchorless:
            delete_nodes.append(name)
    logger.info('Number of nodes in anchorless communities that will be deleted: %d', len(delete_nodes))
    graph.delete_vertices(delete_nodes)


def get_anchor_communities(graph):
    anchor_communities = set()
    anchorless_communities = set(graph.vs['community'])

    for community, node_type in zip(graph.vs['community'], graph.vs['type']):
        if node_type == 'anchor':
            anchor_communities.add(community)
    anchorless_communities -= anchor_communities

    return anchor_communities, anchorless_communities


def knee_filter(graph, community_freq, cutoff):
    delete_nodes = []
    for name, community in zip(graph.vs['name'], graph.vs['community']):
        if community_freq[community] < cutoff:
            delete_nodes.append(name)

    graph.delete_vertices(delete_nodes)


def write_info(graph, header):
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
        

def add_node_attributes(graph, put_anchor_freq, put_barcode_freq):
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


def refine_multi_anchor_communities(graph):
    """Refine Ledien communities for communities with >1 anchors.

    Args:
        graph (iGraph): initial graph with Leiden clusters
    """
    initial_membership = graph.vs['community']
    anchor_counts = defaultdict(int)

    for community, node_type in zip(initial_membership, graph.vs['type']):
        if node_type == 'anchor':
            anchor_counts[community] += 1

    target_communities = {community for community, count in anchor_counts.items() if count > 1}
    if not target_communities:
        logger.info('No multi-anchor Louvain communities found for Leiden refinement.')
        return

    logger.info('Refining %d multi-anchor Louvain communities with Leiden.', len(target_communities))
    refined_membership = list(initial_membership)
    next_label = max(initial_membership) + 1
    refined_count = 0

    for community in sorted(list(target_communities)):
        vertex_indices = [idx for idx, label in enumerate(initial_membership) if label == community]
        if len(vertex_indices) < 2:
            continue

        subgraph = graph.subgraph(sorted(vertex_indices))
        #partition = subgraph.community_multilevel(weights='weight')
        partition = la.find_partition(graph, la.ModularityVertexPartition, seed=seed, weights='weight')

        local_membership = partition.membership
        local_clusters = sorted(set(local_membership))
        if len(local_clusters) <= 1:
            continue

        cluster_to_global = {cluster_id: next_label + offset for offset, cluster_id in enumerate(local_clusters)}
        next_label += len(local_clusters)

        for sub_idx, vertex_idx in enumerate(vertex_indices):
            refined_membership[vertex_idx] = cluster_to_global[local_membership[sub_idx]]

        refined_count += 1

    graph.vs['community'] = refined_membership
    logger.info('Leiden refinement split %d Louvain communities.', refined_count)


def make_graph(put_anchor_freq, put_barcode_freq, knee_scale: float = 1.0, diagnostics: bool = False):
    """Build graph, run Louvain with targeted Leiden refinement, apply knee cutoff, and optionally write diagnostics.

    Args:
        put_anchor_freq (dict):  anchor barcode -> read count
        put_barcode_freq (dict): putative barcode -> read count
        knee_scale (float):      multiplier applied to the knee cutoff (use <1.0 to be less restrictive)
        diagnostics (bool):      when True, write `community_summary.tsv` and extra plots
    """
    logger.info('making the graph')
    edges = []
    weights = []
    nodes = set()
    anchor_nodes = set(put_anchor_freq.keys())
    anchor_anchor_edges = 0

    with open('.graph_edges') as graph_file:
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

    add_node_attributes(graph, put_anchor_freq, put_barcode_freq)

    communities = la.find_partition(graph, la.CPMVertexPartition, seed=seed, weights='weight', resolution_parameter=0.5)
    communities = graph.community_multilevel(weights='weight')
    graph.vs['community'] = communities.membership
    write_info(graph, '\nCOMMUNITY STATS AFTER INITIAL CLUSTERING\n')
    #refine_multi_anchor_communities(graph)
    #write_info(graph, '\nCOMMUNITY STATS AFTER LEDIEN REFINEMENT\n')

	# write graph pre-filtering
    graph.write('graph_edges.ncol', weights='weight')
    graph_data = pd.DataFrame({'name': graph.vs['name'], 'community': graph.vs['community'], 'type': graph.vs['type'], 'freq': graph.vs['freq']})
    graph_data.to_csv('node_attributes.tsv', sep='\t', index=False)

    anchored_communities, anchorless_communities = get_anchor_communities(graph)
    filter_anchorless(graph, anchorless_communities)
    logger.info('Number of communities with anchors: %d', len(anchored_communities))
    logger.info('Number of communities without anchors: %d', len(anchorless_communities))

    # compute per-community total reads
    community_freq = defaultdict(int)
    for community, freq in zip(graph.vs['community'], graph.vs['freq']):
        community_freq[community] += freq

    with open('community_sizes.csv', 'w') as out:
        for name, community, node_type in zip(graph.vs['name'], graph.vs['community'], graph.vs['type']):
            if node_type == 'anchor':
                out.write(f'{name},{community},{community_freq[community]}\n')

    # original knee cutoff (from existing algorithm)
    cutoff = plot_knee(community_freq)
    cutoff = PELT_changepoints(community_freq)
    logger.info('Cutoff threshold from knee plot of filtered community sizes: %d', cutoff)

    # scale cutoff by user-provided multiplier (allows being less/more restrictive)
    #scaled_cutoff = max(1, int(round(cutoff * knee_scale)))
    logger.info('Scaled cutoff (knee_scale=%0.2f): %d', knee_scale, cutoff)

    # diagnostics: write community summary before filtering so you can inspect effect of cutoff
    if diagnostics:
        anchor_count = dict.fromkeys(set(graph.vs['community']), 0)
        putative_count = defaultdict(int)
        for community, node_type in zip(graph.vs['community'], graph.vs['type']):
            if node_type == 'anchor':
                anchor_count[community] += 1
            else:
                putative_count[community] += 1

        df = pd.DataFrame({
            'community': list(community_freq.keys()),
            'total_reads': [community_freq[c] for c in community_freq.keys()],
            'anchor_count': [anchor_count.get(c, 0) for c in community_freq.keys()],
            'putative_count': [putative_count.get(c, 0) for c in community_freq.keys()],
            'passes_scaled_cutoff': [community_freq[c] >= scaled_cutoff for c in community_freq.keys()]
        })
        df.sort_values('total_reads', ascending=False).to_csv('community_summary.tsv', sep='\t', index=False)
        logger.info('Wrote community diagnostics to community_summary.tsv')

    # apply scaled cutoff
    #knee_filter(graph, community_freq, scaled_cutoff)

    #write_info(graph, '\nCOMMUNITY STATS AFTER KNEE FILTERING\n')
    graph_data = pd.DataFrame({'name': graph.vs['name'], 'community': graph.vs['community'], 'type': graph.vs['type'], 'freq': graph.vs['freq']})
    return graph_data

