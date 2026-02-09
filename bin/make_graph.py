from collections import defaultdict, Counter
import igraph as ig
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def find_supported_communities_cutoff(put_community_freq: dict[str, int]):
    """
    Finds realistic communities from the set of communities.

    Args:
        put_community_freq (dict): A dict of all the communities: freq

    Returns:
        int: cutoff to use for community freq
    """
    freqs = np.array(sorted(put_community_freq.values(), reverse=True))
    ranks = np.arange(1, len(freqs) + 1)

    # get line from highest to lowest
    log_ranks = np.log10(ranks)
    log_freqs = np.log10(freqs)
    start_points = np.array([log_ranks[0], log_freqs[0]])
    end_points = np.array([log_ranks[-1], log_freqs[-1]])
    line_vector = end_points - start_points
    line_unit = line_vector / np.linalg.norm(line_vector)
    
    # calculate cutoff
    points = np.column_stack((log_ranks, log_freqs))
    point_vectors = points - start_points
    projection_lengths = np.dot(point_vectors, line_unit)
    projected_points = np.outer(projection_lengths, line_unit) + start_points
    distances = np.linalg.norm(points - projected_points, axis=1)
    knee_index = np.argmax(distances)
    cutoff_freq = freqs[knee_index]

    # plot
    plt.figure(figsize=(6,5))
    plt.plot(ranks, freqs)
    plt.scatter(ranks[knee_index], cutoff_freq, s=80)
    plt.axvline(ranks[knee_index], linestyle='--')
    plt.axhline(cutoff_freq, linestyle='--')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Barcode Rank')
    plt.ylabel('Read Count')
    plt.title('Anchor Barcode Knee Plot')
    plt.tight_layout()
    plt.savefig('anchor_barcode_knee_plot.png')
    plt.close()

    return cutoff_freq


def filter_anchorless(graph, anchorless):
    delete_nodes = []
    for name, community in zip(graph.vs['name'], graph.vs['community']):
        if community in anchorless:
            delete_nodes.append(name)
    print(f'Number of nodes in anchorless communities: {len(delete_nodes)}')
    graph.delete_vertices(delete_nodes)


def get_anchor_communities(graph):
    anchor_communities = set()
    anchorless_communities = set(graph.vs['community'])

    for community, node_type in zip(graph.vs['community'], graph.vs['type']):
        if node_type == 'anchor':
            anchor_communities.add(community)
    anchorless_communities -= anchor_communities

    return anchor_communities, anchorless_communities


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

    print(header)

    print(f'Number of communities: {len(communities)}')
    print(f'Number of zero anchor communities: {no_anchor_communities}')
    print(f'Number of one anchor communities: {one_anchor_communities}')
    print(f'Number of multi anchor communities: {multi_anchor_communities}')

    #print(f'Reads per community distribution:\n{Counter(reads_freq.values())}')  
    print(f'Barcodes per community distribution:\n{Counter(node_count.values())}')   
        

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


def make_graph(put_anchor_freq, put_barcode_freq):
    edges = []
    weights = []
    nodes = set()

    print('making graph')
    with open('.graph_edges') as graph_file:
        for line in graph_file:
            node1, node2, weight = line.strip().split('\t')
            edges.append((node1, node2))
            weights.append(float(weight))
            nodes.add(node1)
            nodes.add(node2)
    
    graph = ig.Graph.TupleList(edges, directed=False)
    graph.es['weight'] = weights

    print('adding attributes')
    add_node_attributes(graph, put_anchor_freq, put_barcode_freq)
    # degree distribution
    print(Counter(graph.vs['degree']))
    #print('remove singletons')
    #remove_putative_singletons(graph)

    print('making community')
    communities = graph.community_leiden(weights='weight')
    graph.vs['community'] = communities.membership

    print('filtering cells')
    anchored_communities, anchorless_communities = get_anchor_communities(graph)
    filter_anchorless(graph, anchorless_communities)
    print(f'Number of communities with anchors: {len(anchored_communities)}')
    print(f'Number of communities without anchors: {len(anchorless_communities)}')

    community_freq = defaultdict(int)
    for community, freq in zip(graph.vs['community'], graph.vs['freq']):
        community_freq[community] += freq
    print(len(community_freq))
    cutoff = find_supported_communities_cutoff(community_freq)
    print(f'Cutoff threshold: {cutoff}')

    write_info(graph, '\nCOMMUNITY STATS\n')

    graph.write('graph_edges.ncol', weights='weight')
    #pd.DataFrame({'name': graph.vs['name'], 'community': graph.vs['community'], 'type': graph.vs['type'], 'freq': graph.vs['freq']}).to_csv('node_attributes.tsv', sep='\t', index=False)
    

        
