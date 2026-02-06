from collections import defaultdict, Counter
import igraph as ig
import math


def remove_putative_singletons(graph):
    deg = graph.degree()

    delete_idx = [
        v.index
        for v, d in zip(graph.vs, deg)
        if d == 0 and v['type'] == 'putative'
    ]

    graph.delete_vertices(delete_idx)


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

    print(f'Reads per community distribution:\n{Counter(reads_freq.values())}')  
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
    graph.vs['log_freq'] = [math.log(x) for x in graph.vs['freq']]

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
    print('remove singletons')
    remove_putative_singletons(graph)
    print('making community')
    communities = graph.community_leiden(weights='weight')
    graph.vs['community'] = communities.membership
    

    write_info(graph, '\nCOMMUNITY STATS\n')
    

        
