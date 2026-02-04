import igraph as ig
import math

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

    
    with open('.graph_edges') as graph_file:
        for line in graph_file:
            node1, node2, weight = line.strip().split('\t')
            edges.append((node1, node2))
            weights.append(float(weight))
            nodes.add(node1)
            nodes.add(node2)
    
    graph = ig.Graph.TupleList(edges, directed=False)
    graph.es['weight'] = weights

    add_node_attributes(graph, put_anchor_freq, put_barcode_freq)

    communities = graph.community_leiden(weights='weight')
    graph.vs['community'] = communities.membership
