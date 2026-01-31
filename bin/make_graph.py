import igraph as ig

def make_graph():
    graph = ig.Graph.Read_Ncol('.graph_edges', directed=False, weights=True)
    clusters = graph.community_leiden(weights='weight', resolution=1)
