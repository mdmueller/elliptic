# Enumerate possible graphs

import networkx as nx
from matplotlib import pyplot as plt

def partitions(n, m):
    # consider all partitions of n of length m (possibly including zeros)
    if m==0:
        return [[]] if n==0 else []
    l = []
    for i in range(n+1):
        for P in partitions(n-i, m-1):
            l.append([i]+P)
    return l

def bipartite_tree_helper(G, i, mu1, mu2):
    # G is a partially constructed bipartite tree, after i left vertices have been dealt with
    # return all possible weighted bipartite tree completions of G
    if i==len(mu1):
        return [G] if nx.is_connected(G) else []

    L = []
    for part in partitions(mu1[i], len(mu2)):
        H = G.copy()
        # part lists the edge weights coming from the ith left vertex
        for w, j in enumerate(part):
            if j!=0: # ignore edges with weight 0
                H.add_edge('l{0}'.format(i), 'r{0}'.format(w), weight=j)
        # keep going unless H has a cycle
        try:
            nx.find_cycle(H)
            continue
        except nx.exception.NetworkXNoCycle:
            L.extend(bipartite_tree_helper(H, i+1, mu1, mu2))
    return L

def weights_add_up(G, mu1, mu2):
    # check if the weights add up for each vertex on the right
    for j in range(len(mu2)):
        weight_from_v = 0 # v is the jth vertex on the right
        for i in range(len(mu1)):
            try:
                weight_from_v += G['l{0}'.format(i)]['r{0}'.format(j)]['weight']
            except KeyError:
                continue
        if weight_from_v != mu2[j]:
            return False
    return True

def bipartite_trees(mu1, mu2):
    # given partitions mu1 and mu2 of d, list potential bipartite weighted trees
    d = sum(mu1)
    assert(d == sum(mu2))
    G = nx.Graph()
    for i in range(len(mu1)):
        G.add_node('l{0}'.format(i), bipartite=0)
        G.nodes['l{0}'.format(i)]['degree'] = mu1[i]
        # 2g - 2 = (-2)d + R => R = 2d + 2g - 2
        G.nodes['l{0}'.format(i)]['R'] = 2*mu1[i]-2 if i!=0 else 2*mu1[i]
        G.nodes['l{0}'.format(i)]['ramif'] = []
    for i in range(len(mu2)):
        G.add_node('r{0}'.format(i), bipartite=1)
        G.nodes['r{0}'.format(i)]['degree'] = mu2[i]
        G.nodes['r{0}'.format(i)]['R'] = 2*mu2[i]-2
        G.nodes['r{0}'.format(i)]['ramif'] = []
    
    trees = bipartite_tree_helper(G, 0, mu1, mu2)
    return [T for T in trees if weights_add_up(T, mu1, mu2)]

def display_bipartite(G, mu1, mu2):
    # draw the bipartite weighted graph G
    L, R = nx.bipartite.sets(G)
    pos = {}
    pos.update((node, (1, index)) for index, node in enumerate(L))
    pos.update((node, (2, index)) for index, node in enumerate(R))
    nx.draw_networkx_nodes(G, pos, node_size=1000)
    nx.draw_networkx_edges(G, pos)
    labels1 = {'l{0}'.format(i): mu1[i] for i in range(len(mu1))}
    labels2 = {'r{0}'.format(i): mu2[i] for i in range(len(mu2))}
    nx.draw_networkx_labels(G, pos, labels=labels1|labels2)
    edge_labels = nx.get_edge_attributes(G, 'weight')
    nx.draw_networkx_edge_labels(G, pos, edge_labels, label_pos=0.2)
    plt.show()

def part(n, k):
    def _part(n, k, pre):
        if n <= 0:
            return []
        if k == 1:
            if n <= pre:
                return [[n]]
            return []
        ret = []
        for i in range(min(pre, n), 0, -1):
            ret += [[i] + sub for sub in _part(n-i, k-1, i)]
        return ret
    return _part(n, k, n)

def Part(n):
    L = []
    for k in range(1, n+1):
        L.extend(part(n, k))
    return L

def place_ramification_helper(T, sigma, side):
    # side = 'r' or 'l'
    # return all placements of sigma ramification on the given side of T
    return 1

def place_ramification(T, sigma0, sigma1, sigma2):
    # T is a weighted bipartite tree
    # sigma_i are partitions; we want these ramification to appear in T
    # Our goal is to allocate ramification across all vertices in T
    L = [T]
    for sigma in sigma0, sigma1, sigma2:
        for side in ['l', 'r']:
            L2 = []
            for T2 in L:
                L2.extend(place_ramification_helper(T2, sigma, side))
    for node in T.nodes():
        deg = T.nodes[node]['degree']
        

def possible_graphs(sigma0, sigma1, sigma2):
    # consider possible graphs with ramification sigma0, sigma1, sigma2
    # this does not consider whether the graph stabilizes appropriately
    d = sum(sigma0)
    assert(d == sum(sigma1) and d == sum(sigma2))
    
    for mu1 in Part(d):
        for mu2 in Part(d):
            # TODO: check if the sigmas are possible just from mu1 and mu2?
            trees = bipartite_trees(mu1, mu2)
            for T in trees:
                place_ramification(T, sigma0, sigma1, sigma2)
                

T = bipartite_trees([3,2],[3,1,1])
G = T[2]
display_bipartite(G, [3,2],[3,1,1])

#possible_graphs([4],[4],[4])

