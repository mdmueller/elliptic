# Enumerate possible graphs

import networkx as nx
import copy
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
        H = copy.deepcopy(G)
        # part lists the edge weights coming from the ith left vertex
        for w, j in enumerate(part):
            if j!=0: # ignore edges with weight 0
                H.add_edge('l{0}'.format(i), 'r{0}'.format(w), weight=j)
                H.nodes['l{0}'.format(i)]['R'] -= (j-1)
                H.nodes['r{0}'.format(w)]['R'] -= (j-1)
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
                #weight_from_v += sum([x['weight'] for x in G['l{0}'.format(i)]['r{0}'.format(j)].values()]) #for multigraph...
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
        G.nodes['l{0}'.format(i)]['ramif'] = [[],[],[]]
        G.nodes['l{0}'.format(i)]['genus'] = 1 if i==0 else 0
    for i in range(len(mu2)):
        G.add_node('r{0}'.format(i), bipartite=1)
        G.nodes['r{0}'.format(i)]['degree'] = mu2[i]
        G.nodes['r{0}'.format(i)]['R'] = 2*mu2[i]-2
        G.nodes['r{0}'.format(i)]['ramif'] = [[],[],[]]
        G.nodes['r{0}'.format(i)]['genus'] = 0

    trees = bipartite_tree_helper(G, 0, mu1, mu2)
    return [T for T in trees if weights_add_up(T, mu1, mu2) and all([T.nodes[node]['R']>=0 for node in T.nodes()])]

def display_bipartite(G, mu1, mu2):
    # draw the bipartite weighted graph G
    L = [node for node in G.nodes() if 'l' in node]
    R = [node for node in G.nodes() if 'r' in node]
    #L, R = nx.bipartite.sets(G)
    pos = {}
    pos.update((node, (1, index)) for index, node in enumerate(L))
    pos.update((node, (2, index)) for index, node in enumerate(R))
    nx.draw_networkx_nodes(G, pos, node_size=4000, node_color='tab:red')
    nx.draw_networkx_edges(G, pos)
    labels1 = {'l{0}'.format(i): 'g={0},d={1}'.format(1 if i==0 else 0, mu1[i]) for i in range(len(mu1))}
    labels2 = {'r{0}'.format(i): 'g=0,d={0}'.format(mu2[i]) for i in range(len(mu2))}
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

def place_ramification_helper(T, sigma, sigmacount, side):
    # side = 'r' or 'l'
    # sigmacount = i if this is sigma_i
    # return all placements of sigma ramification on the given side of T
    if len(sigma)==0:
        return [T]
    ram = sigma[0]
    L = []
    for node in T.nodes():
        if side not in node or T.nodes[node]['R'] < ram-1 or sum(T.nodes[node]['ramif'][sigmacount])+ram>T.nodes[node]['degree']:
            continue
        # try putting ram at this node
        T2 = copy.deepcopy(T)#T.copy()
        T2.nodes[node]['R'] -= (ram-1)
        T2.nodes[node]['ramif'][sigmacount].append(ram)
        #T2.add_edge(node, node, weight=ram)
        L.extend(place_ramification_helper(T2, sigma[1:], sigmacount, side))
    return L

def place_ramification(T, sigmas):
    # T is a weighted bipartite tree
    # sigma_i are partitions; we want these ramification to appear in T
    # Our goal is to allocate ramification across all vertices in T
    L = [T]
    for i, sigma in enumerate(sigmas):
        L2 = []
        for side in ['l', 'r']:
            for T2 in L:
                L2.extend(place_ramification_helper(T2, sigma, i, side))
        L = L2

    return {tuple([tuple([tuple(t) for t in G.nodes[node]['ramif']]) for node in G.nodes()]):G for G in L}.values()

def stabilization(T):
    # produce the stabilization of this graph
    # TODO: adapt this for graphs which have double edges
    T = copy.deepcopy(T)
    while True:
        for node in T.nodes():
            if T.nodes[node]['genus']==1:
                continue
            marked_pts = T.nodes[node]['ramif'][0]
            N = list(T.neighbors(node))
            if len(N) == 1 and len(marked_pts)<=1: # delete node, put its marked points on the neighbor
                neighbor = N[0]
                T.nodes[neighbor]['ramif'][0].extend(T.nodes[node]['ramif'][0])
                T.remove_node(node)
                break
            elif len(N) == 2 and len(marked_pts)==0: # delete node, connect neighbors
                T.add_edge(N[0], N[1])
                T.remove_node(node)
                break
        else:
            break
    return T

def isomorphic(T, G, all_markings=True):
    # check if T looks like G
    # if all_markings, look at all ramification points; otherwise just mu_0
    def matching(n1, n2): # TODO: also check genus matches up?
        return all([n1['ramif'][i] == n2['ramif'][i] for i in range(3 if all_markings else 1)]) and n1['genus'] == n2['genus']
    return nx.is_isomorphic(T, G, node_match=matching)
    

def possible_graphs(sigmas, G):
    # consider possible graphs with ramification sigma_0, sigma_1, sigma_2, ...
    # stabilization should look like G
    d = sum(sigmas[0])
    #n = sum([len(x) for x in sigmas]) # number of marked points
    assert(all([d == sum(sigma) for sigma in sigmas]))
    #assert(sum([sum([x-1 for x in sigma]) for sigma in sigmas])==2*d) # all ramification included
    TOTAL = 0
    seen_graphs = []
    
    for mu1 in Part(d):
        if mu1[0]==1:
            continue # genus 1 component can't be degree 1
        for mu2 in Part(d):
            # TODO: check if the sigmas are possible just from mu1 and mu2?
            trees = bipartite_trees(mu1, mu2)
            for T in trees:
                for T2 in place_ramification(T, sigmas):
                    if not isomorphic(stabilization(T2), G, all_markings=False):
                        continue
                    elif any([isomorphic(T2, Q) for Q in seen_graphs]):
                        continue
                    TOTAL += 1
                    seen_graphs.append(T2)
                    for i in range(len(mu1)):
                        label = 'l{}'.format(i)
                        print(label, T2.nodes[label]['degree'], T2.nodes[label]['ramif'])
                    for i in range(len(mu2)):
                        label = 'r{}'.format(i)
                        print(label, T2.nodes[label]['degree'], T2.nodes[label]['ramif'])
                    display_bipartite(T2, mu1, mu2)
                    print('-'*10)
                    # TODO: display the ramification as legs or loops on the graph...
    return TOTAL

if __name__ == '__main__':
    G = nx.Graph()
    G.add_edge('l0','r0')
    G.nodes['l0']['ramif']=[[2,1],[],[]]
    G.nodes['l0']['genus']=1
    G.nodes['r0']['ramif']=[[1,1],[],[]]
    G.nodes['r0']['genus']=0
    possible_graphs([[2,1,1,1],[4,1],[3,1,1]], G)
