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

def bipartite_tree_helper(G, i, mu1, mu2, double=False):
    # G is a partially constructed bipartite tree, after i left vertices have been dealt with
    # return all possible weighted bipartite tree completions of G
    # if double=True, then give l_0 a double edge with r_0
    if i==len(mu1):
        return [G] if nx.is_connected(G) else []

    L = []
    for part in partitions(mu1[i], len(mu2)+1 if double else len(mu2)):
        H = copy.deepcopy(G)
        if double:
            if part[0]==0 or part[1]==0:
                continue # double edges should not have weight 0
            H.add_edge('l{0}'.format(i), 'r0', weight=part[0])
            H.add_edge('l{0}'.format(i), 'r0', weight=part[1])
            H.nodes['l{0}'.format(i)]['R'] -= (part[0]+part[1]-2)
            H.nodes['r0']['R'] -= (part[0]+part[1]-2)
        # part lists the edge weights coming from the ith left vertex
        for w, j in enumerate(part):
            if double and w<2: # already dealt with
                continue
            if j!=0: # ignore edges with weight 0
                H.add_edge('l{0}'.format(i), 'r{0}'.format(w-1 if double else w), weight=j)
                H.nodes['l{0}'.format(i)]['R'] -= (j-1)
                H.nodes['r{0}'.format(w-1 if double else w)]['R'] -= (j-1)
        # keep going unless H has a cycle (other than possible double edge)
        try:
            H2 = copy.deepcopy(H)
            if double:
                H2.remove_edge('l0','r0')
                H2.remove_edge('l0','r0')
            nx.find_cycle(H2)
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
                weight_from_v += sum([x['weight'] for x in G['l{0}'.format(i)]['r{0}'.format(j)].values()]) #for multigraph...
                #weight_from_v += G['l{0}'.format(i)]['r{0}'.format(j)]['weight']
            except KeyError:
                continue
        if weight_from_v != mu2[j]:
            return False
    return True

def bipartite_trees(mu1, mu2, genus, double=False, num_ramif=10):
    # given partitions mu1 and mu2 of d, list potential bipartite weighted trees of given genus
    # if double=True, then l0 has a double edge
    d = sum(mu1)
    assert(d == sum(mu2))
    G = nx.MultiGraph()
    L = []
    L2 = []
    for i in range(num_ramif):
        L.append([])
        L2.append([])
    for i in range(len(mu1)):
        G.add_node('l{0}'.format(i), bipartite=0)
        G.nodes['l{0}'.format(i)]['degree'] = mu1[i]
        # 2g - 2 = (-2)d + R => R = 2d + 2g - 2
        G.nodes['l{0}'.format(i)]['R'] = 2*mu1[i]-2 if (i!=0 or genus==0) else 2*mu1[i]
        G.nodes['l{0}'.format(i)]['ramif'] = tuple(L)
        # Each list in ...['ramif'] contains tuples (i,r) where i is the index of the point and r is its assigned ramification
        G.nodes['l{0}'.format(i)]['genus'] = genus if i==0 else 0
    for i in range(len(mu2)):
        G.add_node('r{0}'.format(i), bipartite=1)
        G.nodes['r{0}'.format(i)]['degree'] = mu2[i]
        G.nodes['r{0}'.format(i)]['R'] = 2*mu2[i]-2
        G.nodes['r{0}'.format(i)]['ramif'] = tuple(L2)
        G.nodes['r{0}'.format(i)]['genus'] = 0

    trees = bipartite_tree_helper(G, 0, mu1, mu2, double=double)
    return [T for T in trees if weights_add_up(T, mu1, mu2) and all([T.nodes[node]['R']>=0 for node in T.nodes()])]

def display_bipartite(G, mu1, mu2, genus):
    # draw the bipartite weighted graph G
    L = [node for node in G.nodes() if 'l' in node]
    R = [node for node in G.nodes() if 'r' in node]
    #L, R = nx.bipartite.sets(G)
    pos = {}
    pos.update((node, (1, index)) for index, node in enumerate(L))
    pos.update((node, (2, index)) for index, node in enumerate(R))
    nx.draw_networkx_nodes(G, pos, node_size=4000, node_color='tab:red')
    nx.draw_networkx_edges(G, pos)
    labels1 = {'l{0}'.format(i): 'g={0},d={1}'.format(genus if i==0 else 0, mu1[i]) for i in range(len(mu1))}
    labels2 = {'r{0}'.format(i): 'g=0,d={0}'.format(mu2[i]) for i in range(len(mu2))}
    nx.draw_networkx_labels(G, pos, labels=labels1|labels2)
    edge_labels = nx.get_edge_attributes(G, 'weight')
    # deal with multiple edges
    new_labels = {}
    for a,b,c in edge_labels:
        if (a,b) not in new_labels:
            new_labels[(a,b)]=str(edge_labels[a,b,c])
        else:
            new_labels[(a,b)]+=(','+str(edge_labels[a,b,c]))

    nx.draw_networkx_edge_labels(G, pos, new_labels, label_pos=0.2)
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
    if n==0:
        return [[]]
    L = []
    for k in range(1, n+1):
        L.extend(part(n, k))
    return L

def place_ramification_helper(T, node_dict, sigma, sigmacount, pointcount, side, num_fixed):
    # node_dict is a replacement for T.nodes
    # side = 'r' or 'l'
    # sigmacount = i if this is sigma_i
    # pointcount = j if this is the jth point
    # return all placements of sigma ramification on the given side of T
    if len(sigma)==0:
        return [node_dict]
        #return [{node: tuple([tuple(t) for t in node_dict[node]['ramif']]) for node in node_dict}]
    ram = sigma[0]
    L = []
    for node in T.nodes():
        if side not in node or node_dict[node]['R'] < ram-1 or sum([t[1] for t in node_dict[node]['ramif'][sigmacount]])+ram>node_dict[node]['degree']:
            continue
        # try putting ram at this node
        node_dict2 = {x: {y: copy.deepcopy(z) if y=='ramif' else z for y,z in node_dict[x].items()} for x in node_dict}
        node_dict2[node]['R'] -= (ram-1)
        if pointcount <= num_fixed:
            node_dict2[node]['ramif'][sigmacount].append((pointcount, ram))
        else:
            node_dict2[node]['ramif'][sigmacount].append((100, ram))
        L.extend(place_ramification_helper(T, node_dict2, sigma[1:], sigmacount, pointcount+1, side, num_fixed))
        
    return L

def place_ramification(T, sigmas, num_fixed):
    # T is a weighted bipartite tree
    # sigma_i are partitions; we want these ramification to appear in T
    # Our goal is to allocate ramification across all vertices in T
    L = [copy.deepcopy(T.nodes)]
    pointcount = 1
    for i, sigma in enumerate(sigmas):
        L2 = []
        for side in ('l', 'r'):
            for node_dict in L:
                L2.extend(place_ramification_helper(T, node_dict, sigma, i, pointcount, side, num_fixed))
        L = L2
        pointcount += len(sigma)

    # keep only node dicts with distinct node ramification profiles
    node_dicts = {tuple([tuple([tuple(t) for t in node_dict[node]['ramif']]) for node in node_dict]):node_dict for node_dict in L}.values()

    trees = []
    for d in node_dicts:
        T2 = copy.deepcopy(T)
        for node in d:
            for attrib in d[node]:
                T2.nodes[node][attrib] = d[node][attrib]
        trees.append(T2)
    return trees

def stabilizations(T, num_fixed):
    # produce (stabilization of T, all lists [e0, e1, ...] of edges to contract to stabilize T)
    L = []
    stable = True # T is already stable

    for node in T.nodes():
        if T.nodes[node]['genus']>0:
            continue
        marked_pts = sum([len([x for x in L if x[0]<=num_fixed]) for L in T.nodes[node]['ramif']]) # marked points on this node
        edges = [e for e in T.edges if node in e]
        N = list(T.neighbors(node))

        if len(edges) == 1 and N!=[node] and marked_pts<=1: # delete node, put its marked points on the neighbor
            stable = False
            neighbor = N[0]
            T2 = copy.deepcopy(T)
            T2.nodes[neighbor]['ramif'][0].extend(T.nodes[node]['ramif'][0])
            T2.remove_node(node)
            G, L2 = stabilizations(T2, num_fixed)
            L.extend([[T.edges[edges[0]]['ID']]+x for x in L2])
        elif len(edges) == 2 and marked_pts==0:
            stable = False
            # delete either edge
            for j,edge in enumerate(edges):
                n1, n2, _ = edge
                neighbor = n1 if n2==node else n2
                T2 = copy.deepcopy(T)
                T2.nodes[neighbor]['ramif'][0].extend(T.nodes[node]['ramif'][0]) # put marked pts on neighbor, delete node
                T2.remove_node(node)
                other_edge = edges[j-1]
                m1, m2, _ = other_edge
                other_neighbor = m1 if m2==node else m2
                T2.add_edge(other_neighbor, neighbor, weight=T.edges[other_edge]['weight'], ID=T.edges[other_edge]['ID'])
                G, L2 = stabilizations(T2, num_fixed)
                L.extend([[T.edges[edge]['ID']]+x for x in L2])

    if stable:
        return (T, [[]])
    return (G, L)

def stabilization(T, num_fixed):
    # produce the stabilization of this graph and the weights of contracted edges
    #TODO: make weights part work for graphs with double edges...
    T = copy.deepcopy(T)
    for i,edge in enumerate(T.edges):
        T.edges[edge]['ID'] = i
    T2, L = stabilizations(T, num_fixed)
    mult = 0
    for choices in set([frozenset(x) for x in L]):
        factor = 1
        for edge_ID in choices:
            factor *= sum([x['weight'] for x in T.edges.values() if x['ID'] == edge_ID])
        mult += factor
    return (T2, mult)

def isomorphic(T, G, num_fixed, all_markings=True):
    # check if T looks like G
    # if all_markings, look at all ramification points; otherwise just the fixed ones
    #TODO: make an option to ignore order of identifical ramifications; i.e. [[3],[],[]],[[1],[],[]] identical to [[],[3],[]],[[],[1],[]]
    def matching(n1, n2):
        def same_markings(L1, L2):
            return {x for x in L1 if all_markings or x[0]<=num_fixed} == {x for x in L2 if all_markings or x[0]<=num_fixed}
        return all([same_markings(L1,L2) for L1,L2 in zip(n1['ramif'],n2['ramif'])]) and n1['genus'] == n2['genus']
    return nx.is_isomorphic(T, G, node_match=matching)
    
def display_graph(T2, mu1, mu2, genus):
    for i in range(len(mu1)):
        label = 'l{}'.format(i)
        print(label, T2.nodes[label]['degree'], T2.nodes[label]['ramif'])
    for i in range(len(mu2)):
        label = 'r{}'.format(i)
        print(label, T2.nodes[label]['degree'], T2.nodes[label]['ramif'])
    display_bipartite(T2, mu1, mu2, genus)
    print('-'*10)

def possible_graphs_old(sigmas, G, num_fixed='auto', genus=1, double=False, ignore_labels=False): #TODO: remove
    list(possible_graphs(sigmas, G, num_fixed, genus, double, ignore_labels, display=True))

def possible_graphs(sigmas, G, num_fixed='auto', genus=1, double=False, ignore_labels=False, display=False):
    # consider possible graphs with ramification sigma_0, sigma_1, sigma_2, ...
    # stabilization should look like G
    # num_fixed is the # of fixed points, 'auto' means |sigma_0|
    # double=True means include a double edge (genus 0)
    # ignore_labels=True means don't label the ramification points
    d = sum(sigmas[0])
    assert(all([d == sum(sigma) for sigma in sigmas]))
    #assert(sum([sum([x-1 for x in sigma]) for sigma in sigmas])==2*d) # all ramification included
    if num_fixed == 'auto':
        num_fixed = len(sigmas[0])
    TOTAL = 0
    seen_graphs = []

    for mu1 in Part(d):#[[a]+y for a in range(1,d+1) for y in Part(d-a)]:
        if mu1[0]==1 and genus>0:
            continue # positive genus component can't be degree 1

        for mu2 in Part(d):
            # skip if any ramification # is larger than all degrees
            max_deg = max((max(mu1), max(mu2)))
            if max([max(sigma) for sigma in sigmas])>max_deg:
                continue

            smaller_side = min((max(mu1), max(mu2)))
            ramif_one_side = 0 # sum of all ramif that must go on the "larger" side
            for sigma in sigmas:
                if max(sigma)>smaller_side: # must go on the "larger" side
                    ramif_one_side += (sum(sigma)-len(sigma))

            Rright = 2*(1-len(mu2))-2+2*d
            Rleft = 2*(2-len(mu1))-2+2*d

            if smaller_side == max(mu1) and ramif_one_side > Rright: # doesn't fit on "larger" right side
                continue
            elif smaller_side == max(mu2) and ramif_one_side > Rleft: # doesn't fit on "larger" left side
                continue

            trees = bipartite_trees(mu1, mu2, genus, double=double, num_ramif=len(sigmas))

            for T in trees:
                for T2 in place_ramification(T, sigmas, num_fixed):
                    T2.graph['mu1'] = mu1
                    T2.graph['mu2'] = mu2
                    if not isomorphic(stabilization(T2, num_fixed)[0], G, num_fixed, all_markings=False):
                        continue
                    elif any([isomorphic(T2, Q, num_fixed, all_markings=not ignore_labels) for Q in seen_graphs]):
                        continue
                    TOTAL += 1
                    seen_graphs.append(T2)
                    if display:
                        display_graph(T2, mu1, mu2, genus)

                    yield T2

if __name__ == '__main__':
    G = nx.Graph()
    G.add_edge('l0','r0')
    G.nodes['l0']['ramif']=[[(1,2)],[],[],[],[]]
    G.nodes['l0']['genus']=1
    G.nodes['r0']['ramif']=[[(2,2),(3,1)],[],[],[],[]]
    G.nodes['r0']['genus']=0
    possible_graphs_old([[2,2,1],[3,1,1],[3,1,1],[3,1,1]], G, ignore_labels=True)
    '''
    G = nx.MultiGraph()
    G.add_edge('l0','l0')
    G.nodes['l0']['ramif']=[[(1,4),(2,1)],[],[],[],[]]
    G.nodes['l0']['genus']=0
    possible_graphs_old([[4,1],[4,1],[3,1,1],[2,1,1,1],[2,1,1,1]], G, double=True, genus=0) #TODO: why no (2,2) edge?
    '''
