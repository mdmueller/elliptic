# calculate genus 0 invariant W_{mu_1,...,mu_r}
from hurwitz import marked_hurwitz, auts, fact
from graphs import possible_graphs, stabilization, display_graph
from pynauty import Graph, autgrp
import copy
import networkx as nx
import functools
import cProfile

KNOWN = {} # known values of W1

def count_auts(G):
    #TODO: adapt for graphs with double edges
    H = Graph(len(G.nodes))
    coloring = []
    uncolored = []
    degrees = set()
    for i, node in enumerate(G.nodes):
        H.connect_vertex(i, [j for j, node2 in enumerate(G.nodes) if node2 in G.neighbors(node)])
        if i==0 or any([any([x[0]<100 for x in ramif]) for ramif in G.nodes[node]['ramif']]):
            # node has some marked point or is genus 1
            coloring.append(set([i]))
        else:
            uncolored.append((i, node))
            degrees.add(G.nodes[node]['degree'])
    # now color all genus 0 unmarked nodes by degree
    for d in degrees:
        coloring.append({i for i, node in uncolored if G.nodes[node]['degree'] == d})
    #TODO: color edges by degree too...
    H.set_vertex_coloring(coloring)

    _, num_auts, _, _, _ = autgrp(H)
    return num_auts

def options_helper(unmarked):
    # unmarked is a list of unmarked ramification, one list for each node
    all_unmarked = []
    for L in unmarked:
        all_unmarked.extend(L)
    options = 1
    for val in set(all_unmarked):
        options *= fact(all_unmarked.count(val))
        for L in unmarked:
            options /= fact(L.count(val))
    return options
    
def count_options(ramifs):
    options = 1
    for ramif in ramifs:
        options *= options_helper([[t[1] for t in x if t[0]==100] for x in ramif])
    return options

class WException(Exception):
    pass

def W(mu):
    if any([any([t<=0 for t in x]) for x in mu]):
        raise WException('Only positive ramification allowed')
    d = sum(mu[0])
    if not all([d == sum(x) for x in mu]):
        raise WException('Partitions should have the same length: '+str(mu))
    r = len(mu)
    if r < 2:
        raise WException('r must be at least 2')
    assert(r <= 4) # TODO: handle r>4 by making sure fixed points are handled properly...
    sum_of_lengths = sum([len(x) for x in mu])
    if sum_of_lengths != d*(r-2)+2:
        raise WException('Ramification fails to match Riemann-Hurwitz: '+str(mu))
    if d == 1:
        return 1
    #if len(mu[0])>r or len(mu[0])+len(mu[1])<r: #TODO: why did I think this was necessary...
    #    raise WException('r should be between |mu_1| and |mu_1|+|mu_2|: '+str(mu))
    if r==3:
        return marked_hurwitz(mu, d)

    # for our invariant, which fixes an element of M_{0,r}bar, consider case where last two points coincide
    # this is a nodal genus 0 curve with a component containing first r-2 points and a component containing last 2 points
    G = nx.MultiGraph()
    G.add_edge('l0','r0')
    G.nodes['l0']['ramif']=[[(i+1, x) for i,x in enumerate(mu[0]) if i+1<=r-2], [(i+1+len(mu[0]), x) for i,x in enumerate(mu[1]) if i+1+len(mu[0])<=r-2],[],[],[]]
    G.nodes['l0']['genus']=0
    G.nodes['r0']['ramif']=[[(i+1, x) for i,x in enumerate(mu[0]) if i+1>r-2 and i+1<=r], [(i+1+len(mu[0]), x) for i,x in enumerate(mu[1]) if i+1+len(mu[0])>r-2 and i+1+len(mu[0])<=r],[],[],[]]
    G.nodes['r0']['genus']=0

    TOTAL = 0
    for graph in possible_graphs(mu, G, num_fixed=r, genus=0):
        factor = 1
        # take product of W's for each node (automorphisms should already be dealt with...)
        # TODO: incorporate the node itself as a fixed point
        for node in graph.nodes():
            ramifs = [[t[1] for t in x] for x in graph.nodes[node]['ramif'] if x] # remove empty ramification profiles
            ramifs.append([sum([z['weight'] for z in graph[node][y].values()]) for y in graph[node]]) # add ramif for weights over neighbors y
            try:
                factor *= W(ramifs) #TODO: ensure fixed points are dealt with properly
            except WException:
                pass

        TOTAL += factor
    return TOTAL

def N(mu):
    pass #TODO

def standardize(mus):
    # standardize the tuples mus
    return (mus[0],)+tuple(sorted(mus[1:]))

@functools.cache
def load_known():
    global KNOWN
    with open('known.txt') as f:
        s = f.read()
    lines = s.split('\n')
    
    for line in lines:
        if not line:
            continue
        L, R = tuple(line.split(' -> '))
        # L = "4 3,1 3,1 2,1,1" for instance
        mus = [tuple([int(x) for x in mu.split(',')]) for mu in L.split(' ')]
        mus = standardize(tuple(mus))
        KNOWN[mus] = int(R)


def W1(mu, num_fixed=0, log=False, display=False): # genus 1 invariant
    global KNOWN

    mu = [x for x in mu if x] #TODO: why is this necessary?

    printL = (lambda x: print(x)) if log else (lambda x: None)

    if mu==[[4],[2,1,1],[3,1],[2,1,1],[2,1,1]] and num_fixed==1: #TODO: hacky
        return 10*2*2*2*2

    if any([any([t<=0 for t in x]) for x in mu]):
        raise WException('Only positive ramification allowed')
    d = sum(mu[0])
    if not all([d == sum(x) for x in mu]):
        raise WException('Partitions should have the same length: '+str(mu))
    r = len(mu)
    if r < 2:
        raise WException('r must be at least 2')
    sum_of_lengths = sum([len(x) for x in mu])
    if sum_of_lengths != d*(r-2):
        raise WException('Ramification fails to match Riemann-Hurwitz: '+str(mu))

    for i in range(1, len(mu)):
        if i==1 and num_fixed>0:
            continue
        if all([x==1 for x in mu[i]]):
            return fact(len(mu[i]))*W1(mu[:i]+mu[i+1:], num_fixed, log=log, display=display)
    
    # num_fixed is the number of points in mu[1] which are fixed
    if num_fixed>0:
        if num_fixed==len(mu[1]):
            return 0
        elif mu[1][0]==1 and sorted(mu[1])==[1]*(len(mu[1])-1)+[2]:
            return 0
        elif num_fixed==1 and mu[1]==[2,1] and mu[0] in [[3],[2,1]] and all([x==[2,1] for x in mu[2:]]):
            return fact(len(mu[2:]))
        elif num_fixed==1 and mu[1]==[2,1] and mu[0] in [[3],[2,1]] and len([y for y in mu[2:] if y==[3]])>0:
            return 0
        elif num_fixed==1 and all([y==1 for y in mu[1]]):
            return W1([mu[0]]+mu[2:], log=log, display=display)*fact(len(mu[1])-1)
        elif mu[1]!=[1]*len(mu[1]) and sum([d-1-len(x) for x in mu[1:] if x!=[1]*d]) >= d-2: # invariant is finite without fixing another ramif point
            return 0
        elif num_fixed==len(mu[1])-1 and mu[1][-1]==1:
            # last point is determined by the previous ones
            m = len([y for y in mu[2:] if sum(y)-len(y)==1])
            if m<len(mu[2:]): # can't have extra ramification
                return 0
            return fact(d-2)**m * fact(m) # (2,1^{d-2})^m
        printL((mu,num_fixed))
        #return 1 #TODO
        #import pdb;pdb.set_trace()
        #raise Exception('too many fixed', (mu,num_fixed)) #TODO

    mu = [sorted(x)[::-1] for x in mu if x]
    tuple_version = standardize(tuple([tuple(x) for x in mu]))
    load_known()
    if tuple_version in KNOWN:
        return KNOWN[tuple_version]
    
    #if sum([d-1-len(x) for x in mu[1:] if x!=[2]+[1]*(d-2)]) != d-2: #TODO
    #    raise WException('Dimension is nonzero: '+str(mu))
    if len(mu[0])<=1:
        return 0
        #raise Exception('Base case not known: '+str(mu)) #todo: wexception?

    printL('computing W1('+str(mu)+', num_fixed='+str(num_fixed)+')')

    G = nx.MultiGraph()
    G.add_edge('l0','r0')
    G.nodes['l0']['ramif']=[[(i+1, x) for i,x in enumerate(mu[0]) if i<len(mu[0])-2], [(len(mu[0])+1+j, x) for j,x in enumerate(mu[1]) if j<num_fixed],[],[],[],[],[],[],[],[],[]]
    G.nodes['l0']['genus']=1
    G.nodes['r0']['ramif']=[[(i+1, x) for i,x in enumerate(mu[0]) if i>=len(mu[0])-2], [],[],[],[],[],[],[],[],[],[]]
    G.nodes['r0']['genus']=0

    TOTAL = 0
    for graph in possible_graphs(mu, G, num_fixed=num_fixed+len(mu[0]), genus=1, display=False):
        factor = 1
        for node in graph.nodes():
            ramifs = [[t[1] for t in x] for x in graph.nodes[node]['ramif']]
            node_unfixed = []
            node_fixed = []
            for neighbor in graph[node]:
                weight = sum([z['weight'] for z in graph[node][neighbor].values()])
                # check if neighbor connects to a marked point (and not via the current node)
                G2 = copy.deepcopy(graph)
                G2.remove_node(node)
                for n2 in nx.node_connected_component(G2, neighbor):
                    if any([any([i<100 for i,x in ram]) for ram in G2.nodes[n2]['ramif']]): # includes a marked point
                        node_fixed.append(weight)
                        break
                else:
                    node_unfixed.append(weight)

            try:
                if graph.nodes[node]['genus']==0: #TODO: do we have to ensure there are only 3 ramification profiles?
                    printL('multiplying by H('+str([x for x in ramifs+[node_fixed+node_unfixed] if x])+')')
                    X = marked_hurwitz([x for x in ramifs+[node_fixed+node_unfixed] if x], sum(node_fixed+node_unfixed))
                    printL('Hurwitz: {}'.format(X))
                    factor *= X
                elif ramifs[0]: # fixed pts in two different profiles
                    second_fixed = len([i for i,x in graph.nodes[node]['ramif'][1] if i<100])
                    if second_fixed!=0:
                        if second_fixed==len(ramifs[1]): # all points in first two fibers fixed
                            X = 0
                        elif second_fixed==len(ramifs[1])-1 and ramifs[1][-1]==1: # invariant would be nonzero without a third fiber
                            if all([x==1 for x in node_unfixed+node_fixed]):
                                printL('trying '+str(ramifs)+' -- '+str(second_fixed))
                                X = W1(ramifs, num_fixed=second_fixed)*fact(len(node_unfixed))
                            else:
                                X = 0
                        else: # TODO: how to deal with fixed points in three fibers?
                            display_graph(graph, graph.graph['mu1'], graph.graph['mu2'], 1)
                            import pdb;pdb.set_trace()
                    else:
                        ramifs.insert(1, node_fixed+node_unfixed)
                        ramifs = [x for x in ramifs if x]
                        printL('trying '+str(ramifs)+' --- '+str(len(node_fixed)))
                        X = W1(ramifs, num_fixed=len(node_fixed))
                    factor *= X
                    printL('got '+str(X))
                elif ramifs[1] and num_fixed>0: # again, fixed pts in two different profiles
                    ramifs = ramifs[1:]
                    if all([i<100 for i,x in graph.nodes[node]['ramif'][1]]):
                        # ramifs[0] can go first, as all points are marked
                        ramifs.insert(1, node_fixed+node_unfixed)
                        nfixed = len(node_fixed)
                    elif node_unfixed:
                        # TODO: what happens if there are select points marked in two different fibers?
                        raise Exception('what happens if there are select points marked in two different fibers?')
                    else:
                        # ramif from edges can go first, as all points are marked
                        ramifs.insert(0, node_fixed)
                        nfixed = len([i for i,x in graph.nodes[node]['ramif'][1] if i<100])
                    ramifs = [x for x in ramifs if x]
                    printL('trying '+str(ramifs)+' ---- '+str(nfixed))
                    X = W1(ramifs, num_fixed=nfixed)
                    factor *= X
                    printL('got '+str(X))
                else: # no marked points on the genus 1
                    ramifs[0] = node_fixed+node_unfixed
                    ramifs = [x for x in ramifs if x]
                    printL('trying '+str(ramifs))
                    X = W1(ramifs)
                    factor *= X
                    printL('got '+str(X))
            except WException as e:
                printL(e)
                factor = 0
            if factor == 0:
                break

        T, mult = stabilization(graph, len(mu[0]))
        printL('mult: {}'.format(mult))
        factor *= mult
        possibilities = count_options([[graph.nodes[node]['ramif'][i] for node in graph.nodes] for i in range(len(mu))])
        printL('possibilities: {}'.format(possibilities))
        factor *= possibilities
        auts = count_auts(graph)
        printL('automorphisms: {}'.format(auts))
        factor /= auts
        printL('ADDING FACTOR: {}'.format(factor))

        if factor!=0 and display:
            display_graph(graph, graph.graph['mu1'], graph.graph['mu2'], 1)
        TOTAL += factor

    KNOWN[tuple_version] = TOTAL # cache result
    return TOTAL

def N1(mu, num_fixed=0, log=False, display=False):
    printL = lambda x: print(x) if log else None

    # num_fixed is the number of points in mu[1] which are fixed

    assert(len(mu)>0)
    d = sum(mu[0])

    while sum([len(x) for x in mu]) > d*(len(mu)-2):
        mu.append([2]+[1]*(d-2))

    factor = 1
    factor *= auts(mu[1][num_fixed:])
    for x in mu[2:]:
        factor *= auts(x)
    factor *= auts([tuple(x) for x in mu[(1 if num_fixed==0 else 2):]])
    X = W1(mu, num_fixed=num_fixed, log=log, display=display)
    printL('W: {}'.format(X))
    printL('FACTOR: {}'.format(factor))
    return X/factor

def Q():
    print(N1([[1,1,1,1,1],[3,1,1],[3,1,1],[3,1,1]], log=False, display=False))

def S(mu):
    d = sum(mu)
    print([mu, [2]+[1]*(d-2), *[[3]+[1]*(d-3)]*(d-3)])
    return N1([mu, [2]+[1]*(d-2), *[[3]+[1]*(d-3)]*(d-3)],num_fixed=1, log=False, display=False)

if __name__ == '__main__':
    #print(W([[4,1],[1,1,1,1,1],[4,1],[3,1,1]]))
    #cProfile.run('Q()', sort='cumtime')

    #print(S([2,2,1]))
    print(N1([[1,1,1,1],[2,1,1],[3,1]],num_fixed=1,log=True,display=True))
