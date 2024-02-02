# calculate genus 0 invariant W_{mu_1,...,mu_r}
from hurwitz import marked_hurwitz, auts
from graphs import possible_graphs, stabilization
import copy
import networkx as nx
import functools

KNOWN = {} # known values of W1

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
        KNOWN[tuple(mus)] = int(R)


def W1(mu, num_fixed=0, log=False, display=False): # genus 1 invariant
    printL = lambda x: print(x) if log else None
    # num_fixed is the number of points in mu[1] which are fixed
    if num_fixed>0:
        if num_fixed==len(mu[1]):
            return 0
        elif num_fixed==1 and all([y==1 for y in mu[1]]):
            return W1([mu[0]]+mu[2:], log=log)
        printL((mu,num_fixed))
        raise Exception('too many fixed') #TODO
    mu = [sorted(x)[::-1] for x in mu if x]
    tuple_version = tuple([tuple(x) for x in mu])
    load_known()
    if tuple_version in KNOWN:
        return KNOWN[tuple_version]
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
    if sum([d-1-len(x) for x in mu[1:] if x!=[2]+[1]*(d-2)]) != d-2:
        raise WException('Dimension is nonzero: '+str(mu))
    if len(mu[0])<=1:
        raise WException('Base case not known: '+str(mu))
    printL('computing W1('+str(mu)+', num_fixed='+str(num_fixed)+')')

    G = nx.MultiGraph()
    G.add_edge('l0','r0')
    G.nodes['l0']['ramif']=[[(i+1, x) for i,x in enumerate(mu[0]) if i<len(mu[0])-2], [],[],[],[],[],[],[],[],[],[]]
    G.nodes['l0']['genus']=1
    G.nodes['r0']['ramif']=[[(i+1, x) for i,x in enumerate(mu[0]) if i>=len(mu[0])-2], [],[],[],[],[],[],[],[],[],[]]
    G.nodes['r0']['genus']=0

    TOTAL = 0
    for graph in possible_graphs(mu, G, genus=1, display=display):
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
                    if G2.nodes[n2]['ramif'][0]: # includes a marked point
                        node_fixed.append(weight)
                        break
                else:
                    node_unfixed.append(weight)

            try:
                if graph.nodes[node]['genus']==0:
                    printL('multiplying by H('+str([x for x in ramifs+[node_fixed+node_unfixed] if x])+')')
                    X = marked_hurwitz([x for x in ramifs+[node_fixed+node_unfixed] if x], sum(node_fixed+node_unfixed))
                    printL('Hurwitz: {}'.format(X))
                    factor *= X
                elif ramifs[0]: # fixed pts in two different profiles
                    ramifs.insert(1, node_fixed+node_unfixed)
                    ramifs = [x for x in ramifs if x]
                    printL('trying '+str(ramifs)+' ---- '+str(len(node_fixed)))
                    X = W1(ramifs, num_fixed=len(node_fixed))
                    factor *= X
                    printL('got '+str(X))
                else:
                    ramifs[0] = node_fixed+node_unfixed
                    ramifs = [x for x in ramifs if x]
                    printL('trying '+str(ramifs))
                    X = W1(ramifs)
                    factor *= X
                    printL('got '+str(X))
            except WException as e:
                printL(e)
                factor = 0
                pass

        T, mult = stabilization(graph, len(mu[0]))
        printL('mult: {}'.format(mult))
        factor *= mult
        printL('ADDING FACTOR: {}'.format(factor))
        TOTAL += factor

    KNOWN[tuple_version] = TOTAL # cache result
    return TOTAL

def N1(mu, log=False, display=False):
    printL = lambda x: print(x) if log else None

    assert(len(mu)>0)
    d = sum(mu[0])

    while sum([len(x) for x in mu]) > d*(len(mu)-2):
        mu.append([2]+[1]*(d-2))
    
    factor = 1
    for x in mu[1:]:
        factor *= auts(x)
    factor *= auts([tuple(x) for x in mu[1:]])
    X = W1(mu, log=log, display=display)
    printL('W: {}'.format(X))
    printL('FACTOR: {}'.format(factor))
    return X/factor

if __name__ == '__main__':
    #print(W([[4,1],[1,1,1,1,1],[4,1],[3,1,1]]))
    print(N1([[4,1],[4,1],[3,1,1]], log=False, display=False))
