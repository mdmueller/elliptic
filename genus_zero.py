# calculate genus 0 invariant W_{mu_1,...,mu_r}
from hurwitz import marked_hurwitz
from graphs import possible_graphs
import networkx as nx

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
    G = nx.Graph()
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
            ramifs.append([graph[node][y]['weight'] for y in graph[node]]) # add ramif for weights over neighbors y
            try:
                factor *= W(ramifs) #TODO: ensure fixed points are dealt with properly
            except WException:
                pass
        print(factor)
        TOTAL += factor
    return TOTAL

def N(mu):
    pass #TODO

if __name__ == '__main__':
    print(W([[2,1],[1,1,1],[3],[2,1]]))
