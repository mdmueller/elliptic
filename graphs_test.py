import unittest
from graphs import *
from genus_zero import *

class TestGraphs(unittest.TestCase):
    def test_partitions(self):
        self.assertEqual(partitions(4,2), [[0,4],[1,3],[2,2],[3,1],[4,0]])
        self.assertEqual(partitions(1,0), [])

    def test_stabilization(self):
        G = nx.MultiGraph()
        G.add_edge('l0', 'r0', weight=2)
        G.add_edge('l0', 'r1', weight=2)
        G.add_edge('l1', 'r1', weight=1)
        G.nodes['r0']['ramif']=[[(1,3)]]
        G.nodes['r1']['ramif']=[[(2,1),(3,1)]]
        G.nodes['l0']['ramif']=[[]]
        G.nodes['l1']['ramif']=[[]]
        G.nodes['l0']['genus']=1
        G.nodes['l1']['genus']=0
        G.nodes['r0']['genus']=0
        G.nodes['r1']['genus']=0
        T, mult = stabilization(G, 3)
        self.assertEqual(len(T.nodes()), 2)
        self.assertEqual(mult, 2)

    def test_stabilization_two_edges(self):
        G = nx.MultiGraph()
        G.add_edge('l', 'c', weight=2)
        G.add_edge('c', 'r', weight=3)
        G.nodes['l']['ramif']=[[]]
        G.nodes['c']['ramif']=[[]]
        G.nodes['r']['ramif']=[[(1,1),(2,1)]]
        G.nodes['l']['genus']=1
        for node in ['c','r']:
            G.nodes[node]['genus']=0
        T, mult = stabilization(G, 2)
        self.assertEqual(len(T.nodes()), 2)
        self.assertEqual(mult, 5)

    def test_stabilization_zigzag(self):
        G = nx.MultiGraph()
        G.add_edge('l0', 'r0', weight=1)
        G.add_edge('l0', 'r2', weight=1)
        G.add_edge('l1', 'r0', weight=1)
        G.add_edge('l1', 'r1', weight=1)
        G.nodes['l0']['ramif']=[[(1,2)]]
        G.nodes['l1']['ramif']=[[(2,1),(3,1)]]
        G.nodes['r0']['ramif']=[[]]
        G.nodes['r1']['ramif']=[[]]
        G.nodes['r2']['ramif']=[[]]
        G.nodes['l0']['genus']=1
        G.nodes['l1']['genus']=0
        G.nodes['r0']['genus']=0
        G.nodes['r1']['genus']=0
        G.nodes['r2']['genus']=0
        T, mult = stabilization(G, 3)
        self.assertEqual(len(T.nodes()), 2)
        n1, n2 = tuple(T.nodes())
        self.assertEqual(list(T.neighbors(n1)), [n2])
        self.assertEqual(mult, 2)

    def test_bipartite_helper(self):
        pass #TODO

    def test_graphs(self):
        '''
        G = nx.MultiGraph()
        G.add_edge('l0','r0')
        G.nodes['l0']['ramif']=[[(1,4),(2,1)],[],[],[],[]]
        G.nodes['l0']['genus']=0
        G.nodes['r0']['ramif']=[[],[(3,3),(4,1)],[],[],[]]
        G.nodes['r0']['genus']=0
        graphs = possible_graphs([[4,1],[3,1,1],[3,1,1],[2,1,1,1]], G, num_fixed=4, genus=0)
        self.assertEqual(len(list(graphs)), 4) # two graph shapes, the second has 3 different placements
        '''

    def test_genus_zero(self):
        '''self.assertEqual(W([[4,1],[3,1,1],[3,1,1],[2,1,1,1]]), 36)
        self.assertEqual(W([[4,1],[2,1,1,1],[4,1],[2,1,1,1]]), 48)'''

    def test_doubleedge(self):
        '''
        G = nx.MultiGraph()
        G.add_edge('l0','l0')
        G.nodes['l0']['ramif']=[[(1,3)],[],[],[],[]]
        G.nodes['l0']['genus']=0
        graphs = possible_graphs([[3],[3],[2,1],[2,1]], G, double=True, genus=0)
        self.assertEqual(len(list(graphs)), 3)
        '''

    def test_N1(self):
        self.assertEqual(N1([[4],[4]]), 15)
        self.assertEqual(N1([[4],[2,2],[2,2]]), 3)
        for mu in [[3,1],[2,2],[2,1,1],[1,1,1,1]]:
            self.assertEqual(N1([mu,[4]]), 16)
            self.assertEqual(N1([mu,[2,2],[2,2]]), 3*2**(len(mu)-1), msg='failed when mu={}'.format(str(mu)))

if __name__ == '__main__':
    unittest.main()
