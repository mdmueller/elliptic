import unittest
from graphs import *
from genus_zero import *

class TestGraphs(unittest.TestCase):
    def test_partitions(self):
        self.assertEqual(partitions(4,2), [[0,4],[1,3],[2,2],[3,1],[4,0]])
        self.assertEqual(partitions(1,0), [])

    def test_stabilization(self):
        pass #TODO

    def test_bipartite_helper(self):
        pass #TODO

    def test_graphs(self):
        G = nx.MultiGraph()
        G.add_edge('l0','r0')
        G.nodes['l0']['ramif']=[[(1,4),(2,1)],[],[],[],[]]
        G.nodes['l0']['genus']=0
        G.nodes['r0']['ramif']=[[],[(3,3),(4,1)],[],[],[]]
        G.nodes['r0']['genus']=0
        graphs = possible_graphs([[4,1],[3,1,1],[3,1,1],[2,1,1,1]], G, num_fixed=4, genus=0)
        self.assertEqual(len(list(graphs)), 4) # two graph shapes, the second has 3 different placements

    def test_genus_zero(self):
        self.assertEqual(W([[4,1],[3,1,1],[3,1,1],[2,1,1,1]]), 36)
        self.assertEqual(W([[4,1],[2,1,1,1],[4,1],[2,1,1,1]]), 48)

    def test_doubleedge(self):
        G = nx.MultiGraph()
        G.add_edge('l0','l0')
        G.nodes['l0']['ramif']=[[(1,3)],[],[],[],[]]
        G.nodes['l0']['genus']=0
        graphs = possible_graphs([[3],[3],[2,1],[2,1]], G, double=True, genus=0)
        self.assertEqual(len(list(graphs)), 3)

if __name__ == '__main__':
    unittest.main()
