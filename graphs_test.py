import unittest
from graphs import *

class TestGraphs(unittest.TestCase):
    def test_partitions(self):
        self.assertEqual(partitions(4,2), [[0,4],[1,3],[2,2],[3,1],[4,0]])
        self.assertEqual(partitions(1,0), [])

    def test_stabilization(self):
        T = nx.Graph()
        T.add_edge('l0','r0')
        T.add_edge('l0','r1')
        T.nodes['l0']['ramif'] = [[],[3,1],[2,2]]
        T.nodes['r0']['ramif']=[[2,1],[],[]]
        T.nodes['r1']['ramif']=[[1],[],[]]
        T2 = stabilization(T)
        assert(set(T2.nodes()) == {'l0','r0'})
        assert(T2.nodes['l0']['ramif']==[[1],[3,1],[2,2]])
        assert(T2.nodes['r0']['ramif']==[[2,1],[],[]])

    def test_bipartite_helper(self):
        print('hi')

if __name__ == '__main__':
    unittest.main()
