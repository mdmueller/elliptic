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

    def test_genus_zero(self):
        assert(W([[4,1],[3,1,1],[3,1,1],[2,1,1,1]]) == 36)
        assert(W([[4,1],[2,1,1,1],[4,1],[2,1,1,1]]) == 48)

if __name__ == '__main__':
    unittest.main()
