# 2-Approximation for Column Generation

import networkx as nx
from networkx.algorithms.approximation.steinertree import steiner_tree as steiner_tree_2
from ColumnGenerators.MulticastPackingColumnGenerator import MulticastPackingColumnGenerator

class Approx2MulticastPackingColumnGenerator(MulticastPackingColumnGenerator):
    def generate_tree(self, i, prices):
        G = self.instance.graph
        original_weights = nx.get_edge_attributes(G, "weight") 
        nx.set_edge_attributes(G, prices, "weight")
        assert not nx.is_negatively_weighted(G, weight='weight'), print(prices)
        retval = nx.algorithms.approximation.steinertree.steiner_tree(
            G, self.instance.requests[i].multicast_group())
        nx.set_edge_attributes(G, original_weights, "weight")
        return retval