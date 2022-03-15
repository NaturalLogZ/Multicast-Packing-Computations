# 2-Approximation for Column Generation

import networkx as nx
from networkx.algorithms.approximation.steinertree import steiner_tree as steiner_tree_2
from MulticastPackingColumnGenerator import MulticastPackingColumnGenerator

class Approx2MulticastPackingColumnGenerator(MulticastPackingColumnGenerator):
    def generate_tree(self, i, prices):
        original_weights = nx.get_edge_attributes(self.instance.graph, "weight") 
        nx.set_edge_attributes(self.instance.graph, prices, "weight")
        retval = nx.algorithms.approximation.steinertree.steiner_tree(
            self.instance.graph, self.instance.requests[i].multicast_group())
        nx.set_edge_attributes(self.instance.graph, original_weights, "weight")
        return retval