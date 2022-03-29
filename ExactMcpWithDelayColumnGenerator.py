# MCF IP with Delay Constraints for Column Generation

import gurobipy as gp
import networkx as nx
from ExactMulticastPackingColumnGeneratorIP import ExactMulticastPackingColumnGeneratorIP

class ExactMcpWithDelayColumnGenerator(ExactMulticastPackingColumnGeneratorIP):
    def __init__(self, instance, reduced_LP):
        super().__init__(instance, reduced_LP)
        D = nx.DiGraph(instance.graph)
        max_delay = instance.delay
        for i in range(instance.num_requests):
            model = self.GurobiModels[i]
            f = self.variables[i][1]
            
            # Add constraints
            for r in instance.requests[i].recipients:
                model.addConstr(
                    sum(D[u][v]['delay']*f[r][(u,v)] for u,v in D.edges()) <= max_delay,
                    name="Delay constraint for {}-flow".format(r)
                )
            
            # Note we don't add the objective as that is determined by
            # the prices in each iteration.
            model.update()