# Solver for the Integer Program of Multicast Packing Problem

import gurobipy as gp
import networkx as nx

from Solvers.Impls.PureColGenMcpSolver import PureColGenMcpSolver, cost
#from PureColGenMcpSolver import cost

class ColGenIPSolver(PureColGenMcpSolver):
    def __init__(self, instance=None, block_approx=2):
        super().__init__(instance, block_approx)
        
        self.column_generator.varType = gp.GRB.INTEGER
        
        for var in self.reduced_LP.getVars():
            var.setAttr(gp.GRB.Attr.VType, gp.GRB.INTEGER)
        
        self.reduced_LP.update()
        
#         self.column_generator.varType = gp.GRB.INTEGER
#         for var in self.reduced_LP.getVars():
#             print(var)
        
#         for i in range(self.instance.num_requests):
#             for T in self.column_generator.Gurobi_variables[i]:
#                 print(self.column_generator.Gurobi_variables[i][T])
        
            
    def generate_p(self, x, t=0):
        newPriceDict = dict()
        for edge in self.instance.graph.edges():
            e = tuple(sorted(edge))
            if self.f(x)[e] == self.lamb(x):
                newPriceDict[e] = 1
            else:
                newPriceDict[e] = 0
        
        newPriceDict = {e: newPriceDict[e] / sum(list(newPriceDict.values())) for e in newPriceDict}
        self.price[(x,t)] = newPriceDict
                
    def generate_q(self, x, t=0):
        newCostDict = [self.instance.num_edges] * len(self.instance.requests)
        prices = self.p(x)
        for i in range(self.instance.num_requests):
            for T in self.column_generator.Gurobi_variables[i]:
                Tcost = cost(T, prices)
                if Tcost < newCostDict[i]:
                    newCostDict[i] = Tcost
            
        self.multicast_costs[(x,t)] = newCostDict