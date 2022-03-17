# Solver for the Multicast Packing Problem that uses Simplex

import gurobipy as gp
import networkx as nx

from FrozenDict import FrozenDict
from MulticastPackingSolver import MulticastPackingSolver

class PureColGenMcpSolver(MulticastPackingSolver):
    def get_next_solution(self):
        self.reduced_LP.optimize()
        x = [dict() for i in range(self.instance.num_requests)]
        for i in range(self.instance.num_requests):
            for T in self.column_generator.Gurobi_variables[i]:
                value = self.column_generator.Gurobi_variables[i][T].getAttr(gp.GRB.Attr.X)
                if 1 >= value > 0: # Should have better check that var isn't lambda
                    x[i][T] = value
            x[i] = FrozenDict(x[i])
        
        x = tuple(x)
        self.generate_lamb(x)
        self.generate_fVal(x)
        self.generate_p(x)
        self.generate_q(x)
        
        return x
    
    def perform_checks_and_updates(self, x):
        if ( 
            (sum([cost(self.new_trees[i], self.p(x)) 
                  for i in range(self.instance.num_requests)]) >= self.lamb(x))
            or (all([cost(self.new_trees[i], self.p(x)) - self.q(x)[i] >= 0 
                     for i in range(self.instance.num_requests)]))
            or (self.toleranceFunction() <= self.tol/(2+self.tol)) ):
                self.stop_flag = True
                # Could have different stop flags for different stopping criteria
    
    def generate_lamb(self, x):
        self.objVal[x] = self.reduced_LP.getObjective().getValue()
        # Doesn't work if x is not current LP solution
    
    def generate_fVal(self, x):
        self.fVal[x] = dict()
        for e in self.instance.graph.edges():
            self.fVal[x][tuple(sorted(e))] = (
                self.lamb(x) + self.reduced_LP.getConstrByName(
                    "{} congestion".format(sorted(e))).getAttr(gp.GRB.Attr.Slack) 
                # Note: Gurobi signs their slacks stupidly
            )
        
    def generate_p(self, x, t=0):
        self.price[(x,t)] = dict()
        for e in self.instance.graph.edges():
            self.price[(x,t)][tuple(sorted(e))] = self.reduced_LP.getConstrByName(
                "{} congestion".format(sorted(e))).getAttr(gp.GRB.Attr.Pi)
            
    def generate_q(self, x, t=0):
        self.multicast_costs[(x,t)] = [None] * len(self.instance.requests)
        for i in range(self.instance.num_requests):
            self.multicast_costs[(x,t)][i] = self.reduced_LP.getConstrByName(
                "Tree Selection for {}".format(i)).getAttr(gp.GRB.Attr.Pi)
        
# Helper Functions
def cost(G, prices):
    original_weights = nx.get_edge_attributes(G, "weight") 
    nx.set_edge_attributes(G, prices, "weight")
    retval = sum(nx.get_edge_attributes(G, "weight").values())
    nx.set_edge_attributes(G, original_weights, "weight")
    return retval