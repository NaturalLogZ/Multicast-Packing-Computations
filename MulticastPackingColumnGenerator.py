# Abstract Column Generator

from abc import ABC, abstractmethod
import gurobipy as gp
import networkx as nx

class MulticastPackingColumnGenerator(ABC):
    def __init__(self, instance, reduced_LP):
        self.instance = instance
        self.reduced_LP = reduced_LP
        self.Gurobi_variables = [dict() for i in range(self.instance.num_requests)]
        
    @abstractmethod
    def generate_tree(self, i, prices):
        pass
    
    def generate_new_trees(self, prices=1):
        new_trees = list()
        assert nx.is_weighted(self.instance.graph)
        for i in range(self.instance.num_requests):
            new_tree = self.generate_tree(i, prices)
            new_trees.append(new_tree)
            
            coeffs = list()
            constrs = list()
            # the new tree contributes to the constraint for each of its edges
            for e in new_tree.edges():
                coeffs.append(-1)
                constrs.append(self.reduced_LP.getConstrByName(
                    "{} congestion".format(sorted(e))))
            
            # the new tree contributes to the selection constraint 
            # for the corresponding multicast request
            coeffs.append(1)
            constrs.append(self.reduced_LP.getConstrByName(
                "Tree Selection for {}".format(i)))
            self.Gurobi_variables[i][new_tree] = self.reduced_LP.addVar(
                name="x({}, {})".format(i, nx.info(new_tree)), 
                column=gp.Column(coeffs, constrs))
            self.reduced_LP.update()
            
        return new_trees