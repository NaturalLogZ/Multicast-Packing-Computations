# Abstract class for solving Multicast Packing Problems

from abc import ABC, abstractmethod
from math import log

import gurobipy as gp
import networkx as nx
import mpmath as mp
import matplotlib.pyplot as plt

import GlobalConstants
import DebugConstants as db
from MulticastPackingInstance import MulticastPackingInstance
from ColumnGenerators.Impls.Approx2MulticastPackingColumnGenerator import Approx2MulticastPackingColumnGenerator
from ColumnGenerators.Impls.ExactMulticastPackingColumnGeneratorIP import ExactMulticastPackingColumnGeneratorIP
from ColumnGenerators.Impls.ExactMcpWithDelayColumnGenerator import ExactMcpWithDelayColumnGenerator

class MulticastPackingSolver(ABC):
    
    # Constructor
    def __init__(self, instance=None, block_approx=2):
        # Attributes related to problem instance
        if instance is None:
            instance = MulticastPackingInstance(NUM_MULTICAST_REQUESTS, MAX_MULTICAST_SIZE)
        self.instance = instance
        self.reduced_LP = create_LP(self.instance.graph, self.instance.requests)
        
        # Attributes related to how this solver generates columns and solutions
        if block_approx == "Delay":
            block_solver = ExactMcpWithDelayColumnGenerator(self.instance, self.reduced_LP)
        elif block_approx == 1:
            block_solver = ExactMulticastPackingColumnGeneratorIP(self.instance, self.reduced_LP)
        elif block_approx >= 2:
            block_solver = Approx2MulticastPackingColumnGenerator(self.instance, self.reduced_LP)
        self.column_generator = block_solver
        
        # Attributes that track things
        self.tol = GlobalConstants.TOLERANCE
        self.t = None
        self.stop_flag = 0b0
        self.iteration = -1
        self.solution = list()
        self.objVal = dict()
        self.fVal = dict()
        self.price = dict()
        self.multicast_costs = dict()
        self.new_trees = self.column_generator.generate_new_trees()
        
        if db.DEBUG_LEVEL >= db.DEBUG_LEVEL_EXTREME:
            print(reduced_LP.getA().toarray())
            print(reduced_LP.getAttr("RHS"))
            
        if db.DEBUG_LEVEL >= db.DEBUG_LEVEL_FULL:
            G = self.instance.graph
            pos = self.instance.pos
            for i in range(self.instance.num_requests):
                tree = self.new_trees[i]
                plt.figure()
                nx.draw_networkx(G, pos)
                nx.draw_networkx_edges(G, pos, edgelist=tree.edges(), width=5, edge_color="y", style="dashed")
                #print(tree.edges())
            plt.show()
    
    # Abstract methods for generating values of functions needed by the solver
    
    @abstractmethod
    def get_next_solution(self):
        pass
    
    @abstractmethod
    def perform_checks_and_updates(self):
        pass
    
    @abstractmethod
    def generate_lamb(self, x):
        pass
    
    @abstractmethod
    def generate_fVal(self, x):
        pass
     
    @abstractmethod   
    def generate_p(self, x, t):
        pass
    
    @abstractmethod
    def generate_q(self, x, t):
        pass
    
    # Methods for accessing values of functions needed by the solver
    
    def lamb(self, x):
        if x not in self.objVal:
            self.generate_lamb(x)
        return self.objVal[x]
    
    def f(self, x):
        if x not in self.fVal:
            self.generate_fVal(x)
        return self.fVal[x]
    
    def p(self, x, t=None):
        if not t:
            t = self.t
        if (x,t) not in self.price:
            self.generate_p(x,t)
        return self.price[(x,t)]
    
    def q(self, x, t=None):
        if not t:
            t = self.t
        if (x,t) not in self.multicast_costs:
            self.generate_q(x,t)
        return self.multicast_costs[(x,t)]
    
    def theta(self, x, t=None):
        return self.lamb(x)
    
    def phi(self, x, t=None):
        return log(self.lamb(x))
    
    def Phi(self, theta, x, t=None):
        return 
    
    def toleranceFunction(self):
        x = self.solution[self.iteration]
        f_prime = dict()
        for e in self.instance.graph.edges():
            f_prime[tuple(sorted(e))] = mp.mpf(0)
            for i in range(self.instance.num_requests):
                if self.new_trees[i].has_edge(*e):
                    f_prime[tuple(sorted(e))] += 1
    
        pf = dict_innerProd(self.p(x, self.t), self.f(x))
        pf_prime = dict_innerProd(self.p(x, self.t), f_prime)
        return (pf-pf_prime)/(pf+pf_prime)
    
    # Method for printing debug info
    def print_info(self, x, t):
        print(self.iteration)
        print("{:05b}".format(self.stop_flag))
        print("lambda(x): {}".format(self.lamb(x)))
        print("phi_t(x): {}".format(self.phi(x,t)))
        print("tolerance: {}".format(self.toleranceFunction()))
        if db.DEBUG_LEVEL >= db.DEBUG_LEVEL_THEORY_2:
            for e in self.instance.graph.edges():
                e = tuple(sorted(e))
                if self.p(x)[e] > 0.001:
                    print("p_{} = {}".format(e, self.p(x)[e]))
            for i in range(self.instance.num_requests):
                print("ReducedCost_{} = {}".format(i,
                cost(self.new_trees[i], self.p(x)) - self.q(x)[i]))
                print("NewCost_{} = {}".format(i,
                cost(self.new_trees[i], self.p(x))))
        if db.DEBUG_LEVEL >= db.DEBUG_LEVEL_FULL:
            G = self.instance.graph
            pos = self.instance.pos
            for i in range(self.instance.num_requests):
                tree = self.new_trees[i]
                plt.figure()
                plt.title("Tree {} for iteration {}".format(i, self.iteration))
                nx.draw_networkx(G, pos)
                nx.draw_networkx_edges(G, pos, edgelist=tree.edges(), width=5, edge_color="y", style="dashed")
                nx.draw_networkx_edge_labels(tree, pos, self.p(x))
                #print(tree.edges())
            plt.show()
        print("")
    
    # Main Function
    def perform_iteration(self):
        t = self.t
        x = self.get_next_solution()
        self.iteration += 1
        if db.DEBUG_LEVEL > db.DEBUG_LEVEL_FULL:
            print("Iteration {}: num trees = {}, num Gurobi Vars = {}".format(
                self.iteration,
                sum(len(x[i]) for i in range(self.instance.num_requests)),
                sum(len(self.column_generator.Gurobi_variables[i]) for i in range(self.instance.num_requests)),
            ))
        self.solution.append(x)
        self.new_trees = self.column_generator.generate_new_trees(self.p(x, t))
        self.perform_checks_and_updates(x)
        
        if (
            ((db.DEBUG_LEVEL > db.DEBUG_LEVEL_THEORY_0) and
             (self.iteration % round(1/(db.DEBUG_LEVEL-1)) == 0)) or 
            ((db.DEBUG_LEVEL >= db.DEBUG_LEVEL_THEORY_0) and (self.stop_flag))
        ):
            self.print_info(x, t)
        
        
# Helper Functions
def dict_innerProd(x, y):
    assert x.keys() == y.keys()
    return sum(x[key]*y[key] for key in x.keys())

def create_LP(G, multicast_requests):
    nx.set_edge_attributes(G, 1, "weight")
    reduced_LP = gp.Model("Multicast Packing Model - Reduced")
    congestion = reduced_LP.addVar(name="lambda")
    reduced_LP.setObjective(congestion, gp.GRB.MINIMIZE)
    reduced_LP.update()
    # Packing Constraints
    for e in G.edges():
        reduced_LP.addConstr(congestion >= 0, name="{} congestion".format(sorted(e)))
    # Simplex Constraints
    for i in range(len(multicast_requests)):
        reduced_LP.addConstr(0*congestion == 1, name="Tree Selection for {}".format(i))
    reduced_LP.update()
    return reduced_LP

# Helper Functions
def cost(G, prices):
    original_weights = nx.get_edge_attributes(G, "weight") 
    nx.set_edge_attributes(G, prices, "weight")
    retval = sum(nx.get_edge_attributes(G, "weight").values())
    nx.set_edge_attributes(G, original_weights, "weight")
    return retval