# Abstract class for solving Multicast Packing Problems

from abc import ABC, abstractmethod
from math import log

import gurobipy as gp
import networkx as nx
import mpmath as mp

import GlobalConstants
import DebugConstants as db
from MulticastPackingInstance import MulticastPackingInstance
from Approx2MulticastPackingColumnGenerator import Approx2MulticastPackingColumnGenerator
from ExactMulticastPackingColumnGeneratorIP import ExactMulticastPackingColumnGeneratorIP

class MulticastPackingSolver(ABC):
    
    # Constructor
    def __init__(self, instance=None, block_approx=2):
        # Attributes related to problem instance
        if instance is None:
            instance = MulticastPackingInstance(NUM_MULTICAST_REQUESTS, MAX_MULTICAST_SIZE)
        self.instance = instance
        self.reduced_LP = create_LP(self.instance.graph, self.instance.requests)
        
        # Attributes related to how this solver generates columns and solutions
        if block_approx == 1:
            block_solver = ExactMulticastPackingColumnGeneratorIP(self.instance, self.reduced_LP)
        if block_approx >= 2:
            block_solver = ExactMulticastPackingColumnGeneratorIP(self.instance, self.reduced_LP)
        self.column_generator = block_solver
        
        # Attributes that track things
        self.tol = GlobalConstants.TOLERANCE
        self.t = mp.mpf(0)
        self.stop_flag = False
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
    
    def p(self, x, t=0):
        if (x,t) not in self.price:
            self.generate_p(x,t)
        return self.price[(x,t)]
    
    def q(self, x, t=0):
        if (x,t) not in self.multicast_costs:
            self.generate_q(x,t)
        return self.multicast_costs[(x,t)]
    
    def theta(self, x, t=0):
        return self.lamb(x)
    
    def phi(self, x, t=0):
        return log(self.lamb(x))
    
    def Phi(self, theta, x, t=0):
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
    
    # Main Function
    def perform_iteration(self):
        t = self.t
        x = self.get_next_solution()
        self.iteration += 1
        self.solution.append(x)
        self.new_trees = self.column_generator.generate_new_trees(self.p(x, t))
        self.perform_checks_and_updates(x)
        
        
        if (
            (db.DEBUG_LEVEL > db.DEBUG_LEVEL_THEORY_0) and 
                ( (self.iteration % int(10/(db.DEBUG_LEVEL-1)) == 0) or 
                  (self.stop_flag)
                )
        ):
            print(self.iteration)
            print("lambda(x): {}".format(self.lamb(x)))
            print("phi_t(x): {}".format(self.phi(x,t)))
            print("tolerance: {}".format(self.toleranceFunction()))
        
        
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