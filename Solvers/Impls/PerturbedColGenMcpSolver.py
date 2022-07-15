# Solver for the Multicast Packing Problem that uses Simplex, but uses perturbed
# duals for column generation, rather than those found in the tableau

import GlobalConstants
from Solvers.Impls.PureColGenMcpSolver import PureColGenMcpSolver, cost
from Solvers.Impls.JansenZhangMinMaxer import JansenZhangMinMaxer

class PerturbedColGenMcpSolver(PureColGenMcpSolver, JansenZhangMinMaxer):
    def __init__(self, instance=None, block_approx=2):
        super().__init__(instance, block_approx)
        self.t = GlobalConstants.TOLERANCE
        self.true_p = dict()
        
    def perform_checks_and_updates(self, x):
        t = self.t
        if sum([cost(self.new_trees[i], self.true_p[(x,t)]) 
                for i in range(self.instance.num_requests)]) >= self.lamb(x):
            self.stop_flag |= GlobalConstants.STOP_DUALITYMATCH
            
        if all([cost(self.new_trees[i], self.true_p[(x,t)]) - self.q(x)[i] >= 0
                for i in range(self.instance.num_requests)]):
            self.stop_flag |= GlobalConstants.STOP_FLAG_REDCOST
            
        super().perform_checks_and_updates(x)
        
    def generate_p(self, x, t=None):
        if not t:
            t = self.t
        super().generate_p(x)
        self.true_p[(x,t)] = self.p(x)
        JansenZhangMinMaxer.generate_p(self, x, t)
        
    # Method for printing debug info
    def print_info(self, x, t):
        print(self.iteration)
        print("{:05b}".format(self.stop_flag))
        print("lambda(x): {}".format(self.lamb(x)))
        print("phi_t(x): {}".format(self.phi(x,t)))
        print("tolerance: {}".format(self.toleranceFunction()))
        for e in self.instance.graph.edges():
            e = tuple(sorted(e))
            if self.p(x)[e] > 0.001:
                print("p_{} = {}".format(e, self.p(x)[e]))
        for i in range(self.instance.num_requests):
            print("ReducedCost_{} = {}".format(i,
            cost(self.new_trees[i], self.true_p[(x,t)]) - self.q(x)[i]))
        print("")