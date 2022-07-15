# Solver for the Multicast Packing Problem that uses Simplex

from math import log

import mpmath as mp

import GlobalConstants
import DebugConstants as db
from FrozenDict import FrozenDict
from Solvers.MulticastPackingSolver import MulticastPackingSolver, cost, dict_innerProd
#from MulticastPackingSolver import dict_innerProd
   
class JansenZhangMinMaxer(MulticastPackingSolver):
    
    def __init__(self, instance=None, block_approx=2, sigma0=1):
        super().__init__(instance, block_approx)
        self.sigma = mp.mpf(sigma0)
        self.t = self.sigma/6
        self.M = mp.mpf(len(self.instance.graph.edges()))
        self.theta_dict = dict()
        self.phi_dict = dict()
        self.finished_coordination = False
        self.lambda_of_prev_scaling = 0
        self.w = None
        
    def get_next_solution(self):
        if self.iteration == -1:
            x = [dict() for i in range(self.instance.num_requests)]
            step_size = mp.mpf(1)
        else:
            x = self.solution[self.iteration]
            
            f_prime = dict()
            for e in self.instance.graph.edges():
                f_prime[tuple(sorted(e))] = mp.mpf(0)
                for i in range(self.instance.num_requests):
                    if self.new_trees[i].has_edge(*e):
                        f_prime[tuple(sorted(e))] += mp.mpf(1)
    
            pf = dict_innerProd(self.p(x, self.t), self.f(x))
            pf_prime = dict_innerProd(self.p(x, self.t), f_prime)
            step_size = ( (self.t*self.theta(x, self.t)*self.toleranceFunction()) 
                   /(2*self.M*(pf+pf_prime)) )

        
        new_x = [dict() for i in range(self.instance.num_requests)]

        for i in range(self.instance.num_requests):
            for T in x[i]:
                new_x[i][T] = (1-step_size)*x[i][T]
                
            if self.new_trees[i] in x[i]:
                new_x[i][self.new_trees[i]] += step_size
            else:
                new_x[i][self.new_trees[i]] = step_size
            new_x[i] = FrozenDict(new_x[i])
        return tuple(new_x)
    
    def perform_checks_and_updates(self, x):
        t = self.t
        if sum([cost(self.new_trees[i], self.p(x,t)) 
                for i in range(self.instance.num_requests)]) >= self.lamb(x):
            self.stop_flag |= GlobalConstants.STOP_DUALITYMATCH
            
        if all([cost(self.new_trees[i], self.p(x,t)) - self.q(x,t)[i] >= self.tol
                for i in range(self.instance.num_requests)]):
            self.stop_flag |= GlobalConstants.STOP_FLAG_REDCOST
        
        if self.w == None:
            self.sigma = self.sigma/2
            self.t = self.sigma/6
            self.w = mp.mpf((1+self.sigma)/((1+self.sigma/3)*self.M))
            self.lambda_of_prev_scaling = self.lamb(x)
            if db.DEBUG_LEVEL >= db.DEBUG_LEVEL_THEORY_1:
                print("Scaling Phase 0 Over")
            if db.DEBUG_LEVEL >= db.DEBUG_LEVEL_THEORY_2:
                print("sigma = {}".format(self.sigma))
                print("t = {}".format(self.t))
                print("w = {}".format(self.w))
            
            self.new_trees = self.column_generator.generate_new_trees(self.p(x, self.t))
        elif (
                (self.toleranceFunction() <= self.sigma/6) or 
                (self.lamb(x) <= self.w*self.lambda_of_prev_scaling) 
        ):
            if self.sigma <= self.tol:
                self.stop_flag |= GlobalConstants.STOP_FLAG_TOL_MET
            else:
                self.lambda_of_prev_scaling= self.lamb(x)
                self.sigma = self.sigma/2
                self.t = self.sigma/6
                self.w = mp.mpf((1+self.sigma)/(1+2*self.sigma))
                self.new_trees = self.column_generator.generate_new_trees(self.p(x, self.t))
            if db.DEBUG_LEVEL >= db.DEBUG_LEVEL_THEORY_1:
                print("New Scaling Phase")
            if db.DEBUG_LEVEL >= db.DEBUG_LEVEL_THEORY_2:
                print("sigma = {}".format(self.sigma))
                print("t = {}".format(self.t))
                print("w = {}".format(self.w))
            
        if self.iteration >= GlobalConstants.MAX_ITERS:
            self.stop_flag |= GlobalConstants.STOP_FLAG_MAXITER
                
    
    def generate_lamb(self, x):
        self.objVal[x] = max(self.f(x).values())
    
    def generate_fVal(self, x):
        self.fVal[x] = {tuple(sorted(e)): 0 for e in self.instance.graph.edges()}
        for i in range(self.instance.num_requests):
            for T in x[i]:
                for e in T.edges():
                    self.fVal[x][tuple(sorted(e))] += x[i][T]
        
    def generate_p(self, x, t):
        self.price[(x,t)] = dict()
        for e in self.instance.graph.edges():
            self.price[(x,t)][tuple(sorted(e))] = t*self.theta(x,t)/(self.M*(self.theta(x,t) 
                                                                              - self.f(x)[tuple(sorted(e))]))

    def generate_q(self, x, t=None):
        if not t:
            t = self.t
        
        self.multicast_costs[(x,t)] = [
            min(cost(T, self.p(x, t)) for T in x[i].keys()) 
                for i in range(self.instance.num_requests)]
        
    def generate_theta(self, x, t):
        self.theta_dict[(x,t)] = (
            mp.findroot(lambda theta: theta_eq(theta, t, self.M, self.f(x)),
                        (mp.mpf(self.lamb(x)/(1-t/self.M)), mp.mpf(self.lamb(x)/(1-t))), 
                        solver='ridder')
        )
            
    def generate_phi(self, x, t):
        theta = self.theta(x,t)
        phi = mp.mpf(0)
        for e in self.instance.graph.edges():
            phi -= log(theta - self.f(x)[tuple(sorted(e))])
        phi = phi*t/self.M
        phi += log(theta)
        self.phi_dict[(x,t)] = phi
    
    def theta(self, x, t):
        if (x,t) not in self.theta_dict:
            self.generate_theta(x,t)
        return self.theta_dict[(x,t)]
    
    def phi(self, x, t):
        if (x,t) not in self.phi_dict:
            self.generate_phi(x,t)
        return self.phi_dict[(x,t)]
    
# Helper Functions
def theta_eq(theta, t, M, f):
    retval = mp.mpf(0)
    
    if db.DEBUG_LEVEL >= db.DEBUG_LEVEL_EXTREME:
        print("retval: {}".format(type(retval)))
        print("theta: {}: {}".format(type(theta), theta))
    
    for e in f:
        retval += theta/(theta - f[tuple(sorted(e))])
    retval = t*retval/M - 1

    if db.DEBUG_LEVEL >= db.DEBUG_LEVEL_EXTREME:
        print("theta = {}".format(theta))
        print("theta_eq = {}".format(retval))
        

    return retval

def derivative_theta_eq(theta, t, M, f):
    retval = mp.mpf(0)
    for e in f:
        retval += f[tuple(sorted(e))]/((theta - f[tuple(sorted(e))])*(theta - f[tuple(sorted(e))]))
    retval = -t*retval/M

    if db.DEBUG_LEVEL >= db.DEBUG_LEVEL_EXTREME:
        print("d_theta_eq = {}".format(retval))

    return retval