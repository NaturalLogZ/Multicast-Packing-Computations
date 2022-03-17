# Solver for the Multicast Packing Problem that uses Simplex

from math import log

import mpmath as mp

import DebugConstants as db
from FrozenDict import FrozenDict
from MulticastPackingSolver import MulticastPackingSolver
from MulticastPackingSolver import dict_innerProd
   
class JansenZhangMinMaxer(MulticastPackingSolver):
    
    def __init__(self, instance=None, block_solver=None, sigma0=3.6):
        super().__init__(instance, block_solver)
        self.sigma = mp.mpf(sigma0)
        self.t = self.sigma/6
        self.M = mp.mpf(len(self.instance.graph.edges()))
        self.theta_dict = dict()
        self.phi_dict = dict()
        self.finished_coordination = False
        self.lambda_of_prev_scaling = 0
        self.w = None
        
        
    def get_next_solution(self):
        if self.iteration == 0:
            x = [dict() for i in range(self.instance.num_requests)]
            step_size = mp.mpf(1)
        else:
            x = self.solution[self.iteration-1]
            
            f_prime = dict()
            for e in self.instance.graph.edges():
                f_prime[tuple(sorted(e))] = mp.mpf(0)
                for i in range(self.instance.num_requests):
                    if self.new_trees[i].has_edge(*e):
                        f_prime[tuple(sorted(e))] += mp.mpf(1)
    
            pf = dict_innerProd(self.p(x, self.t), self.f(x))
            pf_prime = dict_innerProd(self.p(x, self.t), f_prime)
            step_size = ( (self.t*self.theta(x, self.t)*self.toleranceFunction(1)) 
                   /(2*self.M*(pf+pf_prime)) )

        
        new_x = [dict() for i in range(self.instance.num_requests)]

        for i in range(self.instance.num_requests):
            for T in x[i]:
                new_x[i][T] = (1-step_size)*x[i][T]
                
            if self.new_trees[i] in new_x[i]:
                new_x[i][self.new_trees[i]] += step_size
            else:
                new_x[i][self.new_trees[i]] = step_size
            new_x[i] = FrozenDict(new_x[i])
        return tuple(new_x)
    
    def perform_checks_and_updates(self, x):
        if self.w == None:
            if db.DEBUG_LEVEL >= db.DEBUG_LEVEL_THEORY_2:
                print("Scaling Phase 0 Over")
            self.sigma = self.sigma/2
            self.t = self.sigma/6
            self.w = mp.mpf((1+self.sigma)/((1+self.sigma/3)*self.M))
            self.new_trees = self.column_generator.generate_new_trees(self.p(x, self.t))
        elif ( (self.toleranceFunction() <= self.sigma/6) 
              or (self.lamb(x) <= self.w*self.lambda_of_prev_scaling) ):
            if db.DEBUG_LEVEL >= db.DEBUG_LEVEL_THEORY_2:
                print("New Scaling Phase")
            if self.sigma <= self.tol:
                self.stop_flag = True
            else:
                self.lambda_of_prev_scaling= self.lamb(x)
                self.sigma = self.sigma/2
                self.t = self.sigma/6
                self.w = mp.mpf((1+self.sigma)/(1+2*self.sigma))
                self.new_trees = self.column_generator.generate_new_trees(self.p(x, self.t))
                
    
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
    
    # THIS ISN'T NEEDED BY THIS ALGO. I WOULD LIKE TO ADD THIS FUNCTION BUT CURRENTLY THIS WON'T WORK
    def generate_q(self, x, t=0):
        return
        self.multicast_costs[(x,t)] = [None] * len(self.instance.requests)
        for i in range(self.instance.num_requests):
            self.multicast_costs[(x,t)][i] = self.reduced_LP.getConstrByName(
                "Tree Selection for {}".format(i)).getAttr(gp.GRB.Attr.Pi)
            
    def generate_theta(self, x, t):
        if db.DEBUG_LEVEL >= db.DEBUG_LEVEL_FULL:
            #print(t)
            #print(x)
            #print(self.lamb(x))
            #print(self.f(x))
            print("lambda: {} \t t: {} \t M: {}".format(type(self.lamb(x)), type(t), type(self.M)))
            print("Result: {}".format(type(self.lamb(x)/(1-t/self.M))))
            print(self.lamb(x)/(1-t/self.M))
            #print(theta_eq(self.lamb(x)/(1-t/self.M), t, self.M, self.f(x)))
            print(self.lamb(x)/(1-t))
            #print(theta_eq(self.lamb(x)/(1-t), t, self.M, self.f(x)))
        #self.theta_dict[(x,t)] = scipy.optimize.brentq(
        #    theta_eq, self.lamb(x)/(1-t/self.M), self.lamb(x)/(1-t), args=(t, self.M, self.f(x)))
        self.theta_dict[(x,t)] = (
            mp.findroot(lambda theta: theta_eq(theta, t, self.M, self.f(x)),
                        (mp.mpf(self.lamb(x)/(1-t/self.M)), mp.mpf(self.lamb(x)/(1-t))), 
                        solver='ridder')
        )
        
        
        #self.theta_dict[(x,t)] = scipy.optimize.newton(
        #     theta_eq, self.lamb(x)/(1-t), derivative_theta_eq, args=(t, self.M, self.f(x)))
        
    # THE FOLLOWING IS DEPRECATED AND WILL CAUSE ERRORS
    # Define a function that computes theta_t(x)
    #def theta_via_sympy(self, x, t):
    #    M = self.M
    #    f = self.f
    #    
    #    theta = sympy.Symbol('theta')
    #    expression = 1
    #    for m in constraint_indices:
    #        expression -= (t/M) * theta/(theta - f(x,m))
    #    soln_set = sympy.solvers.solveset(expression, theta, domain=sympy.Interval(self.lamb[x], sympy.S.Infinity, True, True))
    #    #assert len(soln_set == 1)
    #    if db.DEBUG_LEVEL >= db.DEBUG_LEVEL_FULL:
    #        print(soln_set)
    #    for soln in soln_set:
    #        retval = soln
    #    return retval
            
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
    
    if db.DEBUG_LEVEL >= db.DEBUG_LEVEL_FULL:
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