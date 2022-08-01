# Solver for the Multicast Packing Problem that uses Simplex, and performs
# column generation by using heuristics first, and then an exact algo for
# Steiner Tree once the heuristic fails to give an improving direction

import GlobalConstants
from Solvers.Impls.PureColGenMcpSolver import PureColGenMcpSolver, cost
from ColumnGenerators.Impls.ExactMulticastPackingColumnGeneratorIP import ExactMulticastPackingColumnGeneratorIP

class WarmStartColGenMcpSolver(PureColGenMcpSolver):
    def __init__(self, instance=None, block_approx=2):
        super().__init__(instance, 2)
        self.transition_iteration = None
        self.first_stop_flag = 0b00
        
    def perform_checks_and_updates(self, x):
        t = self.t
        super().perform_checks_and_updates(x)
        # Check if we failed to generate trees that will improve the solution
        # and that we haven't already transitioned to the exact column generation algorithm, 
        if (self.stop_flag & ~self.first_stop_flag):
            # Channge the column generation algorithm
            self.column_generator = ExactMulticastPackingColumnGeneratorIP(self.instance, self.reduced_LP)
            # Record the iteration and stop flag
            self.transition_iteration = self.iteration
            self.first_stop_flag = self.stop_flag
             # Reset the stop flag 
            self.stop_flag = 0b0
            
        