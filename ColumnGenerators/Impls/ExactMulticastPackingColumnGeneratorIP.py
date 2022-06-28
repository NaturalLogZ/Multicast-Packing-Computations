# MCF IP for Column Generation

import gurobipy as gp
import networkx as nx
from networkx.algorithms.approximation.steinertree import steiner_tree as steiner_tree_2
from ColumnGenerators.MulticastPackingColumnGenerator import MulticastPackingColumnGenerator

class ExactMulticastPackingColumnGeneratorIP(MulticastPackingColumnGenerator):
    def __init__(self, instance, reduced_LP):
        super().__init__(instance, reduced_LP)
        self.GurobiModels = list()
        self.variables = list()
        D = nx.DiGraph(instance.graph)
        for i in range(instance.num_requests):
            source = instance.requests[i].source
            recipients = instance.requests[i].recipients
            
            model = gp.Model("Steiner Tree IP (MCF) for Request {}".format(i))
            # Add variables
            x = {e : model.addVar(vtype=gp.GRB.BINARY, 
                                  name="{} selection".format(e)) 
                     for e in D.edges()}
            f = {r : {e: model.addVar(vtype=gp.GRB.BINARY, 
                                      name="{} flow to {}".format(e,r)) 
                          for e in D.edges()} 
                     for r in recipients}
            
            # Add constraints
            for r in recipients:
                for v in D.nodes():
                    inflow = sum(f[r][(u,v)] for u in D.predecessors(v))
                    outflow = sum(f[r][(v,w)] for w in D.successors(v))
                    if v == source:
                        netflow = -1
                    elif v == r:
                        netflow = 1
                    else:
                        netflow = 0
                    
                    model.addConstr(inflow - outflow == netflow,
                                    name="{}-flow conservation for {}".format(r,v))
                
                for e in D.edges():
                    model.addConstr(f[r][e] <= x[e],
                                    name="{} availability for {}-flow".format(e,r))
            
            # Note we don't add the objective as that is determined by
            # the prices in each iteration.
            model.update()
            
            self.variables.append((x,f))
            self.GurobiModels.append(model)
    
    
    def generate_tree(self, i, prices):
        model = self.GurobiModels[i]
        x = self.variables[i][0]
        c = dict()
        
        if isinstance(prices, dict):
            for e in x:
                if e in prices:
                    c[e] = prices[e]
                elif tuple(reversed(e)) in prices:
                    c[e] = prices[tuple(reversed(e))]
                else:
                    raise KeyError('Edge not in prices?!')
        else: # should be a constant numeric value
            for e in x:
                c[e] = prices
                    
        
        model.setObjective(sum(c[e]*x[e] for e in x), gp.GRB.MINIMIZE)
        model.update()
        model.optimize()
        status = model.getAttr(gp.GRB.Attr.Status)
        if status != gp.GRB.OPTIMAL:
            print("Column Gen Error")
            print("{} exited with code {}".format(model, status))
            return self.instance.graph.edge_subgraph()
        
        steinerTree = list()
        for e in self.instance.graph.edges():
            if (x[e].X == 1) or (x[tuple(reversed(e))].X == 1):
                steinerTree.append(e)
        
        return self.instance.graph.edge_subgraph(steinerTree)