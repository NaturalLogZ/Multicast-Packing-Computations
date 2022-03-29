# Generating Multicast Packing Problem Instances

# Imports
import networkx as nx
import random

import GlobalConstants


def is_connected(G):
    return nx.is_k_edge_connected(G, 1)

def get_random_connected_graph(n, m):
    degree = 2*m // n
    G = nx.empty_graph()
    while not is_connected(G):
        G = nx.random_regular_graph(degree, n)
    return G

class MulticastRequest:
    def __init__(self, size, graph, source=None, recipients=None):
        if source is None and recipients is None:
            multicast_group = random.sample(list(graph), size)
            source = multicast_group[0]
            recipients = set(multicast_group[1:])
        self.source = source
        self.recipients = recipients
        self.size = len(recipients) + 1
    
    def multicast_group(self):
        return self.recipients.union([self.source])

class MulticastPackingInstance:
    def __init__(self, n=GlobalConstants.NUM_NODES, m=GlobalConstants.NUM_EDGES,
                 num_requests=GlobalConstants.NUM_MULTICAST_REQUESTS, 
                 max_request_size=GlobalConstants.MAX_MULTICAST_SIZE,
                 delay=GlobalConstants.DELAY,
                 graph=None, 
                 requests=None
                ):
        if graph is None:
            graph = get_random_connected_graph(n, m)
            for u,v in graph.edges():
                graph[u][v]['delay'] = round(random.triangular(0, delay/2, delay*m/(n*n)))
                print(graph[u][v]['delay'])
        if requests is None:
            requests = list()
            for i in range(num_requests):
                requests.append(MulticastRequest(max_request_size, graph))
        self.graph = graph
        self.num_edges = len(self.graph.edges())
        self.requests = requests
        self.num_requests = len(requests)
        self.delay = delay