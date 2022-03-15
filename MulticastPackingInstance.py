# Generating Multicast Packing Problem Instances

# Imports
import networkx as nx
import random

import GlobalConstants


def is_connected(G):
    return nx.is_k_edge_connected(G, 1)

def get_random_connected_graph():
    G = nx.empty_graph(GlobalConstants.NUM_NODES)
    while not is_connected(G):
        G = nx.random_regular_graph(GlobalConstants.DEGREE, GlobalConstants.NUM_NODES)
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
    def __init__(self, num_requests, max_request_size, graph=None, requests=None):
        if graph is None:
            graph = get_random_connected_graph()
        if requests is None:
            requests = list()
            for i in range(num_requests):
                k = random.randint(2, max_request_size)
                requests.append(MulticastRequest(k, graph))
        self.graph = graph
        self.num_edges = len(self.graph.edges())
        self.requests = requests
        self.num_requests = len(requests)