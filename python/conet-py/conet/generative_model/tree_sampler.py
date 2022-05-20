from typing import TypeVar, List

import networkx as nx
import random

T = TypeVar('T')


class RandomWalkTreeSampler:
    @staticmethod
    def sample_tree(nodes: List[T]) -> nx.DiGraph:
        """
            Sample random uniform directed tree built from (all) nodes in @nodes.
            Always nodes[0] is assumed to be the root.
        """
        tree = nx.DiGraph()
        tree.add_node(nodes[0])
        current_node = nodes[0]
        while len(nodes) > tree.size() + 1:  # nx.DiGraph.size evaluates to the number of graph edges
            new_node = random.sample(range(0, len(nodes)), 1)[0]
            if nodes[new_node] not in set(tree.nodes):
                tree.add_edge(current_node, nodes[new_node])
            current_node = nodes[new_node]
        return tree
