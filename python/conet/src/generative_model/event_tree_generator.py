from dataclasses import dataclass
from operator import attrgetter
from typing import List

import numpy as np
import networkx as nx

from src.generative_model.cn_sampler import CNSampler
from src.generative_model.event_tree import EventTree
from src.generative_model.node_label import NodeLabel, EventSampler, RootEvent
from src.generative_model.tree_sampler import RandomWalkTreeSampler
from src.generative_model.utils import sample_conditionally, get_logger

logger = get_logger(__name__)


@dataclass
class EventTreeGenerator:
    """
        Generates random EventTree instances such that EventTree.is_valid is true
    """
    no_loci: int
    tree_size: int
    cn_sampler: CNSampler
    max_event_length: int = 100

    def generate_random_tree(self) -> EventTree:
        return sample_conditionally(
            sampler=lambda: self.__random_tree(self.__sample_node_labels()),
            condition=lambda t: EventTree.is_valid(t, self.no_loci, self.cn_sampler.neutral_cn)
        )

    def __sample_node_labels(self) -> List[NodeLabel]:
        events = EventSampler().sample_events(self.no_loci, self.tree_size)
        labels = [NodeLabel(e[0], e[1], self.cn_sampler.sample_non_neutral_CN()) for e in events]
        # Add root...
        labels.insert(0, NodeLabel(RootEvent[0], RootEvent[1], self.cn_sampler.neutral_cn))
        return labels

    def __sample_tree_trunk_size(self) -> int:
        return np.random.randint(int(0.1 * (self.tree_size + 1)), int(0.4 * (self.tree_size + 1))) + 1

    def __random_tree(self, nodes: List[NodeLabel]) -> EventTree:
        logger.info("Starting event tree generation...")
        tree = nx.DiGraph()
        trunk_nodes = nodes[0:self.__sample_tree_trunk_size()]  # these will always contain at least the root
        logger.info(f"Sampled size of trunk equal to {len(trunk_nodes)}.")
        trunk_nodes.sort(key=attrgetter("start", "end"))

        tree.add_node(trunk_nodes[0])
        # Add trunk to the tree
        for i in range(0, len(trunk_nodes) - 1):
            tree.add_edge(trunk_nodes[i], trunk_nodes[i + 1])
        # Add remaining nodes to the tree
        tree.add_edges_from(
            RandomWalkTreeSampler.sample_tree([trunk_nodes[len(trunk_nodes) - 1]] + nodes[len(trunk_nodes):]).edges)
        logger.info("Successfully sampled unvalidated event tree structure.")
        return EventTree(tree=tree, trunk_nodes=trunk_nodes)
