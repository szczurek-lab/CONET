from dataclasses import dataclass
from typing import List, Tuple, Dict

import networkx as nx
import numpy as np

from conet.generative_model.node_label import NodeLabel
from conet.generative_model.utils import get_logger

logger = get_logger(__name__)


@dataclass
class EventTree:
    """
    Represents tree of events in the generative model.

    Each tree has a distinguished section - the trunk (including the root).
    No cells are allowed to be attached to the trunk.
    """
    tree: nx.DiGraph
    trunk_nodes: List[NodeLabel]

    def get_root(self) -> NodeLabel:
        return self.trunk_nodes[0]

    def get_breakpoint_loci(self) -> List[int]:
        indices = list(set([x.start for x in self.tree.nodes] + [x.end for x in self.tree.nodes]))
        indices.sort()
        return indices

    def get_node_counts_and_breakpoints(self, neutral_cn: int, node: NodeLabel, no_loci: int) -> Tuple[
        np.ndarray, np.ndarray]:
        """
            Calculate CN profile and breakpoints of a cell attached to node @node with @no_loci loci in chromosome.
        """
        path = nx.shortest_path(self.tree, self.get_root(), node)
        node_counts = np.zeros(no_loci, dtype=np.float64)
        node_counts.fill(neutral_cn)
        node_brkps = np.zeros(no_loci, dtype=np.float64)
        for ancestor in path:
            if ancestor == self.get_root():
                continue
            for j in range(ancestor.start, ancestor.end):
                node_counts[j] = ancestor.copy_number
            node_brkps[ancestor.start] = 1
            node_brkps[ancestor.end] = 1
        return node_counts, node_brkps

    def get_node_to_counts_and_brkps(self, neutral_cn: int, no_loci: int) -> Tuple[
        Dict[NodeLabel, np.ndarray], Dict[NodeLabel, np.ndarray]]:
        node_to_counts = {}
        node_to_brkps = {}

        for node in self.tree:
            node_counts, node_brkps = self.get_node_counts_and_breakpoints(neutral_cn, node, no_loci)
            node_to_counts[node] = node_counts
            node_to_brkps[node] = node_brkps

        return node_to_counts, node_to_brkps

    @staticmethod
    def is_valid(tree: 'EventTree', no_loci: int, neutral_cn: int) -> bool:
        logger.info("Testing validity of generated tree...")
        for node in tree.tree:
            if node.copy_number == 0 and any(
                    [c.copy_number > 0 and node.overlaps(c) for c in nx.descendants(tree.tree, node)]):
                logger.info("Tree has CN 0 followed by another CN in descendant nodes. Validation failed.")
                return False

            node_counts, node_brkps = tree.get_node_counts_and_breakpoints(neutral_cn, node, no_loci)
            for i in range(0, len(node_brkps)):
                if node_brkps[i] == 1 and node_counts[i] == node_counts[i - 1]:
                    logger.info("Tree has non-identifiable breakpoint loci. Validation failed.")
                    # This condition means that breakpoint locus can't be identified even if one knows correct counts
                    return False
        logger.info("Generated tree is valid.")
        return True
