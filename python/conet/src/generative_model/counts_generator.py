from dataclasses import dataclass
from typing import Generator

import numpy as np
import random
import networkx as nx

from src.generative_model.cn_sampler import CNSampler
from src.generative_model.event_tree import EventTree
from src.generative_model.node_label import NodeLabel


@dataclass
class CountsGenerator:
    cn_sampler: CNSampler
    tree: EventTree

    def generate_data(self, no_loci: int, no_of_cells: int):
        """
              Generates arrays of CN profiles, Corrected Counts, breakpoints
              and attachment for one chromosome consisting of @no_loci bins
              for @no_of_cells cells
        """
        attachment_sampler = self.__node_attachment_sampler(self.tree)

        counts = np.zeros([0, no_loci], dtype=np.float64)
        brkp_matrix = np.zeros([0, no_loci], dtype=np.float64)
        node_to_counts, node_to_brkps = self.tree.get_node_to_counts_and_brkps(self.cn_sampler.neutral_cn, no_loci)
        attachment = [next(attachment_sampler) for _ in range(0, no_of_cells)]

        for cell in range(0, no_of_cells):
            counts = np.vstack([counts, node_to_counts[attachment[cell]]])
            brkp_matrix = np.vstack([brkp_matrix, node_to_brkps[attachment[cell]]])

        return counts, attachment, self.__add_noise_to_counts(counts), brkp_matrix

    def __node_attachment_sampler(self, tree: EventTree) -> Generator[NodeLabel, NodeLabel, None]:
        depths = nx.shortest_path_length(tree.tree, source=tree.get_root())
        depths = dict([(key, np.exp(0.1 * val)) for key, val in depths.items() if key not in tree.trunk_nodes])
        sum_depths = sum(list(depths.values()))
        for key, value in depths.items():
            depths[key] = np.exp(0.1 * depths[key]) / sum_depths

        attachment_probs = list(depths.values())
        while True:
            yield random.choices(list(depths), k=1, weights=attachment_probs)[0]

    def __add_noise_to_counts(self, counts: np.ndarray) -> np.ndarray:
        return np.vectorize(lambda x: self.cn_sampler.sample_corrected_count(x))(counts)
