from dataclasses import dataclass
from typing import Generator, Dict

import numpy as np
import random
import networkx as nx
import pandas as pd

from conet.generative_model.cn_sampler import CNSampler
from conet.generative_model.event_tree import EventTree
from conet.generative_model.node_label import NodeLabel


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
        attachment_sampler = self.__node_attachment_sampler()

        counts = np.zeros([0, no_loci], dtype=np.float64)
        brkp_matrix = np.zeros([0, no_loci], dtype=np.float64)
        node_to_counts, node_to_brkps = self.tree.get_node_to_counts_and_brkps(self.cn_sampler.neutral_cn, no_loci)
        attachment = [next(attachment_sampler) for _ in range(0, no_of_cells)]

        for cell in range(0, no_of_cells):
            counts = np.vstack([counts, node_to_counts[attachment[cell]]])
            brkp_matrix = np.vstack([brkp_matrix, node_to_brkps[attachment[cell]]])

        return counts, attachment, self.__convert_cc_to_CONET_format(self.__add_noise_to_counts(counts)), brkp_matrix

    def __node_attachment_sampler(self) -> Generator[NodeLabel, NodeLabel, None]:
        """
        Generator which yields subsequent random attachment points on the tree.
        """
        probs = self.__create_node_attachment_probabilities()

        while True:
            yield random.choices(list(probs), k=1, weights=list(probs.values()))[0]

    def __create_node_attachment_probabilities(self) -> Dict[NodeLabel, float]:
        """
        Probability of cell being attached to given node is proportional to  e^{0.1 * depth} where depth is node's depth in the tree.
        """
        depths = nx.shortest_path_length(self.tree.tree, source=self.tree.get_root())
        depths = dict([(key, np.exp(0.1 * val)) for key, val in depths.items() if key not in self.tree.trunk_nodes])
        sum_depths = sum(list(depths.values()))
        for key, value in depths.items():
            depths[key] = np.exp(0.1 * depths[key]) / sum_depths
        return depths

    def __add_noise_to_counts(self, counts: np.ndarray) -> np.ndarray:
        return np.vectorize(lambda x: self.cn_sampler.sample_corrected_count(x))(counts)

    def __convert_cc_to_CONET_format(self, cc: np.ndarray) -> pd.DataFrame:
        no_cells = cc.shape[0]
        no_loci = cc.shape[1]
        cc = np.transpose(cc)

        # Add columns required by CONET
        add = np.full([no_loci, 5], 1.0, dtype=np.float64)
        add[:, 1] = range(0, no_loci)# Bin start
        add[:, 2] = range(1, no_loci + 1)# Bin end
        add[:, 4] = 0# Breakpoint markers
        add[self.tree.get_breakpoint_loci(), 4] = 1
        full_counts = np.hstack([add, cc])
        return pd.DataFrame(full_counts, columns=["chr", "start", "end", "width", "candidate_brkp"] + ["cell" + str(i) for i in range(0, no_cells)])

