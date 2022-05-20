from dataclasses import dataclass
from typing import Generator

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
        attachment_sampler = self.__node_attachment_sampler(self.tree)

        counts = np.zeros([0, no_loci], dtype=np.float64)
        brkp_matrix = np.zeros([0, no_loci], dtype=np.float64)
        node_to_counts, node_to_brkps = self.tree.get_node_to_counts_and_brkps(self.cn_sampler.neutral_cn, no_loci)
        attachment = [next(attachment_sampler) for _ in range(0, no_of_cells)]

        for cell in range(0, no_of_cells):
            counts = np.vstack([counts, node_to_counts[attachment[cell]]])
            brkp_matrix = np.vstack([brkp_matrix, node_to_brkps[attachment[cell]]])

        return counts, attachment, self.__create_corrected_counts_df(self.__add_noise_to_counts(counts), self.tree), brkp_matrix

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

    def __create_corrected_counts_df(self, counts: np.ndarray, tree: EventTree) -> pd.DataFrame:
        no_cells = counts.shape[0]
        no_loci = counts.shape[1]
        counts = np.transpose(counts)

        # Add columns required by CONET
        add = np.full([no_loci, 5], 1.0, dtype=np.float64)
        add[:, 1] = range(0, no_loci)# Bin start
        add[:, 2] = range(1, no_loci + 1)# Bin end
        add[:, 4] = 0# Breakpoint markers
        add[tree.get_breakpoint_loci(), 4] = 1
        full_counts = np.hstack([add, counts])
        return pd.DataFrame(full_counts, columns=["chr", "start", "end", "width", "candidate_brkp"] + ["cell" + str(i) for i in range(0, no_cells)])

