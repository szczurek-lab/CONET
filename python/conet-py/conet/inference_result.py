import json
from typing import Tuple, Dict, List

import networkx as nx
import numpy as np

from conet.data_converter.corrected_counts import CorrectedCounts


class InferenceResult:
    def __init__(self, output_path: str, cc: CorrectedCounts):
        if not output_path.endswith('/'):
            output_path = output_path + '/'

        self.__cc = cc
        self.__tree = self.__load_inferred_tree(output_path + "inferred_tree")
        self.__attachment = self.__load_attachment(output_path + "inferred_attachment")
        self.inferred_tree = self.__get_pretty_tree()
        self.attachment = self.__get_pretty_attachment()

    def get_inferred_copy_numbers(self, neutral_cn: int) -> np.ndarray:
        counts = np.full((self.__cc.get_loci_count(), self.__cc.get_cells_count()), neutral_cn)
        bin_bitmap = np.zeros((self.__cc.get_loci_count(), self.__cc.get_cells_count()))
        attachment = self.__attachment.copy()

        for node in nx.traversal.dfs_postorder_nodes(self.__tree, source=(0, 0)):
            if node == (0, 0):
                continue
            attached_cells = [c for c in range(0, len(attachment)) if attachment[c] == node]
            summed_ccs = 0.0
            bin_count = (node[1] - node[0]) * len(attached_cells) - np.sum(
                bin_bitmap[range(node[0], node[1]), :][:, attached_cells])  # non zero elements in bit map
            for cell in attached_cells:
                for bin in range(node[0], node[1]):
                    if bin_bitmap[bin, cell] == 0:
                        summed_ccs += self.__cc.get_locus_corrected_counts(bin)[cell]

            for cell in attached_cells:
                for bin in range(node[0], node[1]):
                    if bin_bitmap[bin, cell] == 0:
                        counts[bin, cell] = round(summed_ccs / bin_count)
                        bin_bitmap[bin, cell] = 1

            parent_node = next(self.__tree.predecessors(node))
            for c in attached_cells:
                attachment[c] = parent_node
        return counts

    def __node_to_pretty(self, node: Tuple[int, int]) -> Dict:
        if node == (0, 0):
            return {}
        return {
            "chr": self.__cc.get_locus_chr(node[0]),
            "bin_start": int(self.__cc.get_locus_bin_start(node[0])),
            "bin_end": int(self.__cc.get_locus_bin_start(node[1]))
        }

    def __get_pretty_attachment(self) -> List[Dict]:
        return [self.__node_to_pretty(n) for n in self.__attachment]

    def __get_pretty_tree(self) -> nx.DiGraph:

        pretty_tree = nx.DiGraph()
        pretty_tree.add_edges_from(
            [(json.dumps(self.__node_to_pretty(e[0])), json.dumps(self.__node_to_pretty(e[1]))) for e in self.__tree.edges]
        )
        return pretty_tree

    def __load_inferred_tree(self, path: str) -> nx.DiGraph:
        def __int_tuple_from_str(s: str) -> Tuple[int, int]:
            s = s.strip().replace('(', '').replace(')', '')
            return int(s.split(',')[0]), int(s.split(',')[1])

        brkp_candidates = self.__cc.get_brkp_candidate_loci_idx()
        with open(path) as f:
            edges = [(__int_tuple_from_str(line.split('-')[0]), __int_tuple_from_str(line.split('-')[1])) for line in
                     f.readlines()]
            edges = [((brkp_candidates[e[0][0]], brkp_candidates[e[0][1]]),
                      (brkp_candidates[e[1][0]], brkp_candidates[e[1][1]])) for e in edges]
            tree = nx.DiGraph()
            tree.add_edges_from(edges)
            return tree

    def __load_attachment(self, path: str) -> List[Tuple[int, int]]:
        brkp_candidates = self.__cc.get_brkp_candidate_loci_idx()
        with open(path) as f:
            return [(brkp_candidates[int(line.split(';')[2])], brkp_candidates[int(line.split(';')[3])]) for line in
                    f.readlines()]

