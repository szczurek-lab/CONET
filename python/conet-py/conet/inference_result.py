import json
from typing import Tuple, Dict, List

import networkx as nx
import numpy
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

    def dump_results_to_dir(self, dir: str, neutral_cn: int) -> None:
        if not dir.endswith('/'):
            dir = dir + '/'
        with open(dir + "inferred_attachment", 'w') as f:
            cells = self.__cc.get_cells_names()
            for i in range(0, len(self.attachment)):
                f.write(";".join(
                    [cells[i], str(i), f"{int(self.attachment[i]['chr'])}_{self.attachment[i]['bin_start']}",
                     f"{int(self.attachment[i]['chr'])}_{self.attachment[i]['bin_end']}\n"]))

        numpy.savetxt(dir + "inferred_counts", X=self.get_inferred_copy_numbers(neutral_cn), delimiter=";")

        def __node_to_str(node: Tuple[int, int]) -> str:
            if node == (0, 0):
                return "(0,0)"
            chr = int(self.__cc.get_locus_chr(node[0]))
            return f"({chr}_{self.__cc.get_locus_bin_start(node[0])},{chr}_{self.__cc.get_locus_bin_start(node[1])})"

        with open(dir + "inferred_tree", 'w') as f:
            for edge in self.__tree.edges:
                f.write(f"{__node_to_str(edge[0])}-{__node_to_str(edge[1])}\n")

    def get_inferred_copy_numbers(self, neutral_cn: int) -> np.ndarray:
        cumulated_attachment = self.__get_cumulated_attachment()
        # Matrix where each bin, cell pair is mapped to integer representing its cluster
        cell_bin_clusters = np.zeros((self.__cc.get_loci_count(), self.__cc.get_cells_count()))

        self.__divide_cell_bin_pairs_into_clusters((0, 0), np.zeros((self.__cc.get_loci_count())), cell_bin_clusters,
                                                   cumulated_attachment)
        regions = list(np.unique(cell_bin_clusters))
        regions.remove(0)  # Cluster with id 0 corresponds to tree's root

        corrected_counts = self.__cc.get_as_numpy()
        counts = np.full((self.__cc.get_loci_count(), self.__cc.get_cells_count()), neutral_cn)

        for r in regions:
            counts[cell_bin_clusters == r] = round(np.median(corrected_counts[cell_bin_clusters == r]))

        # If chromosome ends have been added we want to delete them from the result - this ensures that inferred
        # counts matrix and input CC matrix have the same number of rows
        if self.__cc.added_chr_ends:
            chromosome_ends = [bin for bin in range(0, self.__cc.get_loci_count()) if
                               self.__cc.is_last_locus_in_chr(bin)]
            counts = np.delete(counts, chromosome_ends, axis=0)

        return counts

    def __get_cumulated_attachment(self) -> Dict[Tuple[int, int], List[int]]:
        """
            Create dictionary where every node is mapped to list of cells (represented by their indices) which are
            attached to subtree rooted at the node.
        """
        cum_attach = {}
        for node in nx.traversal.dfs_postorder_nodes(self.__tree, source=(0, 0)):
            cum_attach[node] = [cell for cell in range(0, len(self.__attachment)) if self.__attachment[cell] == node]
            for child in self.__tree.successors(node):
                cum_attach[node].extend(cum_attach[child])
        return cum_attach

    def __divide_cell_bin_pairs_into_clusters(self, node: Tuple[int, int], regions: np.ndarray,
                                              cell_bin_regions: np.ndarray,
                                              cumulated_attachment: Dict[Tuple[int, int], List[int]]) -> None:
        regions_copy = regions.copy()
        self.__update_regions(regions, node, np.max(cell_bin_regions))

        for bin in range(node[0], node[1]):
            cell_bin_regions[bin, cumulated_attachment[node]] = regions[bin]
        for child in self.__tree.successors(node):
            self.__divide_cell_bin_pairs_into_clusters(child, regions, cell_bin_regions, cumulated_attachment)
        regions[0:regions.shape[0]] = regions_copy

    def __update_regions(self, regions: np.ndarray, event: Tuple[int, int], max_region_id: int) -> None:
        region_id_to_new_id = {}
        for bin in range(event[0], event[1]):
            if regions[bin] not in region_id_to_new_id:
                region_id_to_new_id[regions[bin]] = max_region_id + 1
                max_region_id += 1
            regions[bin] = region_id_to_new_id[regions[bin]]

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
            [(json.dumps(self.__node_to_pretty(e[0])), json.dumps(self.__node_to_pretty(e[1]))) for e in
             self.__tree.edges]
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
