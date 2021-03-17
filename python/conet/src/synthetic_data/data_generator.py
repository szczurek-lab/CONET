import numpy as np
import networkx as nx
import random

import conet.src.input_data as cid


class DataGenerator:
    def generate_random_ratios(self, no_loci, tree, ratios_distribution, no_of_cells):
        # dictionary node_label : set of all loci with bps which a cell attached to this node has
        dictionary_nodes_bp_loci = self.__node_bp_loci_dictionary(tree)

        # create empty ratios_matrix
        ratios_matrix, bps_matrix = np.empty([0, no_loci]), np.zeros([0, no_loci], dtype=np.float64)

        # cell index -> node of attachment
        dictionary_cell_node = {}

        attachment_probs = np.random.dirichlet(np.ones(len(tree.nodes) - 1)).transpose()

        for cell in range(0, no_of_cells):
            attachment_node = random.choices(list(tree.nodes)[1:], k=1, weights=list(attachment_probs))[0]
            dictionary_cell_node[cell] = attachment_node
            new_ratios, new_bps = self.__synthetic_data_one_cell_all_loci(no_loci,
                                                                          dictionary_nodes_bp_loci[attachment_node],
                                                                          ratios_distribution)
            ratios_matrix = np.vstack([ratios_matrix, new_ratios])
            bps_matrix = np.vstack([bps_matrix, new_bps])
        return cid.SyntheticInputData(ratios_matrix, bps_matrix, self.__dictionary_cell_node_to_array(dictionary_cell_node),
                                      no_loci)

    def generate_random_tree(self, no_loci, tree_size):
        tree_nodes_list = self.__sample_node_labels(tree_size, no_loci)
        return self.__random_tree(tree_nodes_list)

    def __node_bp_loci_dictionary(self, tree):
        dict_of_node_and_ancestral_bps = dict()
        tree_nodes = tree.nodes() - {(0, 0)}
        for node in tree_nodes:
            temp_bps_list = [node[0], node[1]]
            for el in set(nx.ancestors(tree, node) - {(0, 0)}):
                temp_bps_list.append(el[0])
                temp_bps_list.append(el[1])
            dict_of_node_and_ancestral_bps.update({node: set(temp_bps_list)})
        return dict_of_node_and_ancestral_bps

    def __synthetic_data_one_cell_all_loci(self, no_loci, set_of_loci_with_bps, ratios_distribution):
        ratios = []
        bps = np.zeros([1, no_loci], dtype=np.float64)
        for i in range(0, no_loci):
            if i in set_of_loci_with_bps:
                ratios.append(round(ratios_distribution.sample_breakpoint_ratio(), 3))
                bps[0, i] = 1
            else:
                ratios.append(round(ratios_distribution.sample_non_breakpoint_ratio(), 3))
        return np.array(ratios), bps

    def __dictionary_cell_node_to_array(self, dictionary):
        result = np.zeros((len(dictionary), 3), dtype=int)
        for key, value in dictionary.items():
            result[key, 0] = key
            result[key, 1] = value[0]
            result[key, 2] = value[1]
        return result

    def __sample_node_label(self, labels):
        norming_constant = 0
        for z in labels:
            norming_constant += np.exp(z[0] - z[1])
        while True:
            label = random.sample(labels, 1)[0]
            if np.random.uniform() < np.exp(label[0] - label[1]) / norming_constant:
                return label

    def __sample_node_labels(self, tree_size, no_loci):
        candidates = list(filter(lambda x: x[0] < x[1], [(a, b) for a in range(0, no_loci) for b in range(0, no_loci)]))
        result = list()
        for n in range(0, tree_size):
            label = self.__sample_node_label(candidates)
            candidates.remove(label)
            result.append(label)
        result.insert(0, (0, 0))
        return result

    def __random_tree(self, nodes):
        tree = nx.DiGraph()
        bit_map = np.zeros(len(nodes))
        bit_map[0] = 1
        visited = 1
        current_node = 0
        while len(nodes) > visited:
            new_node = random.sample(range(0, len(nodes)), 1)[0]
            if bit_map[new_node] == 0.0:
                bit_map[new_node] = 1.0
                tree.add_edge(nodes[current_node], nodes[new_node])
                visited += 1
            current_node = new_node
        return tree
