import numpy as np
import random
import networkx as nx
from conet.src.per_bin_model.node_label import *

class TreeGenerator:
    def __init__(self, cn_sampler):
        self.cn_sampler = cn_sampler
        
    def generate_random_tree(self, no_loci, tree_size):
        tree_nodes_list = self.__sample_node_labels(tree_size, no_loci)
        tree, trunk = self.__random_tree(tree_nodes_list)
        while not self.__check_valid(tree, no_loci):
            tree_nodes_list = self.__sample_node_labels(tree_size, no_loci)
            tree, trunk = self.__random_tree(tree_nodes_list)
        return tree, trunk
    
    def __sample_node_label(self, labels):
        norming_constant = 0
        for z in labels:
            norming_constant += np.exp(-0.05 * (z[1] - z[0]))
        while True:
            label = random.sample(labels, 1)[0]
            if np.random.uniform() < np.exp(-0.05 * (label[1] - label[0])) / norming_constant:
                return label

    def __sample_node_labels(self, tree_size, no_loci):
        candidates = list(filter(lambda x: x[0] < x[1] and x[1] - x[0] < 100, [(a, b) for a in range(0, no_loci) for b in range(0, no_loci)]))
        result = list()
        for n in range(0, tree_size):
            label = self.__sample_node_label(candidates)
            candidates.remove(label)
            result.append(NodeLabel(label[0], label[1], self.cn_sampler.sample_CN()))
        result.insert(0, NodeLabel(0, 0, self.cn_sampler.NEUTRAL_CN))
        return result
    
    def __random_tree(self, nodes):
        tree = nx.DiGraph()
        trunk_size = np.random.randint(int(0.1*len(nodes)), int(0.4*len(nodes)))
        current_node = 0
        trunk = [nodes[current_node]]
        for i in range(0, trunk_size):
            for j in range(1, trunk_size):
                if nodes[j].start > nodes[j+1].start or (nodes[j].start == nodes[j+1].start and nodes[j].end > nodes[j+1].end):
                    nodes[j], nodes[j+1] = nodes[j+1], nodes[j]
        
        for i in range(1, trunk_size + 1):
            tree.add_edge(nodes[current_node], nodes[i])
            trunk.append(nodes[i])
            current_node = i
       
        bit_map = np.zeros(len(nodes) - trunk_size)
        bit_map[current_node-trunk_size] = 1
        visited = 1
        while len(nodes) - trunk_size > visited:
            new_node = random.sample(range(trunk_size, len(nodes)), 1)[0]
            if bit_map[new_node - trunk_size] == 0.0:
                bit_map[new_node - trunk_size] = 1.0
                tree.add_edge(nodes[current_node], nodes[new_node])
                visited += 1
            current_node = new_node
        return tree, trunk
    
    def __check_valid(self, tree, no_loci):
        for node in tree:
            if node.copy_number == 0:
                children = [x[1] for x in list(nx.edge_dfs(tree,list(tree.nodes)[0]))]
                for c in children:
                    if c.copy_number > 0 and node.overlaps(c):
                        return False
                    
        counts = np.zeros([1,no_loci], dtype=np.float64)
        counts.fill(self.cn_sampler.NEUTRAL_CN)
        
        brkps = np.zeros([1, no_loci], dtype=np.float64)
        node_to_counts = {}
        node_to_brkps = {}
        
        for node in tree:
            path = nx.shortest_path(tree, list(tree.nodes)[0], node)
            node_counts = np.copy(counts)
            node_brkps = np.copy(brkps)
            for i in range(0, len(path)):
                ancestor = path[i]
                for j in range(ancestor.start, ancestor.end):
                    node_counts[0, j] = ancestor.copy_number           
                    node_brkps[0,ancestor.start] = 1
                    node_brkps[0,ancestor.end] = 1
            node_to_counts[node] = node_counts    
            node_to_brkps[node] = node_brkps
            for i in range(node_to_brkps[node].shape[1]):
                if node_to_brkps[node][0, i] == 1 and node_to_counts[node][0, i] == node_to_counts[node][0, i -1]:
                        return False
        return True
