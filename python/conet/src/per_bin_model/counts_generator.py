import numpy as np
import random
import networkx as nx

class CountsGenerator:
    def __init__(self, cn_sampler):
        self.cn_sampler = cn_sampler
        
    def generate_data(self, no_loci, tree, no_of_cells, trunk):
        counts = np.zeros([0, no_loci], dtype=np.float64)
        brkp_matrix = np.zeros([0, no_loci], dtype=np.float64)
        
        node_to_counts, node_to_brkps = self.__fill_node_to_counts(tree, no_loci)
        depths = nx.shortest_path_length(tree, list(tree.nodes)[0])
        for n in trunk:
            del depths[n]

        for key,value in depths.items():
            depths[key] = np.exp(0.1*depths[key])
            
        sum_depths = sum(list(depths.values()))   
        for key,value in depths.items():
            depths[key] = np.exp(0.1*depths[key]) / sum_depths
        
        attachment_probs = list(depths.values())
        attachment = []
        
        for cell in range(0, no_of_cells):
            attachment_node = random.choices(list(depths), k=1, weights=attachment_probs)[0]
            counts = np.vstack([counts, node_to_counts[attachment_node]])
            brkp_matrix = np.vstack([brkp_matrix, node_to_brkps[attachment_node]])
            attachment.append(attachment_node.get_event())
        
        corrected_counts = self.__add_noise_to_counts(counts)
        diff_matrix = self.__create_diff_matrix(corrected_counts)
        return counts, attachment, corrected_counts, diff_matrix, brkp_matrix
    
    def __fill_node_to_counts(self, tree, no_loci):
        counts = np.zeros([1, no_loci], dtype=np.float64)
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
        return node_to_counts, node_to_brkps
    
    def __add_noise_to_counts(self, counts):
        cc = np.copy(counts)
        for i in range(0, counts.shape[0]):
            for j in range(0, counts.shape[1]):
                cc[i,j] = self.cn_sampler.sample_corrected_count(counts[i,j])
        return cc
    
    def __create_diff_matrix(self, corrected_counts):
        diffs = np.copy(corrected_counts)
        
        for cell in range(0, diffs.shape[0]):
            for loci in range(0, diffs.shape[1]):
                if loci == 0:
                    diffs[cell, loci] = np.abs(corrected_counts[cell, loci] - self.cn_sampler.NEUTRAL_CN)
                else:
                    diffs[cell, loci] = np.abs(corrected_counts[cell, loci] - corrected_counts[cell, loci - 1])
        return diffs
