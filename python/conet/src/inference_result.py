import numpy as np
import networkx as nx

import conet.src.ratios_distribution as rd

class InferenceResult:

    def __init__(self, data_dir):
        self.bp_matrix = self.__read_bps_matrix(data_dir + "inferred_breakpoints")
        self.attachment = self.__read_attachment(data_dir + "inferred_attachment")
        self.tree = self.__read_tree(data_dir + "inferred_tree")
        self.distribution = self.__read_distriubution(data_dir + "inferred_distribution")
	
    def __read_bps_matrix(self, path):
        with open(path, 'r') as f:
            matrix= [[float(x) for x in line.split(';')] for line in f]
            return np.array(matrix)
            
    def __read_attachment(self, path):
        with open(path, 'r') as f:
            attachment = [[x for x in line.split(';')] for line in f]
        for i in range(0, len(attachment)):
            attachment[i][2] = attachment[i][2].replace('\n', "")
        attachment_dir = {}
        for x in attachment:
            attachment_dir[x[0]] = (x[1], x[2])
        return attachment_dir

    def __read_tree(self, path):
        tree = nx.DiGraph()
        with open(path, 'r') as f:
            for line in f:
                parent = self.__node_from_text(line.split("-")[0])
                child = self.__node_from_text(line.split("-")[1])
                tree.add_edge(parent, child)
        return tree

    def __node_from_text(self, text):
        if text == '(0,0)':
            return '0', '0'
        text = text.replace('\n', "").replace('(', '').replace(')', '')
        loci_left = text.split(',')[0]
        loci_right = text.split(',')[1]
        return loci_left, loci_right
    
    def __read_distriubution(self, path):
        sep = ';'
        var_0 = 0
        weights = []
        means = []
        variances = []
        first_line = True
        with open(path, 'r') as f:
            for line in f:
                if first_line:
                    first_line = False
                    var_0 = float(line.split(sep)[1].replace('\n',''))**2
                else:
                    weights.append(float(line.split(sep)[0]))
                    means.append(float(line.split(sep)[1]))
                    variances.append(float(line.split(sep)[0].replace('\n',''))**2)
        return rd.RatiosDistribution(weights, means, variances, var_0)
