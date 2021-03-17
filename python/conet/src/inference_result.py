import numpy as np
import networkx as nx


class InferenceResult:

    def __init__(self, data_dir):
        self.bp_matrix = np.fromfile(data_dir + "inferred_breakpoints", dtype=float, sep=';')
        self.attachment = self.__read_attachment(data_dir + "inferred_attachment")
        self.tree = self.__read_tree(data_dir + "inferred_tree")

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
                child = self.__node_from_text(line.split("-")[0])
                tree.add_edge(parent, child)
        return tree

    def __node_from_text(self, text):
        text = text.replace('\n', "").replace('(', '').replace(')', '')
        loci_left = text.split(',')[0]
        loci_right = text.split(',')[1]
        return loci_left, loci_right
