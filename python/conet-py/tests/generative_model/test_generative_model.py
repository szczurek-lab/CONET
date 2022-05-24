import networkx as nx
import pytest as pytest

from conet.generative_model.cn_sampler import CNSampler
from conet.generative_model.counts_generator import CountsGenerator
from conet.generative_model.event_tree_generator import EventTreeGenerator

CELLS = 200
NO_LOCI = 1000
cn_s = CNSampler.create_default_sampler()
t_gen = EventTreeGenerator(cn_sampler=cn_s, tree_size=10, no_loci=NO_LOCI)


@pytest.mark.parametrize('execution_number', range(5))
def test_counts_breakpoints_generation(execution_number):
    tree = t_gen.generate_random_tree()
    d_gen = CountsGenerator(cn_s, tree)

    counts, attachment, corrected_counts, brkp_matrix = d_gen.generate_data(NO_LOCI, CELLS)

    for cell in range(0, CELLS):
        node = attachment[cell]
        path_to_root = nx.shortest_path(tree.tree, tree.get_root(), node)
        path_to_root.remove(tree.get_root())
        expected_counts = [2 for _ in range(0, NO_LOCI)]

        assert node not in tree.trunk_nodes

        for ancestor in path_to_root:
            for bin in range(ancestor.start, ancestor.end):
                expected_counts[bin] = ancestor.copy_number
            assert brkp_matrix[cell, ancestor.start] == brkp_matrix[cell, ancestor.end] == 1

        assert expected_counts == list(counts[cell, :])
