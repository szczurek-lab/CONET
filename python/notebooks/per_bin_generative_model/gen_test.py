import networkx as nx
import numpy as np

from conet.data_converter.corrected_counts import CorrectedCounts

from conet.data_converter.data_converter import DataConverter
from matplotlib import pyplot as plt
from networkx.drawing.nx_agraph import graphviz_layout

from conet import CONET, CONETParameters, InferenceResult
from conet.generative_model import CNSampler, EventTree, CountsGenerator, EventTreeGenerator

bin_dir = "./"
loci = 1000
tree_size = 10
cells = 260

cn_s = CNSampler.create_default_sampler()

t_gen = EventTreeGenerator(cn_sampler=cn_s, tree_size=tree_size, no_loci=loci)
tree: EventTree = t_gen.generate_random_tree()
d_gen = CountsGenerator(cn_s, tree)
counts, attachment, corrected_counts, brkp_matrix = d_gen.generate_data(loci, cells)

corrected_counts.to_csv(bin_dir + "counts_synthetic", sep=";", index=False)

# Draw simulated tree
plt.figure(3, figsize=(12, 12))
pos = graphviz_layout(tree.tree, prog="dot")
nx.draw(tree.tree, pos=pos, with_labels=True, node_color="grey", node_size=60, verticalalignment="bottom",
        font_size=20, edge_color="grey", font_color="green")
plt.show()

cc = CorrectedCounts(bin_dir + "counts_synthetic", delimiter=';')
data_converter = DataConverter(event_length_normalizer=loci)
data_converter.create_CoNET_input_files(bin_dir, corrected_counts=cc)

conet = CONET(bin_dir + "CONET")
params = CONETParameters(tree_structure_prior_k1=0.5, data_dir=bin_dir, counts_penalty_s1=100000,
                         counts_penalty_s2=100000,
                         param_inf_iters=200000, seed=2167, mixture_size=2, pt_inf_iters=500000,
                         use_event_lengths_in_attachment=False,
                         event_length_penalty_k0=1, output_dir = "./output/")
conet.infer_tree(params)

result = InferenceResult('./output/', cc)

# Draw inferred tree tree
plt.figure(3, figsize=(12, 12))
pos = graphviz_layout(result.inferred_tree, prog="dot")
nx.draw(result.inferred_tree, pos=pos, with_labels=True, node_color="grey", node_size=60, verticalalignment="bottom",
        font_size=20, edge_color="grey", font_color="green")
plt.show()

inferred_counts = np.transpose(result.get_inferred_copy_number(2))
real_counts = counts
