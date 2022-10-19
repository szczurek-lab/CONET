import os

import numpy as np
import random
import json

from conet.data_converter.corrected_counts import CorrectedCounts
from conet.data_converter.data_converter import DataConverter
from conet import CONET, CONETParameters, InferenceResult
from conet.generative_model import CNSampler, EventTree, CountsGenerator, EventTreeGenerator

bin_dir = "./"  # path to dir with CONET executable
loci = 500
tree_size = 10
cells = 260

random.seed(2234)
np.random.seed(234)
ITERS = int(os.environ["TRIES"])


def score():
    # generate event tree and cell data
    cn_s = CNSampler.create_default_sampler()
    t_gen = EventTreeGenerator(cn_sampler=cn_s, tree_size=tree_size, no_loci=loci)
    tree: EventTree = t_gen.generate_random_tree()
    d_gen = CountsGenerator(cn_s, tree)
    counts, attachment, corrected_counts, brkp_matrix = d_gen.generate_data(loci, cells)

    # Convert corrected counts dataframe to input files for CONET
    cc = CorrectedCounts(corrected_counts)
    DataConverter(event_length_normalizer=loci).create_CoNET_input_files(out_path=bin_dir, corrected_counts=cc)

    # Run CONET
    conet = CONET(bin_dir + "CONET")
    params = CONETParameters(tree_structure_prior_k1=0.5, data_dir=bin_dir, counts_penalty_s1=100000,
                             counts_penalty_s2=100000,
                             param_inf_iters=200000, seed=2167, mixture_size=2, pt_inf_iters=100000,
                             use_event_lengths_in_attachment=False,
                             event_length_penalty_k0=1, output_dir="./output/")
    conet.infer_tree(params)

    result = InferenceResult('./output/', cc)

    real_edges = [((e[0].start, e[0].end), (e[1].start, e[1].end)) for e in tree.tree.edges]

    def conv_node(n):
        if n == {}:
            return 0, 0
        return n['bin_start'], n['bin_end']

    inf_edges = [(conv_node(json.loads(e[0])), conv_node(json.loads(e[1]))) for e in result.inferred_tree.edges]

    s = len(set(inf_edges).intersection(set(real_edges))) / tree_size
    print(f"Result of the try {s}")
    return s


if __name__ == "__main__":
    mean_score = sum([score() for _ in range(0, ITERS)]) / ITERS
    if mean_score >= 0.5:
        print(f"Test successful")
    else:
        raise RuntimeError("Test failed")
