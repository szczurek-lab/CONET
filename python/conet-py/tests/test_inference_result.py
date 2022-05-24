import numpy as np
import pandas as pd

from conet import InferenceResult
from conet.data_converter.corrected_counts import CorrectedCounts


def test_counts_inference_synthetic_data():
    corrected_counts_df = pd.read_csv("data/res3/cc", sep=';', header=0, low_memory=False)
    cc = CorrectedCounts(corrected_counts_df)
    result = InferenceResult('./data/res3/', cc)
    inferred_counts = result.get_inferred_copy_numbers(2)

    assert inferred_counts.tolist() == [[2, 2, 2, 2], [2, 2, 2, 2], [2, 2, 9, 9], [1, 1, 4, 4], [2, 2, 2, 2]]

def test_counts_inference_synthetic_data2():
    corrected_counts_df = pd.read_csv("data/res1/cc", sep=';', header=0, low_memory=False)
    cc = CorrectedCounts(corrected_counts_df)
    result = InferenceResult('./data/res1/', cc)
    inferred_counts = result.get_inferred_copy_numbers(2)

    assert inferred_counts.tolist() == [[1, 1, 2, 2], [1, 1, 9, 9], [2, 2, 99, 99], [2, 2, 2, 2], [2, 2, 2, 2]]

def test_counts_inference_biological_data():
    corrected_counts_df = pd.read_csv("data/res2/corr_counts", sep=';', header=0, low_memory=False)
    cc = CorrectedCounts(corrected_counts_df)
    cc.add_chromosome_ends(2, 150000)
    result = InferenceResult('./data/res2/', cc)
    inferred_counts = result.get_inferred_copy_numbers(2)
    for bin in range(0, inferred_counts.shape[0]):
        for cell in range(0, inferred_counts.shape[1]):
            inferred_counts[bin, cell] = min(inferred_counts[bin, cell], 4)

    expected_counts = pd.read_csv("./data/res2/counts", sep=";", header=None).to_numpy()
    print(np.sum(inferred_counts - expected_counts))
    for bin in range(0, inferred_counts.shape[0]):
        for cell in range(0, inferred_counts.shape[1]):
            if inferred_counts[bin, cell] != expected_counts[bin, cell]:
                print(f"Oh {inferred_counts[bin, cell]} vs {expected_counts[bin, cell]}")
    assert (expected_counts == inferred_counts).all()
