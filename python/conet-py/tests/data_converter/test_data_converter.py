import os

import pandas as pd
import pytest

from conet.data_converter.corrected_counts import CorrectedCounts
from conet.data_converter.data_converter import DataConverter


@pytest.fixture(autouse=True)
def test_cleanup():
    yield
    artifacts = ["cell_names", "counts", "counts_squared", "ratios"]
    for a in artifacts:
        if os.path.exists("./data_converter/data/" + a):
            os.remove("./data_converter/data/"+a)


def test_counts_creation():
    cc_df = pd.read_csv("./data_converter/data/data1.csv", sep=';', header=0, low_memory=False)
    cc = CorrectedCounts(cc_df)
    DataConverter(1, 2).create_CoNET_input_files("./data_converter/data/", cc)

    counts = pd.read_csv("./data_converter/data/counts", sep=";", header=None)
    counts_squared = pd.read_csv("./data_converter/data/counts_squared", sep=";", header=None)

    assert counts.shape == (4, 4)
    assert counts_squared.shape == (4, 4)
    assert pd.read_csv("./data_converter/data/expected_summed_counts_data1", sep=";", header=None).equals(counts)
    assert pd.read_csv("./data_converter/data/expected_squared_counts_data1", sep=";", header=None).equals(counts_squared)


