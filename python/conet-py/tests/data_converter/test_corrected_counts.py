import pandas as pd

from conet.data_converter.data_converter import CorrectedCounts



def test_corrected_counts_basic_ops():
    cc = CorrectedCounts(pd.read_csv("./data_converter/data/data1.csv", sep=';', header=0, low_memory=False))

    assert cc.get_loci_count() == 8

    cc.add_chromosome_ends(2, 1)

    assert cc.get_loci_count() == 10
    assert cc.get_cells_count() == 3

    assert cc.get_brkp_candidate_loci_idx() == [1, 4, 5, 6, 7, 9]

    assert cc.get_total_bin_length_between_loci(0, 10) == 249250621 - 1
    assert cc.get_total_bin_length_between_loci(0, 3) == 3