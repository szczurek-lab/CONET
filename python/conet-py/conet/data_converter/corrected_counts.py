from typing import List

import numpy as np
import pandas as pd


__HUMAN_CHR_LENGTHS__ = {
    1: 249250621,
    2: 243199373,
    3: 198022430,
    4: 191154276,
    5: 180915260,
    6: 171115067,
    7: 159138663,
    8: 146364022,
    9: 141213431,
    10: 135534747,
    11: 135006516,
    12: 133851895,
    13: 115169878,
    14: 107349540,
    15: 102531392,
    16: 90354753,
    17: 81195210,
    18: 78077248,
    19: 59128983,
    20: 63025520,
    21: 48129895,
    22: 51304566,
    23: 155270560,
    24: 59373566
}


class CorrectedCounts:
    """
        Thin wrapper over corrected counts matrix.
    """
    CHROMOSOME_COLUMN = 0
    BIN_START_COLUMN = 1
    BIN_END_COLUMN = 2
    BIN_WIDTH_COLUMN = 3
    BRKP_CANDIDATE_LOCUS_COLUMN = 4

    def __init__(self, path: str, delimiter: str = ";"):
        self.corr_counts_df = pd.read_csv(path, sep=delimiter, header=0, low_memory=False)

    def get_loci_count(self):
        return self.corr_counts_df.shape[0]

    def get_cells_count(self) -> int:
        return self.corr_counts_df.shape[1] - 5

    def get_locus_bin_start(self, locus: int) -> int:
        return self.corr_counts_df.iloc[locus][self.BIN_START_COLUMN]

    def is_first_locus_in_chr(self, loc: int) -> bool:
        return loc == 0 or self.get_locus_chr(loc) != self.get_locus_chr(loc - 1)

    def is_last_locus_in_chr(self, loc: int) -> bool:
        return loc == self.get_loci_count() - 1 or self.get_locus_chr(loc) != self.get_locus_chr(loc + 1)

    def get_locus_chr(self, loc: int) -> int:
        return self.corr_counts_df.iloc[loc][CorrectedCounts.CHROMOSOME_COLUMN]

    def get_locus_corrected_counts(self, loc: int):
        return self.corr_counts_df.iloc[loc][5:self.corr_counts_df.shape[1]]

    def get_cells_names(self) -> List[str]:
        return list(self.corr_counts_df.columns.values.tolist())[5:]

    def get_brkp_candidate_loci_idx(self, chromosomes: List[int] = None) -> List[int]:
        if chromosomes is None:
            return [i for i in range(0, self.corr_counts_df.shape[0])
                    if self.corr_counts_df.iloc[i, CorrectedCounts.BRKP_CANDIDATE_LOCUS_COLUMN] == 1]
        return [i for i in range(0, self.corr_counts_df.shape[0])
                if self.corr_counts_df.iloc[i, CorrectedCounts.BRKP_CANDIDATE_LOCUS_COLUMN] == 1 and self.get_locus_chr(
                i) in chromosomes]

    def get_total_bin_length_between_loci(self, loc1: int, loc2: int) -> int:
        if self.is_last_locus_in_chr(loc1):
            return self.corr_counts_df.iloc[loc1][self.BIN_WIDTH_COLUMN]

        return sum([self.corr_counts_df.iloc[i + 1][self.BIN_START_COLUMN] - self.corr_counts_df.iloc[i][self.BIN_START_COLUMN]
                    for i in range(loc1, loc2) if (not self.is_last_locus_in_chr(i)) and self.get_locus_chr(i + 1) == self.get_locus_chr(loc1)])

    def get_counts_between_loci(self, loc1: int, loc2: int):
        return self.corr_counts_df.iloc[loc1:loc2, 5:]

    def add_chromosome_ends(self, neutral_cn: int, end_bin_length: int) -> 'CorrectedCounts':
        i = 0
        while i < self.corr_counts_df.shape[0]:
            if self.is_last_locus_in_chr(i):
                self.corr_counts_df = pd.concat([self.corr_counts_df.iloc[:(i + 1)],
                                                 pd.DataFrame([self.__create_chromosome_end_corrected_counts(self.get_locus_chr(i), neutral_cn, end_bin_length)], columns=self.corr_counts_df.columns),
                                                 self.corr_counts_df.iloc[(i + 1):]
                                                 ], ignore_index=True, axis=0)
                i = i + 2
            else:
                i = i + 1
        return self

    def __create_chromosome_end_corrected_counts(self, chr: int, neutral_cn: int, end_bin_length: int) -> np.ndarray:
        line = np.full(self.corr_counts_df.shape[1], neutral_cn)
        line[self.CHROMOSOME_COLUMN] = chr
        line[self.BIN_START_COLUMN] = __HUMAN_CHR_LENGTHS__[chr]
        line[self.BIN_END_COLUMN] = __HUMAN_CHR_LENGTHS__[chr] + end_bin_length
        line[self.BIN_WIDTH_COLUMN] = end_bin_length
        line[self.BRKP_CANDIDATE_LOCUS_COLUMN] = 1
        return line

    def __is_last_locus_in_chr(self, loc: int) -> bool:
        return loc == self.corr_counts_df.shape[0] - 1 or self.get_locus_chr(loc) != self.get_locus_chr(loc + 1)

