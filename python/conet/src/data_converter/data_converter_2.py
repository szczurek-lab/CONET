from typing import List

import pandas as pd
import numpy as np

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

        return sum([self.corr_counts_df.iloc[i + 1][self.BIN_START_COLUMN] - self.corr_counts_df.iloc[i][self.BIN_START_COLUMN] for i in range(loc1, loc2) if self.get_locus_chr(i) == self.get_locus_chr(loc1)])

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


class DataConverter2:

    def __init__(self,
                 event_length_normalizer: int,  neutral_cn=2.0):

        self.neutral_cn = neutral_cn
        self.event_length_normalizer = event_length_normalizer


    def create_CoNET_input_files(self, out_path: str, corrected_counts: CorrectedCounts, chromosomes=None):
        diffs = self.__create_diff_matrix(corrected_counts, chromosomes)
        counts, squared_counts = self.__create_sum_and_squared_counts_matrices(corrected_counts, chromosomes)

        with open(out_path + "cell_names2", "w") as outfile:
            outfile.write("\n".join(corrected_counts.get_cells_names()))
        np.savetxt(out_path + "ratios2", diffs, delimiter=";", fmt='%.6f')
        np.savetxt(out_path + "counts2", counts, delimiter=";", fmt='%.6f')
        np.savetxt(out_path + "counts_squared2", squared_counts, delimiter=";", fmt='%.6f')

    def __create_diff_matrix(self, corrected_counts: CorrectedCounts, chromosomes: List[int] = None) -> np.ndarray:
        NON_RATIO_DIFF_COLS = {
            "chromosome": 1,
            "bin_start": 0,
            "total_bin_length_to_next_bin": 2
        }

        diffs_columns = corrected_counts.get_cells_count() + len(NON_RATIO_DIFF_COLS)
        brkp_candidates_indices = corrected_counts.get_brkp_candidate_loci_idx(chromosomes)
        diffs = np.zeros((len(brkp_candidates_indices), diffs_columns))

        for i in range(0, len(brkp_candidates_indices)):
            loci_index = brkp_candidates_indices[i]
            loci_chr = corrected_counts.get_locus_chr(loci_index)
            if corrected_counts.is_first_locus_in_chr(loci_index):
                diffs[i, len(NON_RATIO_DIFF_COLS):diffs_columns] = corrected_counts.get_locus_corrected_counts(loci_index) - self.neutral_cn
            else:
                diffs[i, len(NON_RATIO_DIFF_COLS):diffs_columns] = \
                    corrected_counts.get_locus_corrected_counts(loci_index) - \
                    corrected_counts.get_locus_corrected_counts(loci_index - 1)
            diffs[i, 1] = loci_chr
            diffs[i, 0] = corrected_counts.get_locus_bin_start(loci_index)

            if i < len(brkp_candidates_indices) - 1:
                next_loci_index = brkp_candidates_indices[i + 1]
            else:
                next_loci_index = corrected_counts.get_loci_count()
            diffs[i, 2] = corrected_counts.get_total_bin_length_between_loci(loci_index,
                                                                         next_loci_index) / self.event_length_normalizer
        return np.transpose(diffs)

    def __get_sum_of_counts(self, corrected_counts: CorrectedCounts, left_loci, right_loci):
        counts = corrected_counts.get_counts_between_loci(left_loci, right_loci)
        return np.array(counts.sum())

    def __get_square_of_counts(self, corrected_counts: CorrectedCounts, left_loci, right_loci):
        counts = corrected_counts.get_counts_between_loci(left_loci, right_loci)
        squares = np.square(counts)
        return np.array(squares.sum())

    def __create_sum_and_squared_counts_matrices(self, corrected_counts: CorrectedCounts, chromosomes):
        loci_candidate_indices = corrected_counts.get_brkp_candidate_loci_idx(chromosomes)
        counts = np.zeros((len(loci_candidate_indices), corrected_counts.get_cells_count() + 1))
        squared_counts = np.zeros((len(loci_candidate_indices), corrected_counts.get_cells_count() + 1))

        for i in range(0, len(loci_candidate_indices)):
            loci_index = loci_candidate_indices[i]

            if i == len(loci_candidate_indices) - 1 or corrected_counts.is_last_locus_in_chr(loci_index):
                counts[i, :] = np.zeros((corrected_counts.get_cells_count() + 1))
                squared_counts[i, :] = np.zeros((corrected_counts.get_cells_count() + 1))
            else:
                next_loci_index = loci_candidate_indices[i + 1]
                counts[i, 1:] = self.__get_sum_of_counts(corrected_counts, loci_index, next_loci_index)
                squared_counts[i, 1:] = self.__get_square_of_counts(corrected_counts, loci_index, next_loci_index)
                counts[i, 0] = next_loci_index - loci_index
                squared_counts[i, 0] = next_loci_index - loci_index
        return np.transpose(counts), np.transpose(squared_counts)
