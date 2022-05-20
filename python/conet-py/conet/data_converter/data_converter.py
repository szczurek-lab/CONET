from dataclasses import dataclass
from typing import List

import numpy as np

from conet.data_converter.corrected_counts import CorrectedCounts


@dataclass
class DataConverter:
    """
    Converts corrected counts matrix to input files expected by CONET
    """
    event_length_normalizer: int
    neutral_cn: int = 2

    def create_CoNET_input_files(self, out_path: str, corrected_counts: CorrectedCounts, chromosomes: List[int]=None):
        diffs = self.__create_diff_matrix(corrected_counts, chromosomes)
        counts, squared_counts = self.__create_sum_and_squared_counts_matrices(corrected_counts, chromosomes)

        with open(out_path + "cell_names", "w") as outfile:
            outfile.write("\n".join(corrected_counts.get_cells_names()))
        np.savetxt(out_path + "ratios", diffs, delimiter=";", fmt='%.6f')
        np.savetxt(out_path + "counts", counts, delimiter=";", fmt='%.6f')
        np.savetxt(out_path + "counts_squared", squared_counts, delimiter=";", fmt='%.6f')

    def __create_diff_matrix(self, cc: CorrectedCounts, chromosomes: List[int] = None) -> np.ndarray:
        NON_RATIO_DIFF_COLS = {
            "chromosome": 1,
            "bin_start": 0,
            "total_bin_length_to_next_bin": 2
        }

        diffs_columns = cc.get_cells_count() + len(NON_RATIO_DIFF_COLS)
        brkp_candidates_indices = cc.get_brkp_candidate_loci_idx(chromosomes)
        diffs = np.zeros((len(brkp_candidates_indices), diffs_columns))
        brkp_candidates_indices.append(cc.get_loci_count())

        for i in range(0, len(brkp_candidates_indices) - 1):
            loci_index = brkp_candidates_indices[i]
            if cc.is_first_locus_in_chr(loci_index):
                diffs[i, len(NON_RATIO_DIFF_COLS):diffs_columns] = cc.get_locus_corrected_counts(loci_index) - self.neutral_cn
            else:
                diffs[i, len(NON_RATIO_DIFF_COLS):diffs_columns] = cc.get_locus_corrected_counts(loci_index) - cc.get_locus_corrected_counts(loci_index - 1)

            next_loci_index = brkp_candidates_indices[i + 1]

            diffs[i, NON_RATIO_DIFF_COLS["chromosome"]] = cc.get_locus_chr(loci_index)
            diffs[i, NON_RATIO_DIFF_COLS["bin_start"]] = cc.get_locus_bin_start(loci_index)
            diffs[i, NON_RATIO_DIFF_COLS["total_bin_length_to_next_bin"]] = \
                cc.get_total_bin_length_between_loci(loci_index, next_loci_index) / self.event_length_normalizer

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
