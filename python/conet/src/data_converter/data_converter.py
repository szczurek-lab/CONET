import pandas as pd
import numpy as np


class DataConverter:
    human_chr_lengths = {
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

    def __init__(self, corrected_counts_path: str, delimiter: chr, default_bin_length: int,
                 event_length_normalizer: int, add_chromosome_ends=False, neutral_cn=2.0):
        self.corrected_counts = pd.read_csv(corrected_counts_path, sep=delimiter, header=0, low_memory=False)
        self.NON_CELL_COLUMNS = 5
        self.NON_RATIO_DIFFS_COLUMNS = 3
        self.CHROMOSOME_COLUMN = 0
        self.BIN_START_COLUMN = 1
        self.LOCI_COUNT = self.corrected_counts.shape[0]
        self.neutral_cn = neutral_cn
        self.cells = self.corrected_counts.shape[1] - self.NON_CELL_COLUMNS
        self.default_bin_length = default_bin_length
        self.event_length_normalizer = event_length_normalizer
        loci_candidate_indices = []
        for i in range(0, self.corrected_counts.shape[0]):
        	if self.corrected_counts.iloc[i,4] == 1:
        		loci_candidate_indices.append(i)
        self.loci_candidate_indices = loci_candidate_indices
        if add_chromosome_ends:
            self.__add_chromosome_ends()

    # indices must be numbered from 0
    def create_CoNET_input_files(self, out_path: str, add_chr_ends_to_indices=False, chromosomes=None):
        if add_chr_ends_to_indices:
            self.loci_candidate_indices = self.__add_chromosome_ends_to_indices(self.loci_candidate_indices)
        if chromosomes is not None:
            self.loci_candidate_indices = self.__filter_loci_to_chromosomes(self.loci_candidate_indices, chromosomes)

        diffs = self.__create_diff_matrix(self.loci_candidate_indices)
        counts, squared_counts = self.__create_sum_and_squared_counts_matrices(self.loci_candidate_indices)
        
        
        with open(out_path + "cell_names", "w") as outfile:
    	    outfile.write("\n".join(list(self.corrected_counts.columns.values.tolist())[self.NON_CELL_COLUMNS:]))
        np.savetxt(out_path + "ratios", diffs, delimiter=";", fmt='%.6f')
        np.savetxt(out_path + "counts", counts, delimiter=";", fmt='%.6f')
        np.savetxt(out_path + "counts_squared", squared_counts, delimiter=";", fmt='%.6f')

    def __filter_loci_to_chromosomes(self, loci_candidate_indices, chromosomes):
        return [loci for loci in loci_candidate_indices if
                self.corrected_counts.iloc[loci, self.CHROMOSOME_COLUMN] in chromosomes]

    def __get_loci_chromosome(self, loci):
        return self.corrected_counts.iloc[loci][self.CHROMOSOME_COLUMN]

    def __get_bin_length(self, loci):
        if loci == self.LOCI_COUNT - 1 or self.__get_loci_chromosome(loci) != self.__get_loci_chromosome(loci + 1):
            return self.default_bin_length
        else:
            return self.corrected_counts.iloc[loci + 1][self.BIN_START_COLUMN] - self.corrected_counts.iloc[loci][
                self.BIN_START_COLUMN]

    def __create_diff_matrix(self, loci_candidate_indices):
        diffs_columns = self.cells + self.NON_RATIO_DIFFS_COLUMNS
        counts_size = self.corrected_counts.shape[1]
        diffs = np.zeros((len(loci_candidate_indices), diffs_columns))
        for i in range(0, len(loci_candidate_indices)):
            loci_index = loci_candidate_indices[i]
            loci_chr = self.__get_loci_chromosome(loci_index)
            if loci_index == 0 or loci_chr != self.__get_loci_chromosome(loci_index - 1):
                diffs[i, self.NON_RATIO_DIFFS_COLUMNS:diffs_columns] = self.corrected_counts.iloc[loci_index][
                                                                       self.NON_CELL_COLUMNS:counts_size] - self.neutral_cn
            else:
                diffs[i, self.NON_RATIO_DIFFS_COLUMNS:diffs_columns] = \
                    self.corrected_counts.iloc[loci_index][self.NON_CELL_COLUMNS:counts_size] - \
                    self.corrected_counts.iloc[loci_index - 1][self.NON_CELL_COLUMNS:counts_size]
            diffs[i, 1] = loci_chr
            diffs[i, 0] = self.corrected_counts.iloc[loci_index][self.BIN_START_COLUMN]

        for i in range(0, len(loci_candidate_indices)):
            loci_index = loci_candidate_indices[i]
            loci_chr = self.__get_loci_chromosome(loci_index)
            if i < len(loci_candidate_indices) - 1:
                next_loci_index = loci_candidate_indices[i + 1]
            else:
                next_loci_index = self.LOCI_COUNT
            j = loci_index

            bins_lengths = 0
            while j < next_loci_index and self.__get_loci_chromosome(j) == loci_chr:
                bins_lengths += self.__get_bin_length(j)
                j = j + 1
            diffs[i, 2] = bins_lengths / self.event_length_normalizer
        return np.transpose(diffs)

    def __get_sum_of_counts(self, left_loci, right_loci):
        counts = self.corrected_counts.iloc[left_loci:right_loci, self.NON_CELL_COLUMNS:]
        return np.array(counts.sum())

    def __get_square_of_counts(self, left_loci, right_loci):
        counts = self.corrected_counts.iloc[left_loci:right_loci, self.NON_CELL_COLUMNS:]
        squares = np.square(counts)
        return np.array(squares.sum())

    def __create_sum_and_squared_counts_matrices(self, loci_candidate_indices):
        counts = np.zeros((len(loci_candidate_indices), self.cells + 1))
        squared_counts = np.zeros((len(loci_candidate_indices), self.cells + 1))

        for i in range(0, len(loci_candidate_indices)):
            loci_index = loci_candidate_indices[i]
            loci_chr = self.__get_loci_chromosome(loci_index)

            if i == len(loci_candidate_indices) - 1 or loci_chr != self.__get_loci_chromosome(
                    loci_candidate_indices[i + 1]):
                counts[i, :] = np.zeros((self.cells + 1))
                squared_counts[i, :] = np.zeros((self.cells + 1))
            else:
                next_loci_index = loci_candidate_indices[i + 1]
                counts[i, 1:] = self.__get_sum_of_counts(loci_index, next_loci_index)
                squared_counts[i, 1:] = self.__get_square_of_counts(loci_index, next_loci_index)
                counts[i, 0] = next_loci_index - loci_index
                squared_counts[i, 0] = next_loci_index - loci_index
        return np.transpose(counts), np.transpose(squared_counts)

    def __add_chromosome_ends(self):
        i = 0
        while i < self.corrected_counts.shape[0]:
            loci_chr = self.__get_loci_chromosome(i)
            if i == self.corrected_counts.shape[0] - 1 or loci_chr != self.__get_loci_chromosome(i + 1):
                line = np.full(self.corrected_counts.shape[1], self.neutral_cn)
                line[self.CHROMOSOME_COLUMN] = loci_chr
                line[self.BIN_START_COLUMN] = self.human_chr_lengths[loci_chr]
                line[self.BIN_START_COLUMN + 1] = self.human_chr_lengths[loci_chr] + self.default_bin_length
                line[self.BIN_START_COLUMN + 2] = self.default_bin_length
                chunk = self.corrected_counts.iloc[:(i+1)].append(self.corrected_counts.iloc[i], ignore_index=True)
                chunk.iloc[i+1] = line
                self.corrected_counts = chunk.append(self.corrected_counts.iloc[(i+1):], ignore_index=True)
                i = i + 2
            else:
                i = i + 1
        self.LOCI_COUNT = self.corrected_counts.shape[0]

    def __add_chromosome_ends_to_indices(self, loci_indices):
        i = 0
        chr_ends = []
        while i < self.LOCI_COUNT:
            if i == self.LOCI_COUNT - 1 or self.__get_loci_chromosome(i) != self.__get_loci_chromosome(i + 1):
                after_end = [loci for loci in loci_indices if loci >= i]
                before_end = [loci for loci in loci_indices if loci < i]
                after_end = list(map(lambda x: x + 1, after_end))
                loci_indices = before_end + after_end
                chr_ends.append(i)
            i = i +1
        loci_indices = loci_indices + chr_ends
        loci_indices.sort()
        return loci_indices
