import numpy as np


def _check_ratios_data(ratios, names, chromosomes, lengths_between_bins):
    if not isinstance(ratios, np.ndarray) or len(ratios.shape) != 2:
        raise RuntimeError("Ratios should be 2-D numpy array!")
    if len({ratios.shape[1], len(names), len(chromosomes), len(lengths_between_bins)}) > 1:
        print("dupa")
        raise RuntimeError("Length of all other arguments must match number of columns of @ratios!")


class InputData:
    def __init__(self, ratios, names, chromosomes, lengths_between_loci):
        _check_ratios_data(ratios, names, chromosomes, lengths_between_loci)
        self.ratios = ratios
        self.loci_names = names
        self.chromosomes = chromosomes
        self.lengths_between_loci = lengths_between_loci
        self.has_counts = False

    def save(self, data_dir):
        self.__save_ratios(data_dir)
        if self.has_counts:
            self.__save_counts(data_dir)

    def __save_ratios(self, data_dir):
        separator = ";"
        with open(data_dir + "ratios", 'w') as file:
            file.write(separator.join(self.loci_names))
            file.write('\n')
            file.write(separator.join(map(str, self.chromosomes)))
            file.write('\n')
            file.write(separator.join(map(str, self.lengths_between_loci)))
            file.write('\n')

            for i in range(0, self.ratios.shape[0]):
                file.write(separator.join(map(str, list(self.ratios[i, :]))))
                if i < self.ratios.shape[0] - 1:
                    file.write('\n')

    def __save_counts(self, data_dir):
        pass


class SyntheticInputData(InputData):

    def __init__(self, ratios, bps_matrix, cell_to_node, no_loci):
        InputData.__init__(self, ratios, list(map(str, range(0, no_loci))), [1 for i in range(0, no_loci)],
                           range(0, no_loci))
        self.cell_to_node = cell_to_node
        self.bps_matrix = bps_matrix
