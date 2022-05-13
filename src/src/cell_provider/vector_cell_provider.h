#ifndef VECTOR_CELL_PROVIDER_H
#define VECTOR_CELL_PROVIDER_H
#include <vector>
#include <map>

#include "../types.h"

/**
* Container for cell data
* Stores 2D matrix of cell data - [i,j] indices 
* corresponds to i-th loci count of j-th cell
*/
template <class Real_t> class VectorCellProvider {
private:
	const size_t loci_count;
	size_t cell_count{ 0 };
	std::vector<std::vector<Real_t>> cells{ loci_count };
	/**
	* For pair of breakpoints (br1, br2) length of event (br1, br2) 
	* is equal to eventLengths[br1] - eventLengths[br2]
	*/
    
    std::vector<Real_t> betweenBinsLengths;
	/**
	* Contains forst breakpoitn from each chromosome (expcept for the forst one) ad loci_count as a last element;
	*/
	std::vector<size_t> chromosome_markers;


	/* Classes for additional likelihood */
	std::vector<Real_t> counts_scores_lengths;
	std::vector<std::vector<Real_t>> sum_counts;
	std::vector<std::vector<Real_t>> squared_counts;
public:

	VectorCellProvider(size_t loci_count, std::vector<size_t> &chrom_m, std::vector<Real_t> &between) : loci_count{ loci_count }, chromosome_markers {chrom_m},  betweenBinsLengths { between } {
	}

	void post_cell(std::vector<Real_t> &cell) {
		for (size_t i = 0; i < cell.size(); i++) {
			cells[i].push_back(cell[i]);
		}
		cell_count++;
	}

	Real_t get_counts_score_length_of_bin(size_t bin) const {
		return counts_scores_lengths[bin];
	}

	std::vector<std::vector<Real_t>> get_sum_counts() const {
		return sum_counts;
	}

	std::vector<std::vector<Real_t>> get_squared_counts() const {
		return squared_counts;
	}

	void post_counts_dispersion_data(std::vector<Real_t> v1, std::vector<std::vector<Real_t>> v2, std::vector<std::vector<Real_t>> v3) {
		counts_scores_lengths = v1;
		sum_counts = v2;
		squared_counts = v3;
	}

	inline Real_t get_cell_bin_counts(const size_t cell, const size_t bin) const {
		return sum_counts[cell][bin];
	}

	inline Real_t get_cell_bin_squared_counts(const size_t cell, const size_t bin) const {
		return squared_counts[cell][bin];
	}

	std::vector<Real_t> cells_to_vector()
	{
		std::vector<Real_t> result;
		for (auto &el : cells)
		{
			for (auto & cell_ratio : el)
			{
				result.push_back(cell_ratio);
			}
		}
		return result;
	}

	std::vector<size_t> getChromosomeMarkers() {
		return this->chromosome_markers;
	}

	Real_t get_event_length(Event brkp) const {
		if (betweenBinsLengths.empty()) {
			return (Real_t)brkp.second - brkp.first;
		}
        size_t bin = brkp.first;
        Real_t result = 0.0;
        while(bin < brkp.second) {
            result += betweenBinsLengths[bin];
            bin++;
        }
		return result;
	}

	const std::vector<std::vector<Real_t>> &getCellsToLoci() const {
		return cells;
	}

	size_t get_cells_count() const { return cell_count; }

	size_t get_loci_count() const { return loci_count; }

};
#endif // !VECTOR_CELL_PROVIDER_H
