#ifndef VECTOR_CELL_PROVIDER_H
#define VECTOR_CELL_PROVIDER_H
#include <vector>
#include <map>

#include "../types.h"

/**
* Container for CONET input data
* Stores 2D matrix of corrected counts - [i,j] indices 
* corresponds to i-th bin count of j-th cell
* And two matrices for calculation of count dispersion penlty. 
*/
template <class Real_t> class VectorCellProvider {
private:
	const size_t loci_count;
	size_t cell_count{ 0 };
	std::vector<std::vector<Real_t>> corrected_counts{ loci_count };
	/**
	* For pair of breakpoints (br1, br2) length of event (br1, br2) 
	* is equal to eventLengths[br1] - eventLengths[br2]
	*/
    
    std::vector<Real_t> between_bins_lengths;
	/**
	* Contains indices of first breakpoint from each chromosome (expcept for the first one) and loci_count as a last element;
	*/
	std::vector<size_t> chromosome_markers;


	std::vector<Real_t> counts_scores_regions; // Records sizes of segments between consecutive breakpoint candidate loci
	std::vector<std::vector<Real_t>> summed_counts;
	std::vector<std::vector<Real_t>> squared_counts;
public:

	VectorCellProvider(size_t loci_count, std::vector<size_t> chrom_m, std::vector<Real_t> between) : loci_count{ loci_count }, chromosome_markers {chrom_m},  between_bins_lengths { between } {
	}

	void post_cell(std::vector<Real_t> &cell) {
		for (size_t i = 0; i < cell.size(); i++) {
			corrected_counts[i].push_back(cell[i]);
		}
		cell_count++;
	}

	void post_counts_dispersion_data(std::vector<Real_t> &v1, std::vector<std::vector<Real_t>> &v2, std::vector<std::vector<Real_t>> &v3) {
		counts_scores_regions = v1;
		summed_counts = v2;
		squared_counts = v3;
	}

	std::vector<Real_t> get_counts_scores_regions() const {
		return counts_scores_regions;
	}

	std::vector<std::vector<Real_t>> get_summed_counts() const {
		return summed_counts;
	}

	std::vector<std::vector<Real_t>> get_squared_counts() const {
		return squared_counts;
	}


	std::vector<size_t> getChromosomeMarkers() {
		return this->chromosome_markers;
	}

	const std::vector<std::vector<Real_t>> &get_corrected_counts() const {
		return corrected_counts;
	}

	size_t get_cells_count() const { return cell_count; }

	size_t get_loci_count() const { return loci_count; }

	Real_t get_event_length(Event event) const {
        size_t start_locus = get_event_start_locus(event);
        Real_t result = 0.0;
        while(start_locus < get_event_end_locus(event)) {
            result += between_bins_lengths[start_locus];
            start_locus++;
        }
		return result;
	}
};
#endif // !VECTOR_CELL_PROVIDER_H
