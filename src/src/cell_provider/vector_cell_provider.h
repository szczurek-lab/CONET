#ifndef VECTOR_CELL_PROVIDER_H
#define VECTOR_CELL_PROVIDER_H
#include <vector>
#include <map>

#include "../utils/breakpoints/breakpoints.h"

/**
* Container for cell data
* Stores 2D matrix of cell data - [i,j] indices 
* corresponds to i-th log loci count of j-th cell
*/
template <class Real_t> class VectorCellProvider {
private:
	const size_t lociCount;
	size_t cellCount{ 0 };
	std::vector<std::vector<Real_t>> cells{ lociCount };
	/**
	* For pair of breakpoints (br1, br2) length of event (br1, br2) 
	* is equal to eventLengths[br1] - eventLengths[br2]
	*/
    
    std::vector<Real_t> betweenBinsLengths;
	/**
	* Contains forst breakpoitn from each chromosome (expcept for the forst one) ad lociCount as a last element;
	*/
	std::vector<size_t> chromosomeMarkers;

	std::map<size_t, std::string> loci_to_name_map;

	/* Classes for additional likelihood */
	std::vector<Real_t> counts_scores_lengths;
	std::vector<std::vector<Real_t>> sum_counts;
	std::vector<std::vector<Real_t>> squared_counts;
public:

	VectorCellProvider(size_t lociCount) : lociCount{ lociCount }  {
	}

	VectorCellProvider(size_t lociCount, std::vector<Real_t> between, std::map<size_t, std::string> loci_to_name_map) :
		lociCount{ lociCount }, betweenBinsLengths{between}, loci_to_name_map{ loci_to_name_map } {
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



    void setBetweenLengths(std::vector<Real_t> between) {
		this->betweenBinsLengths = between;
	}

	std::vector<size_t> getChromosomeMarkers() {
		return this->chromosomeMarkers;
	}

	void setChromosomeMarkers(std::vector<size_t> chromosomeMarkers) {
		this->chromosomeMarkers = chromosomeMarkers;
	}

	Real_t getEventLength(BreakpointPair brkp) const {
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
	
	void postCell(std::vector<Real_t> cell) {
		for (size_t i = 0; i < cell.size(); i++) {
			cells[i].push_back(cell[i]);
		}
		cellCount++;
	}

	const std::vector<std::vector<Real_t>> &getCellsToLoci() const {
		return cells;
	}

	size_t getCellsCount() const { return cellCount; }

	size_t getLociCount() const { return lociCount; }
};
#endif // !VECTOR_CELL_PROVIDER_H
