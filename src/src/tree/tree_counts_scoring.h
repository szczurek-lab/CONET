#ifndef COUNTS_SCORING_H
#define COUNTS_SCORING_H

#include "pointer_tree.h"
#include "../cell_provider/vector_cell_provider.h"
#include "../parameters/parameters.h"

template <class Real_t> class CountsScoring {
private:
	VectorCellProvider<Real_t> const *cells;
	std::vector<std::vector<Real_t>> sum_counts;
	std::vector<std::vector<Real_t>> squared_counts;
	std::vector<std::vector<bool>> used;
	std::vector<size_t> regions;
	size_t max_region_id{ 0 };
	std::vector<std::vector<size_t>> regions_cache;

	const size_t DEFAULT_CACHE_SIZE = 1000;

public:

    Real_t tree_size = 0.0;
    Real_t all_bins_count = 0.0;
    static bool print;

	void move_cells_to_parent(PointerTree::NodeHandle child, PointerTree::NodeHandle parent, std::map<BreakpointPair, std::set<size_t>> &attachment) {
		if (attachment.find(child->breakpoints) == attachment.end()) {
			return;
		}
		if (attachment.find(parent->breakpoints) == attachment.end()) {
			attachment[parent->breakpoints] = std::set<size_t>();
		}
		attachment[parent->breakpoints].insert(attachment[child->breakpoints].begin(), attachment[child->breakpoints].end());
		attachment.erase(child->breakpoints);
	}

    Real_t calculate_mean_count_score(PointerTree::NodeHandle node, std::map<BreakpointPair, std::set<size_t>> &attachment) {
		if (attachment.find(node->breakpoints) == attachment.end()) {
			return 0.0;
		}
		std::map<size_t, Real_t> region_to_sums;
		std::map<size_t, Real_t> region_to_bins;
		std::map<size_t, Real_t> region_to_squares;

		std::vector<size_t> event_brkps;
		Real_t result = 0.0;
		for (size_t i = node->breakpoints.first; i < node->breakpoints.second; i++) {
			event_brkps.push_back(i);
			region_to_sums[regions[i]] = 0.0;
			region_to_bins[regions[i]] = 0.0;
			region_to_squares[regions[i]] = 0.0;
		}

		for (auto cell : attachment[node->breakpoints]) {
			for (auto brkp : event_brkps) {
				if (!used[cell][brkp]) {
					region_to_sums[regions[brkp]] += sum_counts[cell][brkp];
					region_to_squares[regions[brkp]] += squared_counts[cell][brkp];
					region_to_bins[regions[brkp]] += cells->get_counts_score_length_of_bin(brkp);
                    used[cell][brkp] = true;
				}
			}
		}
		for (const auto &reg : region_to_sums) {
			if (region_to_bins[reg.first] > 0.0) {
				auto bin_count = region_to_bins[reg.first];
				auto square_sum = region_to_squares[reg.first];
				auto counts_sum = region_to_sums[reg.first];
				Real_t mean_count = reg.second / region_to_bins[reg.first];
				result += (square_sum - 2 * mean_count * counts_sum) / (all_bins_count)+mean_count * mean_count * bin_count / all_bins_count;
				if (mean_count >= 1.5 && mean_count < 2.5) {
					result += bin_count / all_bins_count;
				}
			}
		}
		return result;
	}

	void update_regions_id(std::vector<size_t> &region_ids, BreakpointPair event) {
		max_region_id++;
		std::map<size_t, size_t> region_to_id;
		for (size_t i = event.first; i < event.second; i++) {
			if (region_to_id.find(region_ids[i]) == region_to_id.end()) {
				region_to_id[region_ids[i]] = max_region_id;
				max_region_id++;
			}
			region_ids[i] = region_to_id[region_ids[i]];
		}
	}

    Real_t calculate_log_score(PointerTree::NodeHandle node, std::map<BreakpointPair, std::set<size_t>> &attachment, size_t &cache_id) {
		Real_t result = 0.0;
		auto node_cache = cache_id;

		/* Save event partition of regions */
		if (cache_id >= regions_cache.size()) {
			regions_cache.resize(regions_cache.size() + DEFAULT_CACHE_SIZE);
		}
		regions_cache[cache_id] = regions;
		/* Update event partition of regions by @node event */
		update_regions_id(regions, node->breakpoints);

		cache_id++;
		for (auto child : node->children) {
		    result += calculate_log_score(child, attachment, cache_id);
			move_cells_to_parent(child, node, attachment);
		}
		Real_t mean_count = calculate_mean_count_score(node, attachment);

		/* Restore partition induced by parent */
		regions = regions_cache[node_cache];
		return result + mean_count;
	}

    Real_t calculate_log_score_for_root_nodes() {
        tree_size++;
		Real_t counts_sum = 0.0;
		Real_t bin_count = 0.0;
		Real_t square_sum = 0.0;

		for (size_t cell = 0; cell < used.size(); cell++) {
			for (size_t i = 0; i < used[cell].size(); i++) {
				if (!used[cell][i]) {
					counts_sum += sum_counts[cell][i];
					square_sum += squared_counts[cell][i];
					bin_count += cells->get_counts_score_length_of_bin(i);
				}
			}
		}
		Real_t mean_count = 2.0;
		Real_t result = (square_sum - 2 * mean_count * counts_sum) / all_bins_count + mean_count * mean_count * bin_count / all_bins_count;
        if (bin_count == 0) {
            return 0.0;
            tree_size --;
        }
		return result;
	}

	void refresh_used() {
		for (size_t i = 0; i < used.size(); i++) {
			std::fill(used[i].begin(), used[i].end(), false);
		}
		std::fill(regions.begin(), regions.end(), 0);
		max_region_id = 0;
	}
public:
	CountsScoring<Real_t>(VectorCellProvider<Real_t> const *cells, bool counts_score_constant_is_not_zero) : cells{ cells } {
		if (!counts_score_constant_is_not_zero) {
			return;
		}
		this->sum_counts = cells->get_sum_counts();
		this->squared_counts = cells->get_squared_counts();
		std::vector<bool> t;
		regions.resize(cells->getLociCount());
		std::fill(regions.begin(), regions.end(), 0);
		t.resize(squared_counts[0].size());
		std::fill(t.begin(), t.end(), false);
		for (size_t i = 0; i < squared_counts.size(); i++) {
			used.push_back(t);
		}
        for (size_t c = 0; c < cells->getCellsCount(); c++) {
            for (size_t b = 0; b < used[0].size(); b++) {
                all_bins_count += cells->get_counts_score_length_of_bin(b);
            }
        }
		regions_cache.resize(DEFAULT_CACHE_SIZE);
	}

	Real_t calculate_log_score(PointerTree &tree, std::vector<BreakpointPair> &attachment_vec) {
		if (COUNTS_SCORE_CONSTANT == 0.0) {
			return 0.0;
		}
		refresh_used();
        this->tree_size = 0.0;
		size_t cache_id = 0;
		std::map<BreakpointPair, std::set<size_t>> attachment;
		for (size_t cell = 0; cell < attachment_vec.size(); cell++) {
			if (attachment.find(attachment_vec[cell]) == attachment.end()) {
				attachment[attachment_vec[cell]] = std::set<size_t>();
			}
			attachment[attachment_vec[cell]].insert(cell);
		}
		Real_t result = 0.0;
		for (auto node : tree.root->children) {
                result += calculate_log_score(node, attachment, cache_id);
		}
		result +=  calculate_log_score_for_root_nodes();
		if (regions_cache.size() > DEFAULT_CACHE_SIZE) {
			regions_cache.resize(DEFAULT_CACHE_SIZE);
		}
		return  -result;
	}
};

template <class Real_t> bool CountsScoring<Real_t>::print = false;
#endif // !COUNTS_SCORING_H
