#ifndef LIK_CALC___H
#define LIK_CALC___H
#include <algorithm>
#include <numeric>
#include <thread>  
#include <map>

#include "utils/log_sum_accumulator.h"
#include "cell_provider/vector_cell_provider.h"
#include "parameters/parameters.h"
#include "tree/event_tree.h"

template<class Real_t> class LikelihoodCalculatorState {
    private:
        void get_normalized_event_lengths(EventTree* tree, NodeHandle node, Real_t length, std::map<Event, Real_t> &result, Real_t depth)
        {
            if (node != tree->get_root()) {
                length += cells->get_event_length(node->label);
                result[tree->get_node_event(node)] = length / depth;
            } else {
                result[tree->get_node_event(node)] = 0.0;
            } 
            
            for (auto child : tree->get_children(node)) {
                get_normalized_event_lengths(child, length, result, depth + 1);
            }
        }   


    public:
        std::vector<std::vector<Real_t>> breakpoint_likelihoods;
        std::vector<std::vector<Real_t>> no_breakpoint_likelihoods;
        std::vector<LogWeightAccumulator<Real_t>> likelihood_result;
        std::vector<Real_t> root_likelihoods;
        std::vector<Real_t> cell_to_max_attachment_likelihood;
        std::vector<TreeLabel> max_attachment;
        std::map<Event, Real_t> events_lengths;

        LikelihoodCalculatorState<Real>(size_t cells_count, EventTree& tree) {
            likelihood_result.resize(cells_count);
		    root_likelihoods.resize(cells_count);
            cell_to_max_attachment_likelihood.resize(cells_count);
            for (size_t c = 0; c < cells_count; c++) {
                max_attachment.push_back(get_root_label());
            }
        }

        void recalculate_event_lengths(EventTree *tree) {
            get_normalized_event_lengths(tree, tree->get_root(), 0.0, events_lengths, 0.0);
			Real_t sum = std::accumulate(events_lengths.begin(), events_lengths.end(), 0.0, [](const Real_t previous, decltype(*events_lengths.begin()) p) { return previous + std::exp(-p.second); });
			for (auto &el : events_lengths)
			{
				el.second = -el.second - std::log(sum);
			}
        }
};

template <class Real_t> class LikelihoodCalculatorBase {
private: 
    EventTree *tree;
    LikelihoodCalculatorState &state;

    using NodeHandle = EventTree::NodeHandle;

    void calculate_root_likelihood(size_t left, size_t right) {
        for (size_t c = left; c < right; c++) {
            state.root_likelihoods[c] = 0.0;
        }
        for (size_t bin = 0; bin < state.no_breakpoint_likelihoods.size(); bin++) {
            for (size_t c = left; c < right; c++) {
                state.root_likelihoods[c] += state.no_breakpoint_likelihoods[bin][c];
            }
        }
    }

    void extend_likelihood_to_node(NodeHandle node, std::vector<Real_t> &parent_likelihood, size_t left, size_t right) {
        auto breakpoints = state.tree.get_new_breakpoints(node);
        for (auto br : breakpoints) {
            for (size_t c = left; c < right; c++) {
                parent_likelihood[c] += state.breakpoint_likelihoods[br][c] - state.no_breakpoint_likelihoods[br][c];
            }
        }
    }

    void reverse_likelihood_node_extension(NodeHandle node, std::vector<Real_t> &parent_likelihood, size_t left, size_t right) {
        auto breakpoints = state.tree.get_new_breakpoints(node);
        for (auto br : breakpoints) {
            for (size_t c = left; c < right; c++) {
                parent_likelihood[c] += -state.breakpoint_likelihoods[br][c] + state.no_breakpoint_likelihoods[br][c];
            }
        }
    }

    void treeDFS(NodeHandle node, std::vector<Real_t> &parent_likelihood, size_t left, size_t right) {
        extend_likelihood_to_node(node, parent_likelihood, left, right);

        for (size_t c = left; c < right; c++) {
            if (USE_EVENT_LENGTHS_IN_ATTACHMENT)
            {
                state.likelihood_result[c].add(parent_likelihood[c]  + events_lengths[tree->get_node_event(node)]);	
            }
            else {
                state.likelihood_result[c].add(parent_likelihood[c]);
            }
            if (state.cell_to_max_attachment_likelihood[c] < parent_likelihood[c])
            {
                state.cell_to_max_attachment_likelihood[c] = parent_likelihood[c];
                state.max_attachment[c] = state.tree.get_node_label(node);
            }
        }
        for (auto ch : tree->get_children(node)) {
            treeDFS(ch, parent_likelihood, left, right, max_likelihoods, events_lengths);
        }
        reverse_likelihood_node_extension(node, parent_likelihood, left, right);
    }

public:
    void calculate_likelihood(size_t left, size_t right) {
		
        for (size_t s = left; s < right; s++) {
			state.likelihood_result[s].clear();
		}
        
        if (USE_EVENT_LENGTHS_IN_ATTACHMENT) {
            state.recalculate_event_lengths(tree);
        }
        
        calculate_root_likelihood(left, right);

        for (size_t c = left; c < right; c++) {
			state.max_attachment[c] = get_root_label();
            cell_to_max_attachment_likelihood[c] = state.root_likelihoods[c];
		}
		
		for (auto &node : tree->get_children(tree->get_root())) {
			treeDFS(node, state.root_likelihood, left, right);
		}
	}
};

#endif 