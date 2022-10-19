#ifndef LIK_CALC___H
#define LIK_CALC___H
#include <algorithm>
#include <map>
#include <numeric>
#include <thread>

#include "input_data/input_data.h"
#include "parameters/parameters.h"
#include "tree/attachment.h"
#include "tree/event_tree.h"
#include "utils/log_sum_accumulator.h"
#include "utils/matrix.h"

template <class Real_t> class LikelihoodMatrices {
public:
  std::vector<std::vector<Real_t>> breakpoint_likelihoods;
  std::vector<std::vector<Real_t>> no_breakpoint_likelihoods;

  LikelihoodMatrices<Real_t>(size_t bins, size_t cells) {
    this->breakpoint_likelihoods =
        Matrix::create_2d_matrix<Real_t>(bins, cells, 0.0);
    this->no_breakpoint_likelihoods =
        Matrix::create_2d_matrix<Real_t>(bins, cells, 0.0);
  }

  static void swap(LikelihoodMatrices<Real_t> &m1,
                   LikelihoodMatrices<Real_t> &m2) {
    std::swap(m1.breakpoint_likelihoods, m2.breakpoint_likelihoods);
    std::swap(m1.no_breakpoint_likelihoods, m2.no_breakpoint_likelihoods);
  }
};

template <class Real_t> class LikelihoodCalculatorState {
private:
  using NodeHandle = EventTree::NodeHandle;

public:
  std::vector<LogWeightAccumulator<Real_t>> likelihood_result;
  std::vector<Real_t> root_likelihoods;
  std::vector<Real_t> cell_to_max_attachment_likelihood;
  Attachment max_attachment;
  std::map<Event, Real_t> events_lengths;
  Real_t likelihood;

  LikelihoodCalculatorState<Real_t>(size_t cells_count)
      : max_attachment{get_root_label(), cells_count} {
    likelihood_result.resize(cells_count);
    root_likelihoods.resize(cells_count);
    cell_to_max_attachment_likelihood.resize(cells_count);
  }

  static void swap(LikelihoodCalculatorState<Real_t> &s1,
                   LikelihoodCalculatorState<Real_t> &s2) {
    std::swap(s1.likelihood_result, s2.likelihood_result);
    std::swap(s1.root_likelihoods, s2.root_likelihoods);
    std::swap(s1.cell_to_max_attachment_likelihood,
              s2.cell_to_max_attachment_likelihood);
    std::swap(s1.max_attachment, s2.max_attachment);
    std::swap(s1.events_lengths, s2.events_lengths);
    std::swap(s1.likelihood, s2.likelihood);
  }
};

template <class Real_t> class LikelihoodCalculator {
  EventTree &tree;
  LikelihoodCalculatorState<Real_t> &state;
  CONETInputData<Real_t> &cells;
  LikelihoodMatrices<Real_t> &likelihood_matrices;

  using NodeHandle = EventTree::NodeHandle;

  void get_normalized_event_lengths(NodeHandle node, Real_t length,
                                    Real_t depth) {
    state.events_lengths[tree.get_node_event(node)] = 0.0;
    if (node != tree.get_root()) {
      length += cells.get_event_length(node->label);
      state.events_lengths[tree.get_node_event(node)] = length / depth;
    }

    for (auto child : tree.get_children(node)) {
      get_normalized_event_lengths(child, length, depth + 1);
    }
  }

  void recalculate_event_lengths() {
    state.events_lengths.clear();
    get_normalized_event_lengths(tree.get_root(), 0.0, 0.0);
    Real_t sum = std::accumulate(
        state.events_lengths.begin(), state.events_lengths.end(), 0.0,
        [](const Real_t previous, decltype(*state.events_lengths.begin()) p) {
          return previous + std::exp(-p.second);
        });
    for (auto &el : state.events_lengths) {
      el.second = -el.second - std::log(sum);
    }
  }

  void calculate_root_likelihood() {
    for (size_t c = 0; c < cells.get_cells_count(); c++) {
      state.root_likelihoods[c] = 0.0;
    }
    for (size_t bin = 0;
         bin < likelihood_matrices.no_breakpoint_likelihoods.size(); bin++) {
      for (size_t c = 0; c < cells.get_cells_count(); c++) {
        state.root_likelihoods[c] +=
            likelihood_matrices.no_breakpoint_likelihoods[bin][c];
      }
    }
  }

  void extend_likelihood_to_node(NodeHandle node,
                                 std::vector<Real_t> &parent_likelihood) {
    auto breakpoints = tree.get_new_breakpoints(node);
    for (auto br : breakpoints) {
      for (size_t c = 0; c < cells.get_cells_count(); c++) {
        parent_likelihood[c] +=
            likelihood_matrices.breakpoint_likelihoods[br][c] -
            likelihood_matrices.no_breakpoint_likelihoods[br][c];
      }
    }
  }

  void
  reverse_likelihood_node_extension(NodeHandle node,
                                    std::vector<Real_t> &parent_likelihood) {
    auto breakpoints = tree.get_new_breakpoints(node);
    for (auto br : breakpoints) {
      for (size_t c = 0; c < cells.get_cells_count(); c++) {
        parent_likelihood[c] +=
            -likelihood_matrices.breakpoint_likelihoods[br][c] +
            likelihood_matrices.no_breakpoint_likelihoods[br][c];
      }
    }
  }

  void tree_depth_first_likelihood_calculation(
      NodeHandle node, std::vector<Real_t> &parent_likelihood) {
    extend_likelihood_to_node(node, parent_likelihood);

    for (size_t c = 0; c < cells.get_cells_count(); c++) {
      if (USE_EVENT_LENGTHS_IN_ATTACHMENT) {
        state.likelihood_result[c].add(
            parent_likelihood[c] +
            state.events_lengths[tree.get_node_event(node)]);
      } else {
        state.likelihood_result[c].add(parent_likelihood[c]);
      }
      if (state.cell_to_max_attachment_likelihood[c] < parent_likelihood[c]) {
        state.cell_to_max_attachment_likelihood[c] = parent_likelihood[c];
        state.max_attachment.set_attachment(c, tree.get_node_label(node));
      }
    }
    for (auto ch : tree.get_children(node)) {
      tree_depth_first_likelihood_calculation(ch, parent_likelihood);
    }
    reverse_likelihood_node_extension(node, parent_likelihood);
  }

  Real_t sum_cell_likelihoods() {
    Real_t result_ = 0.0;
    if (!USE_EVENT_LENGTHS_IN_ATTACHMENT) {
      std::for_each(state.likelihood_result.rbegin(),
                    state.likelihood_result.rend(),
                    [&](LogWeightAccumulator<Real_t> &acc) {
                      result_ += acc.get_result() -
                                 std::log((Real_t)(tree.get_size() - 1));
                    });
    } else {
      std::for_each(state.likelihood_result.rbegin(),
                    state.likelihood_result.rend(),
                    [&](LogWeightAccumulator<Real_t> &acc) {
                      result_ += acc.get_result();
                    });
    }
    return result_;
  }

public:
  LikelihoodCalculator<Real_t>(EventTree &tree,
                               LikelihoodCalculatorState<Real_t> &state,
                               CONETInputData<Real_t> &cells,
                               LikelihoodMatrices<Real_t> &matrices)
      : tree{tree}, state{state}, cells{cells}, likelihood_matrices{matrices} {}

  Real_t calculate_likelihood() {
    calculate_root_likelihood();

    for (size_t c = 0; c < cells.get_cells_count(); c++) {
      state.likelihood_result[c].clear();
      state.max_attachment.set_attachment(c, get_root_label());
      state.cell_to_max_attachment_likelihood[c] = state.root_likelihoods[c];
    }

    if (USE_EVENT_LENGTHS_IN_ATTACHMENT) {
      recalculate_event_lengths();
    }

    for (auto &node : tree.get_children(tree.get_root())) {
      tree_depth_first_likelihood_calculation(node, state.root_likelihoods);
    }
    state.likelihood = sum_cell_likelihoods();
    return state.likelihood;
  }
};

#endif