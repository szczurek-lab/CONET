#ifndef TREE_SAMPLER_COORDINATOR_H
#define TREE_SAMPLER_COORDINATOR_H
#include <algorithm>
#include <tuple>
#include <vector>

#include "conet_result.h"
#include "tree/attachment.h"
#include "tree/tree_counts_scoring.h"
#include "tree/tree_formatter.h"
#include "tree_mh_steps_executor.h"
#include "utils/logger/logger.h"
#include "utils/utils.h"
/**
 * Class coordinating all components of MH sampling.
 *
 * Tree MH sampling is done in few steps:\n
 * 1) First move type is sampled. \n
 * 2) TreeNodeSampler confirms that sampled move type is possible (for instance
 * <code> DELETE_LEAF</code> move may not be possible \n if the tree has zero
 * size)\n 3) Any required random elements are sampled (for example labels for a
 * new node).\n 4) VertexLabelSampler, EventTree, TreeNodeSampler are called to
 * do their move specific transformations\n 5) If move has not been accepted a
 * move rollback is carried out\n
 */
template <class Real_t> class TreeSamplerCoordinator {
  using NodeHandle = EventTree::NodeHandle;
  using MoveData = typename MHStepsExecutor<Real_t>::MoveData;
  EventTree &tree;
  LikelihoodCoordinator<Real_t> &likelihood_coordinator;
  CountsDispersionPenalty<Real_t> dispersion_penalty_calculator;
  Random<Real_t> random;
  std::map<MoveType, Real_t> move_probabilities;
  Real_t temperature{1.0};
  Real_t tree_count_dispersion_penalty{1.0}; // penalty value of a current tree
  Utils::MaxValueAccumulator<CONETInferenceResult<Real_t>, Real_t>
      best_found_tree;
  MHStepsExecutor<Real_t> mh_step_executor;

  Real_t get_probability_of_reverse_move(MoveType type) {
    switch (type) {
    case ADD_LEAF:
      return move_probabilities[DELETE_LEAF];
    case DELETE_LEAF:
      return move_probabilities[ADD_LEAF];
    case CHANGE_LABEL:
      return move_probabilities[CHANGE_LABEL];
    case SWAP_LABELS:
      return move_probabilities[SWAP_LABELS];
    case PRUNE_REATTACH:
      return move_probabilities[PRUNE_REATTACH];
    case SWAP_SUBTREES:
      return move_probabilities[SWAP_SUBTREES];
    case SWAP_ONE_BREAKPOINT:
    default:
      return move_probabilities[SWAP_ONE_BREAKPOINT];
    }
  }

  void move(MoveType type) {
    recalculate_counts_dispersion_penalty();
    auto before_move_likelihood =
        temperature * likelihood_coordinator.get_likelihood() +
        mh_step_executor.get_log_tree_prior() + tree_count_dispersion_penalty;

    auto move_data = mh_step_executor.execute_move(type);

    auto after_move_likelihood =
        temperature * likelihood_coordinator.calculate_likelihood() +
        mh_step_executor.get_log_tree_prior();
    auto after_move_counts_dispersion_penalty =
        dispersion_penalty_calculator.calculate_log_score(
            tree, likelihood_coordinator.calculate_max_attachment());
    after_move_likelihood += after_move_counts_dispersion_penalty;

    Real_t log_acceptance = after_move_likelihood - before_move_likelihood +
                            move_data.reverse_move_log_kernel -
                            move_data.move_log_kernel +
                            std::log(move_probabilities[type]) -
                            std::log(get_probability_of_reverse_move(type));

    log_debug("Log acceptance ratio: ", log_acceptance, " likelihood before ",
              before_move_likelihood, " likelihood after ",
              after_move_likelihood);

    if (random.log_uniform() <= log_acceptance) {
      likelihood_coordinator.persist_likelihood_calculation_result();
      tree_count_dispersion_penalty = after_move_counts_dispersion_penalty;
      log_debug("Move accepted");
    } else {
      mh_step_executor.rollback_move(type, move_data);
      log_debug("Move rejected");
    }
  }

  MoveType sample_move_type() {
    std::map<size_t, MoveType> type_to_index;
    std::vector<Real_t> weights;
    for (auto x : move_probabilities) {
      type_to_index[weights.size()] = x.first;
      weights.push_back(x.second);
    }
    return type_to_index[random.discrete(weights)];
  }

  void recalculate_counts_dispersion_penalty() {
    tree_count_dispersion_penalty =
        dispersion_penalty_calculator.calculate_log_score(
            tree, likelihood_coordinator.get_max_attachment());
  }

public:
  TreeSamplerCoordinator(EventTree &tree, LikelihoodCoordinator<Real_t> &lC,
                         unsigned int seed, CONETInputData<Real_t> &cells,
                         std::map<MoveType, Real_t> move_probabilities)
      : tree{tree}, likelihood_coordinator{lC},
        dispersion_penalty_calculator{cells}, random{seed},
        move_probabilities{move_probabilities}, mh_step_executor{tree, cells,
                                                                 random} {}

  Real_t get_likelihood_without_priors_and_penalty() {
    return likelihood_coordinator.get_likelihood();
  }

  void set_temperature(Real_t temperature) { this->temperature = temperature; }

  Real_t get_log_tree_prior() { return mh_step_executor.get_log_tree_prior(); }

  Real_t get_total_likelihood() {
    return likelihood_coordinator.get_likelihood() +
           mh_step_executor.get_log_tree_prior() +
           tree_count_dispersion_penalty;
  }

  Real_t get_temperature() const { return this->temperature; }

  size_t get_tree_size() const { return tree.get_size(); }
  Real_t get_current_count_dispersion_penalty() const {
    return tree_count_dispersion_penalty;
  }

  void execute_metropolis_hastings_step() {
    MoveType type = sample_move_type();
    log_debug("Sampled move of type: ", move_type_to_string(type));

    if (mh_step_executor.move_is_possible(type)) {
      move(type);
    }

    auto l = get_total_likelihood();
    best_found_tree.update(
        CONETInferenceResult<Real_t>(
            tree, likelihood_coordinator.get_max_attachment(), l),
        l);
  }

  CONETInferenceResult<Real_t> get_inferred_tree() {
    return best_found_tree.get();
  }
};

#endif // !TREE_SAMPLER_COORDINATOR_H
