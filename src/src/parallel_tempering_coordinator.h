#ifndef PARALLEL_TEMPERING_COORDINATOR_H
#define PARALLEL_TEMPERING_COORDINATOR_H

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>

#include "./adaptive_pt.h"
#include "conet_result.h"
#include "input_data/input_data.h"
#include "likelihood/EM_estimator.h"
#include "likelihood/gaussian_mixture.h"
#include "likelihood/likelihood_data.h"
#include "likelihood_coordinator.h"
#include "moves/move_type.h"
#include "parameters/parameters.h"
#include "tree/tree_formatter.h"
#include "tree/tree_sampler.h"
#include "tree_sampler_coordinator.h"
#include "utils/logger/logger.h"
#include "utils/random.h"
#include "utils/utils.h"

template <class Real_t> class ParallelTemperingCoordinator {
  AdaptivePT<Real_t> adaptive_pt;
  std::vector<Real_t> temperatures;
  CONETInputData<Real_t> &provider;
  Random<Real_t> &random;
  std::vector<EventTree> trees;
  std::vector<std::unique_ptr<TreeSamplerCoordinator<Real_t>>>
      tree_sampling_coordinators;
  std::vector<std::unique_ptr<LikelihoodCoordinator<Real_t>>>
      likelihood_calculators;

  const std::map<MoveType, double> move_probabilities = {
      {DELETE_LEAF, 100.0},       {ADD_LEAF, 30.0},     {PRUNE_REATTACH, 30.0},
      {SWAP_LABELS, 30.0},        {CHANGE_LABEL, 30.0}, {SWAP_SUBTREES, 30.0},
      {SWAP_ONE_BREAKPOINT, 30.0}};

  const size_t INIT_TREE_SIZE = 2;
  const Real_t MIN_COMPONENT_WEIGHT = 0.01;

  EventTree sample_starting_tree_for_chain() {
    log("Sampling initial tree for chain with size ", INIT_TREE_SIZE);
    VertexLabelSampler<Real_t> vertexSet{provider.get_loci_count() - 1,
                                         provider.get_chromosome_end_markers()};
    return sample_tree<Real_t>(INIT_TREE_SIZE, vertexSet, random);
  }

  void prepare_sampling_services(LikelihoodData<Real_t> likelihood) {
    log("Starting preparation of sampling services with ", NUM_REPLICAS, " replicas...");
    for (size_t i = 0; i < NUM_REPLICAS; i++) {
      trees.push_back(sample_starting_tree_for_chain());
    }
    for (size_t i = 0; i < NUM_REPLICAS; i++) {
      likelihood_calculators.push_back(
          std::move(std::make_unique<LikelihoodCoordinator<Real_t>>(
              likelihood, trees[i], provider, random.next_int())));
      tree_sampling_coordinators.push_back(
          std::move(std::make_unique<TreeSamplerCoordinator<Real_t>>(
              trees[i], *likelihood_calculators[i], random.next_int(), provider,
              move_probabilities)));
    }
    log("PID 0 replica will start with temperature ", 1.0);
    temperatures.push_back(1.0);
    for (size_t i = 1; i < NUM_REPLICAS; i++) {
      log("PID ", i, " replica will start with temperature ", temperatures.back() * 0.1);
      temperatures.push_back(temperatures.back() * 0.1);
    }
    for (size_t i = 0; i < NUM_REPLICAS; i++) {
      tree_sampling_coordinators[i]->set_temperature(temperatures[i]);
    }
  }

  LikelihoodData<Real_t>
  estimate_likelihood_parameters(LikelihoodData<Real_t> likelihood,
                                 const size_t iterations) {
    log("Starting parameter MCMC estimation...");
    EventTree tree = sample_starting_tree_for_chain();
    LikelihoodCoordinator<Real_t> calc(likelihood, tree, provider,
                                       random.next_int());
    TreeSamplerCoordinator<Real_t> coordinator(tree, calc, random.next_int(),
                                               provider, move_probabilities);

    for (size_t i = 0; i < iterations; i++) {
      if (i % PARAMETER_RESAMPLING_FREQUENCY == 0) {
        calc.resample_likelihood_parameters(
            coordinator.get_log_tree_prior(),
            coordinator.get_current_count_dispersion_penalty());
      }
      coordinator.execute_metropolis_hastings_step();
    }
    log("Finished parameter estimation");
    auto map_parameters =
        calc.get_map_parameters().remove_components_with_small_weight(
            MIN_COMPONENT_WEIGHT);
    log("Estimated breakpoint distribution: ",
        map_parameters.brkp_likelihood.to_string());
    log("Estimated no-breakpoint distribution: ",
        map_parameters.no_brkp_likelihood.to_string());
    return map_parameters;
  }

  void mcmc_simulation(size_t iterations) {
    for (size_t i = 0; i < iterations / NUMBER_OF_MOVES_BETWEEN_SWAPS; i++) {
      std::vector<std::thread> threads;
      for (size_t th = 0; th < NUM_REPLICAS; th++) {
        threads.emplace_back([this, th] {
          for (size_t i = 0; i < (size_t)NUMBER_OF_MOVES_BETWEEN_SWAPS; i++)
            this->tree_sampling_coordinators[th]
                ->execute_metropolis_hastings_step();
        });
      }
      for (auto &th : threads) {
        th.join();
      }
      swap_step();

      if (VERBOSE && i % 1000 == 0) {
        log("State after ", i * NUMBER_OF_MOVES_BETWEEN_SWAPS, " iterations:");
        log("Tree size: ",
            this->tree_sampling_coordinators[0]->get_tree_size());
        log("Log-likelihood: ",
            this->likelihood_calculators[0]->get_likelihood());
        log("Log-likelihood with penalty: ",
            this->tree_sampling_coordinators[0]->get_total_likelihood());
      }
    }
  }

  void swap_step() {
    if (NUM_REPLICAS == 1) {
      return;
    }
    std::vector<Real_t> states;
    for (size_t i = 0; i < likelihood_calculators.size(); i++) {
      states.push_back(likelihood_calculators[i]->get_likelihood());
    }
    adaptive_pt.update(states);
    this->temperatures = adaptive_pt.get_temperatures();

    for (size_t i = 0; i < likelihood_calculators.size(); i++) {
      tree_sampling_coordinators[i]->set_temperature(temperatures[i]);
    }
    int pid = random.next_int(NUM_REPLICAS - 1);
    auto likelihood_left = tree_sampling_coordinators[pid]
                               ->get_likelihood_without_priors_and_penalty();
    auto likelihood_right = tree_sampling_coordinators[pid + 1]
                                ->get_likelihood_without_priors_and_penalty();
    Real_t swap_acceptance_ratio = (temperatures[pid] - temperatures[pid + 1]) *
                                   (likelihood_right - likelihood_left);
    if (random.log_uniform() <= swap_acceptance_ratio) {
      tree_sampling_coordinators[pid]->set_temperature(temperatures[pid + 1]);
      tree_sampling_coordinators[pid + 1]->set_temperature(temperatures[pid]);
      std::swap(tree_sampling_coordinators[pid],
                tree_sampling_coordinators[pid + 1]);
      std::swap(likelihood_calculators[pid], likelihood_calculators[pid + 1]);
    }
  }

  LikelihoodData<Real_t> prepare_initial_likelihood_parameters() {
    log("Initializing EM estimator...");
    Gauss::EMEstimator<Real_t> EM(
        Utils::flatten<Real_t>(provider.get_corrected_counts()), random);
    log("Starting EM estimation of mixture with ", MIXTURE_SIZE, " components");
    auto result = EM.estimate(MIXTURE_SIZE)
        .remove_components_with_small_weight(MIN_COMPONENT_WEIGHT);
    log("Finished EM estimation");
    return result;
  }

  CONETInferenceResult<Real_t> choose_best_tree_among_replicas() {
    Utils::MaxValueAccumulator<CONETInferenceResult<Real_t>, Real_t> best_tree;
    for (auto &replica : tree_sampling_coordinators) {
      best_tree.update(replica->get_inferred_tree(),
                       replica->get_inferred_tree().likelihood);
    }
    return best_tree.get();
  }

public:
  ParallelTemperingCoordinator(CONETInputData<Real_t> &provider,
                               Random<Real_t> &random)
      : adaptive_pt{NUM_REPLICAS}, provider{provider}, random{random} {}

  CONETInferenceResult<Real_t> simulate(size_t iterations_parameters,
                                        size_t iterations_pt) {
    prepare_sampling_services(estimate_likelihood_parameters(
        prepare_initial_likelihood_parameters(), iterations_parameters));
    mcmc_simulation(iterations_pt);
    return choose_best_tree_among_replicas();
  }
};
#endif // PARALLEL_TEMPERING_COORDINATOR_H
