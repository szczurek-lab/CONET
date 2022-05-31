#ifndef PARALLEL_TEMPERING_COORDINATOR_H
#define PARALLEL_TEMPERING_COORDINATOR_H


#include <vector>
#include <atomic>
#include <condition_variable>
#include <thread>
#include <mutex>
#include <utility>

#include "cell_provider/vector_cell_provider.h"
#include "likelihood/implementations/gaussian_mixture.h"
#include "moves/move_type.h"
#include "tree/tree_sampler.h"
#include "tree/tree_formatter.h"
#include "tree_sampler_coordinator.h"
#include "likelihood/implementations/normal_mixture_likelihood.h"
#include "likelihood_coordinator.h"
#include "utils/random.h"
#include "EM_estimator.h"
#include "utils/logger/logger.h"
#include "./adaptive_pt.h"
#include "parameters/parameters.h"
#include "conet_result.h"
#include "utils/utils.h"

template <class Real_t> class ParallelTemperingCoordinator {
public:
	AdaptivePT<Real_t> adaptive_pt;
	std::vector<Real_t> temperatures;
	VectorCellProvider<Real_t> &provider;
	Random<Real_t> &random;
	std::map<MoveType, double> moveProbability = {
			{DELETE_LEAF, 100.0},
			{ADD_LEAF, 30.0},
			{PRUNE_REATTACH, 30.0},
			{SWAP_LABELS, 30.0},
			{CHANGE_LABEL, 30.0},
			{SWAP_SUBTREES, 30.0},
			{SWAP_ONE_BREAKPOINT, 30.0}
	};

	std::vector<EventTree> trees;
	std::vector<TreeSamplerCoordinator<Real_t>*> tree_sampling_coordinators;
	std::vector<LikelihoodCoordinator<Real_t>*> likelihood_calculators;
	
	EventTree sample_starting_tree_for_chain()
	{
		VertexLabelSampler<double> vertexSet {provider.get_loci_count() - 1, provider.getChromosomeMarkers()};
		return sample_tree<double>(5, vertexSet, random);
	}


	void prepare_sampling_services(LikelihoodData<Real_t> likelihood)
	{
		trees.reserve(THREADS_NUM);
		likelihood_calculators.reserve(THREADS_NUM);
		tree_sampling_coordinators.reserve(THREADS_NUM);
		for (size_t i = 0; i < THREADS_NUM; i++) {
			trees.push_back(sample_starting_tree_for_chain());
		}
		for (size_t i = 0; i < THREADS_NUM; i++) {
			likelihood_calculators.emplace_back(new LikelihoodCoordinator<Real_t>(likelihood, trees[i], &provider, provider.get_loci_count() - 1, random.nextInt()));
		}
		for (size_t i = 0; i < THREADS_NUM; i++) {
			tree_sampling_coordinators.emplace_back(new TreeSamplerCoordinator<Real_t>(trees[i], likelihood_calculators[i], random.nextInt() , provider.get_loci_count() - 1, &provider, moveProbability, i));
		}
		temperatures.push_back(1.0);
		for (size_t i = 1; i < THREADS_NUM; i++) {
			temperatures.push_back(temperatures.back()*0.1);
		}
		for (size_t i = 0; i < THREADS_NUM; i++) {
			tree_sampling_coordinators[i]->set_temperature(temperatures[i]);
		}
	}

	
	LikelihoodData<Real_t> estimateParameters(LikelihoodData<Real_t> likelihood, const size_t iterations)
	{
		log("Starting parameter MCMC estimation...");
		EventTree tree = sample_starting_tree_for_chain();
		LikelihoodCoordinator<Real_t> calc(likelihood, tree, &provider, provider.get_loci_count() - 1, random.nextInt());
		TreeSamplerCoordinator<Real_t> coordinator(tree, &calc, random.nextInt(), provider.get_loci_count() - 1, &provider, moveProbability, false);

		for (size_t i = 0; i < iterations; i++) {
			if (i % PARAMETER_RESAMPLING_FREQUENCY == 0) {
				calc.resampleParameters(coordinator.get_log_tree_prior(), coordinator.tree_count_dispersion_penalty);
			}
			coordinator.execute_metropolis_hastings_step(); 
		}
		log("Finished parameter estimation");
        calc.MAP_parameters.brkp_likelihood.remove_components_with_small_weight(0.01);
		log("Estimated breakpoint distribution:");
        log(calc.MAP_parameters.brkp_likelihood.to_string());
		log("Estimated no-breakpoint distribution:");
        log(calc.MAP_parameters.no_brkp_likelihood.to_string());

		return calc.MAP_parameters;
	}

	void mcmc_simulation(size_t iterations)
	{
		for (size_t i = 0; i < iterations / NUMBER_OF_MOVES_BETWEEN_SWAPS; i++)
		{
			std::vector<std::thread> threads;
			for (size_t th = 0; th < THREADS_NUM; th++) {
				threads.emplace_back([this, th] { for (size_t i = 0; i < (size_t) NUMBER_OF_MOVES_BETWEEN_SWAPS; i++) this->tree_sampling_coordinators[th]->execute_metropolis_hastings_step(); });
			}
			for (auto &th : threads)
			{
				th.join();
			}
			if (THREADS_NUM != 1) swap_step();
            
            if (VERBOSE && i % 1000 == 0) {
                log("State after " , i*NUMBER_OF_MOVES_BETWEEN_SWAPS, " iterations:");
                log("Tree size: ", this->tree_sampling_coordinators[0]->tree.get_size());
                log("Log-likelihood: ", this->likelihood_calculators[0]->get_likelihood());
                log("Log-likelihood with penalty: ",  this->tree_sampling_coordinators[0]->get_total_likelihood());
                
            }
		}
	}

	void swap_step()
	{
		std::vector<Real_t> states;
		for (size_t i = 0; i < likelihood_calculators.size(); i++) {
			states.push_back(likelihood_calculators[i]->get_likelihood());
		}
		adaptive_pt.update(states);
		this->temperatures = adaptive_pt.get_temperatures();

		for (size_t i = 0; i < likelihood_calculators.size(); i ++) {
			tree_sampling_coordinators[i]->set_temperature(temperatures[i]);
		}
		int pid = random.nextInt(THREADS_NUM - 1);
		auto likelihood_left = tree_sampling_coordinators[pid]->get_likelihood_without_priors_and_penalty();
		auto likelihood_right = tree_sampling_coordinators[pid + 1]->get_likelihood_without_priors_and_penalty();
		Real_t swap_acceptance_ratio = (temperatures[pid] - temperatures[pid + 1]) * (likelihood_right - likelihood_left);
		if (random.logUniform() <= swap_acceptance_ratio)
		{
			tree_sampling_coordinators[pid]->set_temperature(temperatures[pid + 1]);
			tree_sampling_coordinators[pid + 1]->set_temperature(temperatures[pid]);
			std::swap(tree_sampling_coordinators[pid], tree_sampling_coordinators[pid + 1]);
			std::swap(likelihood_calculators[pid], likelihood_calculators[pid + 1]);
		}
	}

	LikelihoodData<Real_t> prepare_starting_parameters()
	{

		std::vector<Real_t> means;
		std::vector<Real_t> variances;
		std::vector<Real_t> weights;
		for (int i = 0; i < (int) MIXTURE_SIZE; i++)
		{
			means.push_back(-i);
			variances.push_back(random.uniform());
			weights.push_back(random.uniform());
		}
		Real_t sum = std::accumulate(weights.begin(), weights.end(), 0.0);
		for (size_t i = 0; i < MIXTURE_SIZE; i++) {
			weights[i] = weights[i] / sum;
		}
		EMEstimator < Real_t > EM(Utils::flatten<Real_t>(provider.get_corrected_counts()));
		EM.estimate(means, variances, weights);

		Gauss::Gaussian<Real_t> gaussian(0.0, std::sqrt(variances[0]), random);
		weights[0] = 0;
		sum = std::accumulate(weights.begin(), weights.end(), 0.0);
		for (size_t i = 0; i < MIXTURE_SIZE; i++) {
			weights[i] = weights[i] / sum;
			variances[i] = std::sqrt(variances[i]);
		}
		variances.erase(variances.begin());
		means.erase(means.begin());
		weights.erase(weights.begin());
		Gauss::GaussianMixture<Real_t> mixture(weights, means, variances, random);
		mixture.remove_components_with_small_weight(0.01);
		return LikelihoodData<Real_t>(gaussian, mixture);
	}

	CONETInferenceResult<Real_t> choose_best_tree() {
		auto best = tree_sampling_coordinators[0]->get_inferred_tree();
		for (auto replica : tree_sampling_coordinators) {
			if (replica->get_inferred_tree().likelihood > best.likelihood) {
				best = replica->get_inferred_tree();
			}
		}
		return best;
	}
public:
	ParallelTemperingCoordinator(VectorCellProvider<Real_t> &provider, Random<Real_t> &random): 
		adaptive_pt{ THREADS_NUM }, 
		provider{ provider }, 
		random{ random }
	{}


CONETInferenceResult<Real_t> simulate(size_t iterations_parameters, size_t iterations_pt)
	{
		auto parameters_MAP = estimateParameters(prepare_starting_parameters(), iterations_parameters);
		random.nextInt();// mozesz usuanc pozniej
		prepare_sampling_services(parameters_MAP);
		mcmc_simulation(iterations_pt);
		return choose_best_tree();
	}
	
};

#endif //PARALLEL_TEMPERING_COORDINATOR_H
