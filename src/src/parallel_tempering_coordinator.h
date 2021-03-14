#ifndef PARALLEL_TEMPERING_COORDINATOR_H
#define PARALLEL_TEMPERING_COORDINATOR_H


#include <vector>
#include <atomic>
#include <condition_variable>
#include <thread>
#include <mutex>
#include <utility>

#include "cell_provider/vector_cell_provider.h"
#include "likelihood/implementations/utils/gaussian_mixture.h"
#include "moves/move_type.h"
#include "utils/tree_sampler.h"
#include "tree_sampler_coordinator.h"
#include "likelihood/implementations/normal_mixture_likelihood.h"
#include "likelihood_calculator.h"
#include "utils/random.h"
#include "EM_estimator.h"
#include "utils/logger/logger.h"
#include "utils/adaptive_pt.h"
#include "parameters/parameters.h"

template <class Real_t> class ParallelTemperingCoordinator {
private:
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

	std::vector<PointerTree> trees;
	std::vector<TreeSamplerCoordinator<Real_t>*> tree_sampling_coordinators;
	std::vector<LikelihoodCalculator<Real_t>*> likelihood_calculators;

	std::condition_variable wait_for_mcmc_step;

	std::atomic<bool> finished_swap_step{false };
	std::mutex m_mutex;
	std::atomic<int> finished_mcmc_step{ 0 };
	
	PointerTree sample_tree(unsigned int seed)
	{
		PointerTree tree;
		VertexSet<double> vertexSet {provider.getLociCount() - 1, provider.getChromosomeMarkers(), seed};
		sampleTree<double>(tree, &vertexSet, seed);
		tree.init();
		return tree;
	}


	void prepare_sampling_services(NormalMixtureLikelihood<Real_t> *likelihood)
	{
		trees.reserve(THREADS_NUM);
		likelihood_calculators.reserve(THREADS_NUM);
		tree_sampling_coordinators.reserve(THREADS_NUM);
		for (size_t i = 0; i < THREADS_NUM; i++)
		{
			PointerTree tree = sample_tree(random.nextInt());
			trees.push_back(tree);
		}
		for (size_t i = 0; i < THREADS_NUM; i++)
		{
			trees[i].init();
			likelihood_calculators.emplace_back(new LikelihoodCalculator<Real_t>(likelihood, &trees[i], &provider, provider.getLociCount() - 1, random.nextInt()));
		}
		for (size_t i = 0; i < THREADS_NUM; i++)
		{
			tree_sampling_coordinators.emplace_back(new TreeSamplerCoordinator<Real_t>(&trees[i], likelihood_calculators[i], random.nextInt() , provider.getLociCount() - 1, &provider, moveProbability, i));
		}
		temperatures.push_back(1.0);
		for (size_t i = 1; i < THREADS_NUM; i++)
		{
			temperatures.push_back(temperatures.back()*0.1);
		}
		for (size_t i = 0; i < THREADS_NUM; i++)
		{
			tree_sampling_coordinators[i]->set_temperature(temperatures[i]);
		}
	}

	
	std::tuple<Gauss::GaussianMixture<Real_t>, Gauss::Gaussian<Real_t>, Real_t> estimateParameters(Gauss::GaussianMixture<Real_t> mixture, Gauss::Gaussian<Real_t> no_breakpoint, const size_t iterations)
	{
		log("Starting parameter estimation");
		PointerTree tree = sample_tree(random.nextInt());
		NormalMixtureLikelihood<Real_t> likelihood(no_breakpoint, mixture, random.nextInt());
		LikelihoodCalculator<Real_t> calc(&likelihood, &tree, &provider, provider.getLociCount() - 1, random.nextInt());
		TreeSamplerCoordinator<Real_t> coordinator(&tree, &calc, random.nextInt(), provider.getLociCount() - 1, &provider, moveProbability, false);

		for (size_t i = 0; i < iterations; i++)
		{
			if (i % PARAMETER_RESAMPLING_FREQUENCY == 0)
			{
				calc.resampleParameters(coordinator.getLogTreePrior(), coordinator.tree_count_score);
				coordinator.recount_score();
			}
			coordinator.treeMHStep(); 
		}
		log("Finished parameter estimation");
		auto map = likelihood.get_MAP();
		return std::make_tuple(map.first, map.second, calc.MAP);
	}

	void mcmc_simulation(size_t iterations)
	{
		for (size_t i = 0; i < iterations / NUMBER_OF_MOVES_BETWEEN_SWAPS; i++)
		{
			std::vector<std::thread> threads;
			for (size_t th = 0; th < THREADS_NUM; th++) {
				threads.emplace_back([this, th] { for (size_t i = 0; i < (size_t) NUMBER_OF_MOVES_BETWEEN_SWAPS; i++) this->tree_sampling_coordinators[th]->treeMHStep(); });
			}
			for (auto &th : threads)
			{
				th.join();
			}
			if (THREADS_NUM != 1) swap_step();
		}
	}


	void swap_step()
	{
		std::vector<Real_t> states;
		for (size_t i = 0; i < likelihood_calculators.size(); i++)
		{
			states.push_back(likelihood_calculators[i]->getLikelihood());
		}
		adaptive_pt.update(states);
		this->temperatures = adaptive_pt.get_temperatures();

		for (size_t i = 0; i < likelihood_calculators.size(); i ++)
		{
			tree_sampling_coordinators[i]->set_temperature(temperatures[i]);
		}
		int pid = random.nextInt(THREADS_NUM - 1);
		auto likelihood_left = tree_sampling_coordinators[pid]->get_likelihood();
		auto likelihood_right = tree_sampling_coordinators[pid + 1]->get_likelihood();
		Real_t swap_acceptance_ratio = (temperatures[pid] - temperatures[pid + 1]) * (likelihood_right - likelihood_left);
		if (random.logUniform() <= swap_acceptance_ratio)
		{
			tree_sampling_coordinators[pid]->set_temperature(temperatures[pid + 1]);
			tree_sampling_coordinators[pid + 1]->set_temperature(temperatures[pid]);
			std::swap(tree_sampling_coordinators[pid], tree_sampling_coordinators[pid + 1]);
			std::swap(likelihood_calculators[pid], likelihood_calculators[pid + 1]);
            std::swap(tree_sampling_coordinators[pid]->pid, tree_sampling_coordinators[pid + 1]->pid);
            if (pid == 0) {
                        std::swap(tree_sampling_coordinators[pid]->edge_count, tree_sampling_coordinators[pid + 1]->edge_count);
            }
		}
	}

	std::pair<Gauss::GaussianMixture<Real_t>, Gauss::Gaussian<Real_t>> prepare_starting_parameters()
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
		EMEstimator < Real_t > EM(provider.cells_to_vector());
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
		mixture.prune();
		return std::make_pair(mixture, gaussian);
	}

	std::tuple<PointerTree, std::vector<BreakpointPair>, double> choose_best_tree()
	{
		auto best = tree_sampling_coordinators[0]->getBestTreeData();
		for (size_t i = 1 ; i < tree_sampling_coordinators.size(); i++)
		{
			auto best_local = tree_sampling_coordinators[i]->getBestTreeData();
			if (std::get<2>(best_local) > std::get<2>(best))
			{
				best = best_local;
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


std::tuple<std::tuple<PointerTree, std::vector<BreakpointPair>, Real_t>, std::string, Real_t, std::map<std::pair<BreakpointPair, BreakpointPair>, size_t>> simulate(size_t iterations_parameters, size_t iterations_pt)
	{
		auto parameters = prepare_starting_parameters();
		auto parameters_MAP = estimateParameters(parameters.first, parameters.second, iterations_parameters);
		std::get<0>(parameters_MAP).prune();
		NormalMixtureLikelihood<Real_t> likelihood(std::get<1>(parameters_MAP), std::get<0>(parameters_MAP), random.nextInt());
		prepare_sampling_services(&likelihood);
		mcmc_simulation(iterations_pt);
		auto bestTree = choose_best_tree();
		return std::make_tuple(bestTree, likelihood_calculators[0]->likelihood->toString(), -1 , tree_sampling_coordinators[0]->edge_count);
	}
	
};


















#endif //PARALLEL_TEMPERING_COORDINATOR_H
