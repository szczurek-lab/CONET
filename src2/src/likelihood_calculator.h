#ifndef LIKELIHOOD_CALCULATOR_H
#define LIKELIHOOD_CALCULATOR_H
#include <algorithm>
#include <numeric>
#include <thread>  
#include <map>

#include "utils/log_sum_accumulator.h"
#include "cell_provider/vector_cell_provider.h"
#include "tree/event_tree.h"
#include "likelihood/implementations/normal_mixture_likelihood.h"
#include "utils/random.h"
#include "moves/move_data.h"
#include "moves/move_type.h"
#include "parameters/parameters.h"
#include "tree/tree_counts_scoring.h"

#include "lik_calc.h"
/*
 
*/
template<class Real_t> class LikelihoodCalculator {
public:
	LikelihoodCalculatorState<Real_t> calculator_state; 
	LikelihoodCalculatorState<Real_t> tmp_calculator_state; 

	std::vector<LogWeightAccumulator<Real_t>> likelihood_result;
	std::vector<Real_t> rootLikelihoods;
	NormalMixtureLikelihood<Real_t> *likelihood;
	EventTree *tree;
	VectorCellProvider<Real_t> const *cells;
	using NodeHandle = EventTree::NodeHandle;
	std::vector<std::vector<Real_t>> breakpointLikelihoods;
	std::vector<std::vector<Real_t>> noBreakpointLikelihoods;

	std::vector<std::vector<Real_t>> buffer_breakpointLikelihoods;
	std::vector<std::vector<Real_t>> buffer_noBreakpointLikelihoods;

	std::vector<LogWeightAccumulator<Real_t>> tmp_likelihood_accumulator;
	size_t maxBreakpoint;
	Random<Real_t> random;

	std::vector<Event> attachment;
	std::vector<Event> tmp_attachment;

	Real_t MAP;

	Real_t temperature{ 1.0 };
	CountsDispersionPenalty<Real_t> counts_scoring;

	void set_temperature(Real_t t)
	{
		this->temperature = t;
		
	}

public:
	void swapLikelihoodsBuffers()
	{
		std::swap(breakpointLikelihoods, buffer_breakpointLikelihoods);
		std::swap(noBreakpointLikelihoods, buffer_noBreakpointLikelihoods);
	}
	

	Real_t sumLikelihoods(std::vector<LogWeightAccumulator<Real_t>> &likelihoods) {
		Real_t result_ = 0.0;
		if (!USE_EVENT_LENGTHS_IN_ATTACHMENT) {
			std::for_each(likelihoods.rbegin(), likelihoods.rend(),
				[&](LogWeightAccumulator<Real_t> &acc) { result_ += acc.getResult() - std::log((Real_t)(tree->get_size() - 1)); });
		}
		else {
			std::for_each(likelihoods.rbegin(), likelihoods.rend(),
				[&](LogWeightAccumulator<Real_t> &acc) { result_ += acc.getResult(); });
		}
		return result_;
	}

	void fillLikelihoodMatrices() {
			likelihood->getLogLikelihoodBreakpoint(breakpointLikelihoods, cells->get_corrected_counts());
			likelihood->getLogLikelihoodNoBreakpoint(noBreakpointLikelihoods, cells->get_corrected_counts());
	}

	void init_state_likelihood_matrices() {
		calculator_state.breakpoint_likelihoods.resize(maxBreakpoint + 1);
		calculator_state.no_breakpoint_likelihoods.resize(maxBreakpoint + 1);
		tmp_calculator_state.breakpoint_likelihoods.resize(maxBreakpoint + 1);
		tmp_calculator_state.no_breakpoint_likelihoods.resize(maxBreakpoint + 1);
		
		for (size_t br = 0; br <= maxBreakpoint; br++) {
			calculator_state.breakpoint_likelihoods[br].resize(cells->get_cells_count());
			calculator_state.no_breakpoint_likelihoods[br].resize(cells->get_cells_count());
			tmp.calculator_state.breakpoint_likelihoods[br].resize(cells->get_cells_count());
			tmp_calculator_state.no_breakpoint_likelihoods[br].resize(cells->get_cells_count());
		}
	}

	void clearTmpAccumulator(size_t left, size_t right) {
		for (size_t s = left; s < right; s++) {
			tmp_likelihood_accumulator[s].clear();
		}
	}

	void addAncestorLikelihoods(std::unordered_set<size_t> &ancestorBreakpoints, bool reverse, size_t left, size_t right) {
		Real_t sign = reverse ? -1.0 : 1.0;
		for (auto br : ancestorBreakpoints) {
			for (size_t c = left; c < right; c++) {
				rootLikelihoods[c] += sign* breakpointLikelihoods[br][c];
				rootLikelihoods[c] -= sign* noBreakpointLikelihoods[br][c];
			}
		}
	}

	void updateDataAfterParameterChange() {
		swapLikelihoodsBuffers();
		fillLikelihoodMatrices();
		get_after_move_likelihood();
	}


	void getTreeTmpLikelihoodVanilla(size_t left, size_t right, std::vector<Real_t> &max_likelihood) {
		calculateRootLikelihoods(left, right);

		for (size_t c = left; c < right; c++) {
			tmp_attachment[c] = std::make_pair(0, 0);
			max_likelihood[c] = rootLikelihoods[c];
		}
		
		getSubtreeLikelihood( true, left, right, max_likelihood);
	}
	
public:
	LikelihoodCalculator(NormalMixtureLikelihood<Real_t> *lk, EventTree *tree,
		VectorCellProvider<Real_t> *cells, size_t maxBreakpoint, unsigned int seed) : calculator_state{cells->get_cells_count(), tree}, tmp_calculator_state{cells->get_cells_count(), tree},likelihood{ lk }, tree{ tree }
		, cells{ cells }, maxBreakpoint{ maxBreakpoint }, random{ seed }, attachment{ cells->get_cells_count() }, tmp_attachment{ cells->get_cells_count() },
		counts_scoring{ cells } {
		
		init_state_likelihood_matrices();
		updateDataAfterParameterChange();
		swapLikelihoodWithTmpAccumulator();
		MAP = getLikelihood() + likelihood->getParamsPrior() -100000000.0;
	}

	size_t get_cells_count() {
		return cells->get_cells_count();
	}

	std::vector<Event> &getBestAttachment() {
		return attachment;
	}

	Real_t getLikelihood() {
		return sumLikelihoods(likelihood_result);
	}

	Real_t getTmpLikelihood()
	{
		return sumLikelihoods(tmp_likelihood_accumulator);
	}

	void swapLikelihoodWithTmpAccumulator() {
		std::swap(likelihood_result, tmp_likelihood_accumulator);
		std::swap(attachment, tmp_attachment);
	}

	std::vector<Event> &get_after_move_best_attachment() {
		return tmp_attachment;
	}

	Real_t get_after_move_likelihood() {
		std::vector<Real_t> max_likelihoods;
		for (size_t c = 0; c < cells->get_cells_count(); c++)
		{
			max_likelihoods.push_back(0.0);
		}
		getTreeTmpLikelihoodVanilla(0, cells->get_cells_count(), max_likelihoods);
		return sumLikelihoods(tmp_likelihood_accumulator);
	}
	
	/****** MH sampling of parameters ******/

	void resampleParameters(Real_t log_tree_prior, Real_t tree_count_score) {
		if (DEBUG) {
			logDebug(likelihood->toString());
		}
		auto likelihoodBeforeChange = temperature * getLikelihood();
		auto priorBeforeChange = likelihood->getParamsPrior();
		auto logKernels = likelihood->resampleParameters();
		updateDataAfterParameterChange();
		Attachment tmp_att{tmp_attachment};
		auto counts_score_after_move = counts_scoring.calculate_log_score(*tree, tmp_att);
		auto likelihoodAfterChange =  temperature * getTmpLikelihood();
		auto priorAfterChange = likelihood->getParamsPrior();

		if (!likelihood->valid()) {
			likelihood->rollBackParametersResample();
			swapLikelihoodsBuffers();
			return;
		}

		Real_t acceptanceRatio = likelihoodAfterChange + priorAfterChange  - likelihoodBeforeChange - priorBeforeChange 
			+ logKernels.second - logKernels.first
			+ (counts_score_after_move - tree_count_score);
		if (DEBUG) {
			logDebug("Parameters acceptance ratio equal to ", std::to_string(acceptanceRatio));
			logDebug("Likelihood before change: ", std::to_string(likelihoodBeforeChange));
			logDebug("Likelihood after change: ", std::to_string(likelihoodAfterChange));
			logDebug("Prior after change: ", std::to_string(priorAfterChange));
			logDebug("Prior before change: ", std::to_string(priorBeforeChange));
			logDebug("Log kernels ", std::to_string(logKernels.second), " ", std::to_string(logKernels.first));
		}
		if (random.logUniform() <= acceptanceRatio) {
			if (DEBUG) {
				logDebug("Accepting parameters change");
			}
			likelihood->acceptSampledParameters();
			if (likelihoodAfterChange + priorAfterChange + log_tree_prior + counts_score_after_move > MAP)
			{
				MAP = likelihoodAfterChange + priorAfterChange + log_tree_prior + counts_score_after_move;
				likelihood->update_map();
			}
			swapLikelihoodWithTmpAccumulator();
		}
		else {
			if (DEBUG) {
				logDebug("Rejecting parameters change");
			}
			likelihood->rollBackParametersResample();
			swapLikelihoodsBuffers();
		}
	}

};
#endif // !LIKELIHOOD_CALCULATOR_H

