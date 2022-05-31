#ifndef LIKELIHOOD_COORD_H
#define LIKELIHOOD_COORD_H
#include <algorithm>
#include <numeric>
#include <thread>  
#include <map>

#include "utils/log_sum_accumulator.h"
#include "cell_provider/vector_cell_provider.h"
#include "tree/event_tree.h"
#include "likelihood/implementations/normal_mixture_likelihood.h"
#include "utils/random.h"
#include "moves/move_type.h"
#include "parameters/parameters.h"
#include "tree/tree_counts_scoring.h"
#include "likelihood_calculator.h"


template<class Real_t> class LikelihoodCoordinator {
public:
	LikelihoodCalculatorState<Real_t> calculator_state; 
	LikelihoodCalculatorState<Real_t> tmp_calculator_state; 

	LikelihoodMatrices<Real_t> likelihood_matrices;
	LikelihoodMatrices<Real_t> tmp_likelihood_matrices;
	LikelihoodData<Real_t> likelihood;
	EventTree &tree;
	VectorCellProvider<Real_t> const *cells;
	using NodeHandle = EventTree::NodeHandle;

	Random<Real_t> random;

	Real_t MAP;
	CountsDispersionPenalty<Real_t> counts_scoring;
	size_t step {0};
	LikelihoodData<Real_t> MAP_parameters;

public:
	void swapLikelihoodsBuffers()
	{
		LikelihoodMatrices<Real_t>::swap(likelihood_matrices, tmp_likelihood_matrices);
	}

	void fillLikelihoodMatrices() { 
		likelihood.fill_breakpoint_log_likelihood_matrix(likelihood_matrices.breakpoint_likelihoods, cells->get_corrected_counts());
		likelihood.fill_no_breakpoint_log_likelihood_matrix(likelihood_matrices.no_breakpoint_likelihoods, cells->get_corrected_counts());
	}

	void updateDataAfterParameterChange() {
		fillLikelihoodMatrices();
		calculate_likelihood();
	}

	std::pair<Real_t, Real_t> resampleParametersGibbs() {
		step = (step + 1) % (3*likelihood.brkp_likelihood.number_of_components() + 1);
		if (step == 0) {
			return likelihood.no_brkp_likelihood.resample_standard_deviation();
		}
		else {
			size_t mixture_component = (step - 1)/3;
			if ((step - 1) % 3 == 0) {
				return likelihood.brkp_likelihood.resample_weight(mixture_component);
			} else if ((step - 1) % 3 == 1) {
				return likelihood.brkp_likelihood.resample_sd(mixture_component);
			} 
			return likelihood.brkp_likelihood.resample_mean(mixture_component);
		}
	}


public:
	LikelihoodCoordinator(LikelihoodData<Real_t> lk, EventTree &tree,
		VectorCellProvider<Real_t> *cells, size_t maxBreakpoint, unsigned int seed) : 
					calculator_state{cells->get_cells_count()}, 
					tmp_calculator_state{cells->get_cells_count()}, 
					likelihood_matrices{maxBreakpoint + 1, cells->get_cells_count()},
					tmp_likelihood_matrices{maxBreakpoint + 1, cells->get_cells_count()},
					likelihood{ lk }, tree{ tree }
		, cells{ cells }, random{ seed },
		counts_scoring{ cells }, MAP_parameters{lk.no_brkp_likelihood, lk.brkp_likelihood} {
		
		updateDataAfterParameterChange();
		persist_likelihood_calculation_result();
		MAP = get_likelihood() + likelihood.get_likelihood_parameters_prior() -100000000.0;
	}

	std::vector<TreeLabel> &get_max_attachment() {
		return calculator_state.max_attachment;
	}

	Real_t get_likelihood() {
		return calculator_state.likelihood;
	}

	void persist_likelihood_calculation_result() { 
		LikelihoodCalculatorState<Real_t>::swap(calculator_state, tmp_calculator_state);
	}

	std::vector<TreeLabel> &calculate_max_attachment() {
		return tmp_calculator_state.max_attachment;
	}

	Real_t calculate_likelihood() { 
		LikelihoodCalculator<Real_t> calc{tree, tmp_calculator_state, cells, likelihood_matrices};
		return calc.calculate_likelihood();
	}
	
	/****** MH sampling of parameters ******/

	void resampleParameters(Real_t log_tree_prior, Real_t tree_count_score) {
		auto likelihoodBeforeChange =  get_likelihood();
		auto priorBeforeChange = likelihood.get_likelihood_parameters_prior();
		LikelihoodData<Real_t> previous_parameters = likelihood;
		auto logKernels = resampleParametersGibbs();
		swapLikelihoodsBuffers();
		updateDataAfterParameterChange();
		Attachment tmp_att{tmp_calculator_state.max_attachment};
		auto counts_score_after_move = counts_scoring.calculate_log_score(tree, tmp_att);
		auto likelihoodAfterChange =  tmp_calculator_state.likelihood;
		auto priorAfterChange = likelihood.get_likelihood_parameters_prior();

		if (!likelihood.likelihood_is_valid()) {
			swapLikelihoodsBuffers();
			likelihood = previous_parameters;
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
			if (likelihoodAfterChange + priorAfterChange + log_tree_prior + counts_score_after_move > MAP)
			{
				MAP = likelihoodAfterChange + priorAfterChange + log_tree_prior + counts_score_after_move;
				MAP_parameters = likelihood;
			}
			persist_likelihood_calculation_result();
		}
		else {
			swapLikelihoodsBuffers();
			likelihood = previous_parameters;
			if (DEBUG) {
				logDebug("Rejecting parameters change");
			}
		}
	}

};
#endif // !LIKELIHOOD_COORD_H

