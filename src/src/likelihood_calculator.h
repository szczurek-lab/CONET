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
/*
 
*/
template<class Real_t> class LikelihoodCalculator {
public:
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
	CountsScoring<Real_t> counts_scoring;

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
	
	void calculateRootLikelihoods(size_t left, size_t right) {
		for (size_t c = left; c < right; c++) {
			rootLikelihoods[c] = 0.0;
		}
		for (size_t br = 0; br <= maxBreakpoint; br++) {
			for (size_t c = left; c < right; c++) {
				rootLikelihoods[c] += noBreakpointLikelihoods[br][c];
			}
		}
	}

	void applyDifferencesForNodeLoci(NodeHandle node, std::vector<Real_t> &fatherLikelihood, size_t left, size_t right) {
		auto breakpoints = tree->getNewBreakpoints(node);
		for (auto br : breakpoints) {
			for (size_t c = left; c < right; c++) {
				fatherLikelihood[c] += breakpointLikelihoods[br][c] - noBreakpointLikelihoods[br][c];
			}
		}
	}

	void reverseDifferencesForNodeLoci(NodeHandle node, std::vector<Real_t> &fatherLikelihood,size_t left, size_t right) {
		auto breakpoints = tree->getNewBreakpoints(node);
		for (auto br : breakpoints) {
			for (size_t c = left; c < right; c++) {
				fatherLikelihood[c] += -breakpointLikelihoods[br][c] + noBreakpointLikelihoods[br][c];
			}
		}
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
			likelihood->getLogLikelihoodBreakpoint(breakpointLikelihoods, cells->getCellsToLoci());
			likelihood->getLogLikelihoodNoBreakpoint(noBreakpointLikelihoods, cells->getCellsToLoci());
	}

	void initLikelihoodMatrices() {
		breakpointLikelihoods.resize(maxBreakpoint + 1);
		noBreakpointLikelihoods.resize(maxBreakpoint + 1);
		buffer_breakpointLikelihoods.resize(maxBreakpoint + 1);
		buffer_noBreakpointLikelihoods.resize(maxBreakpoint + 1);
		
		for (size_t br = 0; br <= maxBreakpoint; br++) {
			breakpointLikelihoods[br].resize(cells->getCellsCount());
			noBreakpointLikelihoods[br].resize(cells->getCellsCount());
			buffer_breakpointLikelihoods[br].resize(cells->getCellsCount());
			buffer_noBreakpointLikelihoods[br].resize(cells->getCellsCount());
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

	/** Standard likelihood **/
	void treeDFS(NodeHandle node, std::vector<Real_t> &fatherLikelihood,
		std::vector<LogWeightAccumulator<Real_t>> &result, size_t left, size_t right, std::vector<Real_t> &max_likelihoods,
		std::map<Event, Real_t> &events_lengths) {
		applyDifferencesForNodeLoci(node, fatherLikelihood, left, right);

		for (size_t c = left; c < right; c++) {
			if (USE_EVENT_LENGTHS_IN_ATTACHMENT)
			{
				result[c].add(fatherLikelihood[c]  + events_lengths[tree->get_node_event(node)]);	
			}
			else {
				result[c].add(fatherLikelihood[c]);
			}
			if (max_likelihoods[c] < fatherLikelihood[c])
			{
				max_likelihoods[c] = fatherLikelihood[c];
				tmp_attachment[c] = tree->get_node_event(node);
			}
		}
		auto children = tree->get_children(node);
		for (auto ch : children) {
			treeDFS(ch, fatherLikelihood, result, left, right, max_likelihoods, events_lengths);
		}
		reverseDifferencesForNodeLoci(node, fatherLikelihood, left, right);
	}

	void get_normalized_event_lengths(NodeHandle node, Real_t length, std::map<Event, Real_t> &result, Real_t depth)
	{
		if (node != tree->get_root()) {
			length += cells->getEventLength(node->label);
			result[tree->get_node_event(node)] = length / depth;
		} else {
			result[tree->get_node_event(node)] = 0.0;
		} 
		
		for (auto child : tree->get_children(node)) {
			get_normalized_event_lengths(child, length, result, depth + 1);
		}
	}

	void getSubtreeLikelihood(bool clearAcc, size_t left, size_t right, std::vector<Real_t> &max_likelihoods) {
		std::map<Event, Real_t> events_lengths;
		if (USE_EVENT_LENGTHS_IN_ATTACHMENT)
		{
			get_normalized_event_lengths(tree->get_root(), 0.0, events_lengths, 0.0);
			Real_t sum = std::accumulate(events_lengths.begin(), events_lengths.end(), 0.0,
				[](const Real_t previous, decltype(*events_lengths.begin()) p) { return previous + std::exp(-p.second); });
			for (auto &el : events_lengths)
			{
				el.second = -el.second - std::log(sum);
			}
		}

		if (clearAcc) {
			clearTmpAccumulator(left, right);
		}
		for (auto &node : tree->get_children(tree->get_root())) {
			treeDFS(node, rootLikelihoods, tmp_likelihood_accumulator, left, right, max_likelihoods, events_lengths);
		}
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
		VectorCellProvider<Real_t> *cells, size_t maxBreakpoint, unsigned int seed) : likelihood{ lk }, tree{ tree }
		, cells{ cells }, maxBreakpoint{ maxBreakpoint }, random{ seed }, attachment{ cells->getCellsCount() }, tmp_attachment{ cells->getCellsCount() },
		counts_scoring{ cells , COUNTS_SCORE_CONSTANT_0 != 0.0 && COUNTS_SCORE_CONSTANT_1 != 0.0} {
		likelihood_result.resize(cells->getCellsCount());
		tmp_likelihood_accumulator.resize(cells->getCellsCount());
		rootLikelihoods.resize(cells->getCellsCount());
		initLikelihoodMatrices();
		updateDataAfterParameterChange();
		swapLikelihoodWithTmpAccumulator();
		MAP = getLikelihood() + likelihood->getParamsPrior() -100000000.0;
	}

	size_t getCellsCount() {
		return cells->getCellsCount();
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
		for (size_t c = 0; c < cells->getCellsCount(); c++)
		{
			max_likelihoods.push_back(0.0);
		}
		getTreeTmpLikelihoodVanilla(0, cells->getCellsCount(), max_likelihoods);
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
		auto counts_score_after_move = counts_scoring.calculate_log_score(*tree, tmp_attachment);
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

