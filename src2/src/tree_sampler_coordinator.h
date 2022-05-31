#ifndef TREE_SAMPLER_COORDINATOR_H
#define TREE_SAMPLER_COORDINATOR_H
#include <algorithm>
#include <vector>
#include <tuple>

#include "tree/tree_formatter.h"
#include "utils/logger/logger.h"
#include "tree/tree_counts_scoring.h"
#include "tree/attachment.h"
#include "conet_result.h"
#include "tree_mh_steps_executor.h"
/**
* Class coordinating all components of MH sampling.
* This class is an owner of VertexSetInterface, TreeNodeSampler,
* EventTree, LikelihoodCoordinator and serves as a coordinating service for 
* their tasks.
*
* Tree MH sampling is done in few steps:\n 
* 1) First move type is sampled. \n 
* 2) TreeNodeSampler confirms that sampled move type is possible (for instance <code> DELETE_LEAF</code> move may not be possible \n 
* if the tree has zero size)\n 
* 3) Any required random elements are sampled (for example labels for a new node).\n 
* 4) VertexLabelSampler, EventTree, TreeNodeSampler are ordered to do their move specific transformations\n 
* 5) If move has not been accepted a move rollback is carried out\n
*
* 
*/
template <class Real_t> class TreeSamplerCoordinator {
	using NodeHandle = EventTree::NodeHandle;
	using MoveData = typename MHStepsExecutor<Real_t>::MoveData;
private: 
	
public:
	EventTree &tree;
	LikelihoodCoordinator<Real_t> *likelihoodCalculator;
	VectorCellProvider<Real_t> *cells;
	CountsDispersionPenalty<Real_t> countsScoring;


	Random<Real_t> random;
	std::map<MoveType, Real_t> moveProbability;
	Real_t temperature {1.0};
	Real_t tree_count_dispersion_penalty {1.0};
	CONETInferenceResult<Real_t> best_found_tree;		
	MHStepsExecutor<Real_t> mh_step_executor; 
	

	Real_t get_probability_of_reverse_move(MoveType type) {
		switch (type) {
		case ADD_LEAF:
			return moveProbability[DELETE_LEAF];
		case DELETE_LEAF:
			return moveProbability[ADD_LEAF];
		case CHANGE_LABEL:
			return moveProbability[CHANGE_LABEL];
		case SWAP_LABELS:
			return moveProbability[SWAP_LABELS];
		case PRUNE_REATTACH:
			return moveProbability[PRUNE_REATTACH];
		case SWAP_SUBTREES:
			return moveProbability[SWAP_SUBTREES];
		case SWAP_ONE_BREAKPOINT:
		default:
			return moveProbability[SWAP_ONE_BREAKPOINT];
		}
	}

	void move(MoveType type) {
		recalculate_counts_dispersion_penalty();
		auto beforeMoveLikelihood = temperature * likelihoodCalculator->get_likelihood() + mh_step_executor.get_log_tree_prior() + tree_count_dispersion_penalty;
		auto moveData = mh_step_executor.execute_move(type);

		auto afterMoveLikelihood = temperature * likelihoodCalculator->calculate_likelihood() + mh_step_executor.get_log_tree_prior();
		Attachment after_tmp_att {likelihoodCalculator->calculate_max_attachment()};
		auto after_move_counts_dispersion_penalty = countsScoring.calculate_log_score(tree, after_tmp_att);
		afterMoveLikelihood += after_move_counts_dispersion_penalty;

		Real_t log_acceptance = afterMoveLikelihood - beforeMoveLikelihood + moveData.reverse_move_log_kernel - moveData.move_log_kernel + std::log(moveProbability[type]) - std::log(get_probability_of_reverse_move(type));

		if (DEBUG ) logDebug("Log acceptance ratio: ", log_acceptance, " likelihood before ", beforeMoveLikelihood, " likelihood after ", afterMoveLikelihood);

		if (random.logUniform() <= log_acceptance) {
			likelihoodCalculator->persist_likelihood_calculation_result();
			tree_count_dispersion_penalty = after_move_counts_dispersion_penalty;
			if (DEBUG) logDebug("Move accepted");
		}
		else {
			mh_step_executor.rollback_move(type, moveData);
			if (DEBUG) log("Move rejected");
		}
	}

	MoveType sample_move_type() {
		std::map<size_t, MoveType> type_to_index;
		std::vector<Real_t> weights;
		for (auto x : moveProbability) {
			type_to_index[weights.size()] = x.first;
			weights.push_back(x.second);
		}
		return type_to_index[random.discrete(weights)];
	}

	void recalculate_counts_dispersion_penalty() {
		Attachment at{likelihoodCalculator->get_max_attachment()};
		tree_count_dispersion_penalty = countsScoring.calculate_log_score(tree, at);
	}


public:
	TreeSamplerCoordinator(EventTree &tree, LikelihoodCoordinator<Real_t> *lC, unsigned int seed, size_t maxBrkp, VectorCellProvider<Real_t> *cells, std::map<MoveType, Real_t> moveProbability, int pid): 
		tree{ tree }, 
		likelihoodCalculator{ lC }, 
		cells{ cells },
		countsScoring{ cells },
		random{ seed }, 
		moveProbability{ moveProbability },
		best_found_tree{tree, Attachment{lC->get_max_attachment()}, 0.0},
		mh_step_executor{tree, cells, random}  {
		best_found_tree.likelihood = likelihoodCalculator->get_likelihood() + countsScoring.calculate_log_score(tree, best_found_tree.attachment) + mh_step_executor.get_log_tree_prior();
	}

	Real_t get_likelihood_without_priors_and_penalty() {
		return likelihoodCalculator->get_likelihood();
	}

	void set_temperature(Real_t temperature) {
		this->temperature = temperature;
	}
    
	Real_t get_log_tree_prior() {
		return mh_step_executor.get_log_tree_prior();
	}

    Real_t get_total_likelihood() {
		return likelihoodCalculator->get_likelihood() + mh_step_executor.get_log_tree_prior() + tree_count_dispersion_penalty;
    }

	Real_t get_temperature() const {
		return this->temperature;
	}

	void execute_metropolis_hastings_step() {
		MoveType type = sample_move_type();
		if (DEBUG) logDebug("Sampled move of type: ", moveTypeToString(type));
		
		if (mh_step_executor.move_is_possible(type)) {
			move(type);
		}

		auto l = get_total_likelihood();
		if (l > best_found_tree.likelihood) {
			best_found_tree = CONETInferenceResult<Real_t>(tree, Attachment(likelihoodCalculator->get_max_attachment()), l);
		}
	}

	CONETInferenceResult<Real_t> get_inferred_tree() {
		return best_found_tree;
	}
};

#endif // !TREE_SAMPLER_COORDINATOR_H
