#ifndef TREE_SAMPLER_COORDINATOR_H
#define TREE_SAMPLER_COORDINATOR_H
#include <algorithm>
#include <vector>
#include <tuple>

#include "tree/event_tree.h"
#include "utils/random.h"
#include "moves/move_type.h"
#include "utils/logger/logger.h"
#include "tree/node_set.h"
#include "tree/vertex_label_sampler.h"
#include "moves/move_data.h"
#include "cell_provider/vector_cell_provider.h"
#include "likelihood_calculator.h"
#include "tree/tree_counts_scoring.h"
#include "tree/attachment.h"

/**
* Class coordinating all components of MH sampling.
* This class is an owner of VertexSetInterface, TreeNodeSampler,
* EventTree, LikelihoodCalculator and serves as a coordinating service for 
* their tasks.
*
* Tree MH sampling is done in few steps:\n 
* 1) First move type is sampled. \n 
* 2) TreeNodeSampler confirms that sampled move type is possible (for instance <code> DELETE_LEAF</code> move may not be possible \n 
* if the tree has zero size)\n 
* 3) Any required random elements are sampled (for example labels for a new node).\n 
* 4) VertexLabelSampler, EventTree, TreeNodeSampler are ordered to do their move specific transformations\n 
* 5) If move has not been accepted a move roollback is carried out\n 
*
* 
*/
template <class Real_t> class TreeSamplerCoordinator {
public:
	using NodeHandle = EventTree::NodeHandle;
	EventTree *tree;
	LikelihoodCalculator<Real_t> *likelihoodCalculator;
	VertexLabelSampler<Real_t> *label_sampler;
	TreeNodeSampler<Real_t> *nodeSet;
	VectorCellProvider<Real_t> *cells;
	CountsScoring<Real_t> countsScoring;


	Random<Real_t> random;
	std::map<MoveType, Real_t> moveProbability;
	int pid;
	/*** for statistics gathering ***/
	Real_t bestLikelihood;
       Real_t bestLikelihoodNoPrior;	
    Real_t best_counts_score;
	EventTree bestTree;
	std::vector<Event> bestTreeAttachment;
	size_t moveCount{ 0 };
    std::map<std::pair<Event, Event>, size_t> edge_count;
	Real_t temperature {1.0};

	Real_t tree_count_score = 1.0;
			

	void prune_and_reattach(NodeHandle nodeToPrune, NodeHandle nodeToAttach) {
		auto oldParent = tree->prune_and_reattach(nodeToPrune, nodeToAttach);
		nodeSet->refresh_node_data(oldParent);
		nodeSet->refresh_node_data(nodeToAttach);
	}

	void swapLabels(NodeHandle node1, NodeHandle node2) {
		tree->swapLabels(node1, node2);
	}

	/*
	  Returns handle to parent of @node
	*/
	NodeHandle deleteLeaf(NodeHandle node) {
		nodeSet->detach_leaf(node);
		label_sampler->remove_label(tree->get_node_event(node));
		auto parent = tree->delete_leaf(node);
		nodeSet->refresh_node_data(parent);
		return parent;
	}

	NodeHandle addLeaf(NodeHandle parent, TreeLabel label) {
		NodeHandle newNode = tree->add_leaf(parent, label);
		nodeSet->refresh_node_data(newNode);
		nodeSet->refresh_node_data(parent);
		label_sampler->add_label(label);
		return newNode;
	}

	void change_label(NodeHandle node, std::pair<size_t, size_t> newBrkp) {
		label_sampler->remove_label(tree->get_node_event(node));
		tree->change_label(node, newBrkp);
		label_sampler->add_label(newBrkp);
	}

	void swap_subtrees_non_descendants(NodeHandle root1, NodeHandle root2) {
		tree->swap_subtrees_non_descendants(root1, root2);
	}

	void swap_subtrees_descendants(NodeHandle parent, NodeHandle child, NodeHandle grandChild) {
		tree->swap_subtrees_descendants(parent, child, grandChild);
		nodeSet->refresh_node_data(parent);
		nodeSet->refresh_node_data(grandChild);
		nodeSet->refresh_node_data(child);
	}

	void swapOneBreakpoint(NodeHandle node1, NodeHandle node2, int left, int right) {
		if (!label_sampler->swapOneBreakpointPossible(tree->get_node_event(node1), tree->get_node_event(node2), left, right)) {
			return;
		}
		auto newBrkp = label_sampler->swapOneBreakpoint(tree->get_node_event(node1), tree->get_node_event(node2), left, right);
		tree->change_label(node1, newBrkp.first);
		tree->change_label(node2, newBrkp.second);
	}

	/*
		End to end MCMC moves implementation:
		All of those methods return triple MoveData, Real_t, Real_t
		where MoveData is used to record data needed to revert the move,
		Real_t is kernel value for the move
		the second Real_t is a kernel value for reversed move.
	*/
	std::tuple<MoveData, Real_t, Real_t> deleteLeafMove() {
		MoveData beforeMove;
		Real_t kernelValue = nodeSet->getDeleteLeafKernel();
		auto leaf = nodeSet->sample_leaf(random);
		auto label = tree->get_node_label(leaf);
		auto oldParent = deleteLeaf(leaf);
		beforeMove.oldParent = oldParent;
		beforeMove.old_label = label;
		return std::make_tuple<MoveData, Real_t, Real_t>(std::move(beforeMove), std::move(kernelValue), nodeSet->getAddLeafKernel() + label_sampler->get_sample_label_log_kernel());
	}

	std::tuple<MoveData, Real_t, Real_t> addLeafMove(TreeLabel label) {
		MoveData beforeMove;
		Real_t kernelValue = nodeSet->getAddLeafKernel() + label_sampler->get_sample_label_log_kernel();
		auto parent = nodeSet->sample_node(true, random);
		auto newNode = addLeaf(parent, label);
		beforeMove.node1 = newNode;
		return std::make_tuple<MoveData, Real_t, Real_t>(std::move(beforeMove), std::move(kernelValue), nodeSet->getDeleteLeafKernel());
	}
   
	std::tuple<MoveData, Real_t, Real_t> pruneAndReattachMove() {
		NodeHandle node_to_prune = nodeSet->sample_node(false, random);
		NodeHandle node_to_attach = nodeSet->sample_non_descendant(node_to_prune, random);
		MoveData beforeMove;
		beforeMove.oldParent = tree->get_parent(node_to_prune);
		beforeMove.node1 = node_to_prune;
		prune_and_reattach(node_to_prune, node_to_attach);
		return std::make_tuple<MoveData, Real_t, Real_t>(std::move(beforeMove), 1.0, 1.0);
	}

	std::tuple<MoveData, Real_t, Real_t> swapLabelsMove() {
		MoveData beforeMove;
		auto nodes = nodeSet->sampleTwoNodes(random);
		beforeMove.node1 = nodes.first;
		beforeMove.node2 = nodes.second;
		swapLabels(nodes.first, nodes.second);
		return std::make_tuple<MoveData, Real_t, Real_t>(std::move(beforeMove), 1.0, 1.0);
	}

	std::tuple<MoveData, Real_t, Real_t> changeLabelMove() {
		MoveData beforeMove;
		auto node = nodeSet->sample_node(false, random);
		auto old_label = tree->get_node_label(node);
		auto breakpoints = label_sampler->sample_label(random);
		change_label(node, breakpoints);
		beforeMove.node1 = node;
		beforeMove.old_label = old_label;
		return std::make_tuple<MoveData, Real_t, Real_t>(std::move(beforeMove), 1.0, 1.0);
	}

	std::tuple<MoveData, Real_t, Real_t> swapSubtreesDescendantsMove(NodeHandle parent, NodeHandle child) {
		MoveData beforeMove;
		auto grandChild = nodeSet->sample_descendant(child, random);
		beforeMove.simpleSwap = false;
		beforeMove.oldParent = child;
		beforeMove.node1 = parent;
		beforeMove.node2 = tree->get_parent(child);
		Real_t kernelValue = nodeSet->getSwapSubtreesDescendantsKernel(child);
		swap_subtrees_descendants(parent, child, grandChild);
		return std::make_tuple<MoveData, Real_t, Real_t>(std::move(beforeMove), std::move(kernelValue), nodeSet->getSwapSubtreesDescendantsKernel(parent));
	}

	std::tuple<MoveData, Real_t, Real_t> swapSubtreesNonDescendantsMove(NodeHandle node1, NodeHandle node2) {
		MoveData beforeMove;
		beforeMove.node1 = node1;
		beforeMove.node2 = node2;
		beforeMove.simpleSwap = true;
		swap_subtrees_non_descendants(node1, node2);
		return std::make_tuple<MoveData, Real_t, Real_t>(std::move(beforeMove), 1.0, 1.0);
	}

	std::tuple<MoveData, Real_t, Real_t> swapSubtreesMove() {
		auto nodes = nodeSet->sampleTwoNodes(random);
		int descendants = tree->get_nodes_relation(nodes.first, nodes.second);
		if (descendants != 0) {
			auto parent = descendants == -1 ? nodes.second : nodes.first;
			auto child = descendants == 1 ? nodes.second : nodes.first;
			return swapSubtreesDescendantsMove(parent, child);
		}
		else {
			return swapSubtreesNonDescendantsMove(nodes.first, nodes.second);
		}
	}
	std::tuple<MoveData, Real_t, Real_t> swapOneBreakpointMove() {
		MoveData beforeMove;
		auto nodes = nodeSet->sampleTwoNodes(random);
		beforeMove.node1 = nodes.first;
		beforeMove.node2 = nodes.second;
		beforeMove.old_label = tree->get_node_label(nodes.first);
		swapOneBreakpoint(nodes.first, nodes.second, random.randomIntBit(), random.randomIntBit());
		return std::make_tuple<MoveData, Real_t, Real_t>(std::move(beforeMove), 1.0, 1.0);
	}

	void swapOneBreakpointRollback(MoveData beforeMove) {
		int right = 0;
		int left = 0;
		auto brkp1 = tree->get_node_event(beforeMove.node1);
		if (brkp1.first == beforeMove.old_label.first && brkp1.second == beforeMove.old_label.second) return;
		auto brkp2 = tree->get_node_event(beforeMove.node2);
		if (brkp1.first == beforeMove.old_label.first || brkp1.first == beforeMove.old_label.second) {
			left = 1;
		}
		if (brkp2.second == beforeMove.old_label.first || brkp2.second == beforeMove.old_label.second) {
			right = 1;
		}
		swapOneBreakpoint(beforeMove.node1, beforeMove.node2, left, right);
	}

	void moveRollBack(MoveType type, MoveData &beforeMove) {
		switch (type) {
		case ADD_LEAF:
			deleteLeaf(beforeMove.node1);
			return;
		case DELETE_LEAF:
			addLeaf(beforeMove.oldParent, beforeMove.old_label);
			return;
		case CHANGE_LABEL:
			change_label(beforeMove.node1, beforeMove.old_label);
			return;
		case SWAP_LABELS:
			swapLabels(beforeMove.node1, beforeMove.node2);
			return;
		case PRUNE_REATTACH:
			prune_and_reattach(beforeMove.node1, beforeMove.oldParent);
			return;
		case SWAP_SUBTREES:
			if (beforeMove.simpleSwap) {
				swap_subtrees_non_descendants(beforeMove.node1, beforeMove.node2);
			}
			else {
				swap_subtrees_descendants(beforeMove.oldParent, beforeMove.node1, beforeMove.node2);
			}
			return;
		case SWAP_ONE_BREAKPOINT:
			swapOneBreakpointRollback(beforeMove);
			return;
		}
	}

	std::tuple<MoveData, Real_t, Real_t> doMove(MoveType type) {
		std::pair<size_t, size_t> breakpoints;
		switch (type) {
		case ADD_LEAF:
			return addLeafMove(label_sampler->sample_label(random));
		case DELETE_LEAF:
			return deleteLeafMove();
		case CHANGE_LABEL:
			return changeLabelMove();
		case SWAP_LABELS:
			return swapLabelsMove();
		case PRUNE_REATTACH:
			return pruneAndReattachMove();
		case SWAP_SUBTREES:
			return swapSubtreesMove();
		case SWAP_ONE_BREAKPOINT:
		default:
			return swapOneBreakpointMove();
		}
	}

	Real_t getPriorOfReverseMove(MoveType type) {
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

	Real_t getLogTreePrior() {
		if (!label_sampler->has_free_labels()) {
			return -1000000.0;
		}
		Real_t C = std::log((Real_t)tree->get_size()) - label_sampler->get_sample_label_log_kernel() + nodeSet->getDeleteLeafKernel();
		auto breakpoints = tree->get_all_events();
		Real_t eventsLengths = 0.0;
		std::for_each(breakpoints.begin(), breakpoints.end(), [&eventsLengths, this](Event brkp) { eventsLengths += this->cells->getEventLength(brkp); });
        return -C * (Real_t)tree->get_size() - EVENTS_LENGTH_PENALTY * eventsLengths - DATA_SIZE_PRIOR_CONSTANT * ((Real_t)cells->getCellsCount()) * tree->get_size();
	}	

	void move(MoveType type) {
		if (tree_count_score == 1.0) {
			tree_count_score = countsScoring.calculate_log_score(*tree, likelihoodCalculator->getBestAttachment());;
		}
		auto beforeMoveLikelihood = temperature * likelihoodCalculator->getLikelihood();
		auto beforeMoveTreePrior = getLogTreePrior() ;
		auto moveData = doMove(type);
		auto afterMoveLikelihood = likelihoodCalculator->get_after_move_likelihood();
		afterMoveLikelihood = temperature * afterMoveLikelihood;
		auto after_move_counts_score = countsScoring.calculate_log_score(*tree, likelihoodCalculator->get_after_move_best_attachment());
		if (DEBUG ) {
			logDebug("Likelihood after move: ", afterMoveLikelihood);
			logDebug("Move kernel value: ", std::to_string(std::get<1>(moveData)));
			logDebug("Reversed move kernel value: ", std::to_string(std::get<2>(moveData)));
		}

		Real_t logAcceptance = 
			afterMoveLikelihood - beforeMoveLikelihood 
			+ (after_move_counts_score - tree_count_score)
			- std::get<1>(moveData) + std::get<2>(moveData) 
			- std::log(getPriorOfReverseMove(type)) + std::log(moveProbability[type]) 
			+ getLogTreePrior() - beforeMoveTreePrior;

		if (DEBUG ) {
			logDebug("Log acceptance ratio: ", logAcceptance);
		}
		if (random.logUniform() <= logAcceptance) {
			likelihoodCalculator->swapLikelihoodWithTmpAccumulator();
			tree_count_score = after_move_counts_score;
			if (DEBUG) {
				logDebug("Move accepted");
			}
		}
		else {
			moveRollBack(type, std::get<0>(moveData));
			if (DEBUG) {
				log("Move rejected");
			}
		}
	}

	MoveType sampleMoveType() {
		std::map<size_t, MoveType> typeToIndex;
		std::vector<Real_t> weights;
		for (auto x : moveProbability) {
			typeToIndex[weights.size()] = x.first;
			weights.push_back(x.second);
		}
		return typeToIndex[random.discrete(weights)];
	}

public:
	TreeSamplerCoordinator(EventTree *tree, LikelihoodCalculator<Real_t> *lC, unsigned int seed, size_t maxBrkp, VectorCellProvider<Real_t> *cells, 
		std::map<MoveType, Real_t> moveProbability, int pid): 
tree{ tree }, 
		likelihoodCalculator{ lC }, 
		label_sampler{ new VertexLabelSampler<Real_t>(maxBrkp, cells->getChromosomeMarkers()) },
				cells{ cells },
						countsScoring{ cells , COUNTS_SCORE_CONSTANT_0 !=0.0 && COUNTS_SCORE_CONSTANT_1 != 0.0 },
		random{ seed }, 
		moveProbability{ moveProbability },
		pid{ pid }  {
			
		auto breakpoints = std::move(tree->get_all_events());
		std::for_each(breakpoints.begin(), breakpoints.end(), [this](Event brkp) {this->label_sampler->add_label(brkp); });
		nodeSet = new TreeNodeSampler<Real_t>(*this->tree);
	}
	TreeSamplerCoordinator(TreeSamplerCoordinator &tsc) = default;

	TreeSamplerCoordinator<Real_t>& operator=(TreeSamplerCoordinator &tsc) = default;

	void recount_score() {
		tree_count_score = countsScoring.calculate_log_score(*tree, likelihoodCalculator->getBestAttachment());
	}

	Real_t get_likelihood()
	{
		return likelihoodCalculator->getLikelihood();
	}
	void set_temperature(Real_t temperature)
	{
		this->temperature = temperature;
	}
    
    Real_t get_total_likelihood() {
return likelihoodCalculator->getLikelihood() +getLogTreePrior(false) + tree_count_score;

    }
	Real_t get_temperature() const
	{
		return this->temperature;
	}

	bool moveIsPossible(MoveType type) {
		switch (type) {
		case ADD_LEAF:
			return label_sampler->has_free_labels();
		case CHANGE_LABEL:
			return label_sampler->has_free_labels();
		case DELETE_LEAF:
			return nodeSet->count_leaves() > 0 && tree->get_size() > 2;
		case SWAP_LABELS:
			return tree->get_size() >= 3;
		case PRUNE_REATTACH:
			return tree->get_size() >= 2;
		case SWAP_SUBTREES:
			return tree->get_size() >= 3;
		case SWAP_ONE_BREAKPOINT:
		default:
			return tree->get_size() >= 3;
		}
	}

	void treeMHStep() {
		MoveType type = sampleMoveType();
		if (DEBUG) {
			logDebug("Sampled move of type: ", moveTypeToString(type));
		}
		
		if (moveIsPossible(type)) {
			move(type);
		}
		else if (DEBUG) {
			logDebug("Move of type ", moveTypeToString(type), " is not possibles");
		}
		if (DEBUG) {
			logDebug("Tree size: ", std::to_string(tree->get_size()));
		}

		/** Statistics gathering**/
		if (moveCount == 0) {
			bestTree = *tree;
			bestLikelihoodNoPrior = likelihoodCalculator->getLikelihood();
			bestLikelihood = likelihoodCalculator->getLikelihood() +getLogTreePrior() + tree_count_score;
			bestTreeAttachment = likelihoodCalculator->getBestAttachment();
            best_counts_score = tree_count_score;
		}
		else {
			auto l = likelihoodCalculator->getLikelihood() +getLogTreePrior()+ tree_count_score;
			if (l > bestLikelihood) {
				bestTree = *tree;
				bestLikelihood = l;
				bestLikelihoodNoPrior = likelihoodCalculator->getLikelihood();
				bestTreeAttachment = likelihoodCalculator->getBestAttachment();
                best_counts_score = tree_count_score;
			}
		}
		moveCount++;
	}

	std::tuple<EventTree, std::vector<Event>, Real_t> getBestTreeData() {
		Attachment at{bestTreeAttachment};
	    bestTree.prune_tree(bestTree.get_root(), at);
		return std::make_tuple(bestTree, bestTreeAttachment, bestLikelihood);
	}


};

#endif // !TREE_SAMPLER_COORDINATOR_H
