#ifndef TREE_SAMPLER_COORDINATOR_H
#define TREE_SAMPLER_COORDINATOR_H
#include <algorithm>
#include <vector>
#include <tuple>

#include "tree/pointer_tree.h"
#include "utils/random.h"
#include "moves/move_type.h"
#include "utils/logger/logger.h"
#include "node_set/node_set.h"
#include "vertex_set/vertex_set_factory.h"
#include "moves/move_data.h"
#include "cell_provider/vector_cell_provider.h"
#include "likelihood_calculator.h"
#include "tree/tree_counts_scoring.h"

/**
* Class coordinating all components of MH sampling.
* This class is an owner of VertexSetInterface, NodeSet,
* PointerTree, LikelihoodCalculator and serves as a coordinating service for 
* their tasks.
*
* Tree MH sampling is done in few steps:\n 
* 1) First move type is sampled. \n 
* 2) NodeSet confirms that sampled move type is possible (for instance <code> DELETE_LEAF</code> move may not be possible \n 
* if the tree has zero size)\n 
* 3) Any required random elements are sampled (for example labels for a new node).\n 
* 4) VertexSet, PointerTree, NodeSet are ordered to do their move specific transformations\n 
* 5) If move has not been accepted a move roollback is carried out\n 
*
* 
*/
template <class Real_t> class TreeSamplerCoordinator {
public:
	using NodeHandle = PointerTree::NodeHandle;
	PointerTree *tree;
	LikelihoodCalculator<Real_t> *likelihoodCalculator;
	VertexSet<Real_t> *vertexSet;
	NodeSet<Real_t> *nodeSet;
	VectorCellProvider<Real_t> *cells;
	CountsScoring<Real_t> countsScoring;


	Random<Real_t> random;
	std::map<MoveType, Real_t> moveProbability;
	int pid;
	/*** for statistics gathering ***/
	Real_t bestLikelihood;
       Real_t bestLikelihoodNoPrior;	
    Real_t best_counts_score;
	PointerTree bestTree;
	std::vector<BreakpointPair> bestTreeAttachment;
	size_t moveCount{ 0 };
    std::map<std::pair<BreakpointPair, BreakpointPair>, size_t> edge_count;
	Real_t temperature {1.0};

	Real_t tree_count_score = 1.0;
			
	
	
	bool isLeaf(NodeHandle node) {
		return tree->isLeaf(node);
	}

	std::pair<size_t, size_t> sampleBreakpointsForChangeLabel() {
		return vertexSet->sampleBreakpointsForChangeLabel();
	}

	void pruneAndReattach(NodeHandle nodeToPrune, NodeHandle nodeToAttach) {
		auto oldParent = tree->pruneAndReattach(nodeToPrune, nodeToAttach);
		nodeSet->pruneAndReattachUpdate(nodeToAttach, oldParent);
	}

	void swapLabels(NodeHandle node1, NodeHandle node2) {
		tree->swapLabels(node1, node2);
	}

	/*
	  Returns handle to parent of @node
	*/
	NodeHandle deleteLeaf(NodeHandle node) {
		nodeSet->detachSubtree(node);
		nodeSet->eraseDetachedSubtree(node);
		vertexSet->removeBreakpoints(tree->getNodeBreakpoints(node));
		return tree->deleteLeaf(node);
	}

	NodeHandle addLeaf(NodeHandle parent, size_t leftBrkp, size_t rightBrkp) {
		NodeHandle newNode = tree->addLeaf(parent, leftBrkp, rightBrkp);
		nodeSet->addLeafUpdate(newNode, parent);
		vertexSet->addBreakpoints(std::make_pair(leftBrkp, rightBrkp));
		return newNode;
	}

	void changeLabel(NodeHandle node, std::pair<size_t, size_t> newBrkp) {
		vertexSet->removeBreakpoints(tree->getNodeBreakpoints(node));
		tree->changeLabel(node, newBrkp);
		vertexSet->addBreakpoints(newBrkp);
	}

	void swapSubtreesNonDescendants(NodeHandle root1, NodeHandle root2) {
		tree->swapSubtreesNonDescendants(root1, root2);
	}

	void swapSubtreesDescendants(NodeHandle parent, NodeHandle child, NodeHandle grandChild) {
		nodeSet->detachSubtree(child);
		tree->detachNode(child);

		auto newParent = tree->getParent(parent);
		nodeSet->detachSubtree(parent);
		tree->detachNode(parent);

		nodeSet->attachSubtree(newParent);
		tree->attachNode(child, newParent);

		nodeSet->attachSubtree(grandChild);
		tree->attachNode(parent, grandChild);
		tree->updateSubtreeBreakpoints(child);
	}

	void swapOneBreakpoint(NodeHandle node1, NodeHandle node2, int left, int right) {
		if (!vertexSet->swapOneBreakpointPossible(tree->getNodeBreakpoints(node1), tree->getNodeBreakpoints(node2), left, right)) {
			return;
		}
		auto newBrkp = vertexSet->swapOneBreakpoint(tree->getNodeBreakpoints(node1), tree->getNodeBreakpoints(node2), left, right);
		tree->changeLabel(node1, newBrkp.first);
		tree->changeLabel(node2, newBrkp.second);
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
		auto leaf = nodeSet->sampleLeaf();
		auto breakpoints = tree->getNodeBreakpoints(leaf);
		auto oldParent = deleteLeaf(leaf);
		beforeMove.oldParent = oldParent;
		beforeMove.oldLeftBrkp = breakpoints.first;
		beforeMove.oldRightBrkp = breakpoints.second;
		return std::make_tuple<MoveData, Real_t, Real_t>(std::move(beforeMove), std::move(kernelValue), nodeSet->getAddLeafKernel());
	}

	std::tuple<MoveData, Real_t, Real_t> addLeafMove(size_t leftBrkp, size_t rightBrkp) {
		MoveData beforeMove;
		Real_t kernelValue = nodeSet->getAddLeafKernel();
		auto parent = nodeSet->sampleNode(true);
		auto newNode = addLeaf(parent, leftBrkp, rightBrkp);
		beforeMove.node1 = newNode;
		return std::make_tuple<MoveData, Real_t, Real_t>(std::move(beforeMove), std::move(kernelValue), nodeSet->getDeleteLeafKernel());
	}
    
    NodeHandle nonLeafInsert(size_t leftBrkp, size_t rightBrkp, NodeHandle parent) {
        auto newNode = addLeaf(parent, leftBrkp, rightBrkp);
        std::list<NodeHandle> children = tree->getChildren(parent);
        for (auto child : children) {
            if (child != newNode) {
                tree->detachNode(child);
                tree->attachNode(child, newNode);
            }
        }
        tree->updateSubtreeBreakpoints(parent);
        nodeSet->nodeUpdate(newNode);
		vertexSet->addBreakpoints(std::make_pair(leftBrkp, rightBrkp));
        nodeSet->nodeUpdate(parent);
        return newNode;
    }

    void nonLeafDelete(NodeHandle toDelete) {
        vertexSet->removeBreakpoints(tree->getNodeBreakpoints(toDelete));
        nodeSet->removeNode(toDelete);
        auto parent = tree->getParent(toDelete);
        auto children = tree->getChildren(toDelete);
        for (auto ch : children) {
            tree->detachNode(ch);
            tree->attachNode(ch, parent);
        }
        auto brkp = tree->getNodeBreakpoints(toDelete);
        tree->deleteLeaf(toDelete);
        nodeSet->nodeUpdate(parent);
        tree->updateSubtreeBreakpoints(parent);
        vertexSet->removeBreakpoints(brkp);
    }

    std::tuple<MoveData, Real_t, Real_t> nonLeafInsertMove(size_t leftBrkp, size_t rightBrkp) {
		MoveData beforeMove;
		auto parent = nodeSet->sampleNode(true);
		beforeMove.node1 = nonLeafInsert(leftBrkp, rightBrkp, parent);
		return std::make_tuple<MoveData, Real_t, Real_t>(std::move(beforeMove), 1.0, 1.0);
	}
    
    std::tuple<MoveData, Real_t, Real_t> nonLeafDeleteMove() {
        MoveData beforeMove;
        beforeMove.oldParent = nullptr;
        auto nodes = tree->getNodes();
        std::vector<NodeHandle> oneChildNodes;
        for (auto node: nodes) {
            if (tree->getChildrenCount(node) == 1) oneChildNodes.push_back(node);
        }
        if (oneChildNodes.size() == 0) {
            return std::make_tuple<MoveData, Real_t, Real_t>(std::move(beforeMove), 1.0, 1.0);
        }
        auto chosen_node = oneChildNodes[random.nextInt(oneChildNodes.size())];
        auto child = tree->getChildren(chosen_node).front();
        beforeMove.oldParent = chosen_node;
		beforeMove.oldLeftBrkp = tree->getNodeBreakpoints(child).first;
		beforeMove.oldRightBrkp = tree->getNodeBreakpoints(child).second;
        
        nonLeafDelete(child);
        return std::make_tuple<MoveData, Real_t, Real_t>(std::move(beforeMove), 1.0, 1.0);
    }

	std::tuple<MoveData, Real_t, Real_t> pruneAndReattachMove() {
		auto sample = nodeSet->sampleNodeAndNonDescendant();
		MoveData beforeMove;
		beforeMove.oldParent = tree->getParent(std::get<0>(sample));
		beforeMove.node1 = std::get<0>(sample);
		pruneAndReattach(std::get<0>(sample), std::get<1>(sample));
		return std::make_tuple<MoveData, Real_t, Real_t>(std::move(beforeMove), 1.0, 1.0);
	}

	std::tuple<MoveData, Real_t, Real_t> swapLabelsMove() {
		MoveData beforeMove;
		auto nodes = nodeSet->sampleTwoNodes();
		beforeMove.node1 = nodes.first;
		beforeMove.node2 = nodes.second;
		swapLabels(nodes.first, nodes.second);
		return std::make_tuple<MoveData, Real_t, Real_t>(std::move(beforeMove), 1.0, 1.0);
	}

	std::tuple<MoveData, Real_t, Real_t> changeLabelMove() {
		MoveData beforeMove;
		auto node = nodeSet->sampleNode(false);
		auto oldBreakpoints = tree->getNodeBreakpoints(node);
		auto breakpoints = sampleBreakpointsForChangeLabel();
		changeLabel(node, breakpoints);
		beforeMove.node1 = node;
		beforeMove.oldLeftBrkp = oldBreakpoints.first;
		beforeMove.oldRightBrkp = oldBreakpoints.second;
		return std::make_tuple<MoveData, Real_t, Real_t>(std::move(beforeMove), 1.0, 1.0);
	}

	std::tuple<MoveData, Real_t, Real_t> swapSubtreesDescendantsMove(NodeHandle parent, NodeHandle child) {
		MoveData beforeMove;
		auto grandChild = nodeSet->sampleDescendant(child);
		beforeMove.simpleSwap = false;
		beforeMove.oldParent = child;
		beforeMove.node1 = parent;
		beforeMove.node2 = tree->getParent(child);
		Real_t kernelValue = nodeSet->getSwapSubtreesDescendantsKernel(child);
		swapSubtreesDescendants(parent, child, grandChild);
		return std::make_tuple<MoveData, Real_t, Real_t>(std::move(beforeMove), std::move(kernelValue), nodeSet->getSwapSubtreesDescendantsKernel(parent));
	}

	std::tuple<MoveData, Real_t, Real_t> swapSubtreesNonDescendantsMove(NodeHandle node1, NodeHandle node2) {
		MoveData beforeMove;
		beforeMove.node1 = node1;
		beforeMove.node2 = node2;
		beforeMove.simpleSwap = true;
		swapSubtreesNonDescendants(node1, node2);
		return std::make_tuple<MoveData, Real_t, Real_t>(std::move(beforeMove), 1.0, 1.0);
	}

	std::tuple<MoveData, Real_t, Real_t> swapSubtreesMove() {
		auto nodes = nodeSet->sampleTwoNodes();
		int descendants = tree->areDirectDescendants(nodes.first, nodes.second);
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
		auto nodes = nodeSet->sampleTwoNodes();
		beforeMove.node1 = nodes.first;
		beforeMove.node2 = nodes.second;
		beforeMove.oldLeftBrkp = tree->getNodeBreakpoints(nodes.first).first;
		beforeMove.oldRightBrkp = tree->getNodeBreakpoints(nodes.first).second;
		swapOneBreakpoint(nodes.first, nodes.second, random.randomIntBit(), random.randomIntBit());
		return std::make_tuple<MoveData, Real_t, Real_t>(std::move(beforeMove), 1.0, 1.0);
	}

	void swapOneBreakpointRollback(MoveData beforeMove) {
		int right = 0;
		int left = 0;
		auto brkp1 = tree->getNodeBreakpoints(beforeMove.node1);
		if (brkp1.first == beforeMove.oldLeftBrkp && brkp1.second == beforeMove.oldRightBrkp) return;
		auto brkp2 = tree->getNodeBreakpoints(beforeMove.node2);
		if (brkp1.first == beforeMove.oldLeftBrkp || brkp1.first == beforeMove.oldRightBrkp) {
			left = 1;
		}
		if (brkp2.second == beforeMove.oldLeftBrkp || brkp2.second == beforeMove.oldRightBrkp) {
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
			addLeaf(beforeMove.oldParent, beforeMove.oldLeftBrkp, beforeMove.oldRightBrkp);
			return;
		case CHANGE_LABEL:
			changeLabel(beforeMove.node1, std::make_pair(beforeMove.oldLeftBrkp, beforeMove.oldRightBrkp));
			return;
		case SWAP_LABELS:
			swapLabels(beforeMove.node1, beforeMove.node2);
			return;
		case PRUNE_REATTACH:
			pruneAndReattach(beforeMove.node1, beforeMove.oldParent);
			return;
			return;
		case SWAP_SUBTREES:
			if (beforeMove.simpleSwap) {
				swapSubtreesNonDescendants(beforeMove.node1, beforeMove.node2);
			}
			else {
				swapSubtreesDescendants(beforeMove.oldParent, beforeMove.node1, beforeMove.node2);
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
			breakpoints = vertexSet->sampleBreakpoints();
			return addLeafMove(breakpoints.first, breakpoints.second);
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

	Real_t get_cell_event_length(bool tmp)
	{
		 return likelihoodCalculator->get_attachment_prior(tmp);
	}
	
	Real_t getLogTreePrior() {
		if (!vertexSet->existFreeLabels()) {
			return -1000000.0;
		}
		Real_t C = std::log((Real_t)tree->getSize()) - vertexSet->getSampleLabelKernel() + nodeSet->getDeleteLeafKernel();
		auto breakpoints = tree->getAllBreakpoints();
		Real_t eventsLengths = 0.0;
		std::for_each(breakpoints.begin(), breakpoints.end(), [&eventsLengths, this](BreakpointPair brkp) { eventsLengths += this->cells->getEventLength(brkp); });
        return -C * (Real_t)tree->getSize() - EVENTS_LENGTH_PENALTY * eventsLengths - DATA_SIZE_PRIOR_CONSTANT * ((Real_t)cells->getCellsCount()) * tree->getSize();
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
			+ COUNTS_SCORE_CONSTANT * (after_move_counts_score - tree_count_score)
			- std::get<1>(moveData) + std::get<2>(moveData) 
			- std::log(getPriorOfReverseMove(type)) + std::log(moveProbability[type]) 
			+ getLogTreePrior() - beforeMoveTreePrior;

		if (DEBUG ) {
			logDebug("Log acceptance ratio: ", logAcceptance);
			if (TEST) {
				nodeSet->checkSetIntegrity();
				tree->treeIntegrityCheck();
				vertexSet->checkSetIntegrity(tree->getAllBreakpoints());
			}
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
	TreeSamplerCoordinator(PointerTree *tree, LikelihoodCalculator<Real_t> *lC, unsigned int seed, size_t maxBrkp, VectorCellProvider<Real_t> *cells, 
		std::map<MoveType, Real_t> moveProbability, int pid): 
tree{ tree }, 
		likelihoodCalculator{ lC }, 
		vertexSet{ VertexSetNamespace::create<Real_t>(maxBrkp, cells->getChromosomeMarkers(), seed) },
				cells{ cells },
						countsScoring{ cells , COUNTS_SCORE_CONSTANT !=0.0 },
		random{ seed }, 
		moveProbability{ moveProbability },
		pid{ pid }  {
			
		auto breakpoints = std::move(tree->getAllBreakpoints());
		std::for_each(breakpoints.begin(), breakpoints.end(), [this](BreakpointPair brkp) {this->vertexSet->addBreakpoints(brkp); });
		nodeSet = new NodeSet<Real_t>(this->tree, random, vertexSet);
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
return likelihoodCalculator->getLikelihood() +getLogTreePrior(false) + COUNTS_SCORE_CONSTANT*tree_count_score;

    }
	Real_t get_temperature() const
	{
		return this->temperature;
	}
	
	void treeMHStep() {
		MoveType type = sampleMoveType();
		if (DEBUG) {
			logDebug("Sampled move of type: ", moveTypeToString(type));
		}
		
		if (nodeSet->moveIsPossible(type)) {
			move(type);
		}
		else if (DEBUG) {
			logDebug("Move of type ", moveTypeToString(type), " is not possibles");
		}
		if (DEBUG) {
			logDebug("Tree size: ", std::to_string(tree->getSize()));
			if (TEST) {
				nodeSet->checkSetIntegrity();
				tree->treeIntegrityCheck();
				vertexSet->checkSetIntegrity(tree->getAllBreakpoints());
			}
		}

		/** Statistics gathering**/
		if (moveCount == 0) {
			bestTree = *tree;
			bestLikelihoodNoPrior = likelihoodCalculator->getLikelihood();
			bestLikelihood = likelihoodCalculator->getLikelihood() +getLogTreePrior() + COUNTS_SCORE_CONSTANT*tree_count_score;
			bestTreeAttachment = likelihoodCalculator->getBestAttachment();
            best_counts_score = tree_count_score;
		}
		else {
			auto l = likelihoodCalculator->getLikelihood() +getLogTreePrior()+ COUNTS_SCORE_CONSTANT*tree_count_score;
			if (l > bestLikelihood) {
				bestTree = *tree;
				bestLikelihood = l;
				bestLikelihoodNoPrior = likelihoodCalculator->getLikelihood();
				bestTreeAttachment = likelihoodCalculator->getBestAttachment();
                best_counts_score = tree_count_score;
			}
		}
		moveCount++;
        if (moveCount >= BURNIN) {
            auto attach = likelihoodCalculator->getBestAttachment();
            auto edges = tree->gatherEdges(attach);
            for (auto edge : edges) {
                if (edge_count.find(edge) == edge_count.end()) {
                    edge_count[edge] = 0;
                }
                edge_count[edge]++;
            }
        }
	}

	std::tuple<PointerTree, std::vector<BreakpointPair>, Real_t> getBestTreeData() {
	    bestTree.pruneTree(bestTreeAttachment);
		return std::make_tuple(bestTree, bestTreeAttachment, bestLikelihood);
	}


};

#endif // !TREE_SAMPLER_COORDINATOR_H
