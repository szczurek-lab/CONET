#ifndef NODE_SET_H
#define NODE_SET_H

#include "../tree/pointer_tree.h"
#include "../vertex_set/vertex_set.h"
#include <algorithm>
#include <memory>

template <class Real_t> class NodeSet {
private:
	PointerTree *tree;
	Random<Real_t> &random;
	VertexSet<Real_t> *vertexSet;
	using NodeHandle = PointerTree::NodeHandle;
	using NodeVector = std::vector<NodeHandle>;
	NodeVector leaves;
	NodeVector nodes;

	void eraseFromVector(NodeVector &vec, NodeHandle node) const {
		vec.erase(std::remove_if(vec.begin(), vec.end(), [node](NodeHandle n) { return n == node; }));
	}

	size_t getNodeIndex(NodeVector &vec, NodeHandle node) const {
		auto it = std::find(vec.begin(), vec.end(), node);
		return std::distance(vec.begin(), it);
	}

	NodeHandle sample(NodeVector &vec, bool withRoot) const
	{
		size_t bound = withRoot ? vec.size() + 1 : vec.size();
		size_t node = random.nextInt(bound);
		return node == vec.size() ? tree->getRoot() : vec[node];
	}

	NodeHandle sampleExcludeIndex(NodeVector &vec, size_t index) const
	{
		size_t node = random.nextInt(vec.size() - 1);
		if (node >= index) node++;
		return vec[node];
	}

	bool find(NodeVector &vec, NodeHandle node) {
		return std::find(vec.begin(), vec.end(), node) != vec.end();
	}

public:
	NodeSet(PointerTree *tree, Random<Real_t> &random,VertexSet<Real_t> *vertexSet) : tree{ tree }, random{ random }, vertexSet{ vertexSet } {
		init();
	}

	void init() {
		leaves.clear();
		nodes.clear();
		auto nodes = std::move(tree->getNodes());
		for (auto node : nodes) {
			addNewNode(node);
		}
	}

	void set_tree(PointerTree *tree_) {
		this->tree = tree_;
	}
	NodeHandle sampleNode(bool withRoot) {
		return sample(nodes, withRoot);
	}

	std::tuple<NodeHandle, NodeHandle> sampleNodeAndNonDescendant() {
		auto node = sampleNode(false);
		auto nonDescendants = std::move(tree->getNonDescendants(node));
		auto nonDescendant = nonDescendants[random.nextInt(nonDescendants.size())];
		return std::make_tuple(node, nonDescendant);
	}

	NodeHandle sampleDescendant(NodeHandle node) {
		auto descendants = std::move(tree->getDescendants(node));
		return descendants[random.nextInt(descendants.size())];
	}

	NodeHandle sampleLeaf() {
		return leaves[random.nextInt(leaves.size())];
	}

	std::pair<NodeHandle, NodeHandle> sampleTwoNodes() {
		size_t firstNodeIndex = random.nextInt(nodes.size());
		return std::make_pair(nodes[firstNodeIndex], sampleExcludeIndex(nodes, firstNodeIndex));
	}


	void addNewNode(NodeHandle node) {
		nodes.push_back(node);
		if (tree->isLeaf(node)) leaves.push_back(node);
	}

	void removeNode(NodeHandle node) {
		eraseFromVector(nodes, node);
		if (tree->isLeaf(node)) eraseFromVector(leaves, node);
	}

	void nodeUpdate(const NodeHandle node) {
		if (node == nullptr || node == tree->getRoot()) return;
		if (!find(nodes, node)) nodes.push_back(node);
		if (tree->isLeaf(node) && !find(leaves, node)) leaves.push_back(node);
		if (!tree->isLeaf(node) && find(leaves, node)) eraseFromVector(leaves, node);
	}

	void eraseDetachedSubtree(NodeHandle subtreeRoot) {
		NodeVector subtreeNodes = std::move(tree->getDescendants(subtreeRoot));
		for (auto &node : subtreeNodes) {
			removeNode(node);
		}
	}

	void childLossUpdate(NodeHandle parent) {
		if (parent != tree->getRoot() && tree->isSingleChildParent(parent)) {
			leaves.push_back(parent);
		}
	}

	void childAcquisitionUpdate(const NodeHandle parent) {
		if (tree->isLeaf(parent)) {
			eraseFromVector(leaves, parent);
		}
	}

	void detachSubtree(const NodeHandle subtreeRoot) {
		childLossUpdate(tree->getParent(subtreeRoot));
	}

	void attachSubtree(NodeHandle nodeToAttach) {
		childAcquisitionUpdate(nodeToAttach);
	}

	void pruneAndReattachUpdate( NodeHandle nodeToAttach, NodeHandle oldParent) {
		nodeUpdate(oldParent);
		nodeUpdate(tree->getParent(oldParent));
		nodeUpdate(nodeToAttach);
		nodeUpdate(tree->getParent(nodeToAttach));
	}

	void addSubtreeUpdate(NodeHandle root, NodeHandle parent) {
		nodeUpdate(parent);
		nodeUpdate(tree->getParent(parent));
		NodeVector subtreeNodes = std::move(tree->getDescendants(root));
		std::for_each(subtreeNodes.begin(), subtreeNodes.end(), [this](NodeHandle n) {this->addNewNode(n); });
	}

	void addLeafUpdate(NodeHandle leaf, NodeHandle parent) {
		addNewNode(leaf);
		nodeUpdate(parent);
		nodeUpdate(tree->getParent(parent));
	}

	Real_t getAddLeafKernel() {
		return -std::log((Real_t)nodes.size() + 1) + vertexSet->getSampleLabelKernel();
	}

	Real_t getDeleteLeafKernel() {
		return -std::log((Real_t)leaves.size());
	}

	Real_t getSwapSubtreesDescendantsKernel(NodeHandle node) {
		return  -std::log((Real_t)tree->getDescendants(node).size());
	}

	bool moveIsPossible(MoveType type) {
		switch (type) {
		case ADD_LEAF:
			return vertexSet->existFreeLabels();
		case CHANGE_LABEL:
			return vertexSet->changeLabelIsPossible();
		case DELETE_LEAF:
			return !leaves.empty() && tree->getSize() > 2;
		case SWAP_LABELS:
			return nodes.size() >= 2;
		case PRUNE_REATTACH:
			return !nodes.empty();
		case SWAP_SUBTREES:
			return nodes.size() >= 2;
		case SWAP_ONE_BREAKPOINT:
		default:
			return nodes.size() >= 2;
		}
	}
	//Debug only
	bool checkSetIntegrity() {
		auto nodesReal = tree->getNodes();
		for (auto v : nodesReal) {
			if (std::find(nodes.begin(), nodes.end(), v) == nodes.end()) {
				logDebug("Real node not present in nodes DS, CHECK INTEGRITY FAIL");
				return false;
			}
			if (tree->isLeaf(v) && std::find(leaves.begin(), leaves.end(), v) == leaves.end()) {
				logDebug("Real leaf not present in leaves DS\nCHECK INTEGRITY FAIL\n");
				return false;
			}
		}

		for (auto n : nodes) {
			if (std::find(nodesReal.begin(), nodesReal.end(), n) == nodesReal.end()) {
				logDebug("Tree integrity check fail - nodes contains ghost nodes\n!");
				return false;
			}
		}
		for (auto n : leaves) {
			if (std::find(nodesReal.begin(), nodesReal.end(), n) == nodesReal.end()) {
				logDebug("Tree integrity check fail - leaves contains ghost nodes\n!");
				return false;
			}
			if (!tree->isLeaf(n)) {
				logDebug("Tree integrity check fail - leaf is not!\n");
				return false;
			}
		}
		logDebug("Node set integrity check OK!");
		return true;
	}


};

#endif // !NODE_SET_H

