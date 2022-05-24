#ifndef POINTER_TREE_H
#define POINTER_TREE_H
#include <list>
#include <cassert>
#include <unordered_set>
#include <map>
#include <vector>
#include <string>
#include <set>
#include <algorithm>

#include "../utils/breakpoints/breakpoints.h"
#include "../utils/logger/logger.h"
#include "../cell_provider/vector_cell_provider.h"


class PointerTree {
public:
	class Node {
	public:
		Node *parent;
		std::list<Node*> children;
		BreakpointPair breakpoints;
		std::list<size_t> newBreakpoints;
	public:
		Node() : parent{ nullptr } {}
		~Node() {
			for (auto node : children) {
				delete node;
			}
		}
		friend class PointerTree;
	};
	typedef Node* NodeHandle;

public: 
	size_t size{ 1 };

	PointerTree::Node *root;

	using NodeVector = std::vector<NodeHandle>;

	/*
		Returns -1 if @node1 is a descendant of @node2
		1 if @node2 is a descendant of @node1
		0 otherwise
	*/
	int isDescendant(NodeHandle node1, NodeHandle node2) const {
		auto node = node1->parent;
		while (node != nullptr) {
			if (node2 == node) return -1;
			else node = node->parent;
		}
		node = node2->parent;
		while (node != nullptr) {
			if (node1 == node) return 1;
			else node = node->parent;
		}
		return 0;
	}

	void breakpointsSweep(NodeHandle node, std::unordered_set<size_t> &foundBreakpoints) {
		if (node != root) {
			node->newBreakpoints = std::list<size_t>{ node->breakpoints.first, node->breakpoints.second };
			node->newBreakpoints.remove_if([&foundBreakpoints](size_t br)
			{ return foundBreakpoints.find(br) != foundBreakpoints.end(); });
		}
		for (auto br : node->newBreakpoints) {
			foundBreakpoints.insert(br);
		}
		for (auto ch : node->children) {
			breakpointsSweep(ch, foundBreakpoints);
		}
		for (auto br : node->newBreakpoints) {
			foundBreakpoints.erase(br);
		}
	}

	void updateNewBreakpoints() {
		std::unordered_set<size_t> foundBreakpoints;
		breakpointsSweep(root, foundBreakpoints);
	}

	void gatherAncestorsBreakpoints(Node *node, std::unordered_set<size_t> &breakpoints) const {
		if (node && node != root) {
			breakpoints.insert(node->breakpoints.first);
			breakpoints.insert(node->breakpoints.second);
			gatherAncestorsBreakpoints(node->parent, breakpoints);
		}
	}

	NodeHandle nodeFactory(size_t leftBrkp, size_t rightBrkp) {
		Node *newNode = new Node();
		newNode->breakpoints = std::make_pair(leftBrkp, rightBrkp);
		newNode->newBreakpoints = std::list<size_t>{ leftBrkp, rightBrkp };
		return newNode;
	}

	void collectLeaves(NodeHandle node, NodeVector &leaves) const {
		if (isLeaf(node)) {
			leaves.push_back(node);
		}
		else {
			for (auto &child : node->children) {
				collectLeaves(child, leaves);
			}
		}
	}

	void collectNodes(NodeHandle node, NodeVector &nodes) const {
		if(node != root) nodes.push_back(node);
		for (auto &child : node->children) {
			collectNodes(child, nodes);
		}
	}

	void getNonDescendants(NodeHandle node, NodeHandle excludedChild, NodeVector &found) const {
		for (auto &child : node->children) {
			if (child != excludedChild) {
				collectNodes(child, found);
			}
		}
		found.push_back(node);
		if (node != root) {
			getNonDescendants(node->parent, node, found);
		}
	}

	/*
		@newNode gets all children of @oldParent;
		assumes that @newNode->children.size() == 0 
		@oldParent becomes parent of @node
	*/
	void takeOverChildren(NodeHandle newNode, NodeHandle oldParent) {
		newNode->children = std::move(oldParent->children);
		for (auto &ch : newNode->children) {
			ch->parent = newNode;
		}
		oldParent->children = std::list<Node*>{ newNode };
		newNode->parent = oldParent;
	}
	/*
		@parent takes all children of @child
	*/
	void takeOverChildrenFromChild(NodeHandle parent, NodeHandle child) {
		detachNode(child);
		while (child->children.size() > 0) {
			auto grandChild = child->children.front();
			detachNode(grandChild);
			attachNode(grandChild, parent);
		}
	}

	void collectBreakpoints(NodeHandle node, std::vector<BreakpointPair> &breakpoints) const {
		breakpoints.push_back(node->breakpoints);
		for (auto child : node->children) {
			collectBreakpoints(child, breakpoints);
		}
	}

	//Debug only
	bool checkTreeIntegrity(NodeHandle node, size_t &found) {
		found += node->children.size();
		if (node->children.size() > 100000) return false;
		if (node != root && !Breakpoint::isValid(node->breakpoints)) { 
			logDebug("Reversed breakpoint pair!"); 
			return false; 
		}
		bool result = true;
		for (auto ch : node->children) {
			result = result && checkTreeIntegrity(ch, found);
		}
		return result;
	}

	void copyOtherTree(NodeHandle parent, NodeHandle treeRoot) {
		NodeHandle newNode = nodeFactory(treeRoot->breakpoints.first, treeRoot->breakpoints.second);
		attachNode(newNode, parent);
		for (auto child : treeRoot->children) {
			copyOtherTree(newNode, child);
		}
	}

	void toString(NodeHandle node, std::string &result) {
		for (auto &child : node->children) {
			result += get_node_label(node->breakpoints) + "-" +
				get_node_label(child->breakpoints) + "\n";
		}
		for (auto &child : node->children) {
			toString(child, result);
		}
	}
public:

	void gather_ancestors(NodeHandle node, std::set<NodeHandle> & result)
	{
		if (node == root) return;
		result.insert(node);
		gather_ancestors(node->parent, result);
	}
	PointerTree(const PointerTree &tree) {
		this->size = tree.size;
		root = new Node();
		for (auto child : tree.root->children) {
			copyOtherTree(root, child);
		}
		init();
	}

	PointerTree& operator=(const PointerTree &tree) {
		this->size = tree.size;
		if (root != nullptr) delete root;
		root = new Node();
		for (auto child : tree.root->children) {
			copyOtherTree(root, child);
		}
		init();
		return *this;
	}

	PointerTree() {
		root = new Node();
	}

	~PointerTree() {
		if (root != nullptr) delete root;
	}

	//Debug only
	void treeIntegrityCheck() {
		size_t found = 1;
		bool res = checkTreeIntegrity(root, found);
		if (found != this->size) {
			logDebug("Tree integrity check FAIL ", "sizes do not match ", this->size, " vs ", found);
			return;
		}
		if (!res) {
			logDebug("Tree integrity check FAIL while traversing the tree!");
			return;
		}
		logDebug("Tree integrity check OK!");
	}

	void init() {
		updateNewBreakpoints();
	}

	bool isLeaf(NodeHandle node) const {
		return node->children.size() == 0 && node != root;
	}

	void detachNode(NodeHandle node) {
		node->parent->children.remove(node);
		node->parent = nullptr;
	}

	void attachNode(NodeHandle node, NodeHandle newParent) {
		node->parent = newParent;
		newParent->children.push_back(node);
	}

	void updateSubtreeBreakpoints(NodeHandle subtreeRoot) {
		std::unordered_set<size_t> foundBreakpoints;
		if (subtreeRoot->parent != nullptr) gatherAncestorsBreakpoints(subtreeRoot->parent, foundBreakpoints);
		breakpointsSweep(subtreeRoot, foundBreakpoints);
	}

	NodeHandle postChild(NodeHandle node, size_t leftBrkp, size_t rightBrkp) {
		Node *newNode = nodeFactory(leftBrkp, rightBrkp);
		attachNode(newNode, node);
		updateSubtreeBreakpoints(newNode);
		size++;
		return newNode;
	}

	size_t getSize() const {
		return size;
	}
	
	NodeHandle getRoot() const {
		return root;
	}

	NodeHandle getParent(NodeHandle node) const {
		return node->parent;
	}

	const std::list<NodeHandle> getChildren(const NodeHandle node) const {
		return node->children;
	}

	size_t getChildrenCount(NodeHandle node) const {
		return node->children.size();
	}

	std::list<size_t> const &getNewBreakpoints(const NodeHandle node) const {
		return node->newBreakpoints;
	}

	BreakpointPair getNodeBreakpoints(const NodeHandle node) const {
		return node->breakpoints;
	}
	/***             MCMC moves            ***/


	/*
		Returns handle to new parent of subtree and boolean variable indicating whether
		subtree root was the only child of its parent (if it was, then the former parent becomes leaf)
	*/
	NodeHandle pruneAndReattach(NodeHandle nodeToPrune, NodeHandle nodeToAttach) {
		NodeHandle parent = nodeToPrune->parent;
		detachNode(nodeToPrune);
		attachNode(nodeToPrune, nodeToAttach);
		updateSubtreeBreakpoints(nodeToPrune);
		return parent;
	}

	void swapLabels(NodeHandle node1, NodeHandle node2) {
		std::swap(node1->breakpoints, node2->breakpoints);
		updateSubtreeBreakpoints(node1);
		updateSubtreeBreakpoints(node2);
	}

	NodeHandle deleteLeaf(NodeHandle node) {
		auto parent = node->parent;
		detachNode(node);
		size--;
		delete node;
		return parent;
	}

	NodeHandle addLeaf(NodeHandle parent, size_t leftBrkp, size_t rightBrkp) {
		return postChild(parent, leftBrkp, rightBrkp);
	}

	void changeLabel(NodeHandle node, BreakpointPair breakpoint) {
		node->breakpoints = breakpoint;
		updateSubtreeBreakpoints(node);
	}

	void swapSubtreesNonDescendants(NodeHandle root1, NodeHandle root2) {
		auto parent1 = root1->parent;
		auto parent2 = root2->parent;
		detachNode(root1);
		detachNode(root2);
		attachNode(root1, parent2);
		attachNode(root2, parent1);
		updateSubtreeBreakpoints(root1);
		updateSubtreeBreakpoints(root2);
	}

	void swapSubtreesDescendants(NodeHandle parent, NodeHandle child, NodeHandle grandChild) {
		detachNode(child);
		attachNode(child, parent->parent);
		detachNode(parent);
		attachNode(parent, grandChild);
		updateSubtreeBreakpoints(child);
	}
	/*  Coordinator services */

	NodeHandle getSibling(NodeHandle node) {
		for (auto ch : node->parent->children) {
			if (ch != node) return ch;
		}
		return nullptr;
	}

	NodeHandle getFirstChild(NodeHandle parent) {
		return parent->children.front();
	}

	std::vector<NodeHandle> getLeaves() const {
		std::vector<NodeHandle> leaves;
		collectLeaves(root, leaves);
		return leaves;
	}

	std::vector<NodeHandle> getNodes() const {
		std::vector<NodeHandle> nodes;
		collectNodes(root, nodes);
		return nodes;
	}

	std::unordered_set<size_t> getAncestorsBreakpoints(NodeHandle node) const {
		std::unordered_set<size_t> result;
		if (node != root) gatherAncestorsBreakpoints(node->parent, result);
		return result;
	}

	std::vector<NodeHandle> getNonDescendants(NodeHandle node) const {
		std::vector<NodeHandle> result;
		getNonDescendants(node->parent, node, result);
		return result;
	}

	int areDirectDescendants(NodeHandle node1, NodeHandle node2) const {
		return isDescendant(node1, node2);
	}


	/*
		Returns all nodes from subtree rooted at @node
	*/
	std::vector<NodeHandle> getDescendants(NodeHandle node) const {
		std::vector<NodeHandle> result;
		collectNodes(node, result);
		return result;
	}

	std::vector<BreakpointPair> getAllBreakpoints() const {
		std::vector<BreakpointPair> result;
		for (auto child : root->children) {
			collectBreakpoints(child, result);
		}
		return result;
	}


	bool sortTree(NodeHandle root, std::set<BreakpointPair> &attachment) {
		bool ret = false;
		if (root->children.size() == 1 && root->children.front()->breakpoints < root->breakpoints
			&& attachment.find(root->breakpoints) == attachment.end()
		 && attachment.find(root->children.front()->breakpoints) == attachment.end()) {
			ret = true;
			std::swap(root->breakpoints, root->children.front()->breakpoints);
		}
		for (auto child : root->children) {
			ret = ret || sortTree(child, attachment);
		}
		return ret;
	}

	bool sortTreeNoAttachment(NodeHandle root) {
		bool ret = false;
		if (root->children.size() == 1 && root->children.front()->breakpoints < root->breakpoints) {
			std::swap(root->breakpoints, root->children.front()->breakpoints);
			ret = true;
		}
		for (auto child : root->children) {
			ret = ret || sortTreeNoAttachment(child);
		}
		return ret;
	}
	
	bool isSingleChildParent(NodeHandle node)
	{
		return node->children.size() == 1;
	}
	/**
	* Recursively removes leaves which have no attached cell.
	*/
	void pruneTree(std::vector<BreakpointPair> attachment) {
		std::set<BreakpointPair> attach(attachment.begin(), attachment.end());
		auto leaves = getLeaves();
		bool erased = false;
		do {
			erased = false;
			std::for_each(leaves.begin(), leaves.end(), [this, &attach, &erased](NodeHandle n) { 
				if (attach.find(this->getNodeBreakpoints(n)) == attach.end()) {
					erased = true;
					this->deleteLeaf(n);
				}
			});
			leaves = getLeaves();
		} while (erased);
			while(sortTree(this->root, attach)){};
		this->init();
	}

	void gatherEdges(NodeHandle node, std::set<std::pair<BreakpointPair, BreakpointPair>> &edges) {
		for (auto child : node->children) {
			edges.insert(std::make_pair(node->breakpoints, child->breakpoints));
			gatherEdges(child, edges);
		}
	}

	std::set<std::pair<BreakpointPair, BreakpointPair>> gatherEdges() {
		std::set<std::pair<BreakpointPair, BreakpointPair>> result;
		gatherEdges(this->root, result);
		return result;
	}
    
   void gatherEdges(NodeHandle node, std::set<std::pair<BreakpointPair, BreakpointPair>> &edges, std::set<BreakpointPair> &attachment, std::set<BreakpointPair> &ancestors) {
        ancestors.insert(node->breakpoints);
        bool refresh = false;
        if (node->children.size() == 0 || node->children.size() > 1 || attachment.find(node->children.front()->breakpoints) != attachment.end() || attachment.find(node->breakpoints) != attachment.end()) {
            refresh = true;
        }
        if (refresh && ancestors.size() > 1) {
            std::vector<BreakpointPair> tmp{ancestors.begin(), ancestors.end()};
            for (size_t i = 0; i < tmp.size(); i ++) {
                for (size_t j = 0; j < tmp.size(); j++) {
                    if (i != j) {
                        edges.insert(std::make_pair(tmp[i], tmp[j]));
                    }
                }
            }
        }    
        if (refresh) {
            ancestors.clear();
        }
		for (auto child : node->children) {
			if (refresh) edges.insert(std::make_pair(node->breakpoints, child->breakpoints));
			gatherEdges(child, edges, attachment, ancestors);
		}
	}

    std::set<std::pair<BreakpointPair, BreakpointPair>> gatherEdges(std::vector<BreakpointPair> &attachment) {
        std::set<BreakpointPair> ancestors;
        std::set<BreakpointPair> att{attachment.begin(), attachment.end()};
		std::set<std::pair<BreakpointPair, BreakpointPair>> result;
		gatherEdges(this->root, result, att, ancestors);
		return result;
	}

	int isList() {
		NodeHandle currentNode = root;
		while (currentNode->children.size() == 1) {
			currentNode = currentNode->children.front();
		}
		if (currentNode->children.size() == 0) {
			return 1;
		}
			return 0;
	}

	/*** TMP for statostics grathering ***/
	void gatherEdges(NodeHandle node, std::vector<std::pair< BreakpointPair, BreakpointPair>> &result) {
		for (auto & child : node->children) {
			result.push_back(std::make_pair(node->breakpoints, child->breakpoints));
			gatherEdges(child, result);
		}
	}

	std::vector<std::pair< BreakpointPair, BreakpointPair>> getAllEdges() {
		std::vector<std::pair< BreakpointPair, BreakpointPair>> result;
		for (auto &child : root->children) {
			gatherEdges(child, result);
		}
		return result;
	}

	template<class Real_t> void get_event_lengths(NodeHandle node, Real_t length, std::map<BreakpointPair, Real_t> &result, VectorCellProvider<Real_t> const *cells)
	{
		if (node != root) {
			length += cells->getEventLength(node->breakpoints);
		}
		result[node->breakpoints] = length;
		for (auto child : node -> children)
		{
			get_event_lengths(child, length, result, cells);
		}
	}

	template<class Real_t> void get_normalized_event_lengths(NodeHandle node, Real_t length, std::map<BreakpointPair, Real_t> &result, VectorCellProvider<Real_t> const *cells, Real_t depth)
	{
		if (node != root) {
			length += cells->getEventLength(node->breakpoints);
		}
		if (node == root)
		{
			result[node->breakpoints] = 0.0;
		} else
		{
			result[node->breakpoints] = length / depth;
		}
		for (auto child : node->children)
		{
			get_normalized_event_lengths(child, length, result, cells, depth + 1);
		}
	}

	std::string get_node_label(BreakpointPair node)
	{
		return "(" + std::to_string(node.first) + ","
			+ std::to_string(node.second) + ")";
	}
	std::string toString() {
		std::string result = "";
		for (auto &child : root->children) {
			result += "(0,0)-" + get_node_label(child->breakpoints) + "\n";
		}
		for (auto &child : root->children) {
			toString(child, result);
		}
		return result;
	}

};
#endif // !POINTER_TREE_H
