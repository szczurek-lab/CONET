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

#include "./attachment.h"
#include "../types.h"
#include "../utils/logger/logger.h"

/**
 * @brief 
 * 
 */
class EventTree
{
public:
	class Node
	{
	public:
		Node *parent;
		std::list<Node *> children;
		TreeLabel label;
		/**
		 * Loci which are breakpoints in this node, and are not present in any ancestor nodes
		 */
		std::list<Locus> new_breakpoints;
		Node() : parent{nullptr} {}
		Node(TreeLabel label, std::list<Locus> nb) : parent{nullptr}, children {},  label{label}, new_breakpoints{nb} {}
		~Node()
		{
			for (auto node : this->children)
			{
				delete node;
			}
		}
	};
	typedef Node *NodeHandle;

private:
	// Tree is always rooted
	size_t size{1};
	EventTree::Node *root;
	using NodeVector = std::vector<NodeHandle>;

	NodeHandle create_detached_node(TreeLabel label)
	{
		return new Node(label, get_event_breakpoints(get_event_from_label(label)));
	}

	void gather_ancestor_breakpoints(Node *node, std::unordered_set<Locus> &breakpoints) const
	{
		if (node != root)
		{
			for (Locus l : get_event_breakpoints(get_event_from_label(node->label)))
			{
				breakpoints.insert(l);
			}
			gather_ancestor_breakpoints(node->parent, breakpoints);
		}
	}

	void update_new_breakpoints(NodeHandle node, std::unordered_set<size_t> &breakpoints_from_ancestors)
	{
		if (node != root)
		{
			node->new_breakpoints = get_event_breakpoints(get_event_from_label(node->label));
			node->new_breakpoints.remove_if([&breakpoints_from_ancestors](Locus br)
											{ return breakpoints_from_ancestors.find(br) != breakpoints_from_ancestors.end(); });
		}
		breakpoints_from_ancestors.insert(node->new_breakpoints.begin(), node->new_breakpoints.end());
		
		for (auto ch : node->children)
		{
			update_new_breakpoints(ch, breakpoints_from_ancestors);
		}
		for (auto br : node->new_breakpoints)
		{
			breakpoints_from_ancestors.erase(br);
		}
	}

	void update_new_breakpoints(NodeHandle subtreeRoot)
	{
		std::unordered_set<Locus> found_breakpoints;
		if (subtreeRoot->parent != nullptr)
		{
			gather_ancestor_breakpoints(subtreeRoot->parent, found_breakpoints);
		}
		update_new_breakpoints(subtreeRoot, found_breakpoints);
	}

	void collect_subtree_nodes(NodeHandle node, NodeVector &nodes) const
	{
		nodes.push_back(node);
		for (auto &child : node->children)
		{
			collect_subtree_nodes(child, nodes);
		}
	}

	void detach_node(NodeHandle node)
	{
		node->parent->children.remove(node);
		node->parent = nullptr;
	}

	void attach_node(NodeHandle node, NodeHandle attach_to)
	{
		node->parent = attach_to;
		attach_to->children.push_back(node);
	}

	void copy_subtree(NodeHandle parent, NodeHandle treeRoot)
	{
		NodeHandle newNode = create_detached_node(treeRoot->label);
		attach_node(newNode, parent);
		for (auto child : treeRoot->children)
		{
			copy_subtree(newNode, child);
		}
	}

public:
	EventTree(const EventTree &tree)
	{
		this->size = tree.size;
		root = new Node();
		for (auto child : tree.root->children)
		{
			copy_subtree(root, child);
		}
		update_new_breakpoints(root);
	}

	EventTree &operator=(const EventTree &tree)
	{
		this->size = tree.size;
		if (root != nullptr)
		{
			delete root;
		}
		root = new Node();
		for (auto child : tree.root->children)
		{
			copy_subtree(root, child);
		}
		update_new_breakpoints(root);
		return *this;
	}

	EventTree()
	{
		root = new Node();
	}

	~EventTree()
	{
		if (root != nullptr)
		{
			delete root;
		}
	}

	bool is_leaf(NodeHandle node) const { return node->children.size() == 0 && node != root; }

	size_t get_size() const { return size; }

	NodeHandle get_root() const { return root;}

	NodeHandle get_parent(NodeHandle node) const { return node->parent;	}

	const std::list<NodeHandle> get_children(const NodeHandle node) const { return node->children; }

	std::list<Locus> const &get_new_breakpoints(const NodeHandle node) const { return node->new_breakpoints; }

	/**
	 * @brief Get event represented by node @node
	 */
	Event get_node_event(const NodeHandle node) const { return get_event_from_label(node->label); }

	TreeLabel get_node_label(const NodeHandle node) const { return node->label; }

	/*
		Returns -1 if @node1 is a descendant of @node2
		1 if @node2 is a descendant of @node1
		0 otherwise
	*/
	int get_nodes_relation(NodeHandle node1, NodeHandle node2) const
	{
		auto node = node1->parent;
		while (node != nullptr)
		{
			if (node2 == node)
				return -1;
			else
				node = node->parent;
		}
		node = node2->parent;
		while (node != nullptr)
		{
			if (node1 == node)
				return 1;
			else
				node = node->parent;
		}
		return 0;
	}

	/**
	 * @brief Returns nodes which are not descendants of node @node, root is included.  
	 * 
	 * @return std::vector<NodeHandle> - vector of nodes which are not contained in subtree rooted at @node
	 */
	std::vector<NodeHandle> get_non_descendants(NodeHandle node) const
	{
		std::vector<NodeHandle> all_nodes = get_descendants(root);
		std::list<NodeHandle> non_descendants{all_nodes.begin(), all_nodes.end()};
		for (auto d : get_descendants(node))
		{
			non_descendants.remove(d);
		}
		return std::vector<NodeHandle>(non_descendants.begin(), non_descendants.end());
	}

	/**
	 *	@brief Returns all nodes from subtree rooted at @node. The order is arbitrary
	*/
	std::vector<NodeHandle> get_descendants(NodeHandle node) const
	{
		std::vector<NodeHandle> result;
		collect_subtree_nodes(node, result);
		return result;
	}

	/**
	 * @brief Collect all events from the tree.
	 */
	std::vector<Event> get_all_events() const
	{
		auto nodes = get_descendants(root);
		nodes.erase(std::remove_if(nodes.begin(), nodes.end(), [this](NodeHandle n)
								   { return n == this->root; }));
		std::vector<Event> events;
		std::transform(nodes.begin(), nodes.end(), std::back_inserter(events), [this](NodeHandle n) -> Event
					   { return this->get_node_event(n); });
		return events;
	}

	/**
	 * @brief Recursively removes leaves of subtree rooted at @node which have no attached cell.
	 */
	void prune_tree(NodeHandle node, Attachment &attachment)
	{
		// Node children may be modified while traversing it, that's why we iterate on copy
		auto node_children_copy = node->children; 
		for (auto child : node_children_copy)
		{
			prune_tree(child, attachment);
		}
		if (node->children.empty() && !attachment.has_attached_cells(node->label))
		{
			this->delete_leaf(node);
		}
	}

	/**
	 * This section contains interface for performing MCMC moves on the tree.
	 */


	/**
	 * Adds new child to node @parent with label @label
	 *
	 * @param parent node to which new node will be attached
	 * @param label label for new node
	 * @return NodeHandle - handle to new node
	 */
	NodeHandle add_leaf(NodeHandle parent, TreeLabel label)
	{
		Node *newNode = create_detached_node(label);
		attach_node(newNode, parent);
		update_new_breakpoints(newNode);
		size++;
		return newNode;
	}

	/**
	 * @brief Detaches subtree rooted at @node_to_prune and attaches it to @node_to_attach
	 *
	 * @param node_to_prune - root of subtree which should be reattached
	 * @param node_to_attach - detached subtree will be attached to this node
	 * @return NodeHandle - old parent of reattached node
	 */
	NodeHandle prune_and_reattach(NodeHandle node_to_prune, NodeHandle node_to_attach)
	{
		NodeHandle parent = node_to_prune->parent;
		detach_node(node_to_prune);
		attach_node(node_to_prune, node_to_attach);
		update_new_breakpoints(node_to_prune);
		return parent;
	}

	/**
	 * @brief Swaps labels between two nodes
	 *
	 * @param node1 - will be given label of node @node2
	 * @param node2 - will be given label of node @node1
	 */
	void swapLabels(NodeHandle node1, NodeHandle node2)
	{
		std::swap(node1->label, node2->label);
		update_new_breakpoints(node1);
		update_new_breakpoints(node2);
	}

	/**
	 * @brief Deletes leaf of the tree
	 *
	 * @param node - leaf of the tree which should be deleted.
	 * @return NodeHandle - parent of the deleted leaf
	 */
	NodeHandle delete_leaf(NodeHandle node)
	{
		auto parent = node->parent;
		detach_node(node);
		size--;
		delete node;
		return parent;
	}

	/**
	 * @brief Changes label of node @node 
	 * Whether this will not result in node label duplication is up to the caller.  
	 */
	void change_label(NodeHandle node, Event new_label)
	{
		node->label = new_label;
		update_new_breakpoints(node);
	}

	/**
	 * @brief Swaps subtrees rooted at @root1, @root2
	 * 
	 * It is assumed that @root1 is not a descendant of @root2 and vice versa.
	 */
	void swap_subtrees_non_descendants(NodeHandle root1, NodeHandle root2)
	{
		auto parent1 = root1->parent;
		auto parent2 = root2->parent;
		detach_node(root1);
		detach_node(root2);
		attach_node(root1, parent2);
		attach_node(root2, parent1);
		update_new_breakpoints(root1);
		update_new_breakpoints(root2);
	}

	/**
	 * @brief Performs swap subtree move for descendant nodes
	 * It is assumed that @descendant is a descendant of @parent and @descendant_of_descendant is a descendant of @descendant. 
	 * Subtree rooted at @descendant will be attached to parent of @parent. The remaining subtree rooted at @parent will be attached to @descendant_of_descendant. 
	 */
	void swap_subtrees_descendants(NodeHandle parent, NodeHandle descendant, NodeHandle descendant_of_descendant)
	{
		detach_node(descendant);
		attach_node(descendant, parent->parent);
		detach_node(parent);
		attach_node(parent, descendant_of_descendant);
		update_new_breakpoints(descendant);
	}
};
#endif // !POINTER_TREE_H
