#ifndef NODE_SAMPLER_H
#define NODE_SAMPLER_H

#include "../tree/event_tree.h"
#include "../utils/logger/logger.h"
#include <algorithm>

/**
 * @brief This class is responsible for sampling of EventTree nodes.
 */
template <class Real_t> class TreeNodeSampler {
private:
  using NodeHandle = EventTree::NodeHandle;
  using NodeVector = std::vector<NodeHandle>;

  EventTree &tree;
  NodeVector leaves;
  NodeVector nodes;

  void erase_from_vector(NodeVector &vec, NodeHandle node) const {
    vec.erase(std::remove_if(vec.begin(), vec.end(),
                             [node](NodeHandle n) { return n == node; }));
  }

  bool find(NodeVector &vec, NodeHandle node) {
    return std::find(vec.begin(), vec.end(), node) != vec.end();
  }

  void remove_node(NodeHandle node) {
    erase_from_vector(nodes, node);
    if (tree.is_leaf(node)) {
      erase_from_vector(leaves, node);
    }
  }

public:
  TreeNodeSampler(EventTree &tree) : tree{tree} {
    for (auto node : tree.get_descendants(tree.get_root())) {
      if (node != tree.get_root()) {
        refresh_node_data(node);
      }
    }
  }

  NodeHandle sample_node(bool with_root, Random<Real_t> &random) {
    size_t bound = with_root ? nodes.size() + 1 : nodes.size();
    size_t node = random.next_int(bound);
    return node == nodes.size() ? tree.get_root() : nodes[node];
  }

  NodeHandle sample_leaf(Random<Real_t> &random) {
    return leaves[random.next_int(leaves.size())];
  }

  NodeHandle sample_non_descendant(NodeHandle node, Random<Real_t> &random) {
    auto nonDescendants = std::move(tree.get_non_descendants(node));
    return nonDescendants[random.next_int(nonDescendants.size())];
  }

  NodeHandle sample_descendant(NodeHandle node, Random<Real_t> &random) {
    auto descendants = std::move(tree.get_descendants(node));
    return descendants[random.next_int(descendants.size())];
  }

  std::pair<NodeHandle, NodeHandle> sample_nodes_pair(Random<Real_t> &random) {
    size_t first_node_idx = random.next_int(nodes.size());
    size_t second_node_idx = random.next_int(nodes.size() - 1);
    if (second_node_idx >= first_node_idx) {
      second_node_idx++;
    }
    return std::make_pair(nodes[first_node_idx], nodes[second_node_idx]);
  }

  void delete_leaf(NodeHandle node) {
    erase_from_vector(nodes, node);
    erase_from_vector(leaves, node);
  }

  void refresh_node_data(const NodeHandle node) {
    if (node == nullptr || node == tree.get_root()) {
      return;
    }
    if (!find(nodes, node)) {
      nodes.push_back(node);
    }
    if (tree.is_leaf(node) && !find(leaves, node)) {
      leaves.push_back(node);
    }
    if (!tree.is_leaf(node) && find(leaves, node)) {
      erase_from_vector(leaves, node);
    }
  }

  Real_t get_add_leaf_kernel() { return -std::log((Real_t)nodes.size() + 1); }

  Real_t get_delete_leaf_kernel() { return -std::log((Real_t)leaves.size()); }

  Real_t get_swap_subtrees_descendants_kernel(NodeHandle node) {
    return -std::log((Real_t)tree.get_descendants(node).size());
  }

  size_t count_leaves() { return leaves.size(); }
};

#endif // !NODE_SAMPLER_H
