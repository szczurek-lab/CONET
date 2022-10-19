#ifndef TREE_MH_STEPS_EXECUTOR_H
#define TREE_MH_STEPS_EXECUTOR_
#include <vector>

#include "input_data/input_data.h"
#include "likelihood_coordinator.h"
#include "moves/move_type.h"
#include "tree/event_tree.h"
#include "tree/tree_node_sampler.h"
#include "tree/vertex_label_sampler.h"
#include "utils/random.h"

/**
 * This class is responsible for execution of MCMC moves on the EventTree
 */
template <class Real_t> class MHStepsExecutor {
  using NodeHandle = EventTree::NodeHandle;

public:
  struct MoveData {
    /**
     * This class is used to persist move data so that move reversal is possible
     * only based on this class instance
     */
    std::map<std::string, NodeHandle> nodes;
    TreeLabel label;
    bool boolean_flag;
    Real_t move_log_kernel{1.0};
    Real_t reverse_move_log_kernel{1.0};

    MoveData()
        : label{get_root_label()}, boolean_flag{false}, move_log_kernel{1.0},
          reverse_move_log_kernel{1.0} {}
  };

private:
  EventTree &tree;
  VertexLabelSampler<Real_t> label_sampler;
  TreeNodeSampler<Real_t> node_sampler;
  CONETInputData<Real_t> &cells;
  Random<Real_t> &random;

  void prune_and_reattach(NodeHandle node_to_prune,
                          NodeHandle attachment_node) {
    auto former_node_attachment =
        tree.prune_and_reattach(node_to_prune, attachment_node);
    node_sampler.refresh_node_data(former_node_attachment);
    node_sampler.refresh_node_data(attachment_node);
  }

  void swap_labels(NodeHandle node1, NodeHandle node2) {
    tree.swap_labels(node1, node2);
  }

  NodeHandle delete_leaf(NodeHandle node) {
    node_sampler.delete_leaf(node);
    label_sampler.remove_label(tree.get_node_label(node));
    auto parent = tree.delete_leaf(node);
    node_sampler.refresh_node_data(parent);
    return parent;
  }

  NodeHandle add_leaf(NodeHandle attachment_node, TreeLabel label) {
    NodeHandle leaf = tree.add_leaf(attachment_node, label);
    node_sampler.refresh_node_data(leaf);
    node_sampler.refresh_node_data(attachment_node);
    label_sampler.add_label(label);
    return leaf;
  }

  void change_label(NodeHandle node, TreeLabel new_label) {
    label_sampler.remove_label(tree.get_node_label(node));
    tree.change_label(node, new_label);
    label_sampler.add_label(new_label);
  }

  void swap_subtrees_non_descendants(NodeHandle root1, NodeHandle root2) {
    tree.swap_subtrees_non_descendants(root1, root2);
    node_sampler.refresh_node_data(tree.get_parent(root1));
    node_sampler.refresh_node_data(tree.get_parent(root2));
  }

  void swap_subtrees_descendants(NodeHandle node, NodeHandle descendant,
                                 NodeHandle descendant_of_descendant) {
    auto descendant_attachment = tree.get_parent(descendant);
    tree.swap_subtrees_descendants(node, descendant, descendant_of_descendant);
    node_sampler.refresh_node_data(node);
    node_sampler.refresh_node_data(descendant_of_descendant);
    node_sampler.refresh_node_data(descendant);
    node_sampler.refresh_node_data(tree.get_parent(descendant));
    node_sampler.refresh_node_data(descendant_attachment);
  }

  void swap_breakpoints(NodeHandle node1, NodeHandle node2, int left,
                        int right) {
    if (!label_sampler.can_swap_one_breakpoint(tree.get_node_event(node1),
                                               tree.get_node_event(node2), left,
                                               right)) {
      return;
    }
    auto new_labels = label_sampler.swap_one_breakpoint(
        tree.get_node_event(node1), tree.get_node_event(node2), left, right);
    tree.change_label(node1, new_labels.first);
    tree.change_label(node2, new_labels.second);
  }

  MoveData delete_leaf_move() {
    MoveData move_data;
    move_data.move_log_kernel = node_sampler.get_delete_leaf_kernel();
    auto leaf = node_sampler.sample_leaf(random);
    move_data.label = tree.get_node_label(leaf);
    move_data.nodes["leaf_parent"] = delete_leaf(leaf);
    move_data.reverse_move_log_kernel =
        node_sampler.get_add_leaf_kernel() +
        label_sampler.get_sample_label_log_kernel();
    return move_data;
  }

  MoveData add_leaf_move() {
    MoveData move_data;
    move_data.move_log_kernel = node_sampler.get_add_leaf_kernel() +
                                label_sampler.get_sample_label_log_kernel();
    move_data.nodes["added_leaf"] =
        add_leaf(node_sampler.sample_node(true, random),
                 label_sampler.sample_label(random));
    move_data.reverse_move_log_kernel = node_sampler.get_delete_leaf_kernel();
    return move_data;
  }

  MoveData prune_and_reattach_move() {
    MoveData move_data;
    NodeHandle node_to_prune = node_sampler.sample_node(false, random);
    move_data.nodes["old_subtree_parent"] = tree.get_parent(node_to_prune);
    move_data.nodes["prunned_root"] = node_to_prune;
    prune_and_reattach(node_to_prune, node_sampler.sample_non_descendant(
                                          node_to_prune, random));
    return move_data;
  }

  MoveData swap_labels_move() {
    MoveData move_data;
    auto nodes = node_sampler.sample_nodes_pair(random);
    move_data.nodes["node1"] = nodes.first;
    move_data.nodes["node2"] = nodes.second;
    swap_labels(nodes.first, nodes.second);
    return move_data;
  }

  MoveData change_label_move() {
    MoveData move_data;
    auto node = node_sampler.sample_node(false, random);
    move_data.label = tree.get_node_label(node);
    move_data.nodes["node"] = node;
    change_label(node, label_sampler.sample_label(random));
    return move_data;
  }

  MoveData swap_subtrees_descendant_move(NodeHandle node,
                                         NodeHandle descendant) {
    MoveData move_data;
    move_data.move_log_kernel =
        node_sampler.get_swap_subtrees_descendants_kernel(descendant);

    auto descendant_of_descendant =
        node_sampler.sample_descendant(descendant, random);
    move_data.boolean_flag = false;
    move_data.nodes["middle_node"] = descendant;
    move_data.nodes["first_node"] = node;
    move_data.nodes["parent_of_middle_node"] = tree.get_parent(descendant);

    swap_subtrees_descendants(node, descendant, descendant_of_descendant);
    move_data.reverse_move_log_kernel =
        node_sampler.get_swap_subtrees_descendants_kernel(node);
    return move_data;
  }

  MoveData swap_subtrees_non_descendants_move(NodeHandle node1,
                                              NodeHandle node2) {
    MoveData move_data;
    move_data.nodes["node1"] = node1;
    move_data.nodes["node2"] = node2;
    move_data.boolean_flag = true;
    swap_subtrees_non_descendants(node1, node2);
    return move_data;
  }

  MoveData swap_subtrees_move() {
    auto nodes = node_sampler.sample_nodes_pair(random);
    int descendants = tree.get_nodes_relation(nodes.first, nodes.second);
    if (descendants != 0) {
      auto parent = descendants == -1 ? nodes.second : nodes.first;
      auto child = descendants == 1 ? nodes.second : nodes.first;
      return swap_subtrees_descendant_move(parent, child);
    } else {
      return swap_subtrees_non_descendants_move(nodes.first, nodes.second);
    }
  }
  MoveData swap_breakpoints_move() {
    MoveData move_data;
    auto nodes = node_sampler.sample_nodes_pair(random);
    move_data.nodes["node1"] = nodes.first;
    move_data.nodes["node2"] = nodes.second;
    move_data.label = tree.get_node_label(nodes.first);
    swap_breakpoints(nodes.first, nodes.second, random.random_int_bit(),
                     random.random_int_bit());
    return move_data;
  }

  void swap_breakpoints_rollback(MoveData move_data) {
    int right = 0, left = 0;
    auto label1 = tree.get_node_label(move_data.nodes["node1"]);
    if (label1.first == move_data.label.first &&
        label1.second == move_data.label.second) {
      return;
    }
    auto label2 = tree.get_node_label(move_data.nodes["node2"]);
    if (label1.first == move_data.label.first ||
        label1.first == move_data.label.second) {
      left = 1;
    }
    if (label2.second == move_data.label.first ||
        label2.second == move_data.label.second) {
      right = 1;
    }
    swap_breakpoints(move_data.nodes["node1"], move_data.nodes["node2"], left,
                     right);
  }

  Real_t get_total_events_length() {
    auto events = tree.get_all_events();
    Real_t events_length = 0.0;
    std::for_each(events.begin(), events.end(),
                  [&events_length, this](Event e) {
                    events_length += this->cells.get_event_length(e);
                  });
    return events_length;
  }

public:
  MHStepsExecutor<Real_t>(EventTree &t, CONETInputData<Real_t> &cells,
                          Random<Real_t> &r)
      : tree{t}, label_sampler{cells.get_loci_count() - 1,
                               cells.get_chromosome_end_markers()},
        node_sampler{tree}, cells{cells}, random{r} {
    for (auto event : tree.get_all_events()) {
      label_sampler.add_label(event);
    }
  }

  // reverse move execution
  void rollback_move(MoveType type, MoveData &move_data) {
    switch (type) {
    case ADD_LEAF:
      delete_leaf(move_data.nodes["added_leaf"]);
      return;
    case DELETE_LEAF:
      add_leaf(move_data.nodes["leaf_parent"], move_data.label);
      return;
    case CHANGE_LABEL:
      change_label(move_data.nodes["node"], move_data.label);
      return;
    case SWAP_LABELS:
      swap_labels(move_data.nodes["node1"], move_data.nodes["node2"]);
      return;
    case PRUNE_REATTACH:
      prune_and_reattach(move_data.nodes["prunned_root"],
                         move_data.nodes["old_subtree_parent"]);
      return;
    case SWAP_SUBTREES:
      if (move_data.boolean_flag) {
        swap_subtrees_non_descendants(move_data.nodes["node1"],
                                      move_data.nodes["node2"]);
      } else {
        swap_subtrees_descendants(move_data.nodes["middle_node"],
                                  move_data.nodes["first_node"],
                                  move_data.nodes["parent_of_middle_node"]);
      }
      return;
    case SWAP_ONE_BREAKPOINT:
      swap_breakpoints_rollback(move_data);
      return;
    }
  }

  MoveData execute_move(MoveType type) {
    switch (type) {
    case ADD_LEAF:
      return add_leaf_move();
    case DELETE_LEAF:
      return delete_leaf_move();
    case CHANGE_LABEL:
      return change_label_move();
    case SWAP_LABELS:
      return swap_labels_move();
    case PRUNE_REATTACH:
      return prune_and_reattach_move();
    case SWAP_SUBTREES:
      return swap_subtrees_move();
    case SWAP_ONE_BREAKPOINT:
    default:
      return swap_breakpoints_move();
    }
  }

  bool move_is_possible(MoveType type) {
    switch (type) {
    case ADD_LEAF:
    case CHANGE_LABEL:
      return label_sampler.has_free_labels();
    case DELETE_LEAF:
      return node_sampler.count_leaves() > 0 && tree.get_size() > 2;
    case PRUNE_REATTACH:
      return tree.get_size() >= 2;
    case SWAP_ONE_BREAKPOINT:
    case SWAP_SUBTREES:
    case SWAP_LABELS:
    default:
      return tree.get_size() >= 3;
    }
  }

  Real_t get_log_tree_prior() {
    const Real_t C = std::log((Real_t)tree.get_size()) -
                     label_sampler.get_sample_label_log_kernel() +
                     node_sampler.get_delete_leaf_kernel();
    return -C * (Real_t)tree.get_size() -
           EVENTS_LENGTH_PENALTY * get_total_events_length() -
           DATA_SIZE_PRIOR_CONSTANT * ((Real_t)cells.get_cells_count()) *
               tree.get_size();
  }
};
#endif