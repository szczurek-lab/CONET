#ifndef TREE_SAMPLER_H
#define TREE_SAMPLER_H

#include "../utils/random.h"
#include "./tree_node_sampler.h"
#include "./vertex_label_sampler.h"
#include "event_tree.h"

template <class Real_t>
EventTree sample_tree(const size_t tree_size,
                      VertexLabelSampler<Real_t> &vertex_label_set,
                      Random<Real_t> &random) {
  EventTree tree;
  TreeNodeSampler<Real_t> node_sampler{tree};

  for (size_t i = 0; i < tree_size - 1 && vertex_label_set.has_free_labels();
       i++) {
    auto label = vertex_label_set.sample_label(random);
    auto parent = node_sampler.sample_node(true, random);
    node_sampler.refresh_node_data(tree.add_leaf(parent, label));
    vertex_label_set.add_label(label);
  }

  return tree;
}

#endif