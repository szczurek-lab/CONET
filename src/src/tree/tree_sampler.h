#ifndef TREE_SAMPLER_H
#define TREE_SAMPLER_H


#include "event_tree.h"
#include "./vertex_label_sampler.h"
#include "./node_set.h"
#include "../utils/random.h"

template <class Real_t> EventTree sample_tree(const size_t tree_size, VertexLabelSampler<Real_t> &vertexSet, Random<Real_t> &random) {
	EventTree tree;
	TreeNodeSampler<Real_t> nodeSet{tree};

	for (size_t i = 0; i < tree_size - 1 && vertexSet.has_free_labels(); i++) {
		auto label = vertexSet.sample_label(random);
		auto parent = nodeSet.sample_node(true, random);
		nodeSet.refresh_node_data(tree.add_leaf(parent, label));
		vertexSet.add_label(label);
 	}

	return tree;
}

#endif