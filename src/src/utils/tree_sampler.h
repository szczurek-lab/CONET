#ifndef TREE_SAMPLER_H
#define TREE_SAMPLER_H


#include "../tree/pointer_tree.h"
#include "../vertex_set/vertex_set.h"
#include "../node_set/node_set.h"
#include "random.h"

template <class Real_t> void sampleTree(PointerTree &tree, VertexSet<Real_t> *vertexSet, long seed) {
	Random<Real_t> random(seed);
	NodeSet<Real_t> *nodeSet = new NodeSet<Real_t>(&tree, random, vertexSet);
	size_t size = 4;
	for (size_t i = 0; i < size; i++) {
		auto brkp = vertexSet->sampleBreakpoints();
		auto node = nodeSet->sampleNode(true);
		auto newNode = tree.addLeaf(node, brkp.first, brkp.second);
		nodeSet->addNewNode(newNode);
		vertexSet->addBreakpoints(brkp);
		if (!vertexSet->existFreeLabels()) break;
 	}
	delete nodeSet;
}

#endif