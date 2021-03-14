#ifndef VERTEX_SET_FACTORY_H
#define VERTEX_SET_FACTORY_H

#include "vertex_set.h"

namespace VertexSetNamespace {
	template<class Real_t> VertexSet<Real_t> *create(size_t maxBreakpoint, std::vector<size_t> chromosomeMarkers, long seed) {
			return new VertexSet<Real_t>(maxBreakpoint, chromosomeMarkers, seed);
	}
};
#endif // ! VERTEX_SET_FACTORY_H
