#ifndef BASIC_VERTEX_SET_H
#define BASIC_VERTEX_SET_H
#include <algorithm>

#include "../utils/random.h"
#include "../utils/logger/logger.h"
#include "../utils/pair_container/pair_container.h"
/**
* Implementation of VertexSetInterface for mode in which labels with same breakpoints are allowed.
*
*/
template <class Real_t> class VertexSet {
private:
	size_t maxBreakpoint;
	PairContainer unused_breakpoints;
	std::vector<size_t> chromosomeMarkers;

	Random<Real_t> random;

	size_t getChromosome(size_t brkp) {
		size_t chromosome = 0;
		while (brkp >= chromosomeMarkers[chromosome]) {
			chromosome++;
		}
		return chromosome;
	}

	bool isValid(BreakpointPair brkp) {
		size_t chr1 = getChromosome(brkp.first);
		size_t chr2 = getChromosome(brkp.second);
		return Breakpoint::isValid(brkp)
			&& chr1 == chr2;
	}

	void moveBreakpointToUsedSet(BreakpointPair brkp) {
		unused_breakpoints.erase(brkp);
	}

	void moveBreakpointToUnusedSet(BreakpointPair brkp) {
		unused_breakpoints.insert(brkp);
	}

	void init() {
		for (size_t brkp = 0; brkp <= maxBreakpoint; brkp++) {
			for (size_t brkp2 = brkp + 1; brkp2 <= maxBreakpoint; brkp2++) {
				if (isValid(std::make_pair(brkp, brkp2))) {
					unused_breakpoints.insert(std::make_pair(brkp, brkp2));
				}
			}
		}
	}


public:
	VertexSet(size_t maxBrkp, std::vector<size_t> chromosomeMarkers, unsigned int seed) : 
		maxBreakpoint{ maxBrkp }, 
		unused_breakpoints{ maxBrkp }, 
		chromosomeMarkers{chromosomeMarkers},
		random{ seed } {
		init();
	}

	/* All breakpoints are marked as unused */
	void clear() {
		init();
	}

	void addBreakpoints(BreakpointPair breakpoint) {
		moveBreakpointToUsedSet(breakpoint);
	}

	void removeBreakpoints(BreakpointPair breakpoint) {
		moveBreakpointToUnusedSet(breakpoint);
	}

	Real_t getSampleLabelKernel() {
		return -std::log((Real_t)unused_breakpoints.size());
	}

	BreakpointPair sampleBreakpoints() {
		return unused_breakpoints.getNth(random.nextInt(unused_breakpoints.size()));
	}

	BreakpointPair sampleBreakpointsForChangeLabel() {
		return sampleBreakpoints();
	}

	bool existFreeLabels() {
		return !unused_breakpoints.empty();
	}

	bool changeLabelIsPossible() {
		return !unused_breakpoints.empty();
	}

	bool canAddBreakpoint(BreakpointPair brkp) {
		return unused_breakpoints.find(brkp);
	}

	std::pair<BreakpointPair, BreakpointPair> swapOneBreakpoint(BreakpointPair brkp1, BreakpointPair brkp2, int left, int right) {
		auto newBrkps = Breakpoint::swapBreakpoints(brkp1, brkp2, left, right);
		removeBreakpoints(brkp1);
		removeBreakpoints(brkp2);
		addBreakpoints(newBrkps.first);
		addBreakpoints(newBrkps.second);
		return newBrkps;
	}

	bool swapOneBreakpointPossible(BreakpointPair brkp1, BreakpointPair brkp2, int left, int right) {
		auto newBrkps = Breakpoint::swapBreakpoints(brkp1, brkp2, left, right);
		return isValid(newBrkps.first) && isValid(newBrkps.second)
			&& unused_breakpoints.find(newBrkps.first)
			&& unused_breakpoints.find(newBrkps.second);
	}

	void checkSetIntegrity(std::vector<BreakpointPair> breakpoints) {
		if (breakpoints.size() + unused_breakpoints.size() != ((1 + maxBreakpoint)*maxBreakpoint) / 2) {
			logDebug("Set integrity check FAIL - wrong number of unused breakpoints!");
			return;
		}
		for (size_t brkp = 0; brkp <= maxBreakpoint; brkp++) {
			for (size_t brkp2 = brkp + 1; brkp2 <= maxBreakpoint; brkp2++) {
				auto brkpPair = std::make_pair(brkp, brkp2);
				if (unused_breakpoints.find(brkpPair)) {
					if (!isValid(brkpPair)) {
						logDebug("Set integrity check FAIL - invalid breakpoint is used!");
						return;
					}
					auto it = std::find(breakpoints.begin(), breakpoints.end(), brkpPair);
					if (it != breakpoints.end()) {
						logDebug("Set integrity check FAIL - breakpoint is both used and unused!");
						return;
					}
				}
				else {
					auto it = std::find(breakpoints.begin(), breakpoints.end(), brkpPair);
					if (it == breakpoints.end()) {
						logDebug("Set integrity check FAIL - breakpoint is neither used nor unused!");
						return;
					}
				}
			}
		}
		logDebug("Vertex Set integrity check OK!");
	}
};

#endif // !BASIC_VERTEX_SET_H
