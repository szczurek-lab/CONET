#ifndef BREAKPOINTS_H
#define BREAKPOINTS_H
#include <utility>
#include <cstddef>

using BreakpointPair = std::pair<size_t, size_t>;

namespace Breakpoint {
	std::pair<BreakpointPair, BreakpointPair> swapBreakpoints(BreakpointPair brkp1, BreakpointPair brkp2, int left, int right);

	inline bool isValid(const BreakpointPair brkp) {
		return brkp.first < brkp.second;
	}

	inline bool isRoot(const BreakpointPair brkp) {
		return brkp.first == 0 && brkp.second == 0;
	}
	struct Hash {
		size_t operator() (const BreakpointPair &brkp) const;
	};

	inline bool operator< (const BreakpointPair &lhs, const BreakpointPair &rhs) {
		return lhs.first < rhs.first || (lhs.first == rhs.first && lhs.second < rhs.second);
	}
}

#endif // !BREAKPOINTS_H
