#include "breakpoints.h"

namespace Breakpoint {
	/**
	* Swaps breakpoints between breakpoints pairs.
	* @param left - 0 if the first breakpoint from the @ref brkp1 should be swapped, 1 otherwise
	* @param right - just as for left.
	*/
	std::pair<BreakpointPair, BreakpointPair> swapBreakpoints(BreakpointPair brkp1, BreakpointPair brkp2, int left, int right) {
		BreakpointPair tmpBrkp1 = brkp1;
		if (left == 0) {
			brkp1.first = right == 0 ? brkp2.first : brkp2.second;
		}
		else {
			brkp1.second = right == 0 ? brkp2.first : brkp2.second;
		}
		if (right == 0) {
			brkp2.first = left == 0 ? tmpBrkp1.first : tmpBrkp1.second;
		}
		else {
			brkp2.second = left == 0 ? tmpBrkp1.first : tmpBrkp1.second;
		}
		if (brkp1.first > brkp1.second) {
			size_t tmp = brkp1.first;
			brkp1.first = brkp1.second;
			brkp1.second = tmp;
		}
		if (brkp2.first > brkp2.second) {
			size_t tmp = brkp2.first;
			brkp2.first = brkp2.second;
			brkp2.second = tmp;
		}
		return std::make_pair(brkp1, brkp2);
	}

	size_t Hash::operator() (const BreakpointPair &brkp) const {
		return brkp.first + (brkp.second << 12);
	}
}