#include "types.h"
/**
* Swaps breakpoints between breakpoints pairs.
* @param left - 0 if the first breakpoint from the @ref brkp1 should be swapped, 1 otherwise
* @param right - just as for left.
*/
std::pair<Event, Event> swapBreakpoints(Event brkp1, Event brkp2, int left, int right) {
	Event tmpBrkp1 = brkp1;
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

