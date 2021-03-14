#ifndef PAIR_CONTAINER_H
#define PAIR_CONTAINER_H

#include <vector>
#include "../breakpoints/breakpoints.h"
/**
* Container for storing ordered pairs of positive integers.
* Only pairs <code> (x,y) </code> where <code> x < y </code> are allowed. 
* Allows insertion, removal, lookup and access by index.
*
* Guarantees O(1) insertion, removal, existence check and O(<code>maxBrkp</code>) lookup by index. 
*
* Internally each pair is converted into an index and marked in a bitmap. 
* For fast element retrieval pairs with given first coordinate are counted.
*/
class PairContainer {
private:
	size_t maxBrkp;
	size_t size_{ 0 };
	std::vector<size_t> firstBrkps;
	std::vector<bool> brkpBitMap;

	size_t getBreakpointsIndex(BreakpointPair brkp) const {
		return brkp.first * (maxBrkp + 1) + brkp.second;
	}

	BreakpointPair fromIndex(size_t index) const {
		return std::make_pair(index / (maxBrkp + 1), index % (maxBrkp + 1));
	}

	void init() {
		brkpBitMap.resize((maxBrkp + 1) * (maxBrkp + 1));
		for (size_t i = 0; i < brkpBitMap.size(); i++) {
		  brkpBitMap[i] = false;
		}
		firstBrkps.resize(maxBrkp + 1);
		for (auto &el : firstBrkps) {
			el = 0;
		}
	}

	BreakpointPair getNth(size_t first, size_t index) const {
		size_t i = getBreakpointsIndex(std::make_pair(first, 0));
		size_t found = -1;
		while (found != index) {
			if (brkpBitMap[i] == true) found++;
			i++;
		}
		return fromIndex(i-1);
	}

public:
	/**
	* @param maxBrkp - maximal (inclusive) value for pair coordinates. 
	*/
	PairContainer(size_t maxBrkp) : maxBrkp{ maxBrkp } {
		init();
	}
	
	bool empty() const {
		return this->size_ == 0;
	}

	size_t size() const {
		return this->size_;
	}

	void insert(BreakpointPair brkp) {
		size_++;
		firstBrkps[brkp.first]++;
		brkpBitMap[getBreakpointsIndex(brkp)] = true;
	}

	void erase(BreakpointPair brkp) {
		size_--;
		firstBrkps[brkp.first]--;
		brkpBitMap[getBreakpointsIndex(brkp)] = false;
	}

	bool find(BreakpointPair brkp) {
		return brkpBitMap[getBreakpointsIndex(brkp)];
	}

	BreakpointPair getNth(size_t index) const {
		size_t first = 0;
		size_t cumSum = 0;
		while (cumSum + firstBrkps[first] <= index) {
			cumSum += firstBrkps[first];
			first++;
		}
		auto result = getNth(first, index - cumSum);
		return result;
	}

};

#endif // !PAIR_CONTAINER_H
