#ifndef EVENT_CONTAINER_H
#define EVENT_CONTAINER_H

#include <vector>
#include "../../types.h"
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
class EventContainer {
private:
	size_t max_locus;
	size_t size_{ 0 };
	std::vector<size_t> first_event_loci;
	std::vector<bool> event_bit_map;

	size_t get_event_index(Event ev) const {
		return ev.first * (max_locus + 1) + ev.second;
	}

	Event get_by_index(size_t index) const {
		return std::make_pair(index / (max_locus + 1), index % (max_locus + 1));
	}

	void init() {
		event_bit_map.resize((max_locus + 1) * (max_locus + 1));
		first_event_loci.resize(max_locus + 1);
		std::fill(event_bit_map.begin(), event_bit_map.end(), false);
		std::fill(first_event_loci.begin(), first_event_loci.end(), 0);
	}

	Event get_nth(Locus first_locus, size_t n) const {
		size_t i = get_event_index(std::make_pair(first_locus, 0));
		size_t found = -1;
		while (found != n) {
			if (event_bit_map[i] == true) {
				found++;
			}
			i++;
		}
		return get_by_index(i-1);
	}

public:
	/**
	* @param maxBrkp - maximal (inclusive) value for pair coordinates. 
	*/
	EventContainer(size_t max_locus) : max_locus{ max_locus } {
		init();
	}
	
	bool empty() const {
		return this->size_ == 0;
	}

	size_t size() const {
		return this->size_;
	}

	void insert(Event ev) {
		size_++;
		first_event_loci[ev.first]++;
		event_bit_map[get_event_index(ev)] = true;
	}

	void erase(Event brkp) {
		size_--;
		first_event_loci[brkp.first]--;
		event_bit_map[get_event_index(brkp)] = false;
	}

	bool find(Event brkp) {
		return event_bit_map[get_event_index(brkp)];
	}

	Event get_nth(size_t n) const {
		size_t first = 0;
		size_t cumSum = 0;
		while (cumSum + first_event_loci[first] <= n) {
			cumSum += first_event_loci[first];
			first++;
		}
		auto result = get_nth(first, n - cumSum);
		return result;
	}

};

#endif // !EVENT_CONTAINER_H
