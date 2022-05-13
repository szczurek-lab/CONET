#ifndef BASIC_VERTEX_SET_H
#define BASIC_VERTEX_SET_H
#include <algorithm>

#include "../utils/random.h"
#include "../utils/logger/logger.h"
#include "./utils/event_container.h"
/**
* Implementation of VertexSetInterface for mode in which labels with same breakpoints are allowed.
*
*/
template <class Real_t> class VertexLabelSampler {
private:
	Locus max_loci;
	EventContainer unused_labels;
	std::vector<Locus> chromosome_end_markers;

	size_t get_locus_chromosome(Locus locus) {
		size_t chromosome = 0;
		while (locus >= chromosome_end_markers[chromosome]) {
			chromosome++;
		}
		return chromosome;
	}

	bool is_valid_label(TreeLabel brkp) {
		return is_valid_event(brkp) && get_locus_chromosome(brkp.first) == get_locus_chromosome(brkp.second);
	}

	void init() {
		for (size_t brkp = 0; brkp <= max_loci; brkp++) {
			for (size_t brkp2 = brkp + 1; brkp2 <= max_loci; brkp2++) {
				if (is_valid_label(std::make_pair(brkp, brkp2))) {
					unused_labels.insert(std::make_pair(brkp, brkp2));
				}
			}
		}
	}

public:
	VertexLabelSampler(size_t maxBrkp, std::vector<size_t> chromosomeMarkers) : 
		max_loci{ maxBrkp }, 
		unused_labels{ maxBrkp }, 
		chromosome_end_markers{chromosomeMarkers} {
		init();
	}

	void add_label(TreeLabel brkp) {
		unused_labels.erase(brkp);
	}

	void remove_label(TreeLabel brkp) {
		unused_labels.insert(brkp);
	}

	Real_t get_sample_label_log_kernel() {
		return -std::log((Real_t)unused_labels.size());
	}

	TreeLabel sample_label(Random<Real_t> &random) {
		return unused_labels.get_nth(random.nextInt(unused_labels.size()));
	}

	bool has_free_labels() {
		return !unused_labels.empty();
	}

	std::pair<Event, Event> swapOneBreakpoint(Event brkp1, Event brkp2, int left, int right) {
		auto newBrkps = swapBreakpoints(brkp1, brkp2, left, right);
		remove_label(brkp1);
		remove_label(brkp2);
		add_label(newBrkps.first);
		add_label(newBrkps.second);
		return newBrkps;
	}

	bool swapOneBreakpointPossible(Event brkp1, Event brkp2, int left, int right) {
		auto newBrkps = swapBreakpoints(brkp1, brkp2, left, right);
		return is_valid_event(newBrkps.first) && is_valid_event(newBrkps.second)
			&& unused_labels.find(newBrkps.first)
			&& unused_labels.find(newBrkps.second);
	}
};

#endif // !BASIC_VERTEX_SET_H
