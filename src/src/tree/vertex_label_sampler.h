#ifndef VERTEX_SAMPLER_H
#define VERTEX_SAMPLER_H
#include <algorithm>

#include "../utils/logger/logger.h"
#include "../utils/random.h"
#include "./utils/event_container.h"

/**
 * This class is responsible for sampling labels for new tree vertices
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
    return is_valid_event(brkp) && get_locus_chromosome(brkp.first) ==
                                       get_locus_chromosome(brkp.second);
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
  VertexLabelSampler(size_t max_loci, std::vector<size_t> chr_markers)
      : max_loci{max_loci}, unused_labels{max_loci}, chromosome_end_markers{
                                                         chr_markers} {
    init();
  }

  void add_label(TreeLabel l) { unused_labels.erase(l); }

  void remove_label(TreeLabel l) { unused_labels.insert(l); }

  Real_t get_sample_label_log_kernel() {
    return -std::log((Real_t)unused_labels.size());
  }

  TreeLabel sample_label(Random<Real_t> &random) {
    return unused_labels.get_nth(random.next_int(unused_labels.size()));
  }

  bool has_free_labels() { return !unused_labels.empty(); }

  std::pair<Event, Event> swap_one_breakpoint(Event ev1, Event ev2, int left,
                                              int right) {
    auto new_events = swap_breakpoints(ev1, ev2, left, right);
    remove_label(ev1);
    remove_label(ev2);
    add_label(new_events.first);
    add_label(new_events.second);
    return new_events;
  }

  bool can_swap_one_breakpoint(Event brkp1, Event brkp2, int left, int right) {
    auto newBrkps = swap_breakpoints(brkp1, brkp2, left, right);
    return is_valid_event(newBrkps.first) && is_valid_event(newBrkps.second) &&
           unused_labels.find(newBrkps.first) &&
           unused_labels.find(newBrkps.second);
  }
};

#endif // !VERTEX_SAMPLER_H
