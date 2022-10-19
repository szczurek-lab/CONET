#ifndef TYPES_H
#define TYPES_H
#include <algorithm>
#include <cstddef>
#include <list>
#include <string>
#include <utility>

using Locus = size_t;
using Event = std::pair<Locus, Locus>;
using TreeLabel = Event;

std::pair<Event, Event> swap_breakpoints(Event brkp1, Event brkp2, int left,
                                         int right);

inline Locus get_event_start_locus(Event event) { return event.first; }

inline Locus get_event_end_locus(Event event) { return event.second; }

inline std::list<Locus> get_event_breakpoints(Event event) {
  return std::list<Locus>{event.first, event.second};
}

inline Event get_event_from_label(TreeLabel label) { return label; }

inline TreeLabel get_root_label() { return std::make_pair(0, 0); }

inline std::string label_to_str(TreeLabel label) {
  return "(" + std::to_string(label.first) + "," +
         std::to_string(label.second) + ")";
}

inline TreeLabel label_from_str(std::string l) {
  l.erase(std::remove(l.begin(), l.end(), '('), l.end());
  l.erase(std::remove(l.begin(), l.end(), ')'), l.end());
  return std::make_pair(
      (size_t)std::stoull(l.substr(0, l.find(','))),
      (size_t)std::stoull(l.substr(l.find(',') + 1, l.length())));
}

inline bool is_valid_event(const Event brkp) {
  return brkp.first < brkp.second;
}

inline bool is_root_event(const Event brkp) {
  return brkp.first == 0 && brkp.second == 0;
}

#endif // !TYPES_H
