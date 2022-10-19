#ifndef UTILS_H
#define UTILS_H

#include <algorithm>
#include <optional>
#include <vector>

namespace Utils {
template <class T> std::vector<T> flatten(std::vector<std::vector<T>> vec) {
  std::vector<T> flattened;
  std::for_each(vec.begin(), vec.end(), [&flattened](std::vector<T> &v) {
    std::for_each(v.begin(), v.end(),
                  [&flattened](T e) { flattened.push_back(e); });
  });
  return flattened;
}

template <class T, class Real_t> class MaxValueAccumulator {
  Real_t value;
  std::optional<T> data;

public:
  MaxValueAccumulator() : value{0.0}, data{} {}

  void update(T p, Real_t v) {
    if (!data.has_value() || v > value) {
      data.emplace(p);
      value = v;
    }
  }

  T get() { return data.value(); }
};
} // namespace Utils

#endif