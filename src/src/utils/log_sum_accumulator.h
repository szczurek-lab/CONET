#ifndef LOG_SUM_ACCUMULATOR_H
#define LOG_SUM_ACCUMULATOR_H

#include <cmath>

/**
 *	Iteratively calculates <code>log( exp(w_1) +...+ exp(w_n)) </code> for
 *input <code> w_1,..,w_n </code>
 */
template <class Real_t> class LogWeightAccumulator {
  Real_t max{0};
  Real_t sum{0};
  bool max_set{false};

public:
  void clear() {
    max = 0;
    sum = 0;
    max_set = false;
  }

  void add(const Real_t w) {
    if (!max_set) {
      max = w;
      sum += 1.0;
      max_set = true;
    } else if (w <= max) {
      sum += std::exp(w - max);
    } else {
      sum *= std::exp(max - w);
      sum += 1.0;
      max = w;
    }
  }

  Real_t get_result() const { return std::log(sum) + max; }
};

#endif // !LOG_SUM_ACCUMULATOR_H
