#ifndef RANDOM_H
#define RANDOM_H
#include <cassert>
#include <limits>
#include <random>

/**
 * Encapsulates all random services which may be used by any CONET components.
 */
template <class Real_t> class Random {
  std::mt19937 generator;

public:
  Random(unsigned int seed) : generator{seed} {}
  /**
   * Generates rundom number from <code>{0,.., bound-1}</code> uniformly
   */
  size_t next_int(size_t bound) {
    std::uniform_int_distribution<size_t> dist(0, bound - 1);
    return dist(generator);
  }

  size_t next_int() { return next_int(std::numeric_limits<size_t>::max()); }

  Real_t log_uniform() {
    std::uniform_real_distribution<Real_t> dist(0.0, 1.0);
    return std::log(dist(generator));
  }

  Real_t uniform() {
    std::uniform_real_distribution<Real_t> dist(0.0, 1.0);
    return dist(generator);
  }

  /**
   * Sample from standard gaussian distribution
   */
  Real_t normal() {
    std::normal_distribution<Real_t> dist(0.0, 1.0);
    return dist(generator);
  }

  size_t discrete(std::vector<Real_t> &weights) {
    std::discrete_distribution<size_t> distribution(weights.begin(),
                                                    weights.end());
    return distribution(generator);
  }

  int random_int_bit() { return (int)next_int(2); }
};

#endif // !RANDOM_H
