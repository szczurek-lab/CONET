#ifndef ADAPTIVE_MH_H
#define ADAPTIVE_MH_H

/**
 * Implementation of adaptive algorithm for scaling of random walk proposal's
 * variance.
 */
template <class Real_t> class AdaptiveMH {
  const Real_t epsilon = 0.0001;
  const Real_t initial_variance = 0.0001;
  const Real_t scaling = (1.0) * (1.0);
  const Real_t init_segment = 10;

  Real_t variance{initial_variance};
  Real_t average{0};
  Real_t num_observations{0};

public:
  AdaptiveMH<Real_t>() {}

  AdaptiveMH<Real_t> &operator=(const AdaptiveMH<Real_t> &g) {
    this->variance = g.variance;
    this->average = g.average;
    this->num_observations = g.num_observations;
    return *this;
  }

  Real_t get(Real_t x) {
    const Real_t old_average = average;
    average = average * num_observations / (num_observations + 1) +
              1 / (num_observations + 1) * x;
    if (num_observations == 1) {
      variance = x * x + old_average * old_average - 2 * average * average;
    } else if (num_observations > 0) {
      variance = (num_observations - 1) / (num_observations)*variance +
                 (1 / (num_observations)) *
                     (num_observations * old_average * old_average + x * x -
                      (num_observations + 1) * average * average);
    }
    num_observations++;

    return num_observations < init_segment
               ? initial_variance
               : scaling * variance + scaling * epsilon;
  }
};
#endif // ADAPTIVE_MH_H
