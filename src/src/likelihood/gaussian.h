#ifndef GAUSSIAN_H
#define GAUSSIAN_H
#include <thread>
#include <utility>

#include "../parameters/parameters.h"
#include "../utils/random.h"
#include "adaptive_mh.h"
#include "gaussian_utils.h"

namespace Gauss {
template <class Real_t> class Gaussian {
public:
  Real_t mean;
  Real_t sd;

private:
  Random<Real_t> &random;
  AdaptiveMH<Real_t> adaptive_rw_var_mean;
  AdaptiveMH<Real_t> adaptive_rw_var_variance;

  void fill_log_likelihood_matrix_parallelized(
      std::vector<std::vector<Real_t>> &matrix,
      const std::vector<std::vector<Real_t>> &sample) const {
    std::vector<std::thread> threads;
    size_t rows_per_thread = matrix.size() / THREADS_LIKELIHOOD;
    for (size_t th = 0; th < THREADS_LIKELIHOOD; th++) {
      size_t right = th == THREADS_LIKELIHOOD - 1 ? matrix.size()
                                                  : rows_per_thread * (th + 1);
      threads.emplace_back(
          [&matrix, &sample, th, rows_per_thread, right, this] {
            for (size_t c = rows_per_thread * th; c < right; c++) {
              Gauss::truncated_gaussian_log_likelihood(matrix[c], sample[c],
                                                       this->mean, this->sd);
            }
          });
    }
    for (auto &th : threads) {
      th.join();
    }
  }

public:
  Gaussian(Real_t mean, Real_t sd, Random<Real_t> &random)
      : mean{mean}, sd{sd}, random{random} {}

  Gaussian<Real_t> &operator=(const Gaussian<Real_t> &g) {
    this->mean = g.mean;
    this->sd = g.sd;
    this->random = g.random;
    return *this;
  }

  void fill_log_likelihood_matrix(
      std::vector<std::vector<Real_t>> &matrix,
      const std::vector<std::vector<Real_t>> &sample) const {
    if (THREADS_LIKELIHOOD > 1) {
      fill_log_likelihood_matrix_parallelized(matrix, sample);
    } else {
      for (size_t c = 0; c < matrix.size(); c++) {
        Gauss::truncated_gaussian_log_likelihood(matrix[c], sample[c], mean,
                                                 sd);
      }
    }
  }

  Real_t get_parameters_prior() {
    return Gauss::truncated_gaussian_log_likelihood<Real_t>(mean, 0.0, 1.0) +
           Gauss::truncated_gaussian_log_likelihood<Real_t>(sd, 0.0, 1.0);
  }

  std::pair<Real_t, Real_t> resample_mean() {
    const Real_t step = adaptive_rw_var_mean.get(mean);
    mean += std::sqrt(step) * random.normal();
    return std::make_pair(1.0, 1.0);
  }

  std::pair<Real_t, Real_t> resample_standard_deviation() {
    const Real_t step = adaptive_rw_var_variance.get(sd);
    sd += std::sqrt(step) * random.normal();
    return std::make_pair(1.0, 1.0);
  }

  std::string to_string() {
    std::stringstream stream;
    stream << "mean: " << mean << " sd: " << sd << "\n";
    return stream.str();
  }
};
} // namespace Gauss

#endif
