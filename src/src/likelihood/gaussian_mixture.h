#ifndef GAUSSIAN_MIXTURE_H
#define GAUSSIAN_MIXTURE_H

#include <algorithm>
#include <numeric>
#include <sstream>

#include "../parameters/parameters.h"
#include "../utils/log_sum_accumulator.h"
#include "../utils/matrix.h"
#include "../utils/random.h"
#include "adaptive_mh.h"
#include "gaussian.h"
#include "gaussian_utils.h"

namespace Gauss {
template <class Real_t> class GaussianMixture {
  using vector_r = std::vector<Real_t>;

  vector_r log_weights;
  vector_r log_normalized_weights;
  std::vector<Gauss::Gaussian<Real_t>> components;
  std::vector<AdaptiveMH<Real_t>> rw_step_size_variances;
  Random<Real_t> &random;

  void fill_log_likelihood_matrix(std::vector<Real_t> &matrix,
                                  const std::vector<Real_t> &sample) const {
    std::vector<LogWeightAccumulator<Real_t>> acc(sample.size());

    for (size_t component = 0; component < components.size(); component++) {
      Gauss::truncated_gaussian_log_likelihood(
          matrix, sample, components[component].mean, components[component].sd);
      const Real_t weight = log_normalized_weights[component];
      for (size_t i = 0; i < sample.size(); i++) {
        acc[i].add(matrix[i] + weight);
      }
    }
    for (size_t i = 0; i < sample.size(); i++) {
      matrix[i] = acc[i].get_result();
    }
  }

  void recalculate_log_normalized_weights() {
    std::vector<Real_t> weights = log_weights;
    auto max_weight = std::max_element(log_weights.begin(), log_weights.end());
    for (auto &weight : weights) {
      weight -= *max_weight;
      weight = std::exp(weight);
    }
    Matrix::normalize_1d_matrix_elements(weights);
    for (size_t i = 0; i < components.size(); i++) {
      log_normalized_weights[i] =
          weights[i] <= 0.5 ? std::log(weights[i]) : std::log1p(weights[i] - 1);
    }
  }

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
              this->fill_log_likelihood_matrix(matrix[c], sample[c]);
            }
          });
    }
    for (auto &th : threads) {
      th.join();
    }
  }

  void erase_component(size_t component) {
    components.erase(components.begin() + component);
    log_weights.erase(log_weights.begin() + component);
    log_normalized_weights.erase(log_normalized_weights.begin() + component);
    rw_step_size_variances.erase(rw_step_size_variances.begin() + component);
    recalculate_log_normalized_weights();
  }

public:
  GaussianMixture(vector_r weights, vector_r means, vector_r sds,
                  Random<Real_t> &random)
      : random{random} {
    log_weights.resize(weights.size());
    log_normalized_weights.resize(weights.size());
    rw_step_size_variances.resize(weights.size());
    for (size_t i = 0; i < weights.size(); i++) {
      log_weights[i] = std::log(weights[i]);
      components.push_back(Gauss::Gaussian<Real_t>(means[i], sds[i], random));
    }
    recalculate_log_normalized_weights();
  }

  GaussianMixture<Real_t> &operator=(const GaussianMixture<Real_t> &g) {
    this->components = g.components;
    this->log_weights = g.log_weights;
    this->log_normalized_weights = g.log_normalized_weights;
    this->random = g.random;
    return *this;
  }

  size_t number_of_components() const { return components.size(); }

  std::pair<Real_t, Real_t> resample_component_weight(size_t component) {
    log_weights[component] += std::sqrt(rw_step_size_variances[component].get(
                                  log_weights[component])) *
                              random.normal();
    recalculate_log_normalized_weights();
    return std::make_pair(1.0, 1.0);
  }

  std::vector<Gaussian<Real_t>> get_mixture_components() { return components; }

  std::pair<Real_t, Real_t> resample_component_mean(size_t component) {
    return components[component].resample_mean();
  }

  std::pair<Real_t, Real_t> resample_component_sd(size_t component) {
    return components[component].resample_standard_deviation();
  }

  void remove_components_with_small_weight(const Real_t min_weight) {
    std::vector<Real_t> weights{log_normalized_weights};
    std::transform(log_normalized_weights.begin(), log_normalized_weights.end(),
                   weights.begin(), [](Real_t w) { return std::exp(w); });
    while (true) {
      size_t min_index = std::distance(
          weights.begin(), std::min_element(weights.begin(), weights.end()));
      if (weights[min_index] >= min_weight) {
        break;
      }
      erase_component(min_index);
      weights.erase(weights.begin() + min_index);
    }
  }

  Real_t get_parameters_prior() {
    Real_t result = 0.0;
    for (size_t i = 0; i < components.size(); i++) {
      result += components[i].get_parameters_prior() +
                Gauss::truncated_gaussian_log_likelihood<Real_t>(log_weights[i],
                                                                 0.0, 1.0);
    }
    return result;
  }

  void fill_log_likelihood_matrix(
      std::vector<std::vector<Real_t>> &matrix,
      const std::vector<std::vector<Real_t>> &sample) const {
    if (THREADS_LIKELIHOOD > 1) {
      fill_log_likelihood_matrix_parallelized(matrix, sample);
    } else {
      for (size_t c = 0; c < matrix.size(); c++) {
        fill_log_likelihood_matrix(matrix[c], sample[c]);
      }
    }
  }

  std::string to_string() {
    std::stringstream stream;

    for (size_t i = 0; i < log_weights.size(); i++) {
      stream << "(weight: " << std::exp(log_normalized_weights[i])
             << " mean: " << -components[i].mean << "sd: " << components[i].sd
             << ")\n";
    }
    return stream.str();
  }
};
} // namespace Gauss

#endif // !GAUSSIAN_MIXTURE_H
