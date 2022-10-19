#ifndef EM_ESTIMATOR_H
#define EM_ESTIMATOR_H

#include <numeric>
#include <vector>

#include "../utils/matrix.h"
#include "../utils/random.h"
#include "./likelihood_data.h"
#include "../utils/logger/logger.h"

namespace Gauss {
/**
 * @brief EM estimation of likelihood distribution from corrected counts
 */
template <class Real_t> class EMEstimator {
  std::vector<Real_t> means;
  std::vector<Real_t> variances;
  std::vector<Real_t> weights;
  std::vector<Real_t> data;
  Random<Real_t> &random;

  const size_t EM_ITERS = 4000; // number of iters for the solver
  const Real_t MIN_COMPONENT_VARIANCE =
      0.00001; // minimal variance of mixture component

  std::vector<std::vector<Real_t>> get_component_membership_weights() {
    std::vector<std::vector<Real_t>> result =
        Matrix::create_2d_matrix(means.size(), data.size(), 0.0);
    std::vector<Real_t> densities_sums =
        Matrix::create_1d_matrix(data.size(), 0.0);

    for (size_t component = 0; component < means.size(); component++) {
      for (size_t arg = 0; arg < data.size(); arg++) {
        const Real_t density = gaussian_density(data[arg], means[component],
                                                variances[component]) *
                               weights[component];
        result[component][arg] = density;
        densities_sums[arg] += density;
      }
    }
    for (size_t component = 0; component < means.size(); component++) {
      for (size_t arg = 0; arg < data.size(); arg++) {
        result[component][arg] /= densities_sums[arg];
      }
    }
    return result;
  }

  void re_estimate(std::vector<std::vector<Real_t>> &probabilities) {
    for (size_t component = 0; component < means.size(); component++) {
      const Real_t cluster_probs =
          std::accumulate(probabilities[component].begin(),
                          probabilities[component].end(), 0.0);
      Real_t variance_estimator = 0.0;
      for (size_t arg = 0; arg < data.size(); arg++) {
        variance_estimator += probabilities[component][arg] *
                              (data[arg] - means[component]) *
                              (data[arg] - means[component]);
      }

      variances[component] = variance_estimator / cluster_probs;
      if (variances[component] < MIN_COMPONENT_VARIANCE ||
          std::isnan(variances[component])) {
        variances[component] = MIN_COMPONENT_VARIANCE;
      }
      weights[component] = cluster_probs / data.size();
    }
  }

  void initialize_starting_parameters(const size_t mixture_size) {
    means.clear();
    variances.clear();
    weights.clear();
    for (int i = 0; i < mixture_size; i++) {
      means.push_back(-i);
      variances.push_back(random.uniform());
      weights.push_back(random.uniform());
    }
    Matrix::normalize_1d_matrix_elements<Real_t>(weights);
  }

public:
  EMEstimator<Real_t>(std::vector<Real_t> data, Random<Real_t> &random)
      : data{data}, random{random} {}

  LikelihoodData<Real_t> estimate(const size_t mixture_size) {
    initialize_starting_parameters(mixture_size);

    for (int i = 0; i < EM_ITERS; i++) {
      auto probabilities = get_component_membership_weights();
      re_estimate(probabilities);
    }

    std::transform(variances.begin(), variances.end(), variances.begin(),
                   [](Real_t v) -> Real_t { return std::sqrt(v); });

    Gauss::Gaussian<Real_t> gaussian(0.0, variances[0], random);

    variances.erase(variances.begin());
    means.erase(means.begin());
    weights.erase(weights.begin());

    Matrix::normalize_1d_matrix_elements<Real_t>(weights);
    return LikelihoodData<Real_t>(
        gaussian,
        Gauss::GaussianMixture<Real_t>(weights, means, variances, random));
  }
};
} // namespace Gauss

#endif // EM_ESTIMATOR_H
