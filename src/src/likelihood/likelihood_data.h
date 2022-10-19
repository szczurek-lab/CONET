#ifndef LIKELIHOOD_DATA_H
#define LIKELIHOOD_DATA_H

#include "gaussian.h"
#include "gaussian_mixture.h"
#include "gaussian_utils.h"
#include "../utils/logger/logger.h"
/**
 * @brief Represents matrices of diffs likelihood for breakpoint and
 * no-breakpoint loci
 */
template <class Real_t> class LikelihoodData {
public:
  Gauss::Gaussian<Real_t> no_brkp_likelihood;
  Gauss::GaussianMixture<Real_t> brkp_likelihood;

  LikelihoodData(Gauss::Gaussian<Real_t> noBrkp,
                 Gauss::GaussianMixture<Real_t> mxt)
      : no_brkp_likelihood{noBrkp}, brkp_likelihood{mxt} {}

  void fill_no_breakpoint_log_likelihood_matrix(
      std::vector<std::vector<Real_t>> &matrix,
      const std::vector<std::vector<Real_t>> &corrected_counts) const {
    no_brkp_likelihood.fill_log_likelihood_matrix(matrix, corrected_counts);
  }

  void fill_breakpoint_log_likelihood_matrix(
      std::vector<std::vector<Real_t>> &matrix,
      const std::vector<std::vector<Real_t>> &corrected_counts) const {
    brkp_likelihood.fill_log_likelihood_matrix(matrix, corrected_counts);
  }

  Real_t get_likelihood_parameters_prior() {
    return brkp_likelihood.get_parameters_prior() +
           no_brkp_likelihood.get_parameters_prior() -
           Gauss::truncated_gaussian_log_likelihood<Real_t>(
               no_brkp_likelihood.mean, 0.0, 1.0);
  }

  bool likelihood_is_valid() {
    for (auto &g : brkp_likelihood.get_mixture_components()) {
      if (g.mean >= 0 || g.sd <= 0.0) {
        return false;
      }
    }
    if (no_brkp_likelihood.sd <= 0.0) {
      return false;
    }
    return true;
  }

  LikelihoodData remove_components_with_small_weight(const Real_t min_weight) {
    log("Removing components with weight smaller than ", min_weight);
    log("Mixture before change:\n", brkp_likelihood.to_string());
    brkp_likelihood.remove_components_with_small_weight(min_weight);
    log("Mixture after change:\n", brkp_likelihood.to_string());
    return *this;
  }
};

#endif // !LIKELIHOOD_DATA_H
