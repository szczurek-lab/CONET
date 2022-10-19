#ifndef GAUSSIAN_UTILS_H
#define GAUSSIAN_UTILS_H
#include <cmath>
#include <vector>

namespace Gauss {
template <class Real_t>
Real_t gaussian_density(Real_t arg, Real_t mean, Real_t variance) {
  static const Real_t inv_sqrt_2pi = 0.3989422804014327;
  Real_t a = (arg - mean) * (arg - mean) / variance;
  const Real_t result = inv_sqrt_2pi / std::sqrt(variance) * std::exp(-0.5 * a);
  if (std::isnan(result)) {
    return 0.0;
  }
  return result;
}

template <class Real_t>
Real_t gaussian_CDF(const Real_t x, const Real_t mean,
                    const Real_t sd) // Phi(-∞, x)
{
  const Real_t sqrt2 = 1.414214;
  return 1.0 - std::erfc((x - mean) / (sqrt2 * sd)) / 2;
}

template <class Real_t>
void truncated_gaussian_log_likelihood(std::vector<Real_t> &result,
                                       const std::vector<Real_t> &args,
                                       const Real_t mean, const Real_t sd) {
  const Real_t log_inv_sqrt_2pi = -0.9189385;
  const Real_t log_sd = std::log(sd);
  Real_t correction = std::log(gaussian_CDF<Real_t>(0.0, mean, sd));

  for (size_t i = 0; i < args.size(); i++) {
    Real_t res = -0.5 * (args[i] - mean) * (args[i] - mean) / (sd * sd);
    result[i] = log_inv_sqrt_2pi - log_sd + res - correction;
  }
}

template <class Real_t>
Real_t truncated_gaussian_log_likelihood(const Real_t arg, const Real_t mean,
                                         const Real_t sd) {
  const Real_t log_inv_sqrt_2pi = -0.9189385;
  Real_t res = -0.5 * (arg - mean) * (arg - mean) / (sd * sd);
  Real_t correction = std::log(gaussian_CDF<Real_t>(0.0, mean, sd));
  return log_inv_sqrt_2pi - std::log(sd) + res - correction;
}
} // namespace Gauss
#endif // ! GAUSSIAN_UTILS_H
