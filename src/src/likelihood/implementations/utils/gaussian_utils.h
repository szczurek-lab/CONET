#ifndef  GAUSSIAN_UTILS_H
#define GAUSSIAN_UTILS_H
#include <vector>
#include <cmath>
#include "../../../parameters/parameters.h"

namespace Gauss {
	template <class Real_t> Real_t gaussianCDF(const Real_t x, const Real_t mean, const Real_t sd) // Phi(-∞, x) aka N(x)
	{
		const Real_t sqrt2 = 1.414214;
		return 1.0 - std::erfc((x - mean) / (sqrt2 * sd)) / 2;
	}
	template <class Real_t> void gaussianLogLikelihood(std::vector<Real_t> &resultHolder, const std::vector<Real_t> &args, const Real_t mean, const Real_t sd) {
		const Real_t log_inv_sqrt_2pi = -0.9189385;
		const Real_t logSd = std::log(sd);
		Real_t correction = 0.0;

			correction = std::log(gaussianCDF<Real_t>(0.0, mean, sd));

		for (size_t i = 0; i < args.size(); i++) {
			if (args[i] >= 765432.0) {
				resultHolder[i] = 0.0;
				continue;
			}
			Real_t res = -0.5*(args[i] - mean)*(args[i] - mean) / (sd*sd);
			resultHolder[i] = log_inv_sqrt_2pi - logSd + res - correction;
		}
	}

	template <class Real_t> Real_t gaussianLogLikelihood(const Real_t arg, const Real_t mean, const Real_t sd) {
		const Real_t log_inv_sqrt_2pi = -0.9189385;
		const Real_t logSd = std::log(sd);
		Real_t res = -0.5*(arg - mean)*(arg - mean) / (sd*sd);
		Real_t correction = 0.0;
			correction = std::log(gaussianCDF<Real_t>(0.0, mean, sd));
		return log_inv_sqrt_2pi - logSd + res - correction;
	}

	template <typename T, typename U> std::pair<T, U> operator+(const std::pair<T, U> & l, const std::pair<T, U> & r) {
		return { l.first + r.first,l.second + r.second };
	}
}
#endif // ! GAUSSIAN_UTILS_H
