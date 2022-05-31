#ifndef EM_ESTIMATOR_H
#define  EM_ESTIMATOR_H

#include <vector>
#include <numeric>

#include "parameters/parameters.h"

template <class Real_t> class EMEstimator
{
	std::vector<Real_t> means;
	std::vector<Real_t> variances;
	std::vector<Real_t> weights;
	std::vector<Real_t> data;

	Real_t gaussian_density(Real_t arg, Real_t mean, Real_t variance)
	{
		static const Real_t inv_sqrt_2pi = 0.3989422804014327;
		Real_t a = (arg - mean)*(arg - mean) / variance;
		const Real_t result = inv_sqrt_2pi / std::sqrt(variance) * std::exp(-0.5 * a);
		if (result <= 0.0 || std::isnan(result))
		{
			return 0.0;
		}
		return result;
	}
	std::vector<std::vector<Real_t>> get_posteriors()
	{
		std::vector<std::vector<Real_t>> result(means.size());
		std::vector<Real_t> densities_sums(data.size());
		for (auto &el : densities_sums)
		{
			el = 0.0;
		}
		
		for (size_t i = 0; i < means.size(); i++)
		{
			std::vector<Real_t> densities(data.size());
			for (size_t arg = 0; arg < data.size(); arg ++)
			{
				densities[arg] = gaussian_density(data[arg], means[i], variances[i]) * weights[i];
				densities_sums[arg] += densities[arg];
			}
			result[i] = std::move(densities);
		}
		for (size_t i = 0; i < means.size(); i++)
		{
			for (size_t arg = 0; arg < data.size(); arg++)
			{
				result[i][arg] /= densities_sums[arg];
			}
		}
		return result;
	}

	void re_estimate(std::vector<std::vector<Real_t>> &probabilities) {
		for (size_t i = 0; i < means.size(); i ++) {
			const Real_t cluster_probs = std::accumulate(probabilities[i].begin(), probabilities[i].end(), 0.0);
			Real_t means_data = 0.0;
			Real_t variance_estimator = 0.0;
			for (size_t arg = 0; arg < data.size(); arg++)
			{
				means_data += probabilities[i][arg] * data[arg];
			}
			if ( false && i != 0) {
				means[i] = means_data / cluster_probs;
			}
			for (size_t arg = 0; arg < data.size(); arg++)
			{
				variance_estimator += probabilities[i][arg] * (data[arg] - means[i]) * (data[arg] - means[i]);
			}
			variances[i] = variance_estimator / cluster_probs;
			if (variances[i] < 0.000000001 || std::isnan(variances[i]))
			{
				variances[i] = 0.000000001;
			}
			weights[i] = cluster_probs / data.size();
		}
	}


public:
	EMEstimator<Real_t>(std::vector<Real_t> data): data{data} {}

	/*Real_t get_likelihood()
	{
		Real_t res = 0.0;
		for (size_t arg = 0; arg < data.size(); arg++)
		{
			Real_t l = 0.0;
			for (size_t i = 0; i < means.size(); i++)
			{
				l += gaussian_density(data[arg], means[i], variances[i]) * weights[i];
			}
			res += std::log(l);
		}
		return res;
	}*/
	
	void estimate(std::vector<Real_t> &means_, std::vector<Real_t> &variances_, std::vector<Real_t> &weights_)
	{
		this->means = means_;
		this->variances = variances_;
		this->weights = weights_;
		for (int i = 0; i < 4000; i++)
		{
			auto probabilities = get_posteriors();
			re_estimate(probabilities);
		}
		means_ = means;
		variances_ = variances;
		weights_ = weights;
	}
};











#endif //EM_ESTIMATOR_H
