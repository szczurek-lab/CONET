#ifndef NORMAL_MIXTURE_LIKELIHOOD_H
#define NORMAL_MIXTURE_LIKELIHOOD_H
#include <utility>
#include <fstream>

#include "gaussian_utils.h"
#include "gaussian_mixture.h"
#include "gaussian.h"


template <class Real_t> class LikelihoodData {
public:
	Gauss::Gaussian<Real_t> no_brkp_likelihood;
	Gauss::GaussianMixture<Real_t> brkp_likelihood;
public:

	LikelihoodData(Gauss::Gaussian<Real_t> noBrkp, Gauss::GaussianMixture<Real_t> mxt): no_brkp_likelihood{noBrkp}, brkp_likelihood {mxt}
	 {
	}
	
	void fill_no_breakpoint_log_likelihood_matrix(std::vector<std::vector<Real_t>> &matrix, const std::vector<std::vector<Real_t>> &corrected_counts) const {
		no_brkp_likelihood.fill_matrix_log_likelihood(matrix, corrected_counts);
	}

	void fill_breakpoint_log_likelihood_matrix(std::vector<std::vector<Real_t>> &matrix, const std::vector<std::vector<Real_t>> &corrected_counts) const {
		brkp_likelihood.fill_matrix_log_likelihood(matrix, corrected_counts);
	}
	
	Real_t get_likelihood_parameters_prior() {
		return brkp_likelihood.get_parameters_prior() + no_brkp_likelihood.get_parameters_prior() - Gauss::truncated_gaussian_log_likelihood<Real_t>(no_brkp_likelihood.mean, 0.0, 1.0);
	}

	bool likelihood_is_valid() {
		for (auto &g : brkp_likelihood.gaussians) {
			if (g.mean >= 0  || g.sd <= 0.0)  {
				return false;
			}
		}
		if (no_brkp_likelihood.sd <= 0) {
			return false;
		}
		return true;
	}

	std::string toString() {
		std::string res = "(" + std::to_string(no_brkp_likelihood.mean) + "," + std::to_string(no_brkp_likelihood.sd) + ")\n";
		for (size_t i = 0; i < brkp_likelihood.gaussians.size(); i++) { 
			auto g = brkp_likelihood.gaussians[i];

			res +=  "(" + std::to_string(brkp_likelihood.weights[i]) + " ," + std::to_string(g.mean) + "," + std::to_string(g.sd) + ")\n";
		}
		return res;
	}
};

#endif // !NORMAL_MIXTURE_LIKELIHOOD_H
