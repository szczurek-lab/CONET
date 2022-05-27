#ifndef NORMAL_MIXTURE_LIKELIHOOD_H
#define NORMAL_MIXTURE_LIKELIHOOD_H
#include <utility>
#include <fstream>

#include "utils/gaussian_utils.h"
#include "utils/gaussian_mixture.h"
#include "utils/gaussian.h"

template <class Real_t> class NormalMixtureLikelihood {
public:
	Random<Real_t> random;
	Gauss::Gaussian<Real_t> noBrkpNormal;
	Gauss::GaussianMixture<Real_t> mixture;

	Gauss::Gaussian<Real_t> previousNoBrkpNormal;
	Gauss::GaussianMixture<Real_t> previousMixture;
	size_t step{ 0 };

	Gauss::Gaussian<Real_t> no_brkp_map;
	Gauss::GaussianMixture<Real_t> mixture_map;
public:

	Real_t count = 0;
	Real_t meanSd = 0;
	std::vector < std::vector<Real_t>> params;

	NormalMixtureLikelihood(Gauss::Gaussian<Real_t> noBrkp, Gauss::GaussianMixture<Real_t> mxt, unsigned int seed)
		:random{ seed }, noBrkpNormal{ noBrkp }, mixture{ mxt }, previousNoBrkpNormal{ noBrkp }, previousMixture{ mxt }, no_brkp_map{ noBrkp }, mixture_map{ mxt } {
		params.resize(3);
		for (int i = 0; i < 3; i++) {
			params[i].resize(mxt.size);
		}
	}
	std::pair<Gauss::GaussianMixture<Real_t>, Gauss::Gaussian<Real_t>> get_MAP()
	{
		return std::make_pair(mixture_map, no_brkp_map);
	}
	
	void getLogLikelihoodNoBreakpoint(std::vector<std::vector<Real_t>> &resultHolder, const std::vector<std::vector<Real_t>> &logCounts) const {
		noBrkpNormal.getLogLikelihood(resultHolder, logCounts);
	}

	void getLogLikelihoodBreakpoint(std::vector<std::vector<Real_t>> &resultHolder, const std::vector<std::vector<Real_t>> &logCounts) const {
		mixture.getLogLikelihood(resultHolder, logCounts);
	}
	
	Real_t getParamsPrior() {
		Real_t result = mixture.getParamsPrior() + noBrkpNormal.getParamsPrior();
		result -= Gauss::gaussianLogLikelihood<Real_t>(noBrkpNormal.mean, 0.0, 1.0);
		return result;
	}

	std::pair<Real_t, Real_t> resampleParametersGibbs() {
		step = (step + 1) % (3*mixture.size + 1);
		if (step == 0) {
			return noBrkpNormal.resampleVariance();
		}
		else {
			return mixture.resampleNth((step - 1)/3, (step - 1) % 3);
		}
	}

	bool valid() {
		for (auto &g : mixture.gaussians) {
			if (g.mean >= 0  || g.sd <= 0.0) return false;
		}
		if (noBrkpNormal.sd <= 0) return false;
		return true;
	}

	std::pair<Real_t, Real_t> resampleParameters() {
			return resampleParametersGibbs();
	}

	void rollBackParametersResample() {
		noBrkpNormal = previousNoBrkpNormal;
		mixture = previousMixture;
	}

	void acceptSampledParameters() {
		previousNoBrkpNormal = noBrkpNormal;
		previousMixture = mixture;
	}

	void update_map()
	{
		no_brkp_map = noBrkpNormal;
		mixture_map = mixture;
	}

	std::string toString() {
		std::string res = "(" + std::to_string(noBrkpNormal.mean) + "," + std::to_string(noBrkpNormal.sd) + ")\n";
		for (size_t i = 0; i < mixture.gaussians.size(); i++) { 
			auto g = mixture.gaussians[i];

			res +=  "(" + std::to_string(mixture.weights[i]) + " ," + std::to_string(g.mean) + "," + std::to_string(g.sd) + ")\n";
		}
		return res;
	}
};

#endif // !NORMAL_MIXTURE_LIKELIHOOD_H
