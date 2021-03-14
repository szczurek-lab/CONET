#ifndef GAUSSIAN_H
#define GAUSSIAN_H
#include <utility>
#include <thread>

#include "gaussian_utils.h"
#include "../../../utils/random.h"
#include "../../../parameters/parameters.h"
#include "adaptive_mh.h"

namespace Gauss {
	template <class Real_t> class Gaussian {
	public:
		Real_t mean;
		Real_t sd;
		Random<Real_t> &random;
		AdaptiveMH<Real_t> adaptiveMHMean; 
		AdaptiveMH<Real_t> adaptiveMHVariance; 

		void getLogLikelihoodParallel(std::vector<std::vector<Real_t>> &resultHolder, const std::vector<std::vector<Real_t>> &logCounts) const {
			std::vector<std::thread> threads;
			size_t perThread = resultHolder.size() / THREADS_LIKELIHOOD;
			for (size_t th = 0; th < THREADS_LIKELIHOOD; th++) {
				size_t right = th == THREADS_LIKELIHOOD - 1 ? resultHolder.size() : perThread * (th + 1);
				threads.emplace_back([&resultHolder, &logCounts, th, perThread, right, this] { 
					for (size_t c = perThread * th; c < right; c++) {
						Gauss::gaussianLogLikelihood(resultHolder[c], logCounts[c], this->mean, this->sd);
					}
				 });
			}
			for (auto &th : threads) {
				th.join();
			}
		}

	public:
		Gaussian(Real_t mean, Real_t sd, Random<Real_t> &random) : mean {mean}, sd {sd}, random {random} {}


		void getLogLikelihood(std::vector<std::vector<Real_t>> &resultHolder, const std::vector<std::vector<Real_t>> &logCounts) const {
			if (THREADS_LIKELIHOOD > 1) {
				getLogLikelihoodParallel(resultHolder, logCounts);
			} else {
				for (size_t c = 0; c < resultHolder.size(); c++) {
					Gauss::gaussianLogLikelihood(resultHolder[c], logCounts[c], mean, sd);
				}
			}
		}

		Real_t getParamsPrior() {
			return Gauss::gaussianLogLikelihood<Real_t>(mean, 0.0, 1.0) + Gauss::gaussianLogLikelihood<Real_t>(sd, 0.0, 1.0);
		}

		std::pair<Real_t, Real_t> resample() {
			return resampleVariance() + resampleMean();
		}

		std::pair<Real_t, Real_t> resampleMean() {
			std::pair<Real_t, Real_t> meanResult = std::make_pair(0.0, 0.0);
			const Real_t step = adaptiveMHMean.get(mean);
			mean += std::sqrt(step) * random.normal();
			return meanResult;
		}

		std::pair<Real_t, Real_t> resampleVariance() {
			const Real_t step = adaptiveMHVariance.get(sd);
			sd += std::sqrt(step) * random.normal();
			return std::make_pair(1.0, 1.0);
		}

		Gaussian(const Gaussian<Real_t> &g) : random{g.random} {
			this->mean = g.mean;
			this->sd = g.sd;
			//this->adaptiveMHMean = g.adaptiveMHMean;
			//this->adaptiveMHVariance = g.adaptiveMHVariance;
		}

		Gaussian<Real_t>& operator = (const Gaussian<Real_t> &g) {
			this->mean = g.mean;
			this->sd = g.sd;
			this->random = g.random;
			//this->adaptiveMHMean = g.adaptiveMHMean;
			//this->adaptiveMHVariance = g.adaptiveMHVariance;
			return *this;
		}
	};
}








#endif