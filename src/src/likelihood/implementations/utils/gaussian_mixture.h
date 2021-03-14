#ifndef GAUSSIAN_MIXTURE_H
#define GAUSSIAN_MIXTURE_H

#include <algorithm>
#include <numeric>

#include "gaussian_utils.h"
#include "../../../utils/random.h"
#include "gaussian.h"
#include "../../../parameters/parameters.h"
#include "../../../utils/log_sum_accumulator.h"
#include "adaptive_mh.h"

namespace Gauss {
	template<class Real_t> class GaussianMixture {
		using vector_r = std::vector<Real_t>;
	public:
		size_t size;
		vector_r weights;
		vector_r logWeights;
		vector_r normalizedLogWeights;
		std::vector<Gauss::Gaussian<Real_t>> gaussians;
		std::vector<AdaptiveMH<Real_t>> stepsSizeWeightsRW;
		Random<Real_t> &random;

		void getLogLikelihoodVectorized(std::vector<Real_t> &resultHolder, const std::vector<Real_t> &args) const {
			std::vector<LogWeightAccumulator<Real_t>> acc(args.size());
			for (size_t i = 0; i < size; i++) {
				Gauss::gaussianLogLikelihood(resultHolder, args, gaussians[i].mean, gaussians[i].sd);
				const Real_t weight = normalizedLogWeights[i];
				for (size_t i = 0; i < args.size(); i++) {
					acc[i].add(resultHolder[i] + weight);
				}
			}
			for (size_t i = 0; i < args.size(); i++) {
				resultHolder[i] = acc[i].getResult();
			}
		}

		void weightsResamplingFixup() {
			auto maxWeight = std::max_element(logWeights.begin(), logWeights.end());
			for (auto &weight : weights) {
					weight -= *maxWeight;
					weight = std::exp(weight);
			}
			auto sum = std::accumulate(weights.begin(), weights.end(), 0.0);
			for (auto &weight : weights) {
				weight = weight / sum;
			}
			for (size_t i = 0; i < size; i++) {
				const Real_t weight = weights[i] <= 0.5 ? std::log(weights[i]) : std::log1p(weights[i] - 1);
				normalizedLogWeights[i] = weight;
			}
		}

		void prunningWeightsFixup()
		{
			auto sum = std::accumulate(weights.begin(), weights.end(), 0.0);
			for (auto &weight : weights)
			{
				weight = weight / sum;
			}
			for (size_t i = 0; i < weights.size(); i++)
			{
				logWeights[i] = std::log(weights[i]);
				normalizedLogWeights[i] = std::log(weights[i]);
			}
		}
		/**
		* exclusive 
		*/
		std::pair<Real_t, Real_t> resampleWeights(size_t left, size_t right) {
			for (size_t i = 0; i < weights.size(); i++) {
				if (i >= left && i < right) {
					logWeights[i] = logWeights[i] + std::sqrt(stepsSizeWeightsRW[i].get(logWeights[i]))* random.normal();
				}
				weights[i] = logWeights[i];
			}
			weightsResamplingFixup();
			return std::make_pair(1.0, 1.0);
		}

		void getLogLikelihoodParallel(std::vector<std::vector<Real_t>> &resultHolder, const std::vector<std::vector<Real_t>> &logCounts) const {
			std::vector<std::thread> threads;
			size_t perThread = resultHolder.size() / THREADS_LIKELIHOOD;
			for (size_t th = 0; th < THREADS_LIKELIHOOD; th++) {
				size_t right = th == THREADS_LIKELIHOOD - 1 ? resultHolder.size() : perThread * (th + 1);
				threads.emplace_back([&resultHolder, &logCounts, th, perThread, right, this] {
					for (size_t c = perThread * th; c < right; c++) {
						this->getLogLikelihoodVectorized(resultHolder[c], logCounts[c]);
					}
				});
			}
			for (auto &th : threads) {
				th.join();
			}

		}
	public:
		GaussianMixture(vector_r weights, vector_r means, vector_r sds, Random<Real_t> &random) : random{ random } {


			this->weights = weights;
			this->logWeights = std::vector<Real_t>(weights.size());
			this->stepsSizeWeightsRW.resize(weights.size());
			for (size_t i = 0; i < weights.size(); i++) {
				logWeights[i] = std::log(weights[i]);
			}
			for (size_t i = 0; i < means.size(); i++) {
				gaussians.push_back(Gauss::Gaussian<Real_t>(means[i], sds[i], random));
			}
			size = weights.size();
			for (size_t i = 0; i < size; i++) {
				const Real_t weight = weights[i] <= 0.5 ? std::log(weights[i]) : std::log1p(weights[i] - 1);
				normalizedLogWeights.push_back(weight);
			}
		}

		GaussianMixture(Random<Real_t> &random) : random {random}  {
			size = 0;
		}

		GaussianMixture(size_t size, Random<Real_t> &random) : size{size}, random { random } {
			this->stepsSizeWeightsRW.resize(size);
			for (size_t i = 0; i < size; i++) {
				gaussians.push_back(Gauss::Gaussian<Real_t>(random));
				weights.push_back(std::log(std::sqrt(0.1) * std::abs(random.normal())));
				logWeights.push_back(weights.back());
				normalizedLogWeights.push_back(weights.back());
			}
			weightsResamplingFixup();
		}

		void prune()
		{
			for(size_t i = 0; i < gaussians.size(); i++)
			{
				if (weights[i] <= 0.01)
				{
					gaussians.erase(gaussians.begin() + i);
					weights.erase(weights.begin() + i);
					logWeights.erase(logWeights.begin() + i);
					normalizedLogWeights.erase(normalizedLogWeights.begin() + i);
					stepsSizeWeightsRW.erase(stepsSizeWeightsRW.begin() + i);
					size--;
					prunningWeightsFixup();
					prune();
					return;
				}
			}
		}
		
		Real_t getBiggestMean() {
			Real_t res = gaussians[0].mean;
			for (size_t i = 1; i < gaussians.size(); i++) {
				if (gaussians[i].mean > res) {
					res = gaussians[i].mean;
				}
			}
			return res;

		}

		Real_t getParamsPrior() {
			Real_t result = 0.0;
			std::for_each(gaussians.begin(), gaussians.end(), [&result](Gaussian<Real_t> &n) { result += n.getParamsPrior(); });
			std::for_each(logWeights.begin(), logWeights.end(), [&result](Real_t &w) { result += Gauss::gaussianLogLikelihood<Real_t>(w, 0.0, 1.0); });
			return result;
		}

		void getLogLikelihood(std::vector<std::vector<Real_t>> &resultHolder, const std::vector<std::vector<Real_t>> &logCounts) const {
			if (THREADS_LIKELIHOOD > 1) {
				getLogLikelihoodParallel(resultHolder, logCounts);
			} else {
				for (size_t c = 0; c < resultHolder.size(); c++) {
					getLogLikelihoodVectorized(resultHolder[c], logCounts[c]);
				}
			}
		}

		std::pair<Real_t, Real_t> resampleNth(size_t component, size_t step) {
			if (step == 0) {
				return resampleWeights(component, component + 1);
			}
			else if (step == 1) {
				return gaussians[component].resampleVariance();
			}
			else {
				return gaussians[component].resampleMean();
			}
		}

		std::pair<Real_t, Real_t> resample() {
			std::pair<Real_t, Real_t> result = std::make_pair(0.0, 0.0);
			std::for_each(gaussians.begin(), gaussians.end(), [&result](Gaussian<Real_t> &n) { result = result + n.resample(); });
			result = result + resampleWeights(0, weights.size());
			return result;
		}

		GaussianMixture(const GaussianMixture<Real_t> &g): random{g.random} {
			this->size = g.size;
			this->gaussians = g.gaussians;
			this->weights = g.weights;
			this->logWeights = g.logWeights;
			this->normalizedLogWeights = g.normalizedLogWeights;
			this->stepsSizeWeightsRW = g.stepsSizeWeightsRW;
		}

		GaussianMixture<Real_t>& operator = (const GaussianMixture<Real_t> &g) {
			this->size = g.size;
			this->gaussians = g.gaussians;
			this->weights = g.weights;
			this->logWeights = g.logWeights;
			this->normalizedLogWeights = g.normalizedLogWeights;
			this->random = g.random;
			return *this;
		}

	};
}

#endif // !GAUSSIAN_MIXTURE_H
