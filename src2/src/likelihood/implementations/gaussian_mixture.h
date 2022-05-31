#ifndef GAUSSIAN_MIXTURE_H
#define GAUSSIAN_MIXTURE_H

#include <algorithm>
#include <numeric>
#include <sstream>


#include "gaussian_utils.h"
#include "../../utils/random.h"
#include "gaussian.h"
#include "../../parameters/parameters.h"
#include "../../utils/log_sum_accumulator.h"
#include "adaptive_mh.h"

namespace Gauss {
	template<class Real_t> class GaussianMixture {
		using vector_r = std::vector<Real_t>;
	public:
		vector_r log_weights;
		vector_r log_normalized_weights;
		std::vector<Gauss::Gaussian<Real_t>> gaussians;
		std::vector<AdaptiveMH<Real_t>> stepsSizeWeightsRW;
		Random<Real_t> &random;

		void fill_matrix_log_likelihood(std::vector<Real_t> &matrix, const std::vector<Real_t> &sample) const {
			std::vector<LogWeightAccumulator<Real_t>> acc(sample.size());

			for (size_t component = 0; component < gaussians.size(); component++) {
				Gauss::truncated_gaussian_log_likelihood(matrix, sample, gaussians[component].mean, gaussians[component].sd);
				const Real_t weight = log_normalized_weights[component];
				for (size_t i = 0; i < sample.size(); i++) {
					acc[i].add(matrix[i] + weight);
				}
			}
			for (size_t i = 0; i < sample.size(); i++) {
				matrix[i] = acc[i].getResult();
			}
		}

		void recalculate_log_normalized_weights() {
			std::vector<Real_t> weights = log_weights;
			auto max_weight = std::max_element(log_weights.begin(), log_weights.end());
			for (auto &weight : weights) {
					weight -= *max_weight;
					weight = std::exp(weight);
			}
			auto sum = std::accumulate(weights.begin(), weights.end(), 0.0);
			for (auto &weight : weights) {
				weight = weight / sum;
			}
			for (size_t i = 0; i < gaussians.size(); i++) {
				log_normalized_weights[i] = weights[i] <= 0.5 ? std::log(weights[i]) : std::log1p(weights[i] - 1);
			}
		}

		void fill_matrix_log_likelihood_parallelized(std::vector<std::vector<Real_t>> &matrix, const std::vector<std::vector<Real_t>> &sample) const {
			std::vector<std::thread> threads;
			size_t rows_per_thread = matrix.size() / THREADS_LIKELIHOOD;
			for (size_t th = 0; th < THREADS_LIKELIHOOD; th++) {
				size_t right = th == THREADS_LIKELIHOOD - 1 ? matrix.size() : rows_per_thread * (th + 1);
				threads.emplace_back([&matrix, &sample, th, rows_per_thread, right, this] {
					for (size_t c = rows_per_thread * th; c < right; c++) {
						this->fill_matrix_log_likelihood(matrix[c], sample[c]);
					}
				});
			}
			for (auto &th : threads) {
				th.join();
			}
		}

		void erase_component(size_t component) {
			gaussians.erase(gaussians.begin() + component);
			log_weights.erase(log_weights.begin() + component);
			log_normalized_weights.erase(log_normalized_weights.begin() + component);
			stepsSizeWeightsRW.erase(stepsSizeWeightsRW.begin() + component);
			recalculate_log_normalized_weights();
		}

	public:
		GaussianMixture(vector_r weights, vector_r means, vector_r sds, Random<Real_t> &random) : random{ random } {
			log_weights.resize(weights.size());
			log_normalized_weights.resize(weights.size());
			stepsSizeWeightsRW.resize(weights.size());
			for (size_t i = 0; i < weights.size(); i++) {
				log_weights[i] = std::log(weights[i]);
				gaussians.push_back(Gauss::Gaussian<Real_t>(means[i], sds[i], random));
			}
			recalculate_log_normalized_weights();
		}

		size_t number_of_components() const {
			return gaussians.size();
		}

		std::pair<Real_t, Real_t> resample_weight(size_t component) {
			log_weights[component] += std::sqrt(stepsSizeWeightsRW[component].get(log_weights[component]))* random.normal();
			recalculate_log_normalized_weights();
			return std::make_pair(1.0, 1.0);
		}

		std::pair<Real_t, Real_t> resample_mean(size_t component) {
			return  gaussians[component].resampleMean();
		}

		std::pair<Real_t, Real_t> resample_sd(size_t component) {
			return gaussians[component].resample_standard_deviation();
		}

		void remove_components_with_small_weight(const Real_t min_weight) {
			std::vector<Real_t> weights{log_normalized_weights};
			std::transform(log_normalized_weights.begin(), log_normalized_weights.end(), weights.begin(), [](Real_t w) {return std::exp(w);});
			while (true) {
				size_t min_index = std::distance(weights.begin(), std::min_element(weights.begin(), weights.end()));
				if (weights[min_index] >= min_weight) {
					break;
				}
				erase_component(min_index);
				weights.erase(weights.begin() + min_index);
			}
		}
	

		Real_t get_parameters_prior() {
			Real_t result = 0.0;
			std::for_each(gaussians.begin(), gaussians.end(), [&result](Gaussian<Real_t> &n) { result += n.get_parameters_prior(); });
			std::for_each(log_weights.begin(), log_weights.end(), [&result](Real_t &w) { result += Gauss::truncated_gaussian_log_likelihood<Real_t>(w, 0.0, 1.0); });
			return result;
		}

		void fill_matrix_log_likelihood(std::vector<std::vector<Real_t>> &matrix, const std::vector<std::vector<Real_t>> &sample) const {
			if (THREADS_LIKELIHOOD > 1) {
				fill_matrix_log_likelihood_parallelized(matrix, sample);
			} else {
				for (size_t c = 0; c < matrix.size(); c++) {
					fill_matrix_log_likelihood(matrix[c], sample[c]);
				}
			}
		}

		GaussianMixture<Real_t>& operator =(const GaussianMixture<Real_t> &g) {
			this->gaussians = g.gaussians;
			this->log_weights = g.log_weights;
			this->log_normalized_weights = g.log_normalized_weights;
			this->random = g.random;
			return *this;
		}
        
        std::string to_string() {
            std::stringstream stream;
            
            for (size_t i = 0; i < log_weights.size(); i++) {
                stream << "(weight: " <<std::exp(log_normalized_weights[i]) << " mean: " <<
                        -gaussians[i].mean <<"sd: " <<gaussians[i].sd << ")\n";
            }
            return stream.str();
        }

	};
}

#endif // !GAUSSIAN_MIXTURE_H
