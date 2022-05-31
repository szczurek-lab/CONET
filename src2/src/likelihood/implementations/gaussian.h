#ifndef GAUSSIAN_H
#define GAUSSIAN_H
#include <utility>
#include <thread>

#include "gaussian_utils.h"
#include "../../utils/random.h"
#include "../../parameters/parameters.h"
#include "adaptive_mh.h"

namespace Gauss {
	template <class Real_t> class Gaussian {
	public:
		Real_t mean;
		Real_t sd;
		Random<Real_t> &random;
		AdaptiveMH<Real_t> adaptiveMHMean; 
		AdaptiveMH<Real_t> adaptiveMHVariance; 

		void fill_matrix_log_likelihood_parallelized(std::vector<std::vector<Real_t>> &matrix, const std::vector<std::vector<Real_t>> &sample) const {
			std::vector<std::thread> threads;
			size_t perThread = matrix.size() / THREADS_LIKELIHOOD;
			for (size_t th = 0; th < THREADS_LIKELIHOOD; th++) {
				size_t right = th == THREADS_LIKELIHOOD - 1 ? matrix.size() : perThread * (th + 1);
				threads.emplace_back([&matrix, &sample, th, perThread, right, this] { 
					for (size_t c = perThread * th; c < right; c++) {
						Gauss::truncated_gaussian_log_likelihood(matrix[c], sample[c], this->mean, this->sd);
					}
				 });
			}
			for (auto &th : threads) {
				th.join();
			}
		}

	public:
		Gaussian(Real_t mean, Real_t sd, Random<Real_t> &random) : mean {mean}, sd {sd}, random {random} {}


		void fill_matrix_log_likelihood(std::vector<std::vector<Real_t>> &matrix, const std::vector<std::vector<Real_t>> &sample) const {
			if (THREADS_LIKELIHOOD > 1) {
				fill_matrix_log_likelihood_parallelized(matrix, sample);
			} else {
				for (size_t c = 0; c < matrix.size(); c++) {
					Gauss::truncated_gaussian_log_likelihood(matrix[c], sample[c], mean, sd);
				}
			}
		}

		Real_t get_parameters_prior() {
			return Gauss::truncated_gaussian_log_likelihood<Real_t>(mean, 0.0, 1.0) + Gauss::truncated_gaussian_log_likelihood<Real_t>(sd, 0.0, 1.0);
		}

		std::pair<Real_t, Real_t> resample() {
			return resample_standard_deviation() + resampleMean();
		}

		std::pair<Real_t, Real_t> resampleMean() {
			std::pair<Real_t, Real_t> meanResult = std::make_pair(0.0, 0.0);
			const Real_t step = adaptiveMHMean.get(mean);
			mean += std::sqrt(step) * random.normal();
			return meanResult;
		}

		std::pair<Real_t, Real_t> resample_standard_deviation() {
			const Real_t step = adaptiveMHVariance.get(sd);
			sd += std::sqrt(step) * random.normal();
			return std::make_pair(1.0, 1.0);
		}

		Gaussian(const Gaussian<Real_t> &g) : random{g.random} {
			this->mean = g.mean;
			this->sd = g.sd;
		}

		Gaussian<Real_t>& operator = (const Gaussian<Real_t> &g) {
			this->mean = g.mean;
			this->sd = g.sd;
			this->random = g.random;
			return *this;
		}
        
        std::string to_string() {
             std::stringstream stream;
             stream << "mean: " << mean << " sd: " << sd << "\n";
             return stream.str();
        }
	};
}








#endif
