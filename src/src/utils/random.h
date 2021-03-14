#ifndef RANDOM_REGISTRY_H
#define RANDOM_REGISTRY_H
#include <random>
#include <cassert>
#include <limits>

/**
 * Encapsulates all random services which may be used by any of components of MH sampler.
 */
template <class Real_t> class Random {
	std::mt19937 generator;
	unsigned int seed;
public:
	Random(unsigned int seed) : generator{ seed }, seed{ seed } {}

	/*Random(): generator{2134153}
	{
		
	}*/

	/** 
	* Generates rundom number from <code>{0,.., bound-1}</code> uniformly
	*/
	size_t nextInt(size_t bound) {
		std::uniform_int_distribution<size_t> dist(0, bound - 1);
		return dist(generator);
	}

	size_t nextInt() {
		return nextInt(std::numeric_limits<size_t>::max());
	}

	Real_t logUniform() {
		std::uniform_real_distribution<Real_t> dist(0.0, 1.0);
		return std::log(dist(generator));
	}

	Real_t uniform() {
		std::uniform_real_distribution<Real_t> dist(0.0, 1.0);
		return dist(generator);
	}

	/**
	* Sample from standard gaussian distribution
	*/
	Real_t normal() {
		std::normal_distribution<Real_t> dist(0.0, 1.0);
		return dist(generator);
	}

	size_t discrete(std::vector<Real_t> &weights) {
		std::discrete_distribution<size_t> distribution(weights.begin(), weights.end());
		return distribution(generator);
	}

	int randomIntBit() {
		return (int) nextInt(2);
	}
};

#endif // !RANDOM_REGISTRY_H
