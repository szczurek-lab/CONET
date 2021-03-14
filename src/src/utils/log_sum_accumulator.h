#ifndef LOG_SUM_ACCUMULATOR_H
#define LOG_SUM_ACCUMULATOR_H
#include <cassert>
#include <cmath>
/**
*	Iteratively calculates <code>log( exp(w_1) +...+ exp(w_n)) </code> for input <code> w_1,..,w_n </code>
*/
template<class Real_t> class LogWeightAccumulator {
public:
	Real_t max{ 0 };
	Real_t sum{ 0 };
	bool maxSet{ false };
public:
	void clear() {
		max = 0;
		sum = 0;
		maxSet = false;
	}

	void add(const Real_t logWeight) {
		if (!maxSet) {
			max = logWeight;
			sum += 1.0;
			maxSet = true;
		}
		else if (logWeight <= max) {
			sum += std::exp(logWeight - max);
		}
		else { //logWeight > max && maxSet
			sum *= std::exp(max - logWeight);
			sum += 1.0;
			max = logWeight;
		}
	}

	Real_t getResult() const {
		return std::log(sum) + max;
	}
};


#endif // !LOG_SUM_ACCUMULATOR_H
