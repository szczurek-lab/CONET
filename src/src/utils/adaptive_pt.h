#ifndef ADAPTIVE_PT_H
#define ADAPTIVE_PT_H

#include <cmath>
#include <vector>

template <class Real_t> class AdaptivePT
{
	size_t num_replicas;
	const Real_t min_p = 0.0000000001;
	const Real_t max_p = 100000.0;
	const Real_t alpha = 0.234;
	const Real_t step_size = 0.213;
	const Real_t epsilon_const = 0.75;
	long int step = 0;
	/** of size @num_replicas **/
	std::vector<Real_t> p_constants;

	std::vector<Real_t> get_H_vector(const std::vector<Real_t> &states) const 
	{
		auto temperatures = get_temperatures();
		std::vector<Real_t> result;
		for (size_t i = 0; i < p_constants.size(); i++)
		{
			result.push_back((temperatures[i] - temperatures[i + 1])*(states[i + 1] - states[i]));
		}
		for (auto &el : result)
		{
			if (el >= 0.0 )
			{
				el = -alpha;
			} else
			{
				el = std::exp(el) - alpha;
			}
		}
		return result;
	}


public:
	AdaptivePT(size_t num_replicas): num_replicas{num_replicas}
	{
		Real_t temp = 1.0;
		for (size_t i = 0; i < num_replicas; i++)
		{
			if (i != 0)
			{
				temp = 0.1*temp;
				p_constants.push_back(std::log(-std::log(temp)));
			}
		}
	}
	std::vector<Real_t> get_temperatures() const
	{
		std::vector<Real_t> temperatures;
		temperatures.push_back(1.0);
		for (size_t i = 0; i < p_constants.size(); i++)
		{
			temperatures.push_back(temperatures.back() * std::exp(-std::exp(p_constants[i])));
		}
		return temperatures;
	}



	void update(const std::vector<Real_t> &states)
	{
		auto H_vector = get_H_vector(states);
		step++;
		for (size_t i = 0; i < p_constants.size(); i++)
		{
			p_constants[i] = p_constants[i] + step_size * std::pow((step + 1), -epsilon_const) * H_vector[i];
			if (p_constants[i] <= min_p)
			{
				p_constants[i] = min_p;
			} else if (p_constants[i] >= max_p)
			{
				p_constants[i] = max_p;
			}
		}
	}
	
};





#endif
