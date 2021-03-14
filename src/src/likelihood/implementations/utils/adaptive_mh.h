#ifndef ADAPTIVE_MH_H
#define ADAPTIVE_MH_H


template <class Real_t> class AdaptiveMH {
private:
	const Real_t epsilon = 0.0001;
	const Real_t initialVar = 0.0001;
	const Real_t scaling = (1.0) * (1.0);//2.4
	const Real_t initSegment = 10;
	
	Real_t var{ initialVar };
	Real_t average{ 0 };

	Real_t numberObservations{ 0 };


public:
	AdaptiveMH<Real_t>() {}

	Real_t get()
	{
		return scaling * var + scaling * epsilon;
		
	}
	Real_t get(Real_t x) {
		const Real_t averageOld = average;
		average = average * numberObservations / (numberObservations + 1) + 1 / (numberObservations + 1) * x;
		if (numberObservations == 1) {
			var = x*x + averageOld * averageOld - (2)*average*average;
		} else if (numberObservations > 0) {
			var = (numberObservations - 1) / (numberObservations)* var + (1 / (numberObservations)) * (numberObservations* averageOld * averageOld + x * x - (numberObservations + 1)*average * average);
		}
		numberObservations++;

		if (numberObservations < initSegment) {
			return initialVar;
		}
		else {
			return scaling * var + scaling * epsilon;
		}

	}


	AdaptiveMH<Real_t>(const AdaptiveMH<Real_t> &g) {
		this->var = g.var;
		this->average = g.average;
		this->numberObservations = g.numberObservations;
	}

	AdaptiveMH<Real_t>& operator = (const AdaptiveMH<Real_t> &g) {
		this->var = g.var;
		this->average = g.average;
		this->numberObservations = g.numberObservations;
		return *this;
	}
};








#endif //ADAPTIVE_MH_H
