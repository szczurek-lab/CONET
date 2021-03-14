#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cstddef>
constexpr int MIN_TREE_SIZE_FOR_PARALLEL = 20;

/**
* How many threads shall be used in parallel computations of likelihood.
*/
extern size_t THREADS_NUM;
extern size_t THREADS_LIKELIHOOD;
extern bool USE_EVENT_LENGTHS_IN_ATTACHMENT;
extern double DATA_SIZE_PRIOR_CONSTANT;
extern double COUNTS_SCORE_CONSTANT;
extern double EVENTS_LENGTH_PENALTY;
extern size_t PARAMETER_RESAMPLING_FREQUENCY;
extern size_t NUMBER_OF_MOVES_BETWEEN_SWAPS;
extern size_t MIXTURE_SIZE;
extern long SEED;
extern size_t BURNIN;
#endif // !PARAMETERS_H
