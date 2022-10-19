#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cstddef>

extern size_t NUM_REPLICAS;
extern size_t THREADS_LIKELIHOOD;
extern bool USE_EVENT_LENGTHS_IN_ATTACHMENT;
extern double DATA_SIZE_PRIOR_CONSTANT;
extern double COUNTS_SCORE_CONSTANT_0;
extern double COUNTS_SCORE_CONSTANT_1;
extern double EVENTS_LENGTH_PENALTY;
extern size_t PARAMETER_RESAMPLING_FREQUENCY;
extern size_t NUMBER_OF_MOVES_BETWEEN_SWAPS;
extern size_t MIXTURE_SIZE;
extern long SEED;
extern bool VERBOSE;
extern double NEUTRAL_CN;
#endif // !PARAMETERS_H
