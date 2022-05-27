#ifndef CONET_RESULT_H
#define CONET_RESULT_H

#include "tree/event_tree.h"
#include "types.h"
template <class Real_t> class CONETInferenceResult {
public:
    EventTree tree;
    std::vector<Event> attachment;
};

#endif