#ifndef CONET_RESULT_H
#define CONET_RESULT_H

#include "tree/event_tree.h"
#include "types.h"
#include "tree/attachment.h"
template <class Real_t> class CONETInferenceResult {
public:
    EventTree tree;
    Attachment attachment;
    Real_t likelihood;

    CONETInferenceResult<Real_t>(EventTree t, Attachment a, Real_t likelihood): tree{t}, attachment{a}, likelihood{likelihood} {}
};

#endif