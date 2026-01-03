#ifndef IEXECUTOR_H
#define IEXECUTOR_H
#include "../sam_data_structures/SamAlignments.h"

struct IExecutor {
    virtual ~IExecutor() = default;
    virtual void Process(SamAlignments& batch) = 0;   // run the pipeline work on this batch
    virtual void Finalize() {}                        // optional per-executor finalize
};

#endif //IEXECUTOR_H
