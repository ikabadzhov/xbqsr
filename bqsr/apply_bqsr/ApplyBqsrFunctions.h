#ifndef APPLYBQSRFUNCTIONS_H
#define APPLYBQSRFUNCTIONS_H

#include <memory>
#include "../../device/device_types.h"
#include "ApplyBqsrKernelData.h"


void apply_bqsr_kernel(std::shared_ptr<Batch<AppyBqsrKernelData, AppyBqsrKernelConfiguration>> batch);


#endif //APPLYBQSRFUNCTIONS_H
