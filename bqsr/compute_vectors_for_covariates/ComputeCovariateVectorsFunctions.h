#ifndef GENOMICPIPELINE_COMPUTECOVARIATEVECTORSFUNCTIONS_H
#define GENOMICPIPELINE_COMPUTECOVARIATEVECTORSFUNCTIONS_H

#include "ComputeCovariateVectorsData.h"
#include "../../filters/UtilsFilters.h"
#include "../BQSR.h"

void compute_skip_slice_kernel(std::shared_ptr<Batch<ComputeSkipSliceKernelData, ComputeSkipSliceConfiguration>> batch);

#endif //GENOMICPIPELINE_COMPUTECOVARIATEVECTORSFUNCTIONS_H
