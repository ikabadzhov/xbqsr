#ifndef GENOMICPIPELINE_ALIGNMENTREADERFACTORY_H
#define GENOMICPIPELINE_ALIGNMENTREADERFACTORY_H

#include <memory>
#include "AlignmentIO.h"

class AlignmentReaderFactory {
public:
    AlignmentReaderFactory() = default;
    static std::unique_ptr<AlignmentReader> open(const std::string& filename);
};



#endif //GENOMICPIPELINE_ALIGNMENTREADERFACTORY_H
