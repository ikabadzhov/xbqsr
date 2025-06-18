#ifndef GENOMICPIPELINE_ALIGNMENTWRITERFACTORY_H
#define GENOMICPIPELINE_ALIGNMENTWRITERFACTORY_H

#include <string>
#include <memory>
#include "AlignmentIO.h"

class AlignmentWriterFactory {
public:
    AlignmentWriterFactory() = default;
    static std::unique_ptr<AlignmentWriter> open(const std::string& filename, size_t buffSz = (1<<30));
};


#endif //GENOMICPIPELINE_ALIGNMENTWRITERFACTORY_H
