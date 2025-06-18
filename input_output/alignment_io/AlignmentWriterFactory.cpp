#include "AlignmentWriterFactory.h"
#include "../sam_io/SamWriter.h"
#include "../bam_io/BamWriter.h"
#include "../streams/CompressedOutputStream.h"
#include "../streams/PlainOutputStream.h"

std::unique_ptr<AlignmentWriter> AlignmentWriterFactory::open(const std::string &filename, size_t buffSz) {
    std::string fileExtension = filename.substr(filename.size()-4, 4);
    std::transform(fileExtension.begin(), fileExtension.end(), fileExtension.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    if (fileExtension == ".bam") {
        std::unique_ptr<OutputStream> outputStream = std::make_unique<CompressedOutputStream>(filename);
        return std::make_unique<BamWriter>(std::move(outputStream), buffSz);
    } else {
        std::unique_ptr<OutputStream> outputStream = std::make_unique<PlainOutputStream>(filename);
        return std::make_unique<SamWriter>(std::move(outputStream), buffSz);
    }
}
