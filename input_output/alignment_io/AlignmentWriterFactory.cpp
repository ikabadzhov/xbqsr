#include "AlignmentWriterFactory.h"
#include "../bam_io/BamWriter.h"
#include "../streams/BGZFOutputStream.h"
#include "../streams/CompressedOutputStream.h"
#include "../streams/PlainOutputStream.h"

std::unique_ptr<AlignmentWriter> AlignmentWriterFactory::open(const std::string &filename, size_t buffSz) {
    std::string fileExtension = filename.substr(filename.size()-4, 4);
    std::transform(fileExtension.begin(), fileExtension.end(), fileExtension.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    if (fileExtension == ".bam") {
        std::unique_ptr<OutputStream> outputStream = std::make_unique<BGZFOutputStream>(filename);
        return std::make_unique<BamWriter>(std::move(outputStream));
    } else {
        throw std::runtime_error("File extension not supported");
    }
}
