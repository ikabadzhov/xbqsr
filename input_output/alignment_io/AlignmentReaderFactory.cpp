#include <string_view>

#include "AlignmentReaderFactory.h"
#include "../../zlib.h"
#include "../sam_io/SamReader.h"
#include "../streams/SamStream.h"
#include "../bam_io/BamReader.h"


std::unique_ptr<AlignmentReader> AlignmentReaderFactory::open(const std::string &filename) {
    PlainInputStream fileStream(filename);

    if (SamStream::isBAMFile(fileStream)) {
        std::unique_ptr<InputStream> inputStream = std::make_unique<CompressedInputStream>(filename);
        return std::make_unique<BamReader>(std::move(inputStream));
    }
    throw std::runtime_error("File is not a valid BAM file");
}
