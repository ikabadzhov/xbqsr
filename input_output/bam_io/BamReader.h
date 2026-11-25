#ifndef GENOMICPIPELINE_BAMREADER_H
#define GENOMICPIPELINE_BAMREADER_H


#include "BaiReader.h"
#include "../alignment_io/AlignmentIO.h"
#include "../../zlib.h"
#include "../streams/CompressedInputStream.h"
#include "../streams/BufferReader.h"
#include "../utils/Constant.h"
#include "../parser/IOBatch.h"


class BamReader : public AlignmentReader {

    BufferReader buffer;
    std::string_view getslice(size_t numBytesToRead, size_t &numBytesRead);
    bool header_done_ = false;
    int32_t nextSlice = 0;
    uint64_t slice_id_process = 0;
    uint64_t slice_id_transfer = 0;
    std::vector<Slice> slices_;
    size_t batch_buffer_bytes_ = static_cast<size_t>(DEFAULT_BUFFER_SIZE);

public:
    explicit BamReader(std::unique_ptr<InputStream> inputStream, ssize_t buffSz = 0) : AlignmentReader(std::move(inputStream)),
                                                                                       buffer(static_cast<size_t>(std::min(DEFAULT_BUFFER_SIZE, std::max<ssize_t>(0, buffSz)))) {
        //in->seekg(static_cast<ssize_t>(BAM_MAGIC.size()));
    }

    void ParseHeader(Header& hdr) override;
    void Close() override;
    void ParseAlignmentStreamed(std::shared_ptr<Header> header, std::function<SamAlignments()> get_batch,
                                        std::function<void(SamAlignments&&)> on_batch) override;
    void SeekVirtual(VOff voff_start, VOff voff_end) override;

    void SetSlices(std::vector<Slice>&& slices);
    void SetBatchBufferBytes(size_t bytes) { batch_buffer_bytes_ = bytes; }
};


#endif //GENOMICPIPELINE_BAMREADER_H
