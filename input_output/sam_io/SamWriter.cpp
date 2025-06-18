#include <future>
#include "SamWriter.h"
#include "../formatter/SamHeaderFormatter.h"
#include "../formatter/SamAlignmentFormatter.h"
#include "../parser/IOBatch.h"
#include "../../sycl_parallel_library/scan/scan.h"
#include "../../timer/Timer.h"
#include "../../logger/Logger.h"

void SamWriter::FormatHeader(Header &hdr) {
    std::string headerow;
    SamHeaderFormatter::FormatSamHeader(hdr, headerow);

    size_t startPos = 0;
    size_t bufferSize = buffer.getSize();

    while (startPos < headerow.size()) {
        size_t remainingSize = headerow.size() - startPos;
        size_t chunkSize = std::min(bufferSize, remainingSize);
        std::string partialContent = headerow.substr(startPos, chunkSize);
        startPos += chunkSize;

        buffer.addSlice(partialContent);
        this->out->write(partialContent, -1);
        buffer.flushContent();
    }

    std::string remainedContent;
    buffer.getContent(0, buffer.getLastPos(), remainedContent);
    this->out->write(remainedContent, -1);
    buffer.flushContent();
}

void process_batch(std::shared_ptr<OutputBatchSam> batch, SamAlignments& reads, size_t* lines_offsets, size_t sz_lines_offsets) {
    size_t n_loaded = batch->n_loaded;
    size_t num_elements_to_process = 0;
    if (lines_offsets[sz_lines_offsets-1] - lines_offsets[n_loaded] <= DEFAULT_BUFFER_SIZE) {
        num_elements_to_process = sz_lines_offsets - 1 - n_loaded;
    } else {
        auto itr = std::lower_bound(lines_offsets + n_loaded, lines_offsets + sz_lines_offsets,
                                    lines_offsets[n_loaded] + DEFAULT_BUFFER_SIZE);
        num_elements_to_process = std::distance(lines_offsets + n_loaded, itr);
    }
    size_t from = n_loaded;
    size_t to = n_loaded + num_elements_to_process;
    batch->n_loaded += num_elements_to_process;


    tbb::parallel_for(from, to, [&](size_t idx){
    //for (size_t idx = from; idx < to; idx++) {
        size_t offset = lines_offsets[idx] - lines_offsets[from];
        SamAlignmentFormatter::FormatSamAlignmentInBuffer(reads, idx, batch->buffer.data() + offset);
        //SamAlignmentFormatter::FormatSamAlignment(reads, idx, temporary_line);
        //std::memcpy(batch->buffer.data() + offset, temporary_line.data(), temporary_line.size());
    });
    //}

    batch->n_written_in_buffer = lines_offsets[to] - lines_offsets[from];
    batch->n_processed += (to - from);
}


void load_batch(std::shared_ptr<OutputBatchSam> batch, size_t* lines_offsets, size_t sz_lines_offsets) {

}

//[[clang::optnone]]
void SamWriter::FormatAlignment(SamAlignments &reads) {

    Timer t;
    t.reset();
    sycl::queue q(sycl::gpu_selector_v);
    mem_t<size_t> sam_lines_sz(reads.GetAlnsNum(), q, memory_space_t::memory_space_shared);
    mem_t<size_t> scanned_sam_lines_sz(reads.GetAlnsNum() + 1, q, memory_space_t::memory_space_shared);

    tbb::parallel_for(0UL, reads.GetAlnsNum(), [&](size_t aln_idx){
        sam_lines_sz.data()[aln_idx] = SamAlignmentFormatter::GetSamAlignmentSize(reads, aln_idx);
    });

    Logger::getLogger().log(LOG_INFO, "Time to compute line sizes: " + std::to_string(t.elapsed()));

    t.reset();

    scan(sam_lines_sz.data(), static_cast<int32_t>(sam_lines_sz.size()), scanned_sam_lines_sz.data(), q);
    scanned_sam_lines_sz.data()[scanned_sam_lines_sz.size() - 1] = scanned_sam_lines_sz.data()[scanned_sam_lines_sz.size() - 2] + sam_lines_sz.data()[sam_lines_sz.size()-1];
    Logger::getLogger().log(LOG_INFO, "Time to scan line sizes: " + std::to_string(t.elapsed()));

    std::shared_ptr<OutputBatchSam> batch_transfer = std::make_shared<OutputBatchSam>(DEFAULT_BUFFER_SIZE);
    std::shared_ptr<OutputBatchSam> batch_process = std::make_shared<OutputBatchSam>(DEFAULT_BUFFER_SIZE);

    bool first_batch = true;

    // process until there is no more reads to process
    // except first batch where there is not yet data
    while (first_batch || batch_process->n_processed != reads.GetAlnsNum())
    {
        std::future<void> transferAsync;
        transferAsync = std::async(std::launch::async, [&]() {
            if (batch_transfer->n_written_in_buffer > 0) {
            this->out->write(batch_transfer->buffer, batch_transfer->n_written_in_buffer);
        }
        });
        //store_batch(output, batch_transfer);
        // if (batch_transfer->n_written_in_buffer > 0) {
        //     this->out->write(batch_transfer->buffer, batch_transfer->n_written_in_buffer);
        // }
        auto processAsync = std::async(std::launch::async, process_batch, batch_process, std::ref(reads), scanned_sam_lines_sz.data(), scanned_sam_lines_sz.size()); // does nothing if process_data->n_seqs is 0
        //process_batch(batch_process, reads, scanned_sam_lines_sz.data(), scanned_sam_lines_sz.size());

        transferAsync.wait();

        //auto inputAsync = std::async(std::launch::async, load_batch<kernel_data_t, batch_config_t, input_data_t>, std::ref(input), batch_transfer);

        //load_batch(batch_transfer, scanned_sam_lines_sz.data(), scanned_sam_lines_sz.size());

        //inputAsync.wait();
        processAsync.wait();
        //batch_process->batch_offset += batch_process->_n;

        batch_transfer->n_loaded = batch_process->n_loaded;
        batch_transfer->n_processed = batch_process->n_processed;

        std::swap(batch_process, batch_transfer);

        first_batch = false;
    }

    if (batch_transfer->n_written_in_buffer > 0) {
        this->out->write(batch_transfer->buffer, batch_transfer->n_written_in_buffer);
    }
//    auto transferAsync = std::async(std::launch::async, store_batch<kernel_data_t, batch_config_t, output_data_t>, std::ref(output),
//                                    batch_transfer);
//    transferAsync.wait();


//    for(int idx=0; idx<reads.GetAlnsNum(); idx++) {
//        std::string alignrow;
//        SamAlignmentFormatter::FormatSamAlignment(reads, idx, alignrow);
//        if(alignrow.size() < buffer.available_data()){
//            buffer.addSlice(alignrow);
//        }else {
//            size_t startPos = buffer.getLastPos();
//            size_t bufferSize = buffer.getSize();
//            size_t remainingSize = bufferSize - startPos;
//            buffer.addSlice(alignrow.substr(0, remainingSize));
//            std::string partialContent;
//            buffer.getContent(0, buffer.getLastPos(), partialContent);
//            this->out->write(partialContent, -1);
//            buffer.flushContent();
//
//            if (alignrow.size() - remainingSize < buffer.available_data()){
//                buffer.addSlice(alignrow.substr(remainingSize, alignrow.size()));
//            }else{
//                size_t chunkSize = alignrow.size()-remainingSize;
//                std::string partialContent(&alignrow[remainingSize], chunkSize);
//                this->out->write(partialContent, -1);
//            }
//        }
//    }
//    if(buffer.getLastPos()>0) {
//        std::string remainedContent;
//        buffer.getContent(0, buffer.getLastPos(), remainedContent);
//        this->out->write(remainedContent, -1);
//        buffer.flushContent();
//    }
    Close();
}

void SamWriter::Close() {
    this->out->close();
}

SamWriter::~SamWriter() = default;
