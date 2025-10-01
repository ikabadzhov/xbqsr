#include <future>
#include "SamReader.h"
#include "../parser/SamHeaderParser.h"
#include "../parser/SamAlignmentParser.h"

void SamReader::ParseHeader(Header &hdr) {
    bool isDelimFound = false;
    bool isSameLine = false;
    if (buffer.is_empty()) {
        buffer.fill(in);
    }
    std::string bytes, remaining;
    for (bool first=true; true; first=false) {
        char ch = this->buffer.peek();
        if (ch == -1) {
            buffer.fill(in);
            if (this->buffer.peek() == -1) {
                break;
            }
            ch = this->buffer.peek();
        }
        if (ch != '@' && isSameLine == false)
            break;
        else  isSameLine = true;

        std::string_view slice = this->getslice('\n', isDelimFound);
        if (isDelimFound && remaining.empty()) {
            bytes.append(slice);
            SamHeaderParser::ParseSamHeaderFromSlice(hdr, bytes, first);
            bytes.clear();
            isSameLine = false;
        } else if (isDelimFound && !remaining.empty()) {
            remaining.append(slice);
            SamHeaderParser::ParseSamHeaderFromSlice(hdr, remaining, first);
            remaining.clear();
            isSameLine = false;
        } else if (!isDelimFound && !slice.empty()) { // If this condition is true, then it means that there is a remainder to keep into account
            remaining.append(slice);
        }

        if(buffer.is_empty())
            buffer.fill(in);
    }
}

void parse_header(std::shared_ptr<InputBatchSam> batch, SamAlignments& reads) {
    auto& hdr = reads.GetHdr();
    bool isDelimFound = false;
    bool isSameLine = false;
    if (batch->buffer.is_empty()) {
        return;
    }
    std::string bytes, remaining;
    for (bool first=true; true; first=false) {
        char ch = batch->buffer.peek();
        if (ch == -1) {
            return;
        }
        if (ch != '@' && !isSameLine) {
            batch->header_to_process = false;
            break;
        }
        else  {
            isSameLine = true;
        }

        std::string_view slice = batch->buffer.getslice('\n', isDelimFound);
        if (isDelimFound && remaining.empty()) {
            bytes.append(slice);
            SamHeaderParser::ParseSamHeaderFromSlice(hdr, bytes, first);
            bytes.clear();
            isSameLine = false;
        } else if (isDelimFound && !remaining.empty()) {
            remaining.append(slice);
            SamHeaderParser::ParseSamHeaderFromSlice(hdr, remaining, first);
            remaining.clear();
            isSameLine = false;
        } else if (!isDelimFound && !slice.empty()) { // If this condition is true, then it means that there is a remainder to keep into account
            remaining.append(slice);
        }
    }
}

void parse_alignments(std::shared_ptr<InputBatchSam> batch, SamAlignments& reads) {
    while (batch->buffer.has_data()) {
        if (batch->n_processed == reads.GetAlnsNum()) {
            reads.Resize(reads.GetAlnsNum() + BATCH_SIZE);
        }
        bool delimFound = false;
        std::string_view slice = batch->buffer.getslice('\n', delimFound);

        // If we got nothing from the buffer, exit the loop.
        if (slice.empty() && !delimFound)
            break;

        if (!batch->remaining.empty()) {
            // Append to the already accumulated partial line.
            batch->remaining.append(slice);
            if (delimFound) {
                // Complete line assembled; parse it.
                SamAlignmentParser::ParseSamAlignmentFromSlice(batch->remaining, reads, batch->n_processed);
                batch->n_processed++;
                batch->remaining.clear();
            }
            // Otherwise, if no delimiter was found, keep accumulating and continue the loop.
        } else {
            if (delimFound) {
                // The entire line is contained in the slice.
                SamAlignmentParser::ParseSamAlignmentFromSlice(slice, reads, batch->n_processed);
                batch->n_processed++;
            } else {
                // No newline found: the slice is only a fragment. Accumulate it.
                batch->remaining.append(slice);
            }
        }
    }
}

void load_buffer(std::shared_ptr<InputBatchSam> batch, std::unique_ptr<InputStream>& in) {
    batch->buffer.fill(in);
}


ssize_t SamReader::ParseAlignment(std::vector<SamAlignments> &reads, std::shared_ptr<Header> &header){
    std::shared_ptr<InputBatchSam> batch_transfer = std::make_shared<InputBatchSam>(DEFAULT_BUFFER_SIZE);
    std::shared_ptr<InputBatchSam> batch_process = std::make_shared<InputBatchSam>(DEFAULT_BUFFER_SIZE);

    bool first_batch = true;

    // process until there is no more reads to process
    // except first batch where there is not yet data
    while (first_batch || batch_process->buffer.has_data()) {
        auto processAsync = std::async(std::launch::async,
                                       [=, &reads]() {
                                           if (batch_process->header_to_process) {
                                               parse_header(batch_process, reads[0]);
                                           }
                                           parse_alignments(batch_process, reads[0]);
                                       }
        );
        auto inputAsync = std::async(std::launch::async, load_buffer, batch_transfer, std::ref(in));

        inputAsync.wait();
        processAsync.wait();

        batch_transfer->header_to_process = batch_process->header_to_process;
        batch_transfer->remaining = batch_process->remaining;
        batch_transfer->n_processed = batch_process->n_processed;
        std::swap(batch_process, batch_transfer);
        first_batch = false;
    }

    //reads.Resize(batch_transfer->n_processed);
    //reads.GetHdr().QnameMapper().finalize(reads.QNAME);
    return batch_transfer->n_processed;
}

void SamReader::Close(){
    // this->reader.close();
}

/**
 * Get slice return the lines that are available in the buffer.
 * It does not modify the underlying buffer.
 * The lines are returned as "view" over the current buffer,
 * which means no compies are required BUT they are read-only and temporary:
 * if the buffer is refilled, all the string_views returned will become invalid.
 * It is responsibility of the caller be sure that the buffer is not refilled
 * before all the lines have been returned.
 * */
std::string_view SamReader::getslice(char delim, bool& isDelimFound) {
    return buffer.getslice(delim, isDelimFound);
}
