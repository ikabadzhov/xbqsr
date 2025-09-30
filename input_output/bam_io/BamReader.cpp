#include <future>
#include "BamReader.h"
#include "BamUtilsIO.h"
#include "../../timer/ProfilerIO.h"
#include "../parser/SamHeaderParser.h"
#include "../parser/BamAlignmentParser.h"

/*
 * This function is designed to parse the header utilizing the nested BAMUtilsIO
 * function. Subsequently, upon extracting the information, a validation function
 * is called to verify the accuracy of the quantity of information
 * (not the content itself), ensuring that all header lines have been
 * successfully read.
*/
void BamReader::ParseHeader(Header &hdr) {

    std::string headerData(sizeof(uint32_t), 0);
    uint32_t lengthToRead = 0;
    size_t numBytesRead = in->read(headerData, 0, sizeof(uint32_t));
    assert(numBytesRead == sizeof(uint32_t));
    std::memcpy(&lengthToRead, headerData.data(), sizeof(uint32_t));

    headerData.resize(lengthToRead, 0);
    numBytesRead = in->read(headerData, 0, lengthToRead);
    assert(numBytesRead == lengthToRead);
    size_t pos = 0;

    while (pos < headerData.size()) {
        auto endPos = std::find(headerData.begin() + pos, headerData.end(), '\n');

        // Check if a newline character was found
        if (endPos != headerData.end()) {
            // Extract the line and add it to the lines vector
            std::string line(headerData.begin() + pos, endPos);
            SamHeaderParser::ParseSamHeaderFromSlice(hdr, line);
            // Move pos to the character after the found newline
            pos = std::distance(headerData.begin(), endPos) + 1;
        } else {
            break; // Exit the loop
        }
    }

    headerData.clear();
    //headerData.resize(pos);
    headerData.resize(sizeof(uint32_t), 0);

    uint32_t checkNREF = 0;
    numBytesRead = in->read(headerData, 0, sizeof(uint32_t));
    assert(numBytesRead == sizeof(uint32_t));
    std::memcpy(&checkNREF, headerData.data(), sizeof(uint32_t));
    /**
     *
     * Save these in the header
     * */
    for(int i=0; i<checkNREF; i++) {
        uint32_t tempLengthToRead = 0;
        numBytesRead = in->read(headerData, 0, sizeof(uint32_t));
        assert(numBytesRead == sizeof(uint32_t));
        std::memcpy(&tempLengthToRead, headerData.data(), sizeof(uint32_t));

        std::string tempBuffer(tempLengthToRead, 0);
        numBytesRead = in->read(tempBuffer, 0, tempLengthToRead);

        std::string_view name(tempBuffer.data(), tempBuffer.data() + tempLengthToRead);
        //auto null_terminator_pos = std::find(name.begin(), name.end(), '\0');

        uint32_t tempLengthValue = 0;
        numBytesRead = in->read(headerData, 0, sizeof(uint32_t));
        assert(numBytesRead == sizeof(uint32_t));
        std::memcpy(&tempLengthValue, headerData.data(), sizeof(uint32_t));

        if (!hdr.getSequenceDictionary().getAll().empty()) {
            std::string_view hdrSN = hdr.getSequenceDictionary().getSequenceName(i);
            assert(hdrSN.compare(name.data()) == 0);
        } else {
            SequenceRecord sr(name.data(), tempLengthValue);
            hdr.getSequenceDictionary().add(sr);
        }
    }
}

/*
 * Once the file analysis is complete, there's no requirement to maintain the stream
 * and the gzFile in an open state. Closing them aids in conserving resources.
*/
void BamReader::Close() {
    this->in->close();
}

// Parse BAM header: fills your header structure in `reads`.
// Throws std::runtime_error on invalid/truncated headers.
// Parse BAM header: fills your header structure in `reads`.
// Throws std::runtime_error on invalid/truncated headers.
static void parse_bam_header(std::unique_ptr<InputStream>& in, SamAlignments& reads) {
    Header& hdr = reads.GetHdr();

    // 1) magic "BAM\1"
    {
        std::string magic = read_exact_new(in, 4);
        if (!(magic.size() == 4 && magic[0]=='B' && magic[1]=='A' && magic[2]=='M' && magic[3]==1)) {
            throw std::runtime_error("Not a BAM file: missing BAM\\1 magic");
        }
    }

    // 2) l_text
    uint32_t l_text = 0;
    {
        std::string b = read_exact_new(in, 4);
        l_text = load_le_u32(b.data());
        if (l_text > (64u * 1024u * 1024u)) {
            throw std::runtime_error("BAM header text too large");
        }
    }

    // 3) text (SAM header) — parse line by line like your BamReader::ParseHeader
    if (l_text > 0) {
        std::string headerData = read_exact_new(in, l_text);

        size_t pos = 0;
        while (pos < headerData.size()) {
            const auto* begin = headerData.data();
            const char* lineStart = begin + pos;
            const char* lineEnd   = static_cast<const char*>(
                    memchr(lineStart, '\n', headerData.size() - pos)
            );
            if (!lineEnd) break;

            std::string line(lineStart, lineEnd);
            if (!line.empty()) {
                SamHeaderParser::ParseSamHeaderFromSlice(hdr, line);
            }
            pos = static_cast<size_t>(lineEnd - begin) + 1; // move past '\n'
        }
    }

    // 4) n_ref
    uint32_t n_ref = 0;
    {
        std::string b = read_exact_new(in, 4);
        n_ref = load_le_u32(b.data());
        if (n_ref > (10u * 1000u * 1000u)) {
            throw std::runtime_error("BAM header: unreasonable number of references");
        }
    }

    // 5) references: (l_name [4], name [l_name] NUL-terminated, l_ref [4])
    for (uint32_t i = 0; i < n_ref; ++i) {
        // l_name
        uint32_t l_name = 0;
        {
            std::string b = read_exact_new(in, 4);
            l_name = load_le_u32(b.data());
            if (l_name == 0 || l_name > (1u * 1024u * 1024u)) {
                throw std::runtime_error("BAM header: invalid l_name");
            }
        }

        // name
        std::string name = read_exact_new(in, l_name);
        if (name.empty() || name.back() != '\0') {
            throw std::runtime_error("BAM header: reference name not NUL-terminated");
        }
        name.pop_back(); // drop NUL

        // l_ref
        uint32_t l_ref = 0;
        {
            std::string b = read_exact_new(in, 4);
            l_ref = load_le_u32(b.data());
        }

        // Integrate with your header dictionary like your original code
        if (!hdr.getSequenceDictionary().getAll().empty()) {
            // assert existing header matches BAM dictionary order/content
            std::string_view hdrSN = hdr.getSequenceDictionary().getSequenceName(i);
            // compare name strings (hdrSN is a view)
            if (hdrSN.compare(name) != 0) {
                throw std::runtime_error("BAM header: reference name mismatch vs parsed SAM header");
            }
            // You could also sanity-check length if your dictionary stores it.
        } else {
            SequenceRecord sr(name.c_str(), l_ref);
            hdr.getSequenceDictionary().add(sr);
        }
    }
}

static void load_batch(const std::shared_ptr<InputBatchBam> &bt, std::unique_ptr<InputStream> &in) {
    auto& profiler = ProfilerIO::instance();
    profiler.start(IOStageId::DataLoading);
    bt->buffer.fill(in);
    profiler.stop(IOStageId::DataLoading);
};

// --- Carry state for a single BAM record across batches ---
struct Carry {
    enum Stage : uint8_t { ReadLen, ReadBody } stage = ReadLen;
    int32_t need  = sizeof(int32_t);   // bytes still needed for current stage
    int32_t block = -1;                // block_size after decoding
    std::string bin;                   // accumulation buffer
    Carry() { bin.reserve(4096); }     // small upfront reserve reduces churn
} ;

static void ensure_capacity(SamAlignments& reads, size_t used, size_t GROW) {
    if (used == reads.GetAlnsNum()) {
        reads.Resize(reads.GetAlnsNum() + GROW);  // if you have Reserve(cap), prefer that
    }
};
/*
struct RecordInfo { size_t off; int32_t blk; }; // relative to buf.data()

static void process_batch(const std::shared_ptr<InputBatchBam>& bp,
                          SamAlignments& reads,
                          size_t& alnsRead,
                          Carry& carry)
{
    auto& buf = bp->buffer;

    // 0) Finish any partial record carried from previous batch
    //    (identical to your current logic)
    while (buf.has_data() && carry.stage == Carry::ReadBody && carry.need > 0) {
        const size_t avail = buf.available_data();
        const size_t take = static_cast<size_t>(std::min<int64_t>(carry.need, (int64_t)avail));
        size_t got = 0;
        std::string_view chunk = buf.getslice(take, got);
        if (got == 0) break;
        carry.bin.append(chunk.data(), got);
        carry.need -= static_cast<int32_t>(got);
        if (carry.need == 0) {
            // parse the completed record
            ensure_capacity(reads, alnsRead);
            BamAlignmentParser::ParseBamAlignmentFromSlice(
                std::string_view(carry.bin.data(), carry.bin.size()), reads, alnsRead);
            ++alnsRead;
            carry.bin.clear();
            carry.stage = Carry::ReadLen;
            carry.need  = sizeof(int32_t);
            carry.block = -1;
        }
    }

    // 1) Snapshot the remaining batch bytes (without consuming them)
    const uint8_t* base = reinterpret_cast<const uint8_t*>(buf.data());  // add data()/size() accessors if missing
    size_t         n    = buf.size();
    if (n == 0) return;

    // 2) Prescan: find complete records {offset, block_size}
    std::vector<RecordInfo> recs;
    recs.reserve(1024);
    size_t p = 0;
    while (p + 4 <= n) {
        int32_t block = 0;
        std::memcpy(&block, base + p, 4); // little-endian block_size
        if (block <= 0 || block > (64 * 1024 * 1024))
            throw std::runtime_error("Invalid BAM block_size");
        const size_t rec_bytes = 4u + static_cast<size_t>(block);
        if (p + rec_bytes > n) break;     // incomplete record => becomes tail/carry
        recs.push_back({p, block});
        p += rec_bytes;
    }
    const size_t consumed_bytes = p;      // bytes we will consume from buffer
    const size_t tail_bytes     = n - p;  // keep as carry

    // 3) Pre-size outputs and parallel-parse complete records
    ensure_capacity(reads, alnsRead + recs.size());
    const size_t base_dst = alnsRead;

    // Use your favorite parallel lib (TBB/OpenMP/std::threads). Example: TBB
    tbb::parallel_for(size_t(0), recs.size(), [&](size_t k) {
        const auto off = recs[k].off;
        const auto blk = recs[k].blk;
        // body is after the 4-byte header
        std::string_view body(reinterpret_cast<const char*>(base + off + 4), static_cast<size_t>(blk));
        const size_t dst = base_dst + k;  // unique slot (no locking)
        BamAlignmentParser::ParseBamAlignmentFromSlice(body, reads, dst);
    });

    alnsRead += recs.size();

    // 4) Advance the buffer by the consumed bytes
    buf.consume(consumed_bytes);          // implement this to pop 'consumed_bytes' from the front

    // 5) Move tail into carry for next batch
    if (tail_bytes > 0) {
        // Append tail bytes to carry.bin
        carry.bin.assign(reinterpret_cast<const char*>(base + consumed_bytes), tail_bytes);

        if (tail_bytes < 4) {
            // still in ReadLen stage (need remaining header bytes)
            carry.stage = Carry::ReadLen;
            carry.need  = static_cast<int32_t>(4 - tail_bytes);
            carry.block = -1;
        } else {
            // we have full header of next record; set up to complete its body next time
            int32_t blk = 0;
            std::memcpy(&blk, carry.bin.data(), 4);
            if (blk <= 0 || blk > (64 * 1024 * 1024))
                throw std::runtime_error("Invalid BAM block_size (tail)");
            const int32_t have_body = static_cast<int32_t>(tail_bytes - 4);
            carry.stage = Carry::ReadBody;
            carry.block = blk;
            carry.need  = blk - have_body;    // may be > 0 (incomplete)
        }
    } else {
        // clean carry
        carry.bin.clear();
        carry.stage = Carry::ReadLen;
        carry.need  = 4;
        carry.block = -1;
    }

    bp->n_processed = alnsRead; // cumulative reporting
}*/
//[[clang::optnone]]
static void process_batch(const std::shared_ptr<InputBatchBam>& bp, SamAlignments& reads,
                             size_t& alnsRead, Carry& carry) {
    auto& buf = bp->buffer;
    auto& profiler = ProfilerIO::instance();
    profiler.start(IOStageId::ProcessBatch);

    while (buf.has_data()) {
        // Consume up to 'need' bytes for the current stage
        const size_t avail = buf.available_data();
        if (avail == 0) break;

        const size_t take = static_cast<size_t>(
                std::min<int64_t>(carry.need, static_cast<int64_t>(avail)));
        size_t got = 0;

        profiler.start(IOStageId::GetSlice);
        std::string_view chunk = buf.getslice(take, got);
        profiler.stop(IOStageId::GetSlice);

        if (got == 0) break;
        profiler.start(IOStageId::AppendInBin);
        // Accumulate into carry
        //carry.bin.append(chunk.data(), got);
        carry.need -= static_cast<int32_t>(got);
        profiler.stop(IOStageId::AppendInBin);

        if (carry.need > 0) {
            // Need more bytes → let loader run; we’ll resume after swap
            carry.bin.append(chunk.data(), got);
            break;
        }

        // Stage completed: transition or parse
        if (carry.stage == Carry::ReadLen) {
            // Decode little-endian 32-bit block_size
            int32_t raw_len = 0;
            carry.bin.append(chunk.data(), got);
            std::memcpy(&raw_len, carry.bin.data(), sizeof(raw_len));
            carry.bin.clear();
            const int32_t block = raw_len; // ByteScanner reads LE elsewhere; convert if needed

            // Sanity checks (64MB guard is generous)
            if (block <= 0 || block > (64 * 1024 * 1024)) {
                throw std::runtime_error("Invalid BAM block_size");
            }
            profiler.start(IOStageId::AppendInBin);
            carry.block = block;
            //carry.bin.clear();
            //carry.bin.reserve(static_cast<size_t>(carry.block));
            carry.need  = carry.block;
            carry.stage = Carry::ReadBody;
            profiler.stop(IOStageId::AppendInBin);
            continue;
        }

        // ReadBody complete → parse exactly one record
        profiler.start(IOStageId::EnsureCapacity);
        ensure_capacity(reads, alnsRead, 10*(1<<20));
        profiler.stop(IOStageId::EnsureCapacity);

        if (carry.bin.empty()) {
            BamAlignmentParser::ParseBamAlignmentFromSlice(chunk, reads, alnsRead);
        } else {
            carry.bin.append(chunk.data(), got);
            assert(carry.bin.size() == carry.block);
            BamAlignmentParser::ParseBamAlignmentFromSlice(carry.bin, reads, alnsRead);
            carry.bin.clear();
        }

        ++alnsRead;

        // Reset for next record
        carry.bin.clear();
        carry.stage = Carry::ReadLen;
        carry.need  = sizeof(int32_t);
        carry.block = -1;

        // Loop to see if more records are already present in this batch
    }

    bp->n_processed = alnsRead; // keep your cumulative reporting
    profiler.stop(IOStageId::ProcessBatch);
}

// BamReader::ParseAlignment — optimized, BufferReader untouched
ssize_t BamReader::ParseAlignment(SamAlignments& reads) {

    parse_bam_header(in, reads);

    auto batch_transfer = std::make_shared<InputBatchBam>(DEFAULT_BUFFER_SIZE);
    auto batch_process  = std::make_shared<InputBatchBam>(DEFAULT_BUFFER_SIZE);

    Carry carry;

    load_batch(batch_process, in);


    size_t alnsRead    = 0;
    bool   first_batch = true;

    // Process until both buffer is drained AND no partial body remains
    while (first_batch ||
           batch_process->buffer.has_data()) {

        // Process current batch and load next in parallel
        //auto processF = std::async(std::launch::async,
          //                         [&, batch_process] {
            process_batch(batch_process, reads, alnsRead, carry);
                                   //});

       // auto inputF   = std::async(std::launch::async,
           //                        [&, batch_transfer] {
            load_batch(batch_transfer, in);
           //                        });

        // inputF.get();
        // processF.get();

        // Propagate cumulative counters if you use them elsewhere (kept for parity)
        batch_transfer->header_to_process = batch_process->header_to_process;
        batch_transfer->remaining         = batch_process->remaining;   // unused for BAM, but kept
        batch_transfer->n_processed       = batch_process->n_processed;

        std::swap(batch_process, batch_transfer);
        first_batch = false;
    }

    // Shrink to actual count and finalize
    reads.Resize(alnsRead);
    //reads.GetHdr().QnameMapper().finalize(reads.QNAME);
    auto& profiler = ProfilerIO::instance();
    profiler.report();
    return static_cast<ssize_t>(alnsRead);
}

std::string_view BamReader::getslice(size_t numBytesToRead, size_t &numBytesRead) {
    return buffer.getslice(numBytesToRead, numBytesRead);
}

