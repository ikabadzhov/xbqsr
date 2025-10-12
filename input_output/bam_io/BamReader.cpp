#include <future>
#include <execution>
#include <mutex>
#include "BamReader.h"

#include <oneapi/dpl/functional>
#include <oneapi/tbb/parallel_scan.h>

#include "BamUtilsIO.h"
#include "../../timer/ProfilerIO.h"
#include "../parser/SamHeaderParser.h"
#include "../parser/BamAlignmentParser.h"

void BamReader::Close() {
    this->in->close();
}

// Parse BAM header: fills your header structure in `reads`.
// Throws std::runtime_error on invalid/truncated headers.
// Parse BAM header: fills your header structure in `reads`.
// Throws std::runtime_error on invalid/truncated headers.
static void parse_bam_header(std::unique_ptr<InputStream>& in, Header& hdr) {

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

static void skip_bam_header(std::unique_ptr<InputStream>& in) {
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

    // 3) skip text (SAM header)
    if (l_text > 0) {
        // just read-and-discard
        (void)read_exact_new(in, l_text);
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

    // 5) skip references
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

        // name (NUL terminated string)
        (void)read_exact_new(in, l_name);

        // l_ref
        (void)read_exact_new(in, 4);
    }
}

void BamReader::ParseHeader(Header &hdr) {
    parse_bam_header(in, hdr);
    header_done_ = true;
}

static void load_batch(const std::shared_ptr<InputBatchBam> &bt, std::unique_ptr<InputStream> &in) {
    //auto& profiler = ProfilerIO::instance();
    //profiler.start(IOStageId::DataLoading);
    bt->buffer.fill(in);
    //profiler.stop(IOStageId::DataLoading);
};

// --- Carry state for a single BAM record across batches ---
struct Carry {
    enum Stage : uint8_t { ReadLen, ReadBody } stage = ReadLen;
    int32_t need  = sizeof(int32_t);   // bytes still needed for current stage
    int32_t block = -1;                // block_size after decoding
    std::string bin;                   // accumulation buffer
    std::string bin_tmp;                   // accumulation buffer
    Carry() { bin.reserve(4096); }     // small upfront reserve reduces churn
} ;

static void ensure_capacity(SamAlignments& reads, size_t used, size_t GROW) {
    if (used == reads.GetAlnsNum()) {
        reads.Resize(reads.GetAlnsNum() + GROW);  // if you have Reserve(cap), prefer that
        //reads.Resize(170318378);  // TODO: The iterated resize it's not good!
    }
};

class RecordSpan {
public:
    const char* chunkPtr = nullptr;
    size_t size = 0;
    std::string bin;
    size_t alnIdx = 0;
    RecordSpan(const char* ptr, size_t size, size_t alnIdx, const std::string& bin = {}) : chunkPtr(ptr), size(size), bin(bin), alnIdx(alnIdx) {};
};

//[[clang::optnone]]
static void process_batch(const std::shared_ptr<InputBatchBam>& bp, SamAlignments& reads,
                             size_t& alnsRead, Carry& carry) {
    auto& buf = bp->buffer;
    //auto& profiler = ProfilerIO::instance();
    //profiler.start(IOStageId::ProcessBatch);
    std::vector<RecordSpan> boundaries;
    boundaries.reserve(1024);

    while (buf.has_data()) {
        // Consume up to 'need' bytes for the current stage
        const size_t avail = buf.available_data();
        if (avail == 0) break;

        const size_t take = static_cast<size_t>(std::min<int64_t>(carry.need, static_cast<int64_t>(avail)));
        size_t got = 0;

        //profiler.start(IOStageId::GetSlice);
        std::string_view chunk = buf.getslice(take, got);
        //profiler.stop(IOStageId::GetSlice);

        if (got == 0) break;

        //profiler.start(IOStageId::AppendInBin);
        // Accumulate into carry
        carry.need -= static_cast<int32_t>(got);
        //profiler.stop(IOStageId::AppendInBin);

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
            //profiler.start(IOStageId::AppendInBin);
            carry.block = block;
            carry.need  = carry.block;
            carry.stage = Carry::ReadBody;
            //profiler.stop(IOStageId::AppendInBin);
            continue;
        }

        // ReadBody complete → parse exactly one record
        //profiler.start(IOStageId::EnsureCapacity);
        ensure_capacity(reads, alnsRead, 10*(1<<20));
        //profiler.stop(IOStageId::EnsureCapacity);

        if (carry.bin.empty()) {
            //BamAlignmentParser::ParseBamAlignmentFromSlice(chunk, reads, alnsRead);
            //RecordSpan span;
            boundaries.emplace_back(chunk.data(), chunk.size(), alnsRead);
        } else {
            carry.bin.append(chunk.data(), got);
            RecordSpan rec(nullptr, 0, alnsRead, carry.bin);
            assert(rec.bin.size() == carry.block);
            if (boundaries.size() == boundaries.capacity()) {
                boundaries.reserve(boundaries.capacity() * 2);
            }
            boundaries.push_back(rec);
            //BamAlignmentParser::ParseBamAlignmentFromSlice(carry.bin, reads, alnsRead);
            carry.bin.clear();
        }

        ++alnsRead;

        // Reset for next record
        carry.bin.clear();
        carry.stage = Carry::ReadLen;
        carry.need  = sizeof(int32_t);
        carry.block = -1;
    }

    // If here, it means that the buffer has been processed, it should go to get new bytes.
    // Start the actual rpocessing.
    // Two phases:
    // 1. Parallel Parsing
    //  a) storing of simple field
    //  b) storing of sizes of seq, qual and cigar (variable length fileds)
    // 2. Parallel Scan of sizes
    // 3. Parsing and storing in parallel of the variable length fields
    //profiler.start(IOStageId::Parsing);
    std::vector<VariableLenFieldsVar> fieldsOff(boundaries.size());
    std::vector<VariableLenFieldsVar> fieldsSizes(boundaries.size() + 1);

    //for (size_t recIdx = 0UL; recIdx < boundaries.size(); ++recIdx) {
    tbb::parallel_for (0UL, boundaries.size(), [&boundaries, &reads, &fieldsOff, &fieldsSizes](size_t recIdx) {
        auto begin = boundaries[recIdx].chunkPtr;
        auto end = begin + boundaries[recIdx].size;
        size_t aIdx  = boundaries[recIdx].alnIdx;

        if (begin != nullptr) {
            BamAlignmentParser::ParseBamAlignmentFromSliceSimpleTypes({begin, end}, reads, aIdx, fieldsOff[recIdx], fieldsSizes[recIdx]);
        } else if (begin == nullptr && !boundaries[recIdx].bin.empty()) {
            BamAlignmentParser::ParseBamAlignmentFromSliceSimpleTypes(boundaries[recIdx].bin, reads, aIdx, fieldsOff[recIdx], fieldsSizes[recIdx]);
        }
    });

    // Scan the sizes to know where each thread should write
    std::vector<VariableLenFieldsVar> outputOffsets(fieldsSizes.size());
    std::exclusive_scan(std::execution::par_unseq, fieldsSizes.begin(), fieldsSizes.end(), outputOffsets.begin(), VariableLenFieldsVar(0,0,0,0), [](VariableLenFieldsVar& s1, VariableLenFieldsVar& s2) {
        return VariableLenFieldsVar{s1.name + s2.name, s1.cigar + s2.cigar, s1.seq + s2.seq, s1.qual + s2.qual, s1.aux + s2.aux};
    });

    reads.SEQ.ResizeData(outputOffsets.back().seq);
    reads.QUAL.ResizeData(outputOffsets.back().qual);
    reads.CIGAR.ResizeData(outputOffsets.back().cigar);
    OptionalFields& opt = reads.TAGS;
    opt.tagArena.resize(outputOffsets.back().aux);  // raw arena
    //opt.RG.resize(boundaries.size());
    //opt.tagsBlobOff.assign(N, std::numeric_limits<size_t>::max());
    //opt.tagsBlobSize.assign(N, 0);
    //for (size_t recIdx = 0; recIdx < boundaries.size(); ++recIdx) {

    tbb::parallel_for (0UL, boundaries.size(), [&boundaries, &reads, &fieldsOff, &fieldsSizes, &outputOffsets](size_t recIdx) {

        auto begin = boundaries[recIdx].chunkPtr;
        auto end = begin + boundaries[recIdx].size;
        size_t aIdx  = boundaries[recIdx].alnIdx;

        if (begin != nullptr) {
            BamAlignmentParser::ParseBamAlignmentFromSliceComplexTypes({begin, end}, reads, aIdx, fieldsOff[recIdx], fieldsSizes[recIdx], outputOffsets[recIdx]);
        } else if (begin == nullptr && !boundaries[recIdx].bin.empty()) {
            BamAlignmentParser::ParseBamAlignmentFromSliceComplexTypes(boundaries[recIdx].bin, reads, aIdx, fieldsOff[recIdx], fieldsSizes[recIdx], outputOffsets[recIdx]);
        }
    }
    );
    //profiler.stop(IOStageId::SAMStructUpdates);

    bp->n_processed = alnsRead; // keep your cumulative reporting
    //profiler.stop(IOStageId::ProcessBatch);
}


ssize_t BamReader::ParseAlignment(std::vector<SamAlignments>& reads, std::shared_ptr<Header> &header) {
    reads.clear();
    reads.emplace_back(header);
    if (!header_done_) {
        skip_bam_header(in);
        header_done_ = true;
    }
    auto batch_process  = std::make_shared<InputBatchBam>(DEFAULT_BUFFER_SIZE);
    auto batch_transfer = std::make_shared<InputBatchBam>(DEFAULT_BUFFER_SIZE);

    Carry carry;
    size_t alnsRead    = 0;
    size_t batchIdx    = 0;
    bool   first_batch = true;

    // Preload first batch synchronously
    load_batch(batch_process, in);

    // --- Thread sync primitives ---
    std::mutex m;
    std::condition_variable cv;
    bool ready = false;       // main -> loader
    bool done  = false;       // loader -> main
    bool stop  = false;       // shutdown signal
    std::exception_ptr loader_ex;

    // --- Loader thread ---
    std::thread loader([&] {
        try {
            for (;;) {
                // Wait until main thread signals "ready"
                std::unique_lock lk(m);
                cv.wait(lk, [&]{ return ready || stop; });
                if (stop) break;

                // reset flags
                ready = false;
                lk.unlock();

                // do the work
                load_batch(batch_transfer, in);

                // signal completion
                lk.lock();
                done = true;
                lk.unlock();
                cv.notify_one();
            }
        } catch (...) {
            loader_ex = std::current_exception();
            std::unique_lock lk(m);
            done = true;
            cv.notify_one();
        }
    });

    // --- Main processing loop ---
    while (first_batch || batch_process->buffer.has_data() ||
           (carry.stage == Carry::ReadBody && carry.need > 0)) {

        // Tell loader to load next batch
        {
            std::lock_guard lk(m);
            ready = true;
            done  = false;
        }
        cv.notify_one();

        // Process current batch here
        process_batch(batch_process, reads[batchIdx], alnsRead, carry);

        // Wait for loader to finish
        {
            std::unique_lock lk(m);
            cv.wait(lk, [&]{ return done || loader_ex != nullptr; });
        }
        if (loader_ex) std::rethrow_exception(loader_ex);

        // Transfer counters
        batch_transfer->header_to_process = batch_process->header_to_process;
        batch_transfer->remaining         = batch_process->remaining;
        batch_transfer->n_processed       = batch_process->n_processed;

        // Prepare next iteration
        reads.back().SetAlnsNum(alnsRead);
        reads.emplace_back(header);
        ++batchIdx;
        alnsRead = 0;
        std::swap(batch_process, batch_transfer);
        first_batch = false;
    }

    // --- Shutdown loader thread ---
    {
        std::lock_guard lk(m);
        stop = true;
    }
    cv.notify_one();
    loader.join();

    for (auto& block : reads) {
        block.Resize(block.GetAlnsNum());
    }
    //reads.GetHdr().QnameMapper().finalize(reads.QNAME);
    //auto& profiler = ProfilerIO::instance();
    //profiler.report();
    return static_cast<ssize_t>(alnsRead);

}

std::string_view BamReader::getslice(size_t numBytesToRead, size_t &numBytesRead) {
    return buffer.getslice(numBytesToRead, numBytesRead);
}

