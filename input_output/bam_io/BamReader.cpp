#include <future>
#include <execution>
#include <condition_variable>
#include <stop_token>
#include <mutex>
#include "BamReader.h"
#include "BamUtilsIO.h"
#include "../parser/HeaderParser.h"
#include "../parser/BamAlignmentParser.h"
#include "../../timer/TraceTimeline.h"

void BamReader::Close() {
    this->in->close();
}

// Parse BAM header: fills header structure in reads.
// Throws std::runtime_error on invalid/truncated headers.
// Parse BAM header: fills your header structure in reads.
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
    bt->buffer.fill(in);
};

// --- Carry state for a single BAM record across batches ---
struct Carry {
    enum Stage : uint8_t { ReadLen, ReadBody } stage = ReadLen;
    int32_t need  = sizeof(int32_t);    // bytes still needed for current stage
    int32_t block = -1;                 // block_size after decoding
    std::string bin;                    // accumulation buffer
    int32_t max_block_seen = 0;         // track largest block

    Carry() { bin.reserve(4096); }  // small upfront reserve reduces churn
} ;


class RecordSpan {
public:
    const char* chunkPtr = nullptr;
    size_t size = 0;
    std::string bin;
    size_t alnIdx = 0;
    RecordSpan(const char* ptr, size_t size, size_t alnIdx, const std::string& bin = {}) : chunkPtr(ptr), size(size), bin(bin), alnIdx(alnIdx) {};
};

struct BatchTmpVector {
    std::vector<RecordSpan> boundaries;
    std::vector<VariableLenFieldsVar> fieldsOff;
    std::vector<VariableLenFieldsVar> fieldsSizes;
    std::vector<VariableLenFieldsVar> outputOffsets;
};

//[[clang::optnone]]
[[maybe_unused]]
static void process_batch_parallel(const std::shared_ptr<InputBatchBam>& bp,
                          SamAlignments& reads,
                          size_t& alnsRead,
                          Carry& carry,
                          BatchTmpVector& scratch) {

    auto& buf        = bp->buffer;
    auto& boundaries = scratch.boundaries;

    boundaries.clear();

    // Initial guess: 1024 records per batch; capacity only grows, never shrinks.
    if (boundaries.capacity() < 1024)
        boundaries.reserve(1024);

    while (buf.has_data()) {
        const size_t avail = buf.available_data();
        if (avail == 0) break;

        const size_t take = static_cast<size_t>(
            std::min<int64_t>(carry.need, static_cast<int64_t>(avail))
        );
        size_t got = 0;

        std::string_view chunk = buf.getslice(take, got);
        if (got == 0) break;

        carry.need -= static_cast<int32_t>(got);
        // If we were already accumulating (spanning case) and just finished, append final bytes now.
        if (!carry.bin.empty() && carry.need == 0) {
            carry.bin.append(chunk.data(), got);
        }

        if (carry.need > 0) {
            // Still in ReadLen/ReadBody state: accumulate partial data and bail to loader
            carry.bin.append(chunk.data(), got);
            break;
        }

        if (carry.stage == Carry::ReadLen) {
            const unsigned char* p = carry.bin.empty()
                    ? reinterpret_cast<const unsigned char*>(chunk.data())
                    : reinterpret_cast<const unsigned char*>(carry.bin.data());

            int32_t block_size = p[0] | (p[1] << 8) | (p[2] << 16) | (p[3] << 24);

            if (block_size <= 0 || block_size > (64 * 1024 * 1024)) {
                throw std::runtime_error("Invalid BAM block_size");
            }

            carry.block = block_size;
            carry.need  = block_size;
            carry.stage = Carry::ReadBody;

            // Reserve capacity once when we see a new maximum, but don't reserve every time
            if (block_size > carry.max_block_seen) {
                carry.max_block_seen = block_size;
                carry.bin.reserve(block_size);
            }

            carry.bin.clear();  // Clear for next ReadBody stage
            continue;
        }

        // ReadBody complete: we have exactly one full BAM record
        if (carry.bin.empty()) {
            // Record is contiguous in this chunk
            if (boundaries.size() == boundaries.capacity()) {
                boundaries.reserve(boundaries.capacity() * 2);
            }

            boundaries.emplace_back(chunk.data(), chunk.size(), alnsRead);
        } else {
            // Record spans previous and current chunks → stored in carry.bin
            carry.bin.append(chunk.data(), got);
            RecordSpan rec(nullptr, 0, alnsRead, carry.bin);

            if (boundaries.size() == boundaries.capacity())
                boundaries.reserve(boundaries.capacity() * 2);

            boundaries.push_back(std::move(rec));
            carry.bin.clear();
        }

        ++alnsRead;

        // Prepare for next record
        carry.bin.clear();
        carry.stage = Carry::ReadLen;
        carry.need  = sizeof(int32_t);
        carry.block = -1;
    }

    const size_t nrecs = boundaries.size();
    if (nrecs == 0) {
        bp->n_processed = alnsRead;
        return;
    }
   // reads.EnsureCapacity(nrecs);
    if (nrecs > reads.Capacity()) {
        reads.EnsureCapacity(nrecs);
    }
    reads.Resize(nrecs);

    auto& fieldsOff   = scratch.fieldsOff;
    auto& fieldsSizes = scratch.fieldsSizes;

    // Keep capacity, only grow if needed
    if (fieldsOff.capacity()   < nrecs)   fieldsOff.reserve(nrecs);
    if (fieldsSizes.capacity() < nrecs+1) fieldsSizes.reserve(nrecs+1);

    fieldsOff.resize(nrecs);
    fieldsSizes.resize(nrecs + 1);  // +1 for scan sentinel

    tbb::parallel_for(0UL, nrecs, [&boundaries, &reads, &fieldsOff, &fieldsSizes](size_t recIdx) {
        auto begin = boundaries[recIdx].chunkPtr;
        auto end   = begin + boundaries[recIdx].size;
        size_t aIdx = boundaries[recIdx].alnIdx;

        if (begin != nullptr) {
            BamAlignmentParser::ParseBamAlignmentFromSliceSimpleTypes(
                {begin, end}, reads, aIdx,
                fieldsOff[recIdx], fieldsSizes[recIdx]
            );
        } else if (!boundaries[recIdx].bin.empty()) {
            BamAlignmentParser::ParseBamAlignmentFromSliceSimpleTypes(
                boundaries[recIdx].bin, reads, aIdx,
                fieldsOff[recIdx], fieldsSizes[recIdx]
            );
        }
    });

    auto& outputOffsets = scratch.outputOffsets;
    if (outputOffsets.capacity() < nrecs + 1)
        outputOffsets.reserve(nrecs + 1);

    outputOffsets.resize(nrecs + 1);

    std::exclusive_scan(
        std::execution::par_unseq,
        fieldsSizes.begin(), fieldsSizes.end(),
        outputOffsets.begin(),
        VariableLenFieldsVar(0,0,0,0,0),
        [](const VariableLenFieldsVar& s1, const VariableLenFieldsVar& s2) {
            return VariableLenFieldsVar{
                s1.qname + s2.qname,
                s1.cigar + s2.cigar,
                s1.seq   + s2.seq,
                s1.qual  + s2.qual,
                s1.aux   + s2.aux
            };
        }
    );

    if (outputOffsets.back().seq > reads.SEQ.CapacityData()) {
        reads.SEQ.ReserveData(outputOffsets.back().seq);//.EnsureCapacity(nrecs * (1<<20));
    }
    if (outputOffsets.back().qual > reads.QUAL.CapacityData()) {
        reads.QUAL.ReserveData(outputOffsets.back().qual);//.EnsureCapacity(nrecs * (1<<20));
    }
    if (outputOffsets.back().cigar > reads.CIGAR.CapacityData()) {
        reads.CIGAR.ReserveData(outputOffsets.back().cigar);//.EnsureCapacity(nrecs * (1<<20));
    }
    if (outputOffsets.back().qname > reads.QNAME.CapacityData()) {
        reads.QNAME.ReserveData(outputOffsets.back().qname);//.EnsureCapacity(nrecs * (1<<20));
    }
    reads.SEQ.ResizeData(outputOffsets.back().seq);
    reads.QUAL.ResizeData(outputOffsets.back().qual);
    reads.CIGAR.ResizeData(outputOffsets.back().cigar);
    reads.QNAME.ResizeData(outputOffsets.back().qname);

    OptionalFields& opt = reads.TAGS;

    if (outputOffsets.back().aux > reads.TAGS.CapacityData()) {
        reads.TAGS.ReserveData(outputOffsets.back().aux);//.EnsureCapacity(nrecs * (1<<20));
    }
    opt.tagArena.resize(outputOffsets.back().aux);

    tbb::parallel_for (0UL, nrecs, [&](size_t recIdx) {
        auto begin = boundaries[recIdx].chunkPtr;
        auto end   = begin + boundaries[recIdx].size;
        size_t aIdx = boundaries[recIdx].alnIdx;

        if (begin != nullptr) {
            BamAlignmentParser::ParseBamAlignmentFromSliceComplexTypes(
                {begin, end}, reads, aIdx,
                fieldsOff[recIdx], fieldsSizes[recIdx], outputOffsets[recIdx]
            );
        } else if (!boundaries[recIdx].bin.empty()) {
            BamAlignmentParser::ParseBamAlignmentFromSliceComplexTypes(
                boundaries[recIdx].bin, reads, aIdx,
                fieldsOff[recIdx], fieldsSizes[recIdx], outputOffsets[recIdx]
            );
        }
    });

    bp->n_processed = alnsRead;
}

static void process_batch(const std::shared_ptr<InputBatchBam>& bp,
                          SamAlignments& reads,
                          size_t& alnsRead,
                          Carry& carry,
                          BatchTmpVector& scratch) {

    auto& buf= bp->buffer;
    auto& boundaries = scratch.boundaries;

    boundaries.clear();

    // Initial guess: 1024 records per batch; capacity only grows, never shrinks.
    if (boundaries.capacity() < 1024)
        boundaries.reserve(1024);

    while (buf.has_data()) {
        const size_t avail = buf.available_data();
        if (avail == 0) break;

        const size_t take = static_cast<size_t>(
            std::min<int64_t>(carry.need, static_cast<int64_t>(avail))
        );
        size_t got = 0;

        std::string_view chunk = buf.getslice(take, got);
        if (got == 0) break;

        carry.need -= static_cast<int32_t>(got);
        // If we were already accumulating (spanning case) and just finished, append final bytes now.
        if (!carry.bin.empty() && carry.need == 0) {
            carry.bin.append(chunk.data(), got);
        }

        if (carry.need > 0) {
            // Still in ReadLen/ReadBody state: accumulate partial data and bail to loader
            carry.bin.append(chunk.data(), got);
            break;
        }

        if (carry.stage == Carry::ReadLen) {
            const unsigned char* p = carry.bin.empty()
                    ? reinterpret_cast<const unsigned char*>(chunk.data())
                    : reinterpret_cast<const unsigned char*>(carry.bin.data());

            int32_t block_size = p[0] | (p[1] << 8) | (p[2] << 16) | (p[3] << 24);

            if (block_size <= 0 || block_size > (64 * 1024 * 1024)) {
                throw std::runtime_error("Invalid BAM block_size");
            }

            carry.block = block_size;
            carry.need  = block_size;
            carry.stage = Carry::ReadBody;

            // Reserve capacity once when we see a new maximum, but don't reserve every time
            if (block_size > carry.max_block_seen) {
                carry.max_block_seen = block_size;
                carry.bin.reserve(block_size);
            }

            carry.bin.clear();  // Clear for next ReadBody stage
            continue;
        }

        // ReadBody complete: we have exactly one full BAM record
        if (carry.bin.empty()) {
            // Record is contiguous in this chunk
            if (boundaries.size() == boundaries.capacity()) {
                boundaries.reserve(boundaries.capacity() * 2);
            }

            boundaries.emplace_back(chunk.data(), chunk.size(), alnsRead);
        } else {
            // Record spans previous and current chunks: stored in carry.bin
            carry.bin.append(chunk.data(), got);
            RecordSpan rec(nullptr, 0, alnsRead, carry.bin);

            if (boundaries.size() == boundaries.capacity())
                boundaries.reserve(boundaries.capacity() * 2);

            boundaries.push_back(std::move(rec));
            carry.bin.clear();
        }

        ++alnsRead;

        // Prepare for next record
        carry.bin.clear();
        carry.stage = Carry::ReadLen;
        carry.need  = sizeof(int32_t);
        carry.block = -1;
    }

    const size_t nrecs = boundaries.size();
    if (nrecs == 0) {
        bp->n_processed = alnsRead;
        return;
    }

    if (nrecs > reads.Capacity()) {
        reads.EnsureCapacity(nrecs);
    }
    reads.Resize(nrecs);

    auto& fieldsOff   = scratch.fieldsOff;
    auto& fieldsSizes = scratch.fieldsSizes;

    // capacities
    if (fieldsOff.capacity()   < nrecs)   fieldsOff.reserve(nrecs);
    if (fieldsSizes.capacity() < nrecs)   fieldsSizes.reserve(nrecs);

    fieldsOff.resize(nrecs);
    fieldsSizes.resize(nrecs);

    auto& outputOffsets = scratch.outputOffsets;
    if (outputOffsets.capacity() < nrecs + 1) outputOffsets.reserve(nrecs + 1);
    outputOffsets.resize(nrecs + 1);

    VariableLenFieldsVar tot{0,0,0,0,0};

    for (size_t recIdx = 0; recIdx < nrecs; ++recIdx) {
        outputOffsets[recIdx] = tot;   // Exclusive prefix for this record

        auto begin = boundaries[recIdx].chunkPtr;
        auto end   = begin + boundaries[recIdx].size;
        size_t aIdx = boundaries[recIdx].alnIdx;

        if (begin != nullptr) {
            BamAlignmentParser::ParseBamAlignmentFromSliceSimpleTypes(
                {begin, end}, reads, aIdx,
                fieldsOff[recIdx], fieldsSizes[recIdx]
            );
        } else if (!boundaries[recIdx].bin.empty()) {
            BamAlignmentParser::ParseBamAlignmentFromSliceSimpleTypes(
                boundaries[recIdx].bin, reads, aIdx,
                fieldsOff[recIdx], fieldsSizes[recIdx]
            );
        }

        // advance totals
        tot = VariableLenFieldsVar{
            tot.qname + fieldsSizes[recIdx].qname,
            tot.cigar + fieldsSizes[recIdx].cigar,
            tot.seq   + fieldsSizes[recIdx].seq,
            tot.qual  + fieldsSizes[recIdx].qual,
            tot.aux   + fieldsSizes[recIdx].aux
        };
    }

    outputOffsets[nrecs] = tot; // totals for allocation

    if (outputOffsets.back().seq > reads.SEQ.CapacityData()) {
        reads.SEQ.ReserveData(outputOffsets.back().seq);
    }
    if (outputOffsets.back().qual > reads.QUAL.CapacityData()) {
        reads.QUAL.ReserveData(outputOffsets.back().qual);
    }
    if (outputOffsets.back().cigar > reads.CIGAR.CapacityData()) {
        reads.CIGAR.ReserveData(outputOffsets.back().cigar);
    }
    if (outputOffsets.back().qname > reads.QNAME.CapacityData()) {
        reads.QNAME.ReserveData(outputOffsets.back().qname);
    }
    reads.SEQ.ResizeData(outputOffsets.back().seq);
    reads.QUAL.ResizeData(outputOffsets.back().qual);
    reads.CIGAR.ResizeData(outputOffsets.back().cigar);
    reads.QNAME.ResizeData(outputOffsets.back().qname);

    OptionalFields& opt = reads.TAGS;

    if (outputOffsets.back().aux > reads.TAGS.CapacityData()) {
        reads.TAGS.ReserveData(outputOffsets.back().aux);//.EnsureCapacity(nrecs * (1<<20));
    }
    opt.tagArena.resize(outputOffsets.back().aux);

    for (size_t recIdx = 0; recIdx < nrecs; ++recIdx) {
            auto begin = boundaries[recIdx].chunkPtr;
            auto end   = begin + boundaries[recIdx].size;
            size_t aIdx = boundaries[recIdx].alnIdx;

            if (begin != nullptr) {
                BamAlignmentParser::ParseBamAlignmentFromSliceComplexTypes(
                    {begin, end}, reads, aIdx,
                    fieldsOff[recIdx], fieldsSizes[recIdx], outputOffsets[recIdx]
                );
            } else if (!boundaries[recIdx].bin.empty()) {
                BamAlignmentParser::ParseBamAlignmentFromSliceComplexTypes(
                    boundaries[recIdx].bin, reads, aIdx,
                    fieldsOff[recIdx], fieldsSizes[recIdx], outputOffsets[recIdx]
                );
            }
        }

    bp->n_processed = alnsRead;
}

void BamReader::ParseAlignmentStreamed(std::shared_ptr<Header> header, std::function<SamAlignments()> get_batch,
                                       std::function<void(SamAlignments&&)> on_batch)
{
    int32_t reader_id = slices_[nextSlice].id;
    BQS_TRACE_SET_LANE(trace_tl::Role::Reader, (uint32_t)reader_id);
    BQS_TRACE_SCOPE("bam_reader:ParseAlignmentStreamed_total");

    if (!header_done_) {
        BQS_TRACE_SCOPE("bam_reader:skip_bam_header");
        skip_bam_header(in); header_done_ = true;
    }

    auto batch_process  = std::make_shared<InputBatchBam>(DEFAULT_BUFFER_SIZE);
    auto batch_transfer = std::make_shared<InputBatchBam>(DEFAULT_BUFFER_SIZE);

    Carry carry;
    carry.bin.reserve(4096);
    BatchTmpVector scratch;

    size_t alnsRead    = 0;
    bool   first_batch = true;

    if(nextSlice >= slices_.size()){
        return;
    }
    SeekVirtual(slices_[nextSlice].start, slices_[nextSlice].end);
    slice_id_process = slices_[nextSlice].id;
    nextSlice++;

    // Preload first batch
    {
        BQS_TRACE_SCOPE("bam_reader:load_batch_first");
        load_batch(batch_process, in);
    }

    std::mutex m;
    std::condition_variable_any cv;
    bool ready=false, done=false, stop=false;
    std::exception_ptr loader_ex;

    std::jthread loader([&] (std::stop_token st) {
        BQS_TRACE_SET_LANE(trace_tl::Role::Loader, (uint32_t)reader_id);
        try {
            for (;;) {
                {
                    BQS_TRACE_SCOPE("bam_reader:loader_wait_ready");
                    std::unique_lock lk(m);
                    cv.wait(lk, st, [&]{ return ready || stop; });
                    if (st.stop_requested() || stop) break;
                    ready = false;
                }

                for (;;) {
                    if (nextSlice >= slices_.size())
                        break;

                    // If this batch is empty, we must be starting a new slice
                    if (!batch_transfer->buffer.has_data()) {
                        BQS_TRACE_SCOPE("bam_reader:loader_seek_next_slice");
                        SeekVirtual(slices_[nextSlice].start, slices_[nextSlice].end);
                        slice_id_transfer = slices_[nextSlice].id;
                        ++nextSlice;
                    }

                    BQS_TRACE_SET_LANE(trace_tl::Role::Decompress, reader_id);
                    { BQS_TRACE_SCOPE("bam_reader:loader_load_batch"); load_batch(batch_transfer, in); }
                    BQS_TRACE_SET_LANE(trace_tl::Role::Loader, reader_id);

                    // If we got data → done for this iteration
                    if (batch_transfer->buffer.count() > 0) {
                        break;
                    }

                    // Else: read==0 → slice finished → loop to seek next slice
                }

                {
                    std::lock_guard lk(m);
                    done = true;
                }
                cv.notify_one();
            }
        } catch (...) {
            loader_ex = std::current_exception();
            std::unique_lock lk(m);
            done = true;
            cv.notify_one();
        }
    });

    auto current_batch = get_batch();//(header);
    current_batch.SetAlnsNum(0);

    while (first_batch || batch_process->buffer.has_data()
           || (carry.stage == Carry::ReadBody && carry.need > 0)) {
        { std::lock_guard lk(m); ready = true; done = false; }
        cv.notify_one();

        // Fill current_batch from batch_process
        {
            BQS_TRACE_SCOPE("bam_reader:process_batch");
            process_batch(batch_process, current_batch, alnsRead, carry, scratch);
        }

        // hand off if we produced anything
        if (alnsRead > 0) {
            current_batch.SetAlnsNum(alnsRead);
            current_batch.SetBatchId(slice_id_process);
            on_batch(std::move(current_batch));
            current_batch = get_batch();
            current_batch.SetAlnsNum(0);
            alnsRead = 0;
        }

        // Wait for loader
        {
            BQS_TRACE_SCOPE("bam_reader:wait_loader_done");
            std::unique_lock lk(m);
            cv.wait(lk, [&]{ return done || loader_ex != nullptr; });
        }
        if (loader_ex) std::rethrow_exception(loader_ex);

        // Swap in newly loaded data
        std::swap(batch_process, batch_transfer);
        std::swap(slice_id_process, slice_id_transfer);
        first_batch = false;
    }

    { std::lock_guard lk(m); stop = true; }
    loader.request_stop();
    cv.notify_one();
}

void BamReader::SeekVirtual(VOff voff_start, VOff voff_end) {
    auto* cis = dynamic_cast<CompressedInputStream*>(in.get());
    if (!cis) throw std::runtime_error("SeekVirtual requires CompressedInputStream (BGZF)");
    cis->set_virtual_end(voff_end);
    cis->seek_virtual(voff_start);
    header_done_ = true;
}

std::string_view BamReader::getslice(size_t numBytesToRead, size_t &numBytesRead) {
    return buffer.getslice(numBytesToRead, numBytesRead);
}


void BamReader::SetSlices(std::vector<Slice>&& slices) {
    slices_ = slices;
}