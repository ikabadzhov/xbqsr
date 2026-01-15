// BamBatchFormatter.cpp
#include "BamBatchFormatter.h"

void BamBatchFormatter::FormatHeaderTo(Header& hdr, std::vector<uint8_t>& out) {
    std::call_once(header_once_, [&] {
        std::string sam_text;
        SamHeaderFormatter::FormatSamHeader(hdr, sam_text);

        std::vector<uint8_t> buf;
        buf.reserve(4096 + sam_text.size());

        // magic "BAM\1"
        buf.insert(buf.end(), {'B','A','M','\1'});

        // l_text + text
        push_u32(buf, (uint32_t)sam_text.size());
        buf.insert(buf.end(), sam_text.begin(), sam_text.end());

        // references
        const auto& seqs = hdr.getSequenceDictionary().getAll();
        push_u32(buf, (uint32_t)seqs.size());
        for (const auto& r : seqs) {
            const std::string& name = r.getSequenceName();
            push_u32(buf, (uint32_t)(name.size() + 1));
            buf.insert(buf.end(), name.begin(), name.end());
            buf.push_back('\0');
            push_u32(buf, (uint32_t)r.getSequenceLength());
        }

        out = std::move(buf); // overwrite
        header_written_.store(true, std::memory_order_release);
    });
}

void BamBatchFormatter::FormatAlignmentBatchTo(SamAlignments& reads, std::vector<uint8_t>& out) {
    const size_t total = (size_t)reads.GetAlnsNum();
    out.clear();
    if (total == 0) return;

    std::vector<size_t> sizes(total);
    tbb::parallel_for(size_t(0), total, [&](size_t i) {
        sizes[i] = record_wire_size(reads, i);
    });

    std::vector<size_t> offsets(total);
    std::exclusive_scan(std::execution::par, sizes.begin(), sizes.end(), offsets.begin(), size_t(0));

    const size_t bytes = offsets.back() + sizes.back();
    out.resize(bytes);

    uint8_t* base = out.data();
    tbb::parallel_for(size_t(0), total, [&](size_t i) {
        write_record_into(base + offsets[i], sizes[i], reads, i);
    });
}

void BamBatchFormatter::FormatAlignmentBatchTo(SamAlignments& reads, OutputBatchBam& out) {
    const size_t total = (size_t)reads.GetAlnsNum();
    out.reset();
    if (total == 0) return;

    std::vector<size_t> sizes(total);
    tbb::parallel_for(size_t(0), total, [&](size_t i) {
        sizes[i] = record_wire_size(reads, i);
    });

    std::vector<size_t> offsets(total);
    std::exclusive_scan(std::execution::par, sizes.begin(), sizes.end(), offsets.begin(), size_t(0));

    const size_t bytes = offsets.back() + sizes.back();
    out.ensure_capacity(bytes);

    uint8_t* base = out.data();
    tbb::parallel_for(size_t(0), total, [&](size_t i) {
        write_record_into(base + offsets[i], sizes[i], reads, i);
    });

    out.seal(bytes, total);
}