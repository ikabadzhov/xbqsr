#include <tbb/parallel_for.h>
#include "ComputeCovariateVectorsData.h"
#include "../../bqsr/BQSR.h"
#include "../../timer/Timer.h"
#include "../../logger/Logger.h"
#include "../../sycl_parallel_library/sorting/mergesort.h"


static size_t compute_chunk_sz(const SequencesInput& in, const ComputeSkipSliceKernelData& kernel_data, size_t start_idx, size_t n_elements);


SequencesInput::SequencesInput(mem_t<Interval>& known_sites, mem_t<SpanFieldLong>& known_sites_ptr, mem_t<uint8_t>& reference, mem_t<SpanFieldLong>& reference_ptr, bqsrEntry* quality_table, bqsrEntry* cycle_table, bqsrEntry* context_table, size_t quality_table_sz, size_t cycle_table_sz, size_t context_table_sz, BaseRecalibratorTablesAccessor& table_acc, char* sequence, SpanField* sequence_ptr, char* qualities, SpanField* qualities_ptr, CigarOperation* cigar, SpanField* cigar_ptr, int32_t* pos, int32_t* rnext, int32_t* pnext, int32_t* tlen, int32_t* rname, uint8_t* recalibrated, uint16_t* flags, int32_t* read_group_cov, size_t n_elements) : InputData(n_elements) {
        h_known_sites = known_sites.data();
        h_known_sites_ptr = known_sites_ptr.data();
        h_known_sites_ptr_sz = known_sites_ptr.size();
        h_reference = reference.data();
        h_reference_ptr = reference_ptr.data();
        h_reference_ptr_sz = reference_ptr.size();

        s_quality_table = quality_table;
        s_cycle_table = cycle_table;
        s_context_table = context_table;

        s_quality_table_sz = quality_table_sz;
        s_cycle_table_sz = cycle_table_sz;
        s_context_table_sz = context_table_sz;

        h_sequences = sequence;
        h_sequences_ptr = sequence_ptr;
        h_cigar = cigar;
        h_cigar_ptr = cigar_ptr;
        h_qualities = qualities;
        h_qualities_ptr = qualities_ptr;
        h_rnext = rnext;
        h_pnext = pnext;
        h_tlen = tlen;
        h_pos = pos;
        h_rname = rname;
        h_recalibrated = recalibrated;
        h_flag = flags;
        h_read_group_cov = read_group_cov;
        h_table_acc = table_acc;
        /*h_sz_qual_score_cov = sz_qual_cov;
        h_sz_cycle_cov = sz_cycle_cov;
        h_sz_context_cov = sz_context_cov;*/
}

void ComputeSkipSliceKernelData::swap(ComputeSkipSliceKernelData &rhs) {
    std::swap(_q, rhs._q);
    std::swap(_space, rhs._space);
    //std::swap(max_num_elements, rhs.max_num_elements);
    //std::swap(is_first_copy, rhs.is_first_copy);
    std::swap(d_known_sites, rhs.d_known_sites);
    std::swap(d_known_sites_sz, rhs.d_known_sites_sz);
    std::swap(d_known_sites_ptr, rhs.d_known_sites_ptr);

    std::swap(d_reference, rhs.d_reference);
    std::swap(d_reference_ptr, rhs.d_reference_ptr);
    std::swap(d_reference_sz, rhs.d_reference_sz);

    std::swap(memory_buffer, rhs.memory_buffer);
    std::swap(max_memory, rhs.max_memory);

    std::swap(d_sequences, rhs.d_sequences);
    std::swap(d_sequences_ptr, rhs.d_sequences_ptr);

    std::swap(d_cigar, rhs.d_cigar);
    std::swap(d_clipped_cigar, rhs.d_clipped_cigar);
    std::swap(d_cigar_ptr, rhs.d_cigar_ptr);
    std::swap(d_clipped_cigar_start, rhs.d_clipped_cigar_start);
    std::swap(d_clipped_cigar_len, rhs.d_clipped_cigar_len);

    std::swap(d_qualities, rhs.d_qualities);
    //std::swap(d_qualities_ptr, rhs.d_qualities_ptr);

    std::swap(d_pos, rhs.d_pos);
    std::swap(d_rname, rhs.d_rname);

    std::swap(d_rnext, rhs.d_rnext);
    std::swap(d_pnext, rhs.d_pnext);
    std::swap(d_tlen, rhs.d_tlen);

    std::swap(d_recalibrated, rhs.d_recalibrated);
    std::swap(d_flag, rhs.d_flag);
    std::swap(d_reads_group_covariate, rhs.d_reads_group_covariate);

    std::swap(d_skip_slice, rhs.d_skip_slice);
    //std::swap(d_skip_slice_ptr, rhs.d_skip_slice_ptr);

    std::swap(d_snps, rhs.d_snps);
    //std::swap(d_snps_ptr, rhs.d_snps_ptr);

    std::swap(d_stranded_seq, rhs.d_stranded_seq);
    std::swap(d_stranded_seq_ptr, rhs.d_stranded_seq_ptr);

    std::swap(d_base_context, rhs.d_base_context);
    std::swap(d_base_context_ptr, rhs.d_base_context_ptr);

    std::swap(d_table_acc, rhs.d_table_acc);

    std::swap(d_quality_table_indexes, rhs.d_quality_table_indexes);
    std::swap(d_cycle_table_indexes, rhs.d_cycle_table_indexes);
    std::swap(d_context_table_indexes, rhs.d_context_table_indexes);

    std::swap(d_quality_table_vals, rhs.d_quality_table_vals);
    std::swap(d_cycle_table_vals, rhs.d_cycle_table_vals);
    std::swap(d_context_table_vals, rhs.d_context_table_vals);

    /*
    std::swap(d_segments_quality, rhs.d_segments_quality);
    std::swap(d_segments_cycle, rhs.d_segments_cycle);
    std::swap(d_segments_context, rhs.d_segments_context);*/

    std::swap(prev_batch_max_index, rhs.prev_batch_max_index);
    std::swap(_n_expanded_cigar_ops, rhs._n_expanded_cigar_ops);
}

ComputeSkipSliceKernelData::~ComputeSkipSliceKernelData() {
    if (memory_buffer) {
        sycl::free(memory_buffer, _q);
        memory_buffer = nullptr;
    }
    _q.wait_and_throw();
}

ComputeSkipSliceKernelData::ComputeSkipSliceKernelData(ComputeSkipSliceKernelData &&rhs) noexcept {
    swap(rhs);
}

ComputeSkipSliceKernelData &ComputeSkipSliceKernelData::operator=(ComputeSkipSliceKernelData &&rhs) noexcept {
    swap(rhs);
    return *this;
}

//[[clang::optnone]]
void ComputeSkipSliceKernelData::memcpy_host_to_dev(sycl::queue& transfer_queue, ComputeSkipSliceKernelData& dest, size_t dest_offset, const SequencesInput& src,
                               size_t src_offset, size_t n_elements) {
    Timer t;
    t.reset();
    size_t actual_data_movement = 0;
    try {
        if (dest.is_first_copy) {
            dest.d_known_sites_ptr = src.h_known_sites_ptr; // SHARED MEMORY HOST-DEV
            dest.d_reference_ptr = src.h_reference_ptr;
            dest.d_known_sites = src.h_known_sites;
            dest.d_reference = src.h_reference;

            dest.d_quality_table = src.s_quality_table;
            dest.d_cycle_table = src.s_cycle_table;
            dest.d_context_table = src.s_context_table;

            /*dest.d_sz_qual_score_cov = src.h_sz_qual_score_cov;
            dest.d_sz_cycle_cov = src.h_sz_cycle_cov;
            dest.d_sz_context_cov = src.h_sz_context_cov;
*/
            dest.d_table_acc = src.h_table_acc;
            dest.d_quality_table_sz = src.s_quality_table_sz;
            dest.d_cycle_table_sz = src.s_cycle_table_sz;
            dest.d_context_table_sz = src.s_context_table_sz;

            dest.is_first_copy = false;
            transfer_queue.wait();
        }

        size_t offset = 0;
        if (n_elements > 0) {
            size_t chunk_sz = 0;

            char *base = static_cast<char *>(dest.memory_buffer);
            // 0. span offset
            {
                //transfer_queue.memcpy(dest.d_sequence_span_offset, &src.h_sequences_ptr[src_offset].start, sizeof(src.h_sequences_ptr[src_offset].start));
                dest.d_sequence_span_offset = src.h_sequences_ptr[src_offset].start;
                dest.d_cigar_span_offset = src.h_cigar_ptr[src_offset].start;
            }

            // 1. d_sequences
            {
                size_t chunk_offset = src.h_sequences_ptr[src_offset].start;
                chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                                  src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                                 (src.h_sequences_ptr[src_offset].start);
                dest.d_sequences = alloc_and_copy(transfer_queue, base, offset, src.h_sequences + chunk_offset, chunk_sz);
            }
            // 2. Allocate d_sequences_ptr (an array of SpanField)
            {
                dest.d_sequences_ptr = alloc_and_copy(transfer_queue, base, offset, src.h_sequences_ptr + src_offset, n_elements);
            }
            // 3. Allocate d_cigar (an array of CigarOperation)
            {
                chunk_sz = (src.h_cigar_ptr[src_offset + n_elements - 1].start +
                            src.h_cigar_ptr[src_offset + n_elements - 1].len) - (src.h_cigar_ptr[src_offset].start);
                size_t chunk_offset = src.h_cigar_ptr[src_offset].start;
                dest.d_cigar = alloc_and_copy(transfer_queue, base, offset, src.h_cigar + chunk_offset, chunk_sz);
            }
            // 4. Allocate d_cigar_span
            {
                dest.d_cigar_ptr = alloc_and_copy(transfer_queue, base, offset, src.h_cigar_ptr + src_offset, n_elements);
            }
            // 5. Allocate d_clipped_cigar
            {
                chunk_sz = ((src.h_cigar_ptr[src_offset + n_elements - 1].start +
                            src.h_cigar_ptr[src_offset + n_elements - 1].len) - (src.h_cigar_ptr[src_offset].start)) + (n_elements * PAD);
                dest.d_clipped_cigar = alloc_simple<CigarOperation>(base, offset, chunk_sz);
            }
            // 6 Allocate d_clipped_cigar_start (array of size_t)
            {
                dest.d_clipped_cigar_start =
                    alloc_simple<size_t>(base, offset, n_elements);
            }

            // 7 Allocate d_clipped_cigar_len (array of uint16_t)
            {
                dest.d_clipped_cigar_len =
                    alloc_simple<uint16_t>(base, offset, n_elements);
            }

            // 8. Allocate d_qualities (flat array by spans)
            {
                size_t chunk_offset = src.h_qualities_ptr[src_offset].start;
                chunk_sz = (src.h_qualities_ptr[src_offset + n_elements - 1].start +
                            src.h_qualities_ptr[src_offset + n_elements - 1].len) -
                           src.h_qualities_ptr[src_offset].start;

                dest.d_qualities =
                    alloc_and_copy(transfer_queue, base, offset,
                                   src.h_qualities + chunk_offset,   // src base
                                   chunk_sz);                         // number of chars
            }
            // 9. Allocate d_rnext (array of int32_t)
            {
                dest.d_rnext =
                    alloc_and_copy(transfer_queue, base, offset,
                                   src.h_rnext + src_offset,
                                   n_elements);
            }

            // 10. Allocate d_pnext (array of int32_t)
            {
                dest.d_pnext =
                    alloc_and_copy(transfer_queue, base, offset,
                                   src.h_pnext + src_offset,
                                   n_elements);
            }

            // 11. Allocate d_tlen (array of int32_t)
            {
                dest.d_tlen =
                    alloc_and_copy(transfer_queue, base, offset,
                                   src.h_tlen + src_offset,
                                   n_elements);
            }

            // 12. Allocate d_pos (array of int32_t)
            {
                dest.d_pos =
                    alloc_and_copy(transfer_queue, base, offset,
                                   src.h_pos + src_offset,
                                   n_elements);
            }

            // 13. Allocate d_recalibrated (array of int8_t)
            {
                dest.d_recalibrated = reinterpret_cast<int8_t*>(
                    alloc_and_copy(transfer_queue, base, offset,
                                   src.h_recalibrated + src_offset,
                                   n_elements));
            }

            // 14. Allocate d_rname (array of int32_t)
            {
                dest.d_rname =
                alloc_and_copy(transfer_queue, base, offset,
                               src.h_rname + src_offset,
                               n_elements);
            }

            // 15. Allocate d_flag (array of uint16_t)
            {
                dest.d_flag =
                    alloc_and_copy(transfer_queue, base, offset,
                                   src.h_flag + src_offset,
                                   n_elements);
            }

            // 16. Allocate d_reads_group_covariate (array of int32_t)
            {
                dest.d_reads_group_covariate =
                    alloc_and_copy(transfer_queue, base, offset,
                                   src.h_read_group_cov + src_offset,
                                   n_elements);
            }

            // 17. Allocate d_skip_slice (flat array sized like sequences slice; then fill 0)
            {
                size_t seq_chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                                       src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                                      src.h_sequences_ptr[src_offset].start;

                dest.d_skip_slice =
                    alloc_simple<uint8_t>(base, offset, seq_chunk_sz);

                // zero initialize
                transfer_queue.fill(dest.d_skip_slice, static_cast<uint8_t>(0), seq_chunk_sz);
            }

            // 18. Allocate d_snps (same size as sequences slice; then fill 0)
            {
                size_t seq_chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                                       src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                                      src.h_sequences_ptr[src_offset].start;

                dest.d_snps =
                    alloc_simple<int8_t>(base, offset, seq_chunk_sz);

                transfer_queue.fill(dest.d_snps, static_cast<int8_t>(0), seq_chunk_sz);
            }

            // 19. Allocate d_stranded_seq (flat array by sequences spans)
            {
                size_t chunk_offset = src.h_sequences_ptr[src_offset].start;
                size_t seq_chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                                       src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                                      src.h_sequences_ptr[src_offset].start;

                dest.d_stranded_seq =
                    alloc_and_copy(transfer_queue, base, offset,
                                   src.h_sequences + chunk_offset,
                                   seq_chunk_sz);
            }

            // 20. Allocate d_stranded_seq_ptr (array of SpanField)
            {
                /*dest.d_stranded_seq_ptr =
                    alloc_and_copy(transfer_queue, base, offset,
                                   src.h_sequences_ptr + src_offset,
                                   n_elements);*/
                dest.d_stranded_seq_ptr = alloc_simple<SpanField>(base, offset, n_elements);
            }

            // 21. Allocate d_base_context (flat array sized like sequences slice; no copy now)
            {
                size_t seq_chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                                       src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                                      src.h_sequences_ptr[src_offset].start;

                dest.d_base_context = alloc_simple<int32_t>(base, offset, seq_chunk_sz);
            }

            // 22. Allocate d_base_context_ptr (array of SpanField)
            {
                dest.d_base_context_ptr = alloc_simple<SpanField>(base, offset, n_elements);
            }

            // 23. Allocate d_quality_table_indexes (flat like sequences; no init/fill here)
            {
                size_t seq_chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                                       src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                                      src.h_sequences_ptr[src_offset].start;

                dest.d_quality_table_indexes =
                    alloc_simple<int32_t>(base, offset, seq_chunk_sz);
            }

            // 24. Allocate d_cycle_table_indexes (flat like sequences)
            {
                size_t seq_chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                                       src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                                      src.h_sequences_ptr[src_offset].start;

                dest.d_cycle_table_indexes =
                    alloc_simple<int32_t>(base, offset, seq_chunk_sz);
            }

            // 25. Allocate d_context_table_indexes (flat like sequences)
            {
                size_t seq_chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                                       src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                                      src.h_sequences_ptr[src_offset].start;

                dest.d_context_table_indexes =
                    alloc_simple<int32_t>(base, offset, seq_chunk_sz);
            }

            // 26. Allocate d_quality_table_vals (flat like sequences)
            {
                size_t seq_chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                                       src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                                      src.h_sequences_ptr[src_offset].start;

                dest.d_quality_table_vals =
                    alloc_simple<val_t>(base, offset, seq_chunk_sz);
            }

            // 27. Allocate d_cycle_table_vals (flat like sequences)
            {
                size_t seq_chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                                       src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                                      src.h_sequences_ptr[src_offset].start;

                dest.d_cycle_table_vals =
                    alloc_simple<val_t>(base, offset, seq_chunk_sz);
            }

            // 28. Allocate d_context_table_vals (flat like sequences)
            {
                size_t seq_chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                                       src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                                      src.h_sequences_ptr[src_offset].start;

                dest.d_context_table_vals =
                    alloc_simple<val_t>(base, offset, seq_chunk_sz);
            }

            /*// 29. Allocate d_segments_quality (flat like sequences)
            {
                size_t seq_chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                                       src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                                      src.h_sequences_ptr[src_offset].start;

                dest.d_segments_quality =
                    alloc_simple<segment_t>(base, offset, seq_chunk_sz);
            }

            // 30. Allocate d_segments_cycle (flat like sequences)
            {
                size_t seq_chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                                       src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                                      src.h_sequences_ptr[src_offset].start;

                dest.d_segments_cycle =
                    alloc_simple<segment_t>(base, offset, seq_chunk_sz);
            }

            // 31. Allocate d_segments_context (flat like sequences)
            {
                size_t seq_chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                                       src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                                      src.h_sequences_ptr[src_offset].start;

                dest.d_segments_context =
                    alloc_simple<segment_t>(base, offset, seq_chunk_sz);
            }*/

            // 32. Allocate d_temporary_segments_sequences_offsets (capped)
            {
                size_t tmp_elems = n_elements;
                dest.d_temporary_segments_sequences_offsets =
                    alloc_simple<int32_t>(base, offset, tmp_elems);
            }


            dest.table_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                                        src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                                       src.h_sequences_ptr[src_offset].start;

        }
        if (offset > dest.max_memory) {
            Logger::getLogger().log(LOG_ERROR, "ComputeCovariateVecotrsData: allocated (", offset, ") exceeds max memory (", dest.max_memory, ")\n");
            std::abort();
        }
        transfer_queue.wait_and_throw();

        Logger::getLogger().log(LOG_DEBUG, "    \tActual number of reads in batch: ", dest.max_num_elements);
    } catch (sycl::exception& e) {
        Logger::getLogger().log(LOG_ERROR, "Error during memory copy from host to device: ", e.what());
    }
    Logger::getLogger().log(LOG_DEBUG, "    [Load batch] Time: "+
                                      std::to_string(t.elapsed()) + " sec");
}
//[[clang::optnone]]
void ComputeSkipSliceKernelData::memcpy_dev_to_host(sycl::queue &transfer_queue, TablesOutput &dest, size_t dest_offset,
                                                    ComputeSkipSliceKernelData &src, size_t src_offset,
                                                    size_t n_elements) {

}

//[[clang::optnone]]
size_t bqsr_covariate_vector_kernel::update_batch_num_elements(const SequencesInput& in, const ComputeSkipSliceKernelData& kd, size_t start_idx) {
    if (start_idx >= in._n) return 0;

    const size_t Nmax = in._n - start_idx;

    // 1) Probe: can 1 element fit at all?
    size_t usage = compute_chunk_sz(in, kd, start_idx, 1);
    if (usage > kd.max_memory) {
        return 0;  // caller MUST handle n==0 (no indexing with n-1!)
    }

    // 2) Exponential search for an upper bound that DOESN'T fit
    size_t best = 1;   // last known good (fits)
    size_t high = 1;
    while (true) {
        size_t next = std::min(high * 2, Nmax);
        usage = compute_chunk_sz(in, kd, start_idx, next);
        if (usage <= kd.max_memory) {
            best = next;
            if (next == Nmax) { high = next; break; } // exhausted
            high = next;
        } else {
            // first non-fitting count is our upper bound
            high = next;
            break;
        }
    }

    // 3) Binary search in [best, high] for the max that fits.
    // Invariant: best fits; high may or may not fit.
    size_t low = best;
    while (low < high) {
        size_t mid = low + (high - low + 1) / 2; // upper mid to avoid infinite loop
        usage = compute_chunk_sz(in, kd, start_idx, mid);
        if (usage <= kd.max_memory) {
            low = mid;   // mid fits → move up
        } else {
            high = mid - 1; // mid doesn't fit → move down
        }
    }
    return low;  // >= 1 here
}
//[[clang::optnone]]
size_t compute_chunk_sz(const SequencesInput& in, const ComputeSkipSliceKernelData& kernel_data, size_t start_idx, size_t n_elements) {
    // When we compute this, we need to keep into account two things:
    // + The memory alignment of each pointer
    // + The fact that the flat vectors, such as sequences, qualities etc, may suffer from fragmentation
    //   i.e. few empty chars between two sequences

    size_t offset = 0;

    // 1. d_sequences
    size_t chunk_sz = (in.h_sequences_ptr[start_idx + n_elements - 1].start +
                       in.h_sequences_ptr[start_idx + n_elements - 1].len) -
                      (in.h_sequences_ptr[start_idx].start);
    simulate_alloc<decltype(in.h_sequences[0])>(offset, chunk_sz);

    // 2. d_sequences_ptr (an array of SpanField)
    simulate_alloc<decltype(kernel_data.d_sequences_ptr[0])>(offset, n_elements);

    // 3. d_cigar (an array of CigarOperation)
    chunk_sz = (in.h_cigar_ptr[start_idx + n_elements - 1].start +
                in.h_cigar_ptr[start_idx + n_elements - 1].len) -
               in.h_cigar_ptr[start_idx].start;
    simulate_alloc<decltype(kernel_data.d_cigar[0])>(offset, chunk_sz);

    // Cigar Span
    simulate_alloc<decltype(kernel_data.d_cigar_ptr[0])>(offset, n_elements);

    // Clipped cigar  (capacity = contiguous_ops + PAD * n_elements)
    {
        const size_t first    = in.h_cigar_ptr[start_idx].start;
        const size_t last_end = in.h_cigar_ptr[start_idx + n_elements - 1].start +
                                in.h_cigar_ptr[start_idx + n_elements - 1].len;
        const size_t contig_ops = last_end - first;
        const size_t clipped_cap = contig_ops + PAD * n_elements;
        simulate_alloc<decltype(kernel_data.d_clipped_cigar[0])>(offset, clipped_cap);
    }

    // 4.1 d_clipped_cigar_start (array)
    simulate_alloc<decltype(kernel_data.d_clipped_cigar_start[0])>(offset, n_elements);

    // 4.2 d_clipped_cigar_len (array)
    simulate_alloc<decltype(kernel_data.d_clipped_cigar_len[0])>(offset, n_elements);

    // 5. d_qualities (flat by spans)
    chunk_sz = (in.h_qualities_ptr[start_idx + n_elements - 1].start +
                in.h_qualities_ptr[start_idx + n_elements - 1].len) -
               in.h_qualities_ptr[start_idx].start;
    simulate_alloc<decltype(kernel_data.d_qualities[0])>(offset, chunk_sz);

    // 6. d_rnext / d_pnext / d_tlen (per-read)
    simulate_alloc<decltype(kernel_data.d_rnext[0])>(offset, n_elements);
    simulate_alloc<decltype(kernel_data.d_pnext[0])>(offset, n_elements);
    simulate_alloc<decltype(kernel_data.d_tlen[0])>(offset, n_elements);

    // 7. d_pos (per-read)
    simulate_alloc<decltype(kernel_data.d_pos[0])>(offset, n_elements);

    // 8. d_recalibrated (per-read)
    simulate_alloc<decltype(kernel_data.d_recalibrated[0])>(offset, n_elements);

    // 9. d_rname (per-read)
    simulate_alloc<decltype(kernel_data.d_rname[0])>(offset, n_elements);

    // 10. d_flag (per-read)
    simulate_alloc<decltype(kernel_data.d_flag[0])>(offset, n_elements);

    // 11. d_reads_group_covariate (per-read)
    simulate_alloc<decltype(kernel_data.d_reads_group_covariate[0])>(offset, n_elements);

    // 12. d_skip_slice (per-base like sequences)
    const size_t seq_chunk =
        (in.h_sequences_ptr[start_idx + n_elements - 1].start +
         in.h_sequences_ptr[start_idx + n_elements - 1].len) -
        in.h_sequences_ptr[start_idx].start;
    simulate_alloc<decltype(kernel_data.d_skip_slice[0])>(offset, seq_chunk);

    // 14. d_snps (per-base like sequences)
    simulate_alloc<decltype(kernel_data.d_snps[0])>(offset, seq_chunk);

    // 16. d_stranded_seq (per-base chars, like sequences)
    simulate_alloc<decltype(kernel_data.d_stranded_seq[0])>(offset, seq_chunk);

    // 17. d_stranded_seq_ptr (per-read)
    simulate_alloc<decltype(kernel_data.d_stranded_seq_ptr[0])>(offset, n_elements);

    // 18. d_base_context (per-base like sequences)
    simulate_alloc<decltype(kernel_data.d_base_context[0])>(offset, seq_chunk);

    // 19. d_base_context_ptr (per-read)
    simulate_alloc<decltype(kernel_data.d_base_context_ptr[0])>(offset, n_elements);

    // 20. d_quality_table_indexes (per-base like sequences)
    simulate_alloc<decltype(kernel_data.d_quality_table_indexes[0])>(offset, seq_chunk);

    // 21. d_cycle_table_indexes (per-base like sequences)
    simulate_alloc<decltype(kernel_data.d_cycle_table_indexes[0])>(offset, seq_chunk);

    // 22. d_context_table_indexes (per-base like sequences)
    simulate_alloc<decltype(kernel_data.d_context_table_indexes[0])>(offset, seq_chunk);

    // 23. d_quality_table_vals (per-base like sequences)
    simulate_alloc<decltype(kernel_data.d_quality_table_vals[0])>(offset, seq_chunk);

    // 24. d_cycle_table_vals (per-base like sequences)
    simulate_alloc<decltype(kernel_data.d_cycle_table_vals[0])>(offset, seq_chunk);

    // 25. d_context_table_vals (per-base like sequences)
    simulate_alloc<decltype(kernel_data.d_context_table_vals[0])>(offset, seq_chunk);

    // Temporary segment sequences offsets (per-read; cap if you need)
    simulate_alloc<decltype(kernel_data.d_temporary_segments_sequences_offsets[0])>(offset, n_elements);

    // Segments (per-base like sequences)
    /*simulate_alloc<decltype(kernel_data.d_segments_quality[0])>(offset, seq_chunk);
    simulate_alloc<decltype(kernel_data.d_segments_cycle[0])>(offset, seq_chunk);
    simulate_alloc<decltype(kernel_data.d_segments_context[0])>(offset, seq_chunk);*/
    return offset;
}