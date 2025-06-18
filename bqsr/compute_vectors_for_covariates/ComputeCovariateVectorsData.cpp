#include <tbb/parallel_for.h>
#include "ComputeCovariateVectorsData.h"
#include "../../bqsr/BQSR.h"
#include "../../timer/Timer.h"
#include "../../logger/Logger.h"
#include "../../sycl_parallel_library/sorting/mergesort.h"

static size_t compute_chunk_sz(const SequencesInput& in, const ComputeSkipSliceKernelData& kernel_data, size_t start_idx, size_t n_elements);


SequencesInput::SequencesInput(mem_t<Interval>& known_sites, mem_t<SpanFieldLong>& known_sites_ptr, mem_t<uint8_t>& reference, mem_t<SpanFieldLong>& reference_ptr, bqsrEntry* quality_table, bqsrEntry* cycle_table, bqsrEntry* context_table, int32_t sz_qual_cov, int32_t sz_cycle_cov, int32_t sz_context_cov, char* sequence, SpanField* sequence_ptr, char* qualities, SpanField* qualities_ptr, CigarOperation* cigar, size_t* cigar_start, uint16_t* cigar_len, int64_t* pos, int32_t* rname, uint8_t* recalibrated, uint16_t* flags, int32_t* read_group_cov, size_t n_elements) : InputData(n_elements) {
        h_known_sites = known_sites.data();
        h_known_sites_ptr = known_sites_ptr.data();
        h_known_sites_ptr_sz = known_sites_ptr.size();
        h_reference = reference.data();
        h_reference_ptr = reference_ptr.data();
        h_reference_ptr_sz = reference_ptr.size();

        s_quality_table = quality_table;
        s_cycle_table = cycle_table;
        s_context_table = context_table;

        h_sequences = sequence;
        h_sequences_ptr = sequence_ptr;
        h_cigar = cigar;
        h_cigar_start = cigar_start;
        h_cigar_len = cigar_len;
        h_qualities = qualities;
        h_qualities_ptr = qualities_ptr;
        h_pos = pos;
        h_rname = rname;
        h_recalibrated = recalibrated;
        h_flag = flags;
        h_read_group_cov = read_group_cov;
        h_sz_qual_score_cov = sz_qual_cov;
        h_sz_cycle_cov = sz_cycle_cov;
        h_sz_context_cov = sz_context_cov;
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
    std::swap(d_cigar_start, rhs.d_cigar_start);
    std::swap(d_cigar_len, rhs.d_cigar_len);

    std::swap(d_qualities, rhs.d_qualities);
    //std::swap(d_qualities_ptr, rhs.d_qualities_ptr);

    std::swap(d_pos, rhs.d_pos);
    std::swap(d_rname, rhs.d_rname);
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

    std::swap(d_sz_qual_score_cov, rhs.d_sz_qual_score_cov);
    std::swap(d_sz_cycle_cov, rhs.d_sz_cycle_cov);
    std::swap(d_sz_context_cov, rhs.d_sz_context_cov);

    std::swap(d_quality_table_indexes, rhs.d_quality_table_indexes);
    std::swap(d_cycle_table_indexes, rhs.d_cycle_table_indexes);
    std::swap(d_context_table_indexes, rhs.d_context_table_indexes);

    std::swap(d_quality_table_vals, rhs.d_quality_table_vals);
    std::swap(d_cycle_table_vals, rhs.d_cycle_table_vals);
    std::swap(d_context_table_vals, rhs.d_context_table_vals);

    std::swap(d_segments_quality, rhs.d_segments_quality);
    std::swap(d_segments_cycle, rhs.d_segments_cycle);
    std::swap(d_segments_context, rhs.d_segments_context);

//    std::swap(num_segments_quality, rhs.num_segments_quality);
//    std::swap(num_segments_cycle, rhs.num_segments_cycle);
//    std::swap(num_segments_context, rhs.num_segments_context);

    std::swap(d_read_index, rhs.d_read_index);
    std::swap(d_i_offset, rhs.d_i_offset);
    std::swap(d_j_offset, rhs.d_j_offset);
    //std::swap(d_table_indexes_ptr, rhs.d_table_indexes_ptr);

    std::swap(prev_batch_max_index, rhs.prev_batch_max_index);
    std::swap(_n_expanded_cigar_ops, rhs._n_expanded_cigar_ops);
//    std::swap(d_index_offset, rhs.d_index_offset);
}

ComputeSkipSliceKernelData::~ComputeSkipSliceKernelData() {

//    if (d_known_sites) {
//        sycl::free(d_known_sites, _q);
//        d_known_sites = nullptr;
//    }
//    if (d_known_sites_ptr) {
//        sycl::free(d_known_sites_ptr, _q);
//        d_known_sites_ptr = nullptr;
//    }
//    if (d_reference) {
//        sycl::free(d_reference, _q);
//        d_reference = nullptr;
//    }
//    if (d_reference_ptr) {
//        sycl::free(d_reference_ptr, _q);
//        d_reference_ptr = nullptr;
//    }
//    if (d_index_offset) {
//        sycl::free(d_index_offset, _q);
//        d_index_offset = nullptr;
//    }
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


void ComputeSkipSliceKernelData::memcpy_host_to_dev(sycl::queue& transfer_queue, ComputeSkipSliceKernelData& dest, size_t dest_offset, const SequencesInput& src,
                               size_t src_offset, size_t n_elements) {
    Timer t;
    t.reset();
    size_t actual_data_movement = 0;
    try {
        if (dest.is_first_copy) {
            dest.d_known_sites_ptr = src.h_known_sites_ptr; // SHARED MEMMORY HOST-DEV
            dest.d_reference_ptr = src.h_reference_ptr;
            dest.d_known_sites = src.h_known_sites;
            dest.d_reference = src.h_reference;

            dest.d_quality_table = src.s_quality_table;
            dest.d_cycle_table = src.s_cycle_table;
            dest.d_context_table = src.s_context_table;

            dest.d_sz_qual_score_cov = src.h_sz_qual_score_cov;
            dest.d_sz_cycle_cov = src.h_sz_cycle_cov;
            dest.d_sz_context_cov = src.h_sz_context_cov;
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
                dest.d_cigar_span_offset = src.h_cigar_start[src_offset];
            }

            // 1. d_sequences
            {
                dest.d_sequences = reinterpret_cast<char *>(base + offset);
                chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                                   src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                                  (src.h_sequences_ptr[src_offset].start);
                size_t chunk_offset = src.h_sequences_ptr[src_offset].start;
                size_t sequences_data_size = sizeof(src.h_sequences[0]) * chunk_sz;
                transfer_queue.memcpy(dest.d_sequences, src.h_sequences + chunk_offset, sequences_data_size).wait();
                offset += sequences_data_size;
                actual_data_movement += sequences_data_size;
            }
            // 2. Allocate d_sequences_ptr (an array of SpanField)
            //    Align the pointer properly.
            {
                size_t alignment = alignof(decltype(dest.d_sequences_ptr[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_sequences_ptr = reinterpret_cast<SpanField *>(base + offset);
                size_t sequences_ptr_size = sizeof(dest.d_sequences_ptr[0]) * n_elements;
                transfer_queue.memcpy(dest.d_sequences_ptr, src.h_sequences_ptr + src_offset,
                                      sequences_ptr_size);
                actual_data_movement += sequences_ptr_size;

                offset += sequences_ptr_size;
            }
            // 3. Allocate d_cigar (an array of CigarOperation)
            {
                size_t alignment = alignof(decltype(dest.d_cigar[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_cigar = reinterpret_cast<CigarOperation *>(base + offset);
                chunk_sz = (src.h_cigar_start[src_offset + n_elements - 1] +
                            src.h_cigar_len[src_offset + n_elements - 1]) - (src.h_cigar_start[src_offset]);
                size_t chunk_offset = src.h_cigar_start[src_offset];
                size_t cigar_data_size = sizeof(src.h_cigar[0]) * chunk_sz;
                transfer_queue.memcpy(dest.d_cigar, src.h_cigar + chunk_offset, cigar_data_size);
                actual_data_movement += cigar_data_size;
                offset += cigar_data_size;
            }

            // 4.1 Allocate d_cigar_start (an array of SpanField)
            {
                size_t alignment = alignof(decltype(dest.d_cigar_start[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_cigar_start = reinterpret_cast<size_t*>(base + offset);
                size_t cigar_ptr_size = sizeof(dest.d_cigar_start[0]) * n_elements;
                transfer_queue.memcpy(dest.d_cigar_start, src.h_cigar_start + src_offset, cigar_ptr_size);
                actual_data_movement += cigar_ptr_size;
                offset += cigar_ptr_size;
            }
            // 4.2 Allocate d_cigar_len (an array of SpanField)
            {
                size_t alignment = alignof(decltype(dest.d_cigar_len[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_cigar_len = reinterpret_cast<uint16_t*>(base + offset);
                size_t cigar_ptr_size = sizeof(dest.d_cigar_len[0]) * n_elements;
                transfer_queue.memcpy(dest.d_cigar_len, src.h_cigar_len + src_offset, cigar_ptr_size);
                actual_data_movement += cigar_ptr_size;
                offset += cigar_ptr_size;
            }

            // 5. Allocate d_qualities (an array of char)
            {
                size_t alignment = alignof(decltype(dest.d_qualities[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_qualities = reinterpret_cast<char *>(base + offset);
                chunk_sz = (src.h_qualities_ptr[src_offset + n_elements - 1].start +
                            src.h_qualities_ptr[src_offset + n_elements - 1].len) -
                           (src.h_qualities_ptr[src_offset].start);
                size_t chunk_offset = src.h_qualities_ptr[src_offset].start;
                size_t qualities_data_size = sizeof(src.h_qualities[0]) * chunk_sz;
                transfer_queue.memcpy(dest.d_qualities, src.h_qualities + chunk_offset, qualities_data_size);
                actual_data_movement += qualities_data_size;
                offset += qualities_data_size;
            }

            // 6. Allocate d_qualities_ptr (an array of SpanField)
            //            {
            //                size_t alignment = alignof(decltype(dest.d_qualities_ptr[0]));
            //                size_t padding = (alignment - (offset % alignment)) % alignment;
            //                offset += padding;
            //                dest.d_qualities_ptr = reinterpret_cast<SpanField *>(base + offset);
            //                size_t qualities_ptr_size = sizeof(dest.d_qualities_ptr[0]) * n_elements;
            //                transfer_queue.memcpy(dest.d_qualities_ptr, src.h_qualities_ptr + src_offset, qualities_ptr_size);
            //
            //                offset += qualities_ptr_size;
            //            }

            // 7. Allocate d_pos (an array of int64_t)
            {
                size_t alignment = alignof(decltype(dest.d_pos[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_pos = reinterpret_cast<int64_t *>(base + offset);
                size_t pos_data_size = sizeof(dest.d_pos[0]) * n_elements;
                transfer_queue.memcpy(dest.d_pos, src.h_pos + src_offset, pos_data_size);
                actual_data_movement += pos_data_size;
                offset += pos_data_size;
            }
            // 8. Allocate d_recalibrated (an array of int8_t)
            {
                size_t alignment = alignof(decltype(dest.d_recalibrated[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_recalibrated = reinterpret_cast<int8_t *>(base + offset);
                size_t recalibrated_data_size = sizeof(dest.d_recalibrated[0]) * n_elements;
                transfer_queue.memcpy(dest.d_recalibrated, src.h_recalibrated + src_offset,
                                      recalibrated_data_size);
                actual_data_movement += recalibrated_data_size;
                offset += recalibrated_data_size;
            }
            // 9. Allocate d_rname (an array of int8_t)
            {
                size_t alignment = alignof(decltype(dest.d_rname[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_rname = reinterpret_cast<int32_t*>(base + offset);
                size_t rname_data_size = sizeof(dest.d_rname[0]) * n_elements;
                transfer_queue.memcpy(dest.d_rname, src.h_rname + src_offset, rname_data_size);
                actual_data_movement += rname_data_size;
                offset += rname_data_size;
            }
            // 10. Allocate d_flag (an array of int8_t)
            {
                size_t alignment = alignof(decltype(dest.d_flag[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_flag = reinterpret_cast<uint16_t *>(base + offset);
                size_t flag_data_size = sizeof(dest.d_flag[0]) * n_elements;
                transfer_queue.memcpy(dest.d_flag, src.h_flag + src_offset, flag_data_size);
                actual_data_movement += flag_data_size;
                offset += flag_data_size;
            }
            // 11. Allocate d_reads_group_covariate (an array of int8_t)
            {
                size_t alignment = alignof(decltype(dest.d_reads_group_covariate[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_reads_group_covariate = reinterpret_cast<int32_t *>(base + offset);
                size_t reads_group_covariate_data_size = sizeof(dest.d_reads_group_covariate[0]) * n_elements;
                transfer_queue.memcpy(dest.d_reads_group_covariate, src.h_read_group_cov + src_offset,
                                      reads_group_covariate_data_size);
                actual_data_movement += reads_group_covariate_data_size;
                offset += reads_group_covariate_data_size;
            }
            // 12. Allocate d_skip_slice
            {
                size_t alignment = alignof(decltype(dest.d_skip_slice[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_skip_slice = reinterpret_cast<uint8_t *>(base + offset);
                chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                            src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                           src.h_sequences_ptr[src_offset].start;
                size_t skip_slice_size = sizeof(dest.d_skip_slice[0]) * chunk_sz;
                transfer_queue.fill(dest.d_skip_slice, static_cast<uint8_t>(0), chunk_sz);
                offset += skip_slice_size;
            }

            // 13. Allocate d_skip_slice_ptr (an array of SpanField)
            //            {
            //                size_t alignment = alignof(decltype(dest.d_skip_slice_ptr[0]));
            //                size_t padding = (alignment - (offset % alignment)) % alignment;
            //                offset += padding;
            //                dest.d_skip_slice_ptr = reinterpret_cast<SpanField *>(base + offset);
            //                size_t skip_slice_ptr_size = sizeof(dest.d_skip_slice_ptr[0]) * n_elements;
            //                transfer_queue.memcpy(dest.d_skip_slice_ptr, src.h_sequences_ptr + src_offset, skip_slice_ptr_size);
            //                offset += skip_slice_ptr_size;
            //            }

            // 14. Allocate d_snps
            {
                size_t alignment = alignof(decltype(dest.d_snps[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_snps = reinterpret_cast<int8_t *>(base + offset);
                chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                            src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                           src.h_sequences_ptr[src_offset].start;
                size_t snps_size = sizeof(dest.d_snps[0]) * chunk_sz;
                transfer_queue.fill(dest.d_snps, static_cast<int8_t>(0), chunk_sz);
                offset += snps_size;
            }
            // 15. Allocate d_snps_ptr (an array of SpanField)
            //            {
            //                size_t alignment = alignof(decltype(dest.d_snps_ptr[0]));
            //                size_t padding = (alignment - (offset % alignment)) % alignment;
            //                offset += padding;
            //                dest.d_snps_ptr = reinterpret_cast<SpanField *>(base + offset);
            //                size_t snps_ptr_size = sizeof(dest.d_snps_ptr[0]) * n_elements;
            //                transfer_queue.memcpy(dest.d_snps_ptr, src.h_sequences_ptr + src_offset, snps_ptr_size);
            //                offset += snps_ptr_size;
            //            }

            // 16. Allocate d_stranded_seq (an array of char)
            {
                size_t alignment = alignof(decltype(dest.d_stranded_seq[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_stranded_seq = reinterpret_cast<char *>(base + offset);
                chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                            src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                           src.h_sequences_ptr[src_offset].start;
                size_t stranded_seq_size = sizeof(dest.d_stranded_seq[0]) * chunk_sz;
                size_t chunk_offset = src.h_sequences_ptr[src_offset].start;
                transfer_queue.memcpy(dest.d_stranded_seq, src.h_sequences + chunk_offset, stranded_seq_size);
                actual_data_movement += stranded_seq_size;
                offset += stranded_seq_size;
            }

            // 17. Allocate d_stranded_seq_ptr (an array of SpanField)
            {
                size_t alignment = alignof(decltype(dest.d_stranded_seq_ptr[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_stranded_seq_ptr = reinterpret_cast<SpanField*>(base + offset);
                size_t stranded_seq_ptr_size = sizeof(dest.d_stranded_seq_ptr[0]) * n_elements;
                transfer_queue.memcpy(dest.d_stranded_seq_ptr, src.h_sequences_ptr + src_offset,
                                      stranded_seq_ptr_size);
                actual_data_movement += stranded_seq_ptr_size;
                offset += stranded_seq_ptr_size;
            }

            // 18. Allocate d_base_context (an array of int32_t)
            {
                size_t alignment = alignof(decltype(dest.d_base_context[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_base_context = reinterpret_cast<int32_t *>(base + offset);
                chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                            src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                           src.h_sequences_ptr[src_offset].start;
                size_t base_context_size = sizeof(dest.d_base_context[0]) * chunk_sz;
                offset += base_context_size;
            }

            // 19. Allocate d_base_context_ptr (an array of SpanField)
            {
                size_t alignment = alignof(decltype(dest.d_base_context_ptr[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_base_context_ptr = reinterpret_cast<SpanField *>(base + offset);
                size_t base_context_ptr_size = sizeof(dest.d_base_context_ptr[0]) * n_elements;
                transfer_queue.memcpy(dest.d_base_context_ptr, src.h_sequences_ptr + src_offset, base_context_ptr_size);
                actual_data_movement += base_context_ptr_size;
                offset += base_context_ptr_size;
            }

            // 20. Allocate d_quality_table_indexes
            {
                size_t alignment = alignof(decltype(dest.d_quality_table_indexes[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_quality_table_indexes = reinterpret_cast<int32_t *>(base + offset);
                chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                            src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                           src.h_sequences_ptr[src_offset].start;
                size_t quality_table_indexes_size = sizeof(dest.d_quality_table_indexes[0]) * chunk_sz;
                //transfer_queue.fill(dest.d_quality_table_indexes, static_cast<int32_t>(INVALID_INDEX), chunk_sz);
                offset += quality_table_indexes_size;
            }

            // 21. Allocate d_cycle_table_indexes
            {
                size_t alignment = alignof(decltype(dest.d_cycle_table_indexes[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_cycle_table_indexes = reinterpret_cast<int32_t *>(base + offset);
                chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                            src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                           src.h_sequences_ptr[src_offset].start;
                size_t cycle_table_indexes_size = sizeof(dest.d_cycle_table_indexes[0]) * chunk_sz;
                offset += cycle_table_indexes_size;
            }

            // 22. Allocate d_context_table_indexes
            {
                size_t alignment = alignof(decltype(dest.d_context_table_indexes[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_context_table_indexes = reinterpret_cast<int32_t *>(base + offset);
                chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                            src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                           src.h_sequences_ptr[src_offset].start;
                size_t context_table_indexes_size = sizeof(dest.d_context_table_indexes[0]) * chunk_sz;
                offset += context_table_indexes_size;
            }

            // 23. Allocate d_quality_table_vals
            {
                size_t alignment = alignof(decltype(dest.d_quality_table_vals[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_quality_table_vals = reinterpret_cast<val_t*>(base + offset);
                chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                            src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                           src.h_sequences_ptr[src_offset].start;
                size_t quality_table_vals_size = sizeof(dest.d_quality_table_vals[0]) * chunk_sz;
                offset += quality_table_vals_size;
            }

            // 24. Allocate d_cycle_table_vals
            {
                size_t alignment = alignof(decltype(dest.d_cycle_table_vals[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_cycle_table_vals = reinterpret_cast<val_t*>(base + offset);
                chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                            src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                           src.h_sequences_ptr[src_offset].start;
                size_t cycle_table_vals_size = sizeof(dest.d_cycle_table_vals[0]) * chunk_sz;
                offset += cycle_table_vals_size;
            }

            // 25. Allocate d_context_table_vals
            {
                size_t alignment = alignof(decltype(dest.d_context_table_vals[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_context_table_vals = reinterpret_cast<val_t*>(base + offset);
                chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                            src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                           src.h_sequences_ptr[src_offset].start;
                size_t context_table_vals_size = sizeof(dest.d_context_table_vals[0]) * chunk_sz;
                offset += context_table_vals_size;
            }

            // 26. Allocate d_segments_quality
            {
                size_t alignment = alignof(decltype(dest.d_segments_quality[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_segments_quality = reinterpret_cast<segment_t*>(base + offset);
                chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                            src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                           src.h_sequences_ptr[src_offset].start;
                size_t segments_size = sizeof(dest.d_segments_quality[0]) * chunk_sz;
                offset += segments_size;
            }
            // 27. Allocate d_segments_cycle
            {
                size_t alignment = alignof(decltype(dest.d_segments_cycle[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_segments_cycle = reinterpret_cast<segment_t*>(base + offset);
                chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                            src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                           src.h_sequences_ptr[src_offset].start;
                size_t segments_size = sizeof(dest.d_segments_cycle[0]) * chunk_sz;
                offset += segments_size;
            }
            // 28. Allocate d_segments_context
            {
                size_t alignment = alignof(decltype(dest.d_segments_context[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_segments_context = reinterpret_cast<segment_t *>(base + offset);
                chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                            src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                           src.h_sequences_ptr[src_offset].start;
                size_t segments_size = sizeof(dest.d_segments_context[0]) * chunk_sz;
                offset += segments_size;
            }
            // 28. Allocate d_temporary_keys, capped at 2 BILLIONS
            {
                size_t alignment = alignof(decltype(dest.d_temporary_keys[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_temporary_keys = reinterpret_cast<segment_t *>(base + offset);
                size_t tmp_chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                                       src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                                      src.h_sequences_ptr[src_offset].start;
                chunk_sz = std::min(tmp_chunk_sz, syclgal_mergesort::MAX_NUM_ELEMENTS);
                size_t segments_size = sizeof(dest.d_temporary_keys[0]) * chunk_sz;
                offset += segments_size;
            }
            // 28. Allocate d_temporary_segments_offsets, capped at 2 BILLIONS
            {
                size_t alignment = alignof(decltype(dest.d_temporary_segments_sequences_offsets[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_temporary_segments_sequences_offsets = reinterpret_cast<int32_t*>(base + offset);
                size_t tmp_chunk_sz = n_elements;
                chunk_sz = std::min(tmp_chunk_sz, syclgal_mergesort::MAX_NUM_ELEMENTS);
                size_t segments_size = sizeof(dest.d_temporary_segments_sequences_offsets[0]) * chunk_sz;
                offset += segments_size;
            }
            // 28. Allocate d_temporary_vals, capped at 2 BILLIONS
            {
                size_t alignment = alignof(decltype(dest.d_temporary_vals[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_temporary_vals = reinterpret_cast<segment_t *>(base + offset);
                size_t tmp_chunk_sz = (src.h_sequences_ptr[src_offset + n_elements - 1].start +
                            src.h_sequences_ptr[src_offset + n_elements - 1].len) -
                           src.h_sequences_ptr[src_offset].start;
                chunk_sz = std::min(tmp_chunk_sz, syclgal_mergesort::MAX_NUM_ELEMENTS);
                size_t segments_size = sizeof(dest.d_temporary_vals[0]) * chunk_sz;
                offset += segments_size;
            }
            dest.table_sz = chunk_sz;

            // other temporary arrays
            {
                size_t alignment = alignof(decltype(dest.d_read_index[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_read_index = reinterpret_cast<int32_t*>(base + offset);
                chunk_sz = (src.h_cigar_start[src_offset + n_elements - 1] +
                            src.h_cigar_len[src_offset + n_elements - 1] + 1) - (src.h_cigar_start[src_offset]);
                size_t cigar_data_size = sizeof(dest.d_read_index[0]) * chunk_sz;
                offset += cigar_data_size;
                dest._n_expanded_cigar_ops = chunk_sz;
            }
            {

                size_t alignment = alignof(decltype(dest.d_i_offset[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_i_offset = reinterpret_cast<int32_t*>(base + offset);
                size_t i_off_size = sizeof(dest.d_i_offset[0]) * chunk_sz;
                offset += i_off_size;
            }
            {
                size_t alignment = alignof(decltype(dest.d_j_offset[0]));
                size_t padding = (alignment - (offset % alignment)) % alignment;
                offset += padding;
                dest.d_j_offset = reinterpret_cast<int32_t*>(base + offset);
                size_t j_off_size = sizeof(dest.d_j_offset[0]) * chunk_sz;
                offset += j_off_size;
            }

        }
        if (offset >= dest.max_memory) {
            exit(-1);
        }
        dest.max_num_elements = bqsr_covariate_vector_kernel::update_batch_num_elements(src, dest, src_offset + n_elements);
        transfer_queue.wait_and_throw();
        Logger::getLogger().log(LOG_DEBUG, "    \tActual data movement: ", actual_data_movement, " bytes out of ", offset, " bytes allocated");
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
    // std::vector<int32_t> tmp_vals(src.table_sz, 0);
    // std::vector<SpanField> tmp_span(src.max_num_elements);
    // transfer_queue.memcpy(tmp_vals.data(), src.d_quality_table_indexes, sizeof(src.d_quality_table_indexes[0])*src.table_sz);
    // transfer_queue.memcpy(tmp_span.data(), src.d_sequences_ptr, sizeof(src.d_sequences_ptr[0])*src.max_num_elements);
    // transfer_queue.wait_and_throw();

    /*
    for (size_t i = 0; i < src.max_num_elements; i++) {
        for (size_t j = 0; j < tmp_span[i].len; j++) {
            if (tmp_vals[tmp_span[i].start + j] != 1) {
                std::cout << "ERROR" << std::endl;
            }
        }
    }*/
    // for (size_t j = 0; j < src.max_num_elements; j++) {
    //     for (size_t i = tmp_span[j].start; i < tmp_span[j].start + tmp_span[j].len; i++) {
    //         if (tmp_vals[i] != 1) {
    //         //if (tmp_vals[i] != j) {
    //             std::cout << "ERROR" << std::endl;
    //         }
    //     }
    // }
//    Timer t;
//    //Logger::getLogger().log(LOG_INFO, "Updating tables");
//    size_t num_segments = src.num_segments_quality;
//    auto d_quality_table_indexes = src.d_quality_table_indexes;
//    auto d_quality_table_vals = src.d_quality_table_vals;
//    auto d_segments_quality = src.d_segments_quality;
//
//    auto d_quality_table = dest.d_quality_table;
//
//    if (num_segments > 1) {
//
//        transfer_queue.parallel_for(sycl::range<1>(num_segments - 1), [=](sycl::id<1> segment_idx) {
//
//            int32_t covariate_idx = d_quality_table_indexes[segment_idx];
//            int32_t err_val = d_quality_table_vals[segment_idx];
//            int32_t segment_sz = d_segments_quality[segment_idx + 1] - d_segments_quality[segment_idx];
//
//            if (covariate_idx == INVALID_INDEX) {
//                return;
//            }
//            if (d_quality_table[covariate_idx].Observations == -1) {
//                d_quality_table[covariate_idx].Observations = 0;
//            }
//            d_quality_table[covariate_idx].Mismatches += err_val;
//            d_quality_table[covariate_idx].Observations += segment_sz;
//        });
//    }
//
//    num_segments = src.num_segments_cycle;
//    auto d_cycle_table_indexes = src.d_cycle_table_indexes;
//    auto d_cycle_table_vals = src.d_cycle_table_vals;
//    auto d_segments_cycle = src.d_segments_cycle;
//
//    auto d_cycle_table = dest.d_cycle_table;
//    if (num_segments > 1) {
//        transfer_queue.parallel_for(sycl::range<1>(num_segments - 1), [=](sycl::id<1> segment_idx) {
//
//            int32_t covariate_idx = d_cycle_table_indexes[segment_idx];
//            int32_t err_val = d_cycle_table_vals[segment_idx];
//            int32_t segment_sz = d_segments_cycle[segment_idx + 1] - d_segments_cycle[segment_idx];
//
//            if (covariate_idx == INVALID_INDEX) {
//                return;
//            }
//            if (d_cycle_table[covariate_idx].Observations == -1) {
//                d_cycle_table[covariate_idx].Observations = 0;
//            }
//            d_cycle_table[covariate_idx].Mismatches += err_val;
//            d_cycle_table[covariate_idx].Observations += segment_sz;
//        });
//    }
//
//    num_segments = src.num_segments_context;
//    auto d_context_table_indexes = src.d_context_table_indexes;
//    auto d_context_table_vals = src.d_context_table_vals;
//    auto d_segments_context = src.d_segments_context;
//
//    auto d_context_table = dest.d_context_table;
//    if (num_segments > 1) {
//
//        transfer_queue.parallel_for(sycl::range<1>(num_segments - 1), [=](sycl::id<1> segment_idx) {
//
//            int32_t covariate_idx = d_context_table_indexes[segment_idx];
//            int32_t err_val = d_context_table_vals[segment_idx];
//            int32_t segment_sz = d_segments_context[segment_idx + 1] - d_segments_context[segment_idx];
//
//            if (covariate_idx == INVALID_INDEX) {
//                return;
//            }
//            if (d_context_table[covariate_idx].Observations == -1) {
//                d_context_table[covariate_idx].Observations = 0;
//            }
//            d_context_table[covariate_idx].Mismatches += err_val;
//            d_context_table[covariate_idx].Observations += segment_sz;
//        });
//    }
//    transfer_queue.wait_and_throw();
    //Logger::getLogger().log(LOG_INFO, "    [UpdateTables kernels] Time: " + std::to_string(t.elapsed()) + " sec");
    //dest.table_offset += src.table_sz;
}


size_t bqsr_covariate_vector_kernel::update_batch_num_elements(const SequencesInput& in, const ComputeSkipSliceKernelData& kernel_data, size_t start_idx) {
    size_t n_elements = 1;
    size_t idx = start_idx + n_elements;
    //size_t base = 0;

    if (idx > in._n) {
        return 0;
    }
    size_t offset = compute_chunk_sz(in, kernel_data, start_idx, n_elements);

    while (offset < kernel_data.max_memory && idx < in._n) {
        n_elements = min(n_elements * 2, in._n-start_idx);
        offset = compute_chunk_sz(in, kernel_data, start_idx, n_elements);
        idx = start_idx + n_elements;
    }

    // At this point, n_elements is an upper bound that may exceed memory.
    // Set up binary search between [low, high], where:
    //   low = n_elements/2 (last count that fit) and high = n_elements.
    size_t low = n_elements / 2;
    size_t high = n_elements;
    size_t best = low;  // best so far that fits.
    size_t usage = 0;

    while (low <= high) {
        size_t mid = low + (high - low) / 2;
        if (start_idx + mid > in._n) {  // avoid going out of bounds
            if (mid==0) continue;
            high = mid - 1;
            continue;
        }
        usage = compute_chunk_sz(in, kernel_data, start_idx, mid);
        if (usage < kernel_data.max_memory) {
            best = mid;    // mid fits within the memory
            low = mid + 1; // try for more elements
        } else {
            if (mid==0) continue;
            high = mid - 1;
        }
    }
    usage = compute_chunk_sz(in, kernel_data, start_idx, best);
    n_elements = best;
    return n_elements;
}

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
    size_t sequences_data_size = sizeof(in.h_sequences[0]) * chunk_sz;
    offset += sequences_data_size;

    // 2. d_sequences_ptr (an array of SpanField)
    //    Align the pointer properly.
    {
        size_t alignment = alignof(decltype(kernel_data.d_sequences_ptr[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        size_t sequences_ptr_size = sizeof(kernel_data.d_sequences_ptr[0]) * n_elements;
        offset += sequences_ptr_size;
    }
    // 3. d_cigar (an array of CigarOperation)
    {
        size_t alignment = alignof(decltype(kernel_data.d_cigar[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        chunk_sz = (in.h_cigar_start[start_idx + n_elements - 1] +
                    in.h_cigar_len[start_idx + n_elements - 1]) - (in.h_cigar_start[start_idx]);
        size_t cigar_data_size = sizeof(in.h_cigar[0]) * chunk_sz;
        offset += cigar_data_size;
    }

    // 4.1 d_cigar_ptr (an array of SpanField)
    {
        size_t alignment = alignof(decltype(kernel_data.d_cigar_start[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        size_t cigar_ptr_size = sizeof(kernel_data.d_cigar_start[0]) * n_elements;
        offset += cigar_ptr_size;
    }
    // 4.2 d_cigar_ptr (an array of SpanField)
    {
        size_t alignment = alignof(decltype(kernel_data.d_cigar_len[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        size_t cigar_ptr_size = sizeof(kernel_data.d_cigar_len[0]) * n_elements;
        offset += cigar_ptr_size;
    }

    // 5. d_qualities (an array of char)
    {
        size_t alignment = alignof(decltype(kernel_data.d_qualities[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        chunk_sz = (in.h_qualities_ptr[start_idx + n_elements - 1].start +
                    in.h_qualities_ptr[start_idx + n_elements - 1].len) -
                   (in.h_qualities_ptr[start_idx].start);
        size_t qualities_data_size = sizeof(in.h_qualities[0]) * chunk_sz;
        offset += qualities_data_size;
    }

    // 6. d_qualities_ptr (an array of SpanField)
//    {
//        size_t alignment = alignof(decltype(kernel_data.d_qualities_ptr[0]));
//        size_t padding = (alignment - (offset % alignment)) % alignment;
//        offset += padding;
//        size_t qualities_ptr_size = sizeof(kernel_data.d_qualities_ptr[0]) * n_elements;
//        offset += qualities_ptr_size;
//    }

    // 7. d_pos (an array of int64_t)
    {
        size_t alignment = alignof(decltype(kernel_data.d_pos[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        size_t pos_data_size = sizeof(kernel_data.d_pos[0]) * n_elements;
        offset += pos_data_size;
    }
    // 8. d_recalibrated (an array of int8_t)
    {
        size_t alignment = alignof(decltype(kernel_data.d_recalibrated[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        size_t recalibrated_data_size = sizeof(kernel_data.d_recalibrated[0]) * n_elements;
        offset += recalibrated_data_size;
    }
    // 9. d_rname (an array of int8_t)
    {
        size_t alignment = alignof(decltype(kernel_data.d_rname[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        size_t rname_data_size = sizeof(kernel_data.d_rname[0]) * n_elements;
        offset += rname_data_size;
    }
    // 10. d_flag (an array of int8_t)
    {
        size_t alignment = alignof(decltype(kernel_data.d_flag[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        size_t flag_data_size = sizeof(kernel_data.d_flag[0]) * n_elements;
        offset += flag_data_size;
    }
    // 11. d_reads_group_covariate (an array of int8_t)
    {
        size_t alignment = alignof(decltype(kernel_data.d_reads_group_covariate[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        size_t reads_group_covariate_data_size = sizeof(kernel_data.d_reads_group_covariate[0]) * n_elements;
        offset += reads_group_covariate_data_size;
    }
    // 12. d_skip_slice
    {
        size_t alignment = alignof(decltype(kernel_data.d_skip_slice[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        chunk_sz = (in.h_sequences_ptr[start_idx + n_elements - 1].start +
                    in.h_sequences_ptr[start_idx + n_elements - 1].len) -
                   in.h_sequences_ptr[start_idx].start;
        size_t skip_slice_size = sizeof(kernel_data.d_skip_slice[0]) * chunk_sz;
        offset += skip_slice_size;
    }

    // 13. d_skip_slice_ptr (an array of SpanField)
//    {
//        size_t alignment = alignof(decltype(kernel_data.d_skip_slice_ptr[0]));
//        size_t padding = (alignment - (offset % alignment)) % alignment;
//        offset += padding;
//        size_t skip_slice_ptr_size = sizeof(kernel_data.d_skip_slice_ptr[0]) * n_elements;
//        offset += skip_slice_ptr_size;
//    }

    // 14. d_snps
    {
        size_t alignment = alignof(decltype(kernel_data.d_snps[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        chunk_sz = (in.h_sequences_ptr[start_idx + n_elements - 1].start +
                    in.h_sequences_ptr[start_idx + n_elements - 1].len) -
                   in.h_sequences_ptr[start_idx].start;
        size_t snps_size = sizeof(kernel_data.d_snps[0]) * chunk_sz;
        offset += snps_size;
    }

    // 15. d_snps_ptr (an array of SpanField)
//    {
//        size_t alignment = alignof(decltype(kernel_data.d_snps_ptr[0]));
//        size_t padding = (alignment - (offset % alignment)) % alignment;
//        offset += padding;
//        size_t snps_ptr_size = sizeof(kernel_data.d_snps_ptr[0]) * n_elements;
//        offset += snps_ptr_size;
//    }

    // 16. d_stranded_seq (an array of char)
    {
        size_t alignment = alignof(decltype(kernel_data.d_stranded_seq[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        chunk_sz = (in.h_sequences_ptr[start_idx + n_elements - 1].start +
                    in.h_sequences_ptr[start_idx + n_elements - 1].len) -
                   in.h_sequences_ptr[start_idx].start;
        size_t stranded_seq_size = sizeof(kernel_data.d_stranded_seq[0]) * chunk_sz;
        offset += stranded_seq_size;
    }

    // 17. d_stranded_seq_ptr (an array of SpanField)
    {
        size_t alignment = alignof(decltype(kernel_data.d_stranded_seq_ptr[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        size_t stranded_seq_ptr_size = sizeof(kernel_data.d_stranded_seq_ptr[0]) * n_elements;
        offset += stranded_seq_ptr_size;
    }

    // 18. d_base_context (an array of int32_t)
    {
        size_t alignment = alignof(decltype(kernel_data.d_base_context[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        chunk_sz = (in.h_sequences_ptr[start_idx + n_elements - 1].start +
                    in.h_sequences_ptr[start_idx + n_elements - 1].len) -
                   in.h_sequences_ptr[start_idx].start;
        size_t base_context_size = sizeof(kernel_data.d_base_context[0]) * chunk_sz;
        offset += base_context_size;
    }

    // 19. d_base_context_ptr (an array of SpanField)
    {
        size_t alignment = alignof(decltype(kernel_data.d_base_context_ptr[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        size_t base_context_ptr_size = sizeof(kernel_data.d_base_context_ptr[0]) * n_elements;
        offset += base_context_ptr_size;
    }

    // 20. d_quality_table_indexes
    {
        size_t alignment = alignof(decltype(kernel_data.d_quality_table_indexes[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        chunk_sz = (in.h_sequences_ptr[start_idx + n_elements - 1].start +
                    in.h_sequences_ptr[start_idx + n_elements - 1].len) -
                   in.h_sequences_ptr[start_idx].start;
        size_t quality_table_indexes_size = sizeof(kernel_data.d_quality_table_indexes[0]) * chunk_sz;
        offset += quality_table_indexes_size;
        //extra_table_indexes_sz = quality_table_indexes_size;
    }

    // 21. d_cycle_table_indexes
    {
        size_t alignment = alignof(decltype(kernel_data.d_cycle_table_indexes[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        chunk_sz = (in.h_sequences_ptr[start_idx + n_elements - 1].start +
                    in.h_sequences_ptr[start_idx + n_elements - 1].len) -
                   in.h_sequences_ptr[start_idx].start;
        size_t cycle_table_indexes_size = sizeof(kernel_data.d_cycle_table_indexes[0]) * chunk_sz;
        offset += cycle_table_indexes_size;
    }

    // 22. d_context_table_indexes
    {
        size_t alignment = alignof(decltype(kernel_data.d_context_table_indexes[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        chunk_sz = (in.h_sequences_ptr[start_idx + n_elements - 1].start +
                    in.h_sequences_ptr[start_idx + n_elements - 1].len) -
                   in.h_sequences_ptr[start_idx].start;
        size_t context_table_indexes_size = sizeof(kernel_data.d_context_table_indexes[0]) * chunk_sz;
        offset += context_table_indexes_size;
    }

    // 23. d_quality_table_vals
    {
        size_t alignment = alignof(decltype(kernel_data.d_quality_table_vals[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        chunk_sz = (in.h_sequences_ptr[start_idx + n_elements - 1].start +
                    in.h_sequences_ptr[start_idx + n_elements - 1].len) -
                   in.h_sequences_ptr[start_idx].start;
        size_t quality_table_vals_size = sizeof(kernel_data.d_quality_table_vals[0]) * chunk_sz;
        offset += quality_table_vals_size;
        //extra_table_values_sz = quality_table_vals_size;
    }

    // 24. d_cycle_table_vals
    {
        size_t alignment = alignof(decltype(kernel_data.d_cycle_table_vals[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        chunk_sz = (in.h_sequences_ptr[start_idx + n_elements - 1].start +
                    in.h_sequences_ptr[start_idx + n_elements - 1].len) -
                   in.h_sequences_ptr[start_idx].start;
        size_t cycle_table_vals_size = sizeof(kernel_data.d_cycle_table_vals[0]) * chunk_sz;
        offset += cycle_table_vals_size;
    }

    // 25. d_context_table_vals
    //size_t table_indexes_sz = 0;
    {
        size_t alignment = alignof(decltype(kernel_data.d_context_table_vals[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        chunk_sz = (in.h_sequences_ptr[start_idx + n_elements - 1].start +
                    in.h_sequences_ptr[start_idx + n_elements - 1].len) -
                   in.h_sequences_ptr[start_idx].start;
        size_t context_table_vals_size = sizeof(kernel_data.d_context_table_vals[0]) * chunk_sz;
        offset += context_table_vals_size;
    }
    // add the memory for temporary segment sequences offsets
    {
        size_t alignment = alignof(decltype(kernel_data.d_temporary_segments_sequences_offsets[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        chunk_sz = n_elements;
        size_t segments_size = sizeof(kernel_data.d_temporary_segments_sequences_offsets[0]) * chunk_sz;
        offset += segments_size;
    }
    // add the memory for segments 1
    {
        size_t alignment = alignof(decltype(kernel_data.d_segments_quality[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        chunk_sz = (in.h_sequences_ptr[start_idx + n_elements - 1].start +
                    in.h_sequences_ptr[start_idx + n_elements - 1].len) -
                   in.h_sequences_ptr[start_idx].start;
        size_t segments_size = sizeof(kernel_data.d_segments_quality[0]) * chunk_sz;
        offset += segments_size;
    }
    //  // add the memory for segments 2
    {
        size_t alignment = alignof(decltype(kernel_data.d_segments_cycle[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        chunk_sz = (in.h_sequences_ptr[start_idx + n_elements - 1].start +
                    in.h_sequences_ptr[start_idx + n_elements - 1].len) -
                   in.h_sequences_ptr[start_idx].start;
        size_t segments_size = sizeof(kernel_data.d_segments_cycle[0]) * chunk_sz;
        offset += segments_size;
    }
    // add the memory for segments 3
    {
        size_t alignment = alignof(decltype(kernel_data.d_segments_context[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        chunk_sz = (in.h_sequences_ptr[start_idx + n_elements - 1].start +
                    in.h_sequences_ptr[start_idx + n_elements - 1].len) -
                   in.h_sequences_ptr[start_idx].start;
        size_t segments_size = sizeof(kernel_data.d_segments_context[0]) * chunk_sz;
        offset += segments_size;
    }

    // add the memory for the temporary buffer keys in sorting
    {
        size_t alignment = alignof(decltype(kernel_data.d_temporary_keys[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        chunk_sz = (in.h_sequences_ptr[start_idx + n_elements - 1].start +
                    in.h_sequences_ptr[start_idx + n_elements - 1].len) -
                   in.h_sequences_ptr[start_idx].start;
        size_t segments_size = sizeof(kernel_data.d_temporary_keys[0]) * chunk_sz;
        offset += segments_size;
    }
    // add the memory for the temporary buffer keys in sorting
    {
        size_t alignment = alignof(decltype(kernel_data.d_temporary_vals[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        chunk_sz = (in.h_sequences_ptr[start_idx + n_elements - 1].start +
                    in.h_sequences_ptr[start_idx + n_elements - 1].len) -
                   in.h_sequences_ptr[start_idx].start;
        size_t segments_size = sizeof(kernel_data.d_temporary_vals[0]) * chunk_sz;
        offset += segments_size;
    }

    // d_cigar_segments
    {
        size_t alignment = alignof(decltype(kernel_data.d_read_index[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;

        size_t cigar_data_size = sizeof(kernel_data.d_read_index[0]) * n_elements;
        offset += cigar_data_size;
    }
    {
        chunk_sz = (in.h_cigar_start[start_idx + n_elements - 1] +
                    in.h_cigar_len[start_idx + n_elements - 1]) - (in.h_cigar_start[start_idx]);
        size_t alignment = alignof(decltype(kernel_data.d_i_offset[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        size_t i_offset_size = sizeof(kernel_data.d_i_offset[0]) * chunk_sz;
        offset += i_offset_size;
    }
    {
        size_t alignment = alignof(decltype(kernel_data.d_j_offset[0]));
        size_t padding = (alignment - (offset % alignment)) % alignment;
        offset += padding;
        size_t j_offset_size = sizeof(kernel_data.d_j_offset[0]) * chunk_sz;
        offset += j_offset_size;
    }

    return offset;
}