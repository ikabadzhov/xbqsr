#ifndef GENOMICPIPELINE_COMPUTECOVARIATEVECTORSDATA_H
#define GENOMICPIPELINE_COMPUTECOVARIATEVECTORSDATA_H

#include <sycl/sycl.hpp>
#include <cstdint>
#include "../../sycl_parallel_library/memory/memory.h"
#include "../../device/device_types.h"
#include "../../sam_data_structures/sam_fields/Cigar.h"
#include "../../intervals/Interval.h"
#include "../TablesBQSR.h"

typedef int32_t val_t;
typedef int32_t segment_t;

const size_t PAD = 8;

class TablesOutput : public OutputData {
public:
//    int32_t* h_quality_table_indexes;
//    int32_t* h_cycle_table_indexes;
//    int32_t* h_context_table_indexes;
//    val_t* h_quality_table_vals;
//    val_t* h_cycle_table_vals;
//    val_t* h_context_table_vals;

    explicit TablesOutput() = default;
};


class SequencesInput : public InputData {
public:
    Interval* h_known_sites;
    SpanFieldLong* h_known_sites_ptr;
    size_t h_known_sites_ptr_sz;

    uint8_t* h_reference;
    SpanFieldLong* h_reference_ptr;
    size_t h_reference_ptr_sz;

    char* h_sequences;
    SpanField* h_sequences_ptr;

    CigarOperation* h_cigar;
    SpanField* h_cigar_ptr;

    size_t* h_cigar_start;
    uint16_t* h_cigar_len;

    char* h_qualities;
    SpanField* h_qualities_ptr;

    int32_t* h_rnext = nullptr;
    int32_t* h_pnext = nullptr;
    int32_t* h_tlen = nullptr;

    int32_t* h_pos;
    uint16_t* h_flag;
    int32_t* h_read_group_cov;
    int32_t* h_rname;
    uint8_t* h_recalibrated;
    int32_t h_sz_qual_score_cov;
    int32_t h_sz_cycle_cov;
    int32_t h_sz_context_cov;

    bqsrEntry* s_quality_table;
    size_t s_quality_table_sz;
    bqsrEntry* s_cycle_table;
    size_t s_cycle_table_sz;
    bqsrEntry* s_context_table;
    size_t s_context_table_sz;


    explicit SequencesInput(mem_t<Interval>& known_sites, mem_t<SpanFieldLong>& known_sites_ptr, mem_t<uint8_t>& reference, mem_t<SpanFieldLong>& reference_ptr, bqsrEntry* quality_table, bqsrEntry* cycle_table, bqsrEntry* context_table, size_t quality_table_sz, size_t cycle_table_sz, size_t context_table_sz, int32_t sz_qual_cov, int32_t sz_cycle_cov, int32_t sz_context_cov, char* sequence, SpanField* sequence_ptr, char* qualities, SpanField* qualities_ptr, CigarOperation* cigar, SpanField* cigar_ptr, int32_t* pos, int32_t* rnext, int32_t* pnext, int32_t* tlen, int32_t* rname, uint8_t* recalibrated, uint16_t* flags, int32_t* read_group_cov, size_t n_elements);
 };


class ComputeSkipSliceConfiguration {
public:
    size_t max_memory = 0;
    size_t interval_memory = 0;
    size_t interval_ptr_memory = 0;
    size_t reference_memory = 0;
    size_t reference_ptr_memory = 0;
    size_t max_n_elements_first_batch = 0;
};



class ComputeSkipSliceKernelData {
public:
    sycl::queue _q;
    memory_space_t _space = memory_space_device;
    size_t max_num_elements {};
    bool is_first_copy = true;
    Interval* d_known_sites = nullptr;
    size_t d_known_sites_sz = 0;
    SpanFieldLong* d_known_sites_ptr = nullptr;

    uint8_t* d_reference = nullptr;
    SpanFieldLong* d_reference_ptr = nullptr;
    size_t d_reference_sz = 0;

    size_t d_sequence_span_offset = 0;
    size_t d_cigar_span_offset = 0;

    char* memory_buffer = nullptr;
    size_t max_memory = 0;

    char* d_sequences = nullptr;
    SpanField* d_sequences_ptr = nullptr;

    CigarOperation* d_cigar = nullptr;
    CigarOperation* d_clipped_cigar = nullptr;
    SpanField* d_cigar_ptr = nullptr;
    size_t* d_clipped_cigar_start = nullptr;
    uint16_t* d_clipped_cigar_len = nullptr;

    char* d_qualities = nullptr;

    int32_t* d_pos = nullptr;
    int32_t* d_rname = nullptr;
    int32_t* d_rnext = nullptr;
    int32_t* d_pnext = nullptr;
    int32_t* d_tlen = nullptr;

    int8_t* d_recalibrated = nullptr;
    uint16_t* d_flag = nullptr;
    int32_t* d_reads_group_covariate = nullptr;

    uint8_t* d_skip_slice = nullptr;
    int8_t* d_snps = nullptr;
    char* d_stranded_seq = nullptr;
    SpanField* d_stranded_seq_ptr = nullptr;
    int32_t* d_base_context = nullptr;
    SpanField* d_base_context_ptr = nullptr;

    int32_t d_sz_qual_score_cov = 0;
    int32_t d_sz_cycle_cov = 0;
    int32_t d_sz_context_cov = 0;

    int32_t* d_quality_table_indexes = nullptr;
    int32_t* d_cycle_table_indexes = nullptr;
    int32_t* d_context_table_indexes = nullptr;

    val_t* d_quality_table_vals = nullptr;
    val_t* d_cycle_table_vals = nullptr;
    val_t* d_context_table_vals = nullptr;

    // segment_t* d_segments_quality = nullptr;
    // segment_t* d_segments_cycle = nullptr;
    // segment_t* d_segments_context = nullptr;

//    key_t* d_temporary_keys = nullptr;
//    val_t* d_temporary_vals = nullptr;
    int32_t* d_temporary_segments_sequences_offsets = nullptr;

    size_t table_sz = 0;
    int32_t prev_batch_max_index = 0;

    bqsrEntry* d_quality_table = nullptr;
    bqsrEntry* d_cycle_table = nullptr;
    bqsrEntry* d_context_table = nullptr;

    size_t d_quality_table_sz = 0;
    size_t d_cycle_table_sz = 0;
    size_t d_context_table_sz = 0;

    size_t _n_expanded_cigar_ops = 0;

public:
    ComputeSkipSliceKernelData() = default;

    ComputeSkipSliceKernelData(sycl::queue &q, ComputeSkipSliceConfiguration &config,
                                                           memory_space_t space = memory_space_device) : _q(q), _space(space) {
        //assert(config.max_memory > (config.interval_memory + config.reference_memory));

        size_t max_memory = config.max_memory;// - config.interval_memory - config.reference_memory;

        if (space != memory_space_t::none) {

            memory_buffer = sycl::malloc_device<char>(max_memory, q);

            if (!memory_buffer) {
                std::cerr << "Error: Memory allocation failed for one or more pointers." << std::endl;
                exit(-1);
            }
        }
        this->max_memory = max_memory;
        this->max_num_elements = config.max_n_elements_first_batch;

    }

    void swap(ComputeSkipSliceKernelData& rhs);

    // Destructor to manage dynamic memory if ownership is taken
    ~ ComputeSkipSliceKernelData();

    // Disable copy semantics to prevent accidental copying
    ComputeSkipSliceKernelData(const ComputeSkipSliceKernelData&) = delete;
    ComputeSkipSliceKernelData& operator=(const ComputeSkipSliceKernelData&) = delete;

    // Enable move semantics
    ComputeSkipSliceKernelData(ComputeSkipSliceKernelData&& rhs) noexcept;

    ComputeSkipSliceKernelData& operator=(ComputeSkipSliceKernelData&& rhs) noexcept;

    /**
     * @brief Generic memcpy function to copy specified fields from src to dest.
     *
     * @tparam Members Variadic template parameters representing member pointers.
     * @param dest Destination DeviceData object.
     * @param src Source DeviceData object.
     * @param src_from Starting index in the source.
     * @param dest_to Starting index in the destination.
     * @param n Number of elements to copy.
     * @param members Member pointers specifying which fields to copy.
     */
    static void memcpy_host_to_dev(sycl::queue& transfer_queue, ComputeSkipSliceKernelData& dest, size_t dest_offset, const SequencesInput& src,
                                   size_t src_offset, size_t n_elements);

    static void memcpy_dev_to_host(sycl::queue& transfer_queue, TablesOutput& dest, size_t dest_offset, ComputeSkipSliceKernelData& src,
                                   size_t src_offset, size_t n_elements);

    template <typename InputT>
    static size_t compute_max_elements(const InputT& in,
                                      const ComputeSkipSliceKernelData& kd,
                                      size_t start_idx);
};

namespace bqsr_covariate_vector_kernel {
    size_t update_batch_num_elements(const SequencesInput &in, const ComputeSkipSliceKernelData &kd,
                                     size_t start_idx);
}

template<typename InputT>
size_t ComputeSkipSliceKernelData::compute_max_elements(const InputT &in, const ComputeSkipSliceKernelData &kd,
    size_t start_idx) {
    return bqsr_covariate_vector_kernel::update_batch_num_elements(in, kd, start_idx);
}



#endif //GENOMICPIPELINE_COMPUTECOVARIATEVECTORSDATA_H
