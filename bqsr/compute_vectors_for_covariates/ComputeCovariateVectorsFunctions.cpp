#include "ComputeCovariateVectorsFunctions.h"

#include <boost/iostreams/categories.hpp>
#include "../UtilsBQSR.h"
#include "../../sycl_parallel_library/sorting/mergesort.h"
#include "../../timer/Timer.h"
#include "../../sycl_parallel_library/segmented_reduce/segreduce.h"
#include "../../sycl_parallel_library/intervalsmove/intervalmove.h"
#include "ComputeCovariateVectorsData.h"
#include "../../logger/Logger.h"
#include "../../timer/Profiler.h"

// Quality
template <typename IndexT, typename ValueT>
inline void aggregate_quality(
    sycl::queue& q,
    const IndexT* d_idx, const ValueT* d_val, size_t N,
    bqsrEntry* d_table, size_t table_bins
) {
    //aggregate_dense_atomic(q, d_idx, d_val, N, d_table, table_bins);
    aggregate_quality_localWG(q, d_idx, d_val, N, d_table, table_bins);
}

// Cycle
template <typename IndexT, typename ValueT>
inline void aggregate_cycle(
    sycl::queue& q,
    const IndexT* d_idx, const ValueT* d_val, size_t N,
    bqsrEntry* d_table, size_t table_bins
) {
    aggregate_dense_atomic(q, d_idx, d_val, N, d_table, table_bins);
    //aggregate_quality_globalWG(q, d_idx, d_val, N, d_table, table_bins);
}

// Context (sparse index space; still fine with atomics)
template <typename IndexT, typename ValueT>
inline void aggregate_context(
    sycl::queue& q,
    const IndexT* d_idx, const ValueT* d_val, size_t N,
    bqsrEntry* d_table, size_t table_bins
) {
    aggregate_dense_atomic(q, d_idx, d_val, N, d_table, table_bins);
}

template <typename IndexT, typename ValueT>
void aggregate_dense_atomic(
    sycl::queue& q,
    const IndexT* d_indexes,     // size N
    const ValueT* d_vals,        // size N (err flags or counts)
    size_t N,
    bqsrEntry* d_table,         // size = table_bins
    size_t table_bins            // for optional bounds check
) {
    const size_t wg = 256;                         // good default
    const size_t g  = div_up(std::max<size_t>(N, wg), wg) * wg;

    q.submit([&](sycl::handler& cgh){
        cgh.parallel_for(
            sycl::nd_range<1>(sycl::range<1>(g), sycl::range<1>(wg)),
            [=](sycl::nd_item<1> it) {
                const size_t stride = it.get_global_range(0);
                size_t i = it.get_global_id(0);
                for (; i < N; i += stride) {
                    IndexT idx = d_indexes[i];
                    if (idx < 0) continue;          // skip invalid
                    if (static_cast<size_t>(idx) >= table_bins) continue; // guard

                    // Relaxed atomics are fine here.
                    auto obs_ref = sycl::atomic_ref<int64_t,
                        sycl::memory_order::relaxed,
                        sycl::memory_scope::device,
                        sycl::access::address_space::global_space>(d_table[idx].Observations);

                    auto mis_ref = sycl::atomic_ref<int64_t,
                        sycl::memory_order::relaxed,
                        sycl::memory_scope::device,
                        sycl::access::address_space::global_space>(d_table[idx].Mismatches);

                    mis_ref.fetch_add(static_cast<int64_t>(d_vals[i]));
                    obs_ref.fetch_add(1);

                }
            }
        );
    });
}

// N.B. d_idx in [0,93] or -1 for "invalid"; d_val is the error (0/1 or any int to add).
// d_table has exactly 94 entries and is already zeroed (or you intend to accumulate into it).
inline void aggregate_quality_localWG(
    sycl::queue& q,
    const int32_t* d_idx,           // [N]
    const int32_t* d_val,           // [N]
    size_t N,
    bqsrEntry* d_table,
    size_t table_size// [94]
) {
    const size_t wg = 256; // good default; tune if needed
    // make global size a multiple of wg and >= N (via grid-stride loop)
    const size_t g  = ((std::max(N, wg) + wg - 1) / wg) * wg;

    q.submit([&](sycl::handler& cgh) {
        // local histograms (Observations and Mismatches)
        sycl::local_accessor<int32_t, 1> local_obs (sycl::range<1>(table_size), cgh);
        sycl::local_accessor<int32_t, 1> local_mis (sycl::range<1>(table_size), cgh);

        cgh.parallel_for(
            sycl::nd_range<1>(sycl::range<1>(g), sycl::range<1>(wg)),
            [=](sycl::nd_item<1> it) {
                const int  lid    = it.get_local_id(0);
                const int  lsize  = it.get_local_range(0);
                const size_t gsize = it.get_global_range(0);
                size_t i = it.get_global_id(0);

                // --- Zero local histograms in parallel
                for (int b = lid; b < table_size; b += lsize) {
                    local_obs[b] = 0;
                    local_mis[b] = 0;
                }
                it.barrier(sycl::access::fence_space::local_space);

                // --- Grid-stride loop over all elements
                for (; i < N; i += gsize) {
                    const int idx = d_idx[i];
                    if ((unsigned)idx < (unsigned)table_size) {
                        // local-space relaxed atomics to reduce bank conflicts
                        auto a_obs = sycl::atomic_ref<
                            int32_t,
                            sycl::memory_order::relaxed,
                            sycl::memory_scope::work_group,
                            sycl::access::address_space::local_space>(local_obs[idx]);
                        auto a_mis = sycl::atomic_ref<
                            int32_t,
                            sycl::memory_order::relaxed,
                            sycl::memory_scope::work_group,
                            sycl::access::address_space::local_space>(local_mis[idx]);
                        a_mis.fetch_add(d_val[i]);
                        a_obs.fetch_add(1);
                    }
                }

                it.barrier(sycl::access::fence_space::local_space);

                // --- Flush: first BINS threads (or strided) push to global with one atomic each
                for (int b = lid; b < table_size; b += lsize) {
                    const int obs = local_obs[b];
                    const int mis = local_mis[b];
                    if (obs == 0 && mis == 0) continue;

                    // Atomics on global table entries
                    auto g_obs = sycl::atomic_ref<
                        int64_t,
                        sycl::memory_order::relaxed,
                        sycl::memory_scope::device,
                        sycl::access::address_space::global_space>(d_table[b].Observations);
                    auto g_mis = sycl::atomic_ref<
                        int64_t,
                        sycl::memory_order::relaxed,
                        sycl::memory_scope::device,
                        sycl::access::address_space::global_space>(d_table[b].Mismatches);
                        g_obs.fetch_add(obs);
                        g_mis.fetch_add(mis);
                }
            }
        );
    });
}

inline void aggregate_quality_globalWG(
    sycl::queue& q,
    const int32_t* d_idx,           // [N]
    const int32_t* d_val,           // [N]
    size_t N,
    bqsrEntry* d_table,
    size_t table_size // [94]
) {
    const size_t wg = 256; // good default; tune if needed
    // make global size a multiple of wg and >= N (via grid-stride loop)
    const size_t g  = ((std::max(N, wg) + wg - 1) / wg) * wg;
    mem_t<int32_t> obs(g*table_size, q);
    mem_t<int32_t> mis(g*table_size, q);
    auto* local_obs = obs.data();
    auto* local_mis = mis.data();

    q.submit([&](sycl::handler& cgh) {
        // local histograms (Observations and Mismatches)
        //sycl::local_accessor<int32_t, 1> local_obs (sycl::range<1>(table_size), cgh);
        //sycl::local_accessor<int32_t, 1> local_mis (sycl::range<1>(table_size), cgh);

        cgh.parallel_for(
            sycl::nd_range<1>(sycl::range<1>(g), sycl::range<1>(wg)),
            [=](sycl::nd_item<1> it) {
                const int  lid    = it.get_local_id(0);
                const int  lsize  = it.get_local_range(0);
                const size_t gsize = it.get_global_range(0);
                size_t i = it.get_global_id(0);
                const size_t group  = it.get_group(0);

                const size_t base = group * table_size; // this WG's private histogram

                // --- Zero local histograms in parallel
                for (int b = lid; b < table_size; b += lsize) {
                    local_obs[base + b] = 0;
                    local_mis[base + b] = 0;
                }
                it.barrier(sycl::access::fence_space::local_space);

                // --- Grid-stride loop over all elements
                for (; i < N; i += gsize) {
                    const int idx = d_idx[i];
                    if ((unsigned)idx < (unsigned)table_size) {
                        // local-space relaxed atomics to reduce bank conflicts
                        auto a_obs = sycl::atomic_ref<
                            int32_t,
                            sycl::memory_order::relaxed,
                            sycl::memory_scope::work_group,
                            sycl::access::address_space::global_space>(local_obs[base + idx]);
                        auto a_mis = sycl::atomic_ref<
                            int32_t,
                            sycl::memory_order::relaxed,
                            sycl::memory_scope::work_group,
                            sycl::access::address_space::global_space>(local_mis[base + idx]);
                        a_mis.fetch_add(d_val[i]);
                        a_obs.fetch_add(1);
                    }
                }

                it.barrier(sycl::access::fence_space::local_space);

                // --- Flush: first BINS threads (or strided) push to global with one atomic each
                for (int b = lid; b < table_size; b += lsize) {
                    const int obs = local_obs[base + b];
                    const int mis = local_mis[base + b];
                    if (obs == 0 && mis == 0) continue;

                    // Atomics on global table entries
                    auto g_obs = sycl::atomic_ref<
                        int64_t,
                        sycl::memory_order::relaxed,
                        sycl::memory_scope::device,
                        sycl::access::address_space::global_space>(d_table[b].Observations);
                    auto g_mis = sycl::atomic_ref<
                        int64_t,
                        sycl::memory_order::relaxed,
                        sycl::memory_scope::device,
                        sycl::access::address_space::global_space>(d_table[b].Mismatches);
                        g_obs.fetch_add(obs);
                        g_mis.fetch_add(mis);
                }
            }
        );
    });
}

//[[clang::optnone]]
void compute_skip_slice_kernel(std::shared_ptr<Batch<ComputeSkipSliceKernelData, ComputeSkipSliceConfiguration>> batch) {
    if (batch->_n <= 0) {
        return;
    }
    auto& profile = ProfilerBQSR::instance();
    Timer t;
    auto d_recalibrated = batch->_d_data.d_recalibrated;
    auto d_sequences = batch->_d_data.d_sequences;
    auto d_sequences_ptr = batch->_d_data.d_sequences_ptr;
    auto d_qualities = batch->_d_data.d_qualities;
    auto d_cigar = batch->_d_data.d_cigar;
    auto d_cigar_ptr = batch->_d_data.d_cigar_ptr;
    auto d_clipped_cigar = batch->_d_data.d_clipped_cigar;
    auto d_clipped_cigar_start = batch->_d_data.d_clipped_cigar_start;
    auto d_clipped_cigar_len = batch->_d_data.d_clipped_cigar_len;
    auto d_pos = batch->_d_data.d_pos;
    auto d_rname = batch->_d_data.d_rname;
    auto d_rnext = batch->_d_data.d_rnext;
    auto d_pnext = batch->_d_data.d_pnext;
    auto d_tlen = batch->_d_data.d_tlen;
    auto d_known_site = batch->_d_data.d_known_sites;
    auto d_known_site_ptr = batch->_d_data.d_known_sites_ptr;
    auto d_reference = batch->_d_data.d_reference;
    auto d_reference_ptr = batch->_d_data.d_reference_ptr;

    auto d_flag = batch->_d_data.d_flag;
    auto d_reads_group_covariate = batch->_d_data.d_reads_group_covariate;
    auto d_skip_slice = batch->_d_data.d_skip_slice;
    auto d_snps = batch->_d_data.d_snps;
    auto d_stranded_seq = batch->_d_data.d_stranded_seq;
    auto d_stranded_seq_ptr = batch->_d_data.d_stranded_seq_ptr;
    auto d_base_context = batch->_d_data.d_base_context;
    auto d_base_context_ptr = batch->_d_data.d_base_context_ptr;

    auto d_table_accessor = batch->_d_data.d_table_acc;

    auto d_quality_table_indexes = batch->_d_data.d_quality_table_indexes;
    auto d_cycle_table_indexes = batch->_d_data.d_cycle_table_indexes;
    auto d_context_table_indexes = batch->_d_data.d_context_table_indexes;
    auto d_quality_table_vals = batch->_d_data.d_quality_table_vals;
    auto d_cycle_table_vals = batch->_d_data.d_cycle_table_vals;
    auto d_context_table_vals = batch->_d_data.d_context_table_vals;

    /*auto d_segments_quality = batch->_d_data.d_segments_quality;
    auto d_segments_cycle = batch->_d_data.d_segments_cycle;
    auto d_segments_context = batch->_d_data.d_segments_context;*/

    auto d_quality_table = batch->_d_data.d_quality_table;
    auto d_cycle_table = batch->_d_data.d_cycle_table;
    auto d_context_table = batch->_d_data.d_context_table;

    auto d_quality_table_sz = batch->_d_data.d_quality_table_sz;
    auto d_cycle_table_sz = batch->_d_data.d_cycle_table_sz;
    auto d_context_table_sz = batch->_d_data.d_context_table_sz;

    auto d_temporary_segments_sequences_offsets = batch->_d_data.d_temporary_segments_sequences_offsets;

    size_t offset_sequences = batch->_d_data.d_sequence_span_offset;
    size_t offset_cigar = batch->_d_data.d_cigar_span_offset;

    assert(batch->_d_data.table_sz < 2000000000);

    profile.start(StageId::ClipSequences);
    t.reset();
    batch->_q_process.parallel_for<class SetOffset>(sycl::range<1>(batch->_n), [=](sycl::id<1> idx) {
        auto base = d_cigar_ptr[idx].start - offset_cigar;
        d_cigar_ptr[idx].start -= offset_cigar;

        d_sequences_ptr[idx].start -= offset_sequences;
        d_stranded_seq_ptr[idx].start -= offset_sequences;
        d_base_context_ptr[idx].start -= offset_sequences;
        d_clipped_cigar_start[idx]    = base + (PAD * (size_t)idx);
        d_clipped_cigar_len[idx]      = d_cigar_ptr[idx].len;
    });
    batch->_q_process.wait_and_throw();
    Logger::getLogger().log(LOG_DEBUG, "    [SetOffset kernel] Time: " + std::to_string(t.elapsed()) + " sec");

    // Copy the source cigar ops into their clipped windows
    // (Use your favorite copy kernel; shown here as a simple loop kernel.)
    batch->_q_process.parallel_for<class CopyCigarToClipped>(
        sycl::range<1>(batch->_n), [=](sycl::id<1> id) {
            size_t i = id[0];
            size_t src_s = d_cigar_ptr[i].start;
            size_t L     = d_cigar_ptr[i].len;
            size_t dst_s = d_clipped_cigar_start[i];

            for (size_t k = 0; k < L; ++k) {
                d_clipped_cigar[dst_s + k] = d_cigar[src_s + k];
            }
        });
    batch->_q_process.wait_and_throw();

    batch->_q_process.parallel_for<class Clipping>(sycl::range<1>(batch->_n), [=](sycl::id<1> idx) {
        if (!d_recalibrated[idx]) {
            return;
        }
        SpanField clip_span;
        clip_span.start = d_clipped_cigar_start[idx];
        clip_span.len   = d_clipped_cigar_len[idx];  // initialized to original L

        deviceHardClipAdaptorSequence(d_clipped_cigar, clip_span, d_pos[idx], d_sequences, d_sequences_ptr[idx],
                                      d_qualities, d_sequences_ptr[idx], idx, d_rname[idx],
                                      d_flag[idx],
                                      d_rnext[idx], d_tlen[idx],
                                      d_pnext[idx]);

        if (d_sequences_ptr[idx].len == 0) {
            d_clipped_cigar_start[idx] = clip_span.start;
            d_clipped_cigar_len[idx]   = clip_span.len;
            return;
        }
        deviceHardClipSoftClippedBases(d_clipped_cigar, clip_span, d_pos[idx], d_sequences, d_sequences_ptr[idx],
                                      d_qualities, d_sequences_ptr[idx], d_rname[idx],
                                   d_flag[idx]);

        d_clipped_cigar_start[idx] = clip_span.start;
        d_clipped_cigar_len[idx]   = clip_span.len;
        d_stranded_seq_ptr[idx] = d_sequences_ptr[idx];
        d_base_context_ptr[idx] = d_sequences_ptr[idx];

    });
    batch->_q_process.wait_and_throw();
    Logger::getLogger().log(LOG_DEBUG, "    [Clip Sequences kernel] Time: " + std::to_string(t.elapsed()) + " sec");

    profile.stop(StageId::ClipSequences);

    t.reset();
    profile.start(StageId::ComputeSkipVector);
    Timer t2;
    t2.reset();
    batch->_q_process.parallel_for<class ComputeSkipSliceVector>(sycl::range<1>(batch->_n), [=](sycl::id<1> idx) {

        if (!d_recalibrated[idx]) {
            return;
        }
        if (d_sequences_ptr[idx].len == 0) {
            return;
        }

        int64_t sftStart = 0;
        {
            auto cigar = d_clipped_cigar + d_clipped_cigar_start[idx];
            auto start = d_pos[idx];
            for (size_t i = 0; i < d_clipped_cigar_len[idx]; i++) {
                CigarOperation op = cigar[i];
                if (op.Operation == 'S') {
                    start -= op.Length;
                } else if (op.Operation != 'H') {
                    break;
                }
            }
            sftStart = static_cast<int32_t>(start);
        }

        /*
         * Reduce by segment
         * Segments: cigars
         * Return: length per string
         */
        int32_t length = 0;
        for (size_t i = 0; i < d_clipped_cigar_len[idx]; i++) {
            CigarOperation op = d_clipped_cigar[d_clipped_cigar_start[idx] + i];
            // We can comment out this
            if (CigarOperatorConsumesReferenceBases[static_cast<uint8_t>(op.Operation)] > 0) {
                length += static_cast<int32_t>(CigarOperatorConsumesReferenceBases[static_cast<uint8_t>(op.Operation)] * op.Length);
            }
        }

        int64_t end = d_pos[idx] + length - 1;

        int64_t sftEnd = end;
        {
            bool saw_nonHS = false;
            auto cigar = d_clipped_cigar + d_clipped_cigar_start[idx];

            for (int64_t i = static_cast<int64_t>(d_clipped_cigar_len[idx]) - 1; i >= 0; i--) {
                 CigarOperation op = cigar[i];
                 if (op.Operation == 'S') {
                     sftEnd += static_cast<int32_t>(op.Length);
                 } else if (op.Operation != 'H') {
                     //return softEnd;
                     saw_nonHS = true;
                     break;
                 }
             }
            if (!saw_nonHS) {
                sftEnd = end;
            }
        }
        if (d_known_site_ptr[d_rname[idx]].len > 0) {
            SpanField cigar_ptr(d_clipped_cigar_start[idx], d_clipped_cigar_len[idx]);
            calculateSkipSlice(d_sequences, d_sequences_ptr[idx], d_qualities, d_sequences_ptr[idx]/*d_qualities_ptr[idx]*/,
                               d_clipped_cigar, cigar_ptr, d_pos[idx], sftStart, sftEnd,
                               d_known_site, d_known_site_ptr[d_rname[idx]],
                               d_skip_slice + d_sequences_ptr[idx].start/*d_skip_slice_ptr[idx].start*/);
        }
    });
    batch->_q_process.wait_and_throw();
    Logger::getLogger().log(LOG_DEBUG, "        [ComputeSkipSliceVector kernel] Time: " + std::to_string(t2.elapsed()) + " sec");
    profile.stop(StageId::ComputeSkipVector);

    profile.start(StageId::ComputeSNPSVector);
    t2.reset();
    {
        batch->_q_process.parallel_for<class ComputeSNPSVector>(sycl::range<1>(batch->_n), [=](sycl::id<1> idx) {

            if (!d_recalibrated[idx]) {
                return;
            }
            if (d_sequences_ptr[idx].len == 0) {
                return;
            }

            computeSnpEvents(d_sequences, d_sequences_ptr[idx], d_pos[idx],
                             d_clipped_cigar, d_clipped_cigar_start[idx], d_clipped_cigar_len[idx], d_reference, d_reference_ptr[d_rname[idx]],
                            d_snps + d_sequences_ptr[idx].start/*d_snps_ptr[idx].start*/);
    });
        batch->_q_process.wait_and_throw();
    }

    Logger::getLogger().log(LOG_DEBUG, "        [ComputeSNPSVector kernel] Time: " + std::to_string(t2.elapsed()) + " sec");
    profile.stop(StageId::ComputeSNPSVector);

    profile.start(StageId::ComputeContextVector);
    t2.reset();

    batch->_q_process.parallel_for<class ComputeContextKeysVector>(sycl::range<1>(batch->_n), [=](sycl::id<1> idx) {

        if (!d_recalibrated[idx]) {
            return;
        }
        if (d_sequences_ptr[idx].len == 0) {
            return;
        }

        bool is_reversed = SamAlignments::IsReversed(d_flag[idx]);
        computeStrandedClippedSeq(d_sequences, d_sequences_ptr[idx], d_qualities, d_sequences_ptr[idx]/*d_qualities_ptr[idx]*/,
                                  is_reversed, d_stranded_seq, d_stranded_seq_ptr[idx]);
        computeBaseContextCovariate(is_reversed, d_stranded_seq, d_stranded_seq_ptr[idx],
                                    d_base_context + d_base_context_ptr[idx].start,
                                    d_base_context_ptr[idx].len);
    });
    batch->_q_process.wait_and_throw();

    Logger::getLogger().log(LOG_DEBUG, "        [ComputeContextKeysVector kernel] Time: " + std::to_string(t2.elapsed()) + " sec");
    profile.stop(StageId::ComputeContextVector);

    Logger::getLogger().log(LOG_DEBUG, "    [ComputeCovariateVectors kernel] Time: " + std::to_string(t.elapsed()) + " sec");

    profile.start(StageId::UpdateIndexes);
    profile.start(StageId::Preprocessing);

    t.reset();
    {
        /*enum { nt = 256, vt = 1, vt0 = 5, nv = nt*vt };
        uint32_t num_ctas = div_up(batch->_d_data.table_sz, size_t(nv));

        batch->_q_process.parallel_for<class SetDefaultValues>(
                sycl::nd_range<1>(sycl::range<1>(num_ctas * nt),
                sycl::range<1>(nt)),
                [=, count = batch->_d_data.table_sz](sycl::nd_item<1> item) {

                    int tid = item.get_local_id(); // Current work-item index
                    int group_id = item.get_group(0); // Workgroup index
                    size_t global_base = group_id * vt * nt; // Global offset for this workgroup

                    strided_iterate<nt, vt>([&](int i, int j) {
                        d_quality_table_indexes[global_base + j] = INVALID_INDEX;
                    }, tid, count);

                    strided_iterate<nt, vt>([&](int i, int j) {
                        d_cycle_table_indexes[global_base + j] = INVALID_INDEX;
                    }, tid, count);

                    strided_iterate<nt, vt>([&](int i, int j) {
                        d_context_table_indexes[global_base + j] = INVALID_INDEX;
                    }, tid, count);

                    strided_iterate<nt, vt>([&](int i, int j) {
                        d_quality_table_vals[global_base + j] = 0;
                    }, tid, count);

                    strided_iterate<nt, vt>([&](int i, int j) {
                        d_cycle_table_vals[global_base + j] = 0;
                    }, tid, count);

                    strided_iterate<nt, vt>([&](int i, int j) {
                        d_context_table_vals[global_base + j] = 0;
                    }, tid, count);

                    strided_iterate<nt, vt>([&](int i, int j) {
                        d_segments_quality[global_base + j] = 0;
                    }, tid, count);

                    strided_iterate<nt, vt>([&](int i, int j) {
                        d_segments_cycle[global_base + j] = 0;
                    }, tid, count);

                    strided_iterate<nt, vt>([&](int i, int j) {
                        d_segments_context[global_base + j] = 0;
                    }, tid, count);
        });
    }*/

        const size_t n = batch->_d_data.table_sz;

        // 0-inits are free of kernels:
        batch->_q_process.fill(d_quality_table_vals,  static_cast<val_t>(0), n);
        batch->_q_process.fill(d_cycle_table_vals,    static_cast<val_t>(0), n);
        batch->_q_process.fill(d_context_table_vals,  static_cast<val_t>(0), n);
        /*batch->_q_process.memset(d_segments_quality,    0, n * sizeof(*d_segments_quality));
        batch->_q_process.memset(d_segments_cycle,      0, n * sizeof(*d_segments_cycle));
        batch->_q_process.memset(d_segments_context,    0, n * sizeof(*d_segments_context));
*/
        // If INVALID_INDEX == -1 and index type is 2's-complement int, this is valid and fastest:
        batch->_q_process.fill(d_quality_table_indexes, static_cast<int32_t>(INVALID_INDEX), n);
        batch->_q_process.fill(d_cycle_table_indexes,   static_cast<int32_t>(INVALID_INDEX), n);
        batch->_q_process.fill(d_context_table_indexes, static_cast<int32_t>(INVALID_INDEX), n);
    }

    batch->_q_process.wait_and_throw();
    Logger::getLogger().log(LOG_DEBUG, "    [SetDefaultValues kernel] Time: " + std::to_string(t.elapsed()) + " sec");

    t.reset();

    batch->_q_process.parallel_for<class CopySequencesOffset>(sycl::range<1>(batch->_n), [=](sycl::id<1> idx) {
        d_temporary_segments_sequences_offsets[idx] = d_sequences_ptr[idx].start - d_sequences_ptr[0].start; // It can also be:  -d_sequences_ptr->start[0]
    });
    batch->_q_process.wait_and_throw();
    Logger::getLogger().log(LOG_DEBUG, "    [CopyingSequencesOffsets kernel] Time: " + std::to_string(t.elapsed()) + " sec");

    t.reset();
    profile.stop(StageId::Preprocessing);
    profile.start(StageId::EmitIndexes);

       {
        enum { nt = 128, vt = 7, vt0 = 5, nv = nt*vt };
        typedef cta_load_balance_t<nt, vt> load_balance_t;

        union SharedMemory {
            load_balance_t::storage_t lbs;
        };


        uint32_t num_ctas = div_up(batch->_d_data.table_sz + batch->_n, size_t(nv));
        mem_t<int> mp = load_balance_partitions(static_cast<int64_t>(batch->_d_data.table_sz), d_temporary_segments_sequences_offsets, batch->_n, nv, batch->_d_data._q);

        auto mp_data = mp.data();

        batch->_q_process.submit([&](sycl::handler& cgh) {
            sycl::local_accessor<SharedMemory, 1> shared(sycl::range<1>(1), cgh);

            cgh.parallel_for<class ComputeCovariateTableIndexes>(
                sycl::nd_range<1>(sycl::range<1>(num_ctas * nt),
                sycl::range<1>(nt)),
                [=, count=batch->_d_data.table_sz, num_reads=(batch->_n)](sycl::nd_item<1> item) {

                    int tid = item.get_local_id(); // Current work-item index
                    int cta = item.get_group(0); // Current work-group index

                    auto lbs = cta_load_balance_t<nt, vt>().load_balance(count, d_temporary_segments_sequences_offsets, num_reads, tid, cta, mp_data, shared[0].lbs, item);

                    strided_iterate<nt, vt, vt>([&](int k, int j) {
                        int output_i = lbs.merge_range.a_begin + j;
                        int idx = lbs.segments[k];
                        int i = lbs.ranks[k];
                        d_quality_table_indexes[output_i]++;
                        if (!d_recalibrated[idx]) {
                            return;
                        }

                        if (d_sequences_ptr[idx].len == 0) {
                            return;
                        }

                        size_t l = d_sequences_ptr[idx].len;
                        int32_t readGrpCov = d_reads_group_covariate[idx];

                        auto [cycleFactor, cycleIncrement] = prepareCycleCovariates(d_sequences_ptr[idx].len, d_flag[idx]);

                        if (i >= l) {
                            return;
                        }

                        if (d_skip_slice[d_sequences_ptr[idx].start + i] > 0) {
                            return;
                        }

                        int32_t ok = simpleBaseToBaseIndexTable[static_cast<uint8_t>(d_sequences[output_i])];

                        if (ok < 0) {
                            // continue;
                            return;
                        }
                        char qual = d_qualities[output_i];

                        if (qual < MIN_INTERESTING_QUAL) {
                            //    continue;
                            return;
                        }

                        int32_t errVal = d_snps[output_i];

                        bqsrTableKey key1(qual, 0, readGrpCov);
                        auto idx1 = static_cast<int32_t>(d_table_accessor.getQualityScoresIdx(key1));

                        d_quality_table_indexes[output_i] = idx1;
                        d_quality_table_vals[output_i] = errVal;

                        auto baseCycleCovariate = computeBaseCycleCovariate(cycleFactor, cycleIncrement, i, bqsr::defaults::MAX_CYCLE);

                        bqsrTableKey key2(qual, baseCycleCovariate, readGrpCov);
                        auto idx2 = static_cast<int32_t>(d_table_accessor.getCycleCovIdx(key2));

                        d_cycle_table_indexes[output_i] = idx2;
                        d_cycle_table_vals[output_i] = errVal;


                        if (d_base_context_ptr[idx].len > 0 && d_base_context[output_i] >= 0) {
                            bqsrTableKey key3(qual, d_base_context[output_i], readGrpCov);
                            auto idx3 = static_cast<int32_t>(d_table_accessor.getContextCovIdx(key3));

                            d_context_table_indexes[output_i] = idx3;
                            d_context_table_vals[output_i] = errVal;
                        }
                    }, tid, lbs.merge_range.a_count());
            });
        });
    }
    batch->_q_process.wait_and_throw();
    Logger::getLogger().log(LOG_DEBUG, "    [ComputeCovariateTableIndexes kernel] Time: " + std::to_string(t.elapsed()) + " sec");

    profile.stop(StageId::EmitIndexes);

    profile.start(StageId::AggregateQual);
    aggregate_quality(batch->_q_process,
                  d_quality_table_indexes, d_quality_table_vals, batch->_d_data.table_sz,
                  d_quality_table, d_quality_table_sz);
    batch->_q_process.wait_and_throw();
    profile.stop(StageId::AggregateQual);

    profile.start(StageId::AggregateCycle);
    aggregate_cycle(batch->_q_process,
                    d_cycle_table_indexes, d_cycle_table_vals, batch->_d_data.table_sz,
                    d_cycle_table, d_cycle_table_sz);
    batch->_q_process.wait_and_throw();
    profile.stop(StageId::AggregateCycle);

    profile.start(StageId::AggregateContext);
    aggregate_context(batch->_q_process,
                      d_context_table_indexes, d_context_table_vals, batch->_d_data.table_sz,
                      d_context_table, d_context_table_sz);

    batch->_q_process.wait_and_throw();
    profile.stop(StageId::AggregateContext);

    Logger::getLogger().log(LOG_DEBUG, "    [Tables updates] Time: " + std::to_string(t.elapsed()) + " sec");
    profile.stop(StageId::UpdateIndexes);
}
