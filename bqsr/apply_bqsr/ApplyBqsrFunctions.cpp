#include "ApplyBqsrFunctions.h"
#include "../../sycl_parallel_library/sorting/mergesort.h"
#include "../../timer/Timer.h"
#include "../../logger/Logger.h"
#include "../../sycl_parallel_library/segmented_reduce/segreduce.h"
#include "../../sycl_parallel_library/intervalsmove/intervalmove.h"
#include "../TablesBQSR.h"
#include "../BQSR.h"


const int32_t globalQualityScorePrior = -1;

// //[[clang::optnone]]
 void apply_bqsr_kernel(std::shared_ptr<Batch<AppyBqsrKernelData, AppyBqsrKernelConfiguration>> batch) {
    if (batch->_n <= 0) {
        return;
    }

    Timer t;

    auto d_combined_bqsr_table = batch->_d_data.d_combined_bqsr_table;
    auto d_quantized_quality_scores = batch->_d_data.d_quantized_quality_scores;
    auto d_static_quantized_scores = batch->_d_data.d_static_quantized_scores;
//    auto d_empirical_read_group_entrys = batch->_d_data.d_empirical_read_group_entries;
    auto static_quantized_scores_size = batch->_d_data.static_quantized_scores_size;
    auto d_sequences = batch->_d_data.d_sequences;
    auto d_sequences_ptr = batch->_d_data.d_sequences_ptr;
    auto d_qualities = batch->_d_data.d_qualities;

    auto d_flag = batch->_d_data.d_flag;
    auto d_reads_group_covariate = batch->_d_data.d_reads_group_covariate;

    auto d_stranded_seq = batch->_d_data.d_stranded_seq;
    auto d_stranded_seq_ptr = batch->_d_data.d_stranded_seq_ptr;
    auto d_base_context = batch->_d_data.d_base_context;
    auto d_base_context_ptr = batch->_d_data.d_base_context_ptr;

    auto d_temporary_segments_sequences_offsets = batch->_d_data.d_temporary_segments_sequences_offsets;

    auto table_acc = batch->_d_data.d_table_acc;

    auto d_quality_table = batch->_d_data.d_quality_table;
    auto d_cycle_table = batch->_d_data.d_cycle_table;
    auto d_context_table = batch->_d_data.d_context_table;

    auto d_recalibrated_score_cache = batch->_d_data.d_recalibrated_score_cache;

    size_t offset_sequences = batch->_d_data.sequence_span_offset;

    t.reset();
    batch->_q_process.parallel_for<class ApplySetOffset>(sycl::range<1>(batch->_n), [=](sycl::id<1> idx) {
        d_sequences_ptr[idx].start -= offset_sequences;
        d_stranded_seq_ptr[idx].start -= offset_sequences;
        d_base_context_ptr[idx].start -= offset_sequences;
    });

    batch->_q_process.wait_and_throw();
    Logger::getLogger().log(LOG_DEBUG, "    [SetOffset kernel] Time: " + std::to_string(t.elapsed()) + " sec");

    t.reset();
    batch->_q_process.parallel_for<class ApplyComputeContextKeysVector>(sycl::range<1>(batch->_n), [=](sycl::id<1> idx) {
        auto readGrpIndex = d_reads_group_covariate[idx];
        if (!d_combined_bqsr_table[readGrpIndex].isValid) {
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
    Logger::getLogger().log(LOG_DEBUG, "    [ComputeBaseContextCovariate kernel] Time: " + std::to_string(t.elapsed()) + " sec");

    t.reset();
    batch->_q_process.parallel_for<class ApplyCopySequencesOffset>(sycl::range<1>(batch->_n), [=](sycl::id<1> idx) {
       d_temporary_segments_sequences_offsets[idx] = d_sequences_ptr[idx].start - d_sequences_ptr[0].start; // It can also be:  -d_sequences_ptr->start[0]
    });

    batch->_q_process.wait_and_throw();
    Logger::getLogger().log(LOG_DEBUG, "    [Extract sequences segments kernel] Time: " + std::to_string(t.elapsed()) + " sec");

    t.reset();

    {
        enum { nt = 128, vt = 7, vt0 = 5, nv = nt*vt };
        typedef cta_load_balance_t<nt, vt> load_balance_t;

        union SharedMemory {
            load_balance_t::storage_t lbs;
        };

        uint32_t num_ctas = div_up(batch->_d_data.chunk_sz + batch->_n, size_t(nv));
        mem_t<int> mp = load_balance_partitions(static_cast<int64_t>(batch->_d_data.chunk_sz), d_temporary_segments_sequences_offsets, batch->_n, nv, batch->_d_data._q);

        auto mp_data = mp.data();

        batch->_q_process.submit([&](sycl::handler& cgh) {
            sycl::local_accessor<SharedMemory, 1> shared(sycl::range<1>(1), cgh);

            cgh.parallel_for<class ApplyComputeCovariateTableIndexes>(
                sycl::nd_range<1>(sycl::range<1>(num_ctas * nt),
                sycl::range<1>(nt)),
                [=, count=batch->_d_data.chunk_sz, num_reads=(batch->_n), GLOBAL_QUALITY_SCORE_PRIOR = globalQualityScorePrior](sycl::nd_item<1> item) {

                    int tid = item.get_local_id(); // Current work-item index
                    int cta = item.get_group(0); // Current work-group index

                    auto lbs = cta_load_balance_t<nt, vt>().load_balance(count, d_temporary_segments_sequences_offsets, num_reads, tid, cta, mp_data, shared[0].lbs, item);

                    strided_iterate<nt, vt, vt>([&](int k, int j) {
                        //int output_i = lbs.merge_range.a_begin + j;
                        int idx = lbs.segments[k];
                        int i = lbs.ranks[k];

                        auto readGrpCov = d_reads_group_covariate[idx];
                        auto empiricalReadGroupEntry = d_combined_bqsr_table[readGrpCov];

                        double epsilon;
                        if (GLOBAL_QUALITY_SCORE_PRIOR > 0.0) {
                            epsilon = GLOBAL_QUALITY_SCORE_PRIOR;
                        } else {
                            epsilon = empiricalReadGroupEntry.reportedQuality;
                        }

                        applyKey key{};
                        key.readGroupCovariate = readGrpCov;

                        int32_t base_i = d_sequences_ptr[idx].start + i;

                        uint8_t qual = d_qualities[base_i];

                        if (qual >= MIN_INTERESTING_QUAL) {

                            auto [cycleFactor, cycleIncrement] = prepareCycleCovariates(d_sequences_ptr[idx].len, d_flag[idx]);

                            key.qual = qual;
                            key.cycleCovariate = computeBaseCycleCovariate(cycleFactor, cycleIncrement, i, bqsr::defaults::MAX_CYCLE);
                            key.contextCovariate = d_base_context[d_base_context_ptr[idx].start + i];

                            size_t buffer_idx = 0;

                            // Precompute strides (compiler should hoist this if invariant, but explicit is better)
                            const int32_t strideCtx   = (int32_t)table_acc.nCycleCov;
                            const int32_t strideQual  = (int32_t)(table_acc.nContextCov + 1) * strideCtx;
                            const int32_t strideRG    = (int32_t)table_acc.nQualScore * strideQual;
                            const int32_t maxCycleVal = (int32_t)bqsr::defaults::MAX_CYCLE;

                            int32_t cyc = key.cycleCovariate;
                            int32_t cyc_idx = (cyc < 0) ? (-cyc + maxCycleVal - 1) : cyc;

                            buffer_idx = (size_t)key.readGroupCovariate * strideRG +
                                         (size_t)key.qual * strideQual +
                                         (size_t)(key.contextCovariate + 1) * strideCtx +
                                         (size_t)cyc_idx;

                            sycl::atomic_ref<int, sycl::memory_order::acq_rel, sycl::memory_scope::device, sycl::access::address_space::global_space> atom_cached_value(d_recalibrated_score_cache[buffer_idx]);
                            int old = atom_cached_value.load();
                            if (old > 0) {
                                d_qualities[base_i] = static_cast<uint8_t>(old);
                                return;
                            }

                            bqsrTableKey bqKey;
                            bqKey.ReadGroup = readGrpCov;
                            bqKey.Qual = qual;
                            auto idx1 = table_acc.getQualityScoresIdx(bqKey);
                            auto empiricalBqEntry = d_quality_table[idx1];

                            bqsrTableKey cycleKey;
                            cycleKey.Covariate = key.cycleCovariate;
                            cycleKey.ReadGroup = readGrpCov;
                            cycleKey.Qual = qual;

                            auto idx2 = static_cast<int32_t>(table_acc.getCycleCovIdx(cycleKey));
                            auto empiricalCycleEntry = d_cycle_table[idx2];

                            bqsrTableKey contextKey;
                            contextKey.Covariate = key.contextCovariate;
                            contextKey.ReadGroup = readGrpCov;
                            contextKey.Qual = qual;

                            bqsrEntry empiricalContextEntry(-1,0);

                            auto idx3 = static_cast<int32_t>(table_acc.getContextCovIdx(contextKey));

                            if (contextKey.Covariate >= 0) {
                                empiricalContextEntry = d_context_table[idx3];
                            }

                            double recalibratedQual = estimateHierarchicalBayesianQuality(epsilon, empiricalReadGroupEntry,
                                                                                          empiricalBqEntry,
                                                                                          empiricalContextEntry,
                                                                                          empiricalCycleEntry);

                            uint8_t recalibratedQualityScore = d_quantized_quality_scores[std::max(1,
                                                                                            std::min(fastRound(recalibratedQual),
                                                                                                   maxRecalibratedQualScore))];
                            if (static_quantized_scores_size > 0) {
                                recalibratedQualityScore = d_static_quantized_scores[recalibratedQualityScore];
                            }


                            if(key.cycleCovariate < 0) {
                                buffer_idx = key.readGroupCovariate * table_acc.nQualScore * (table_acc.nContextCov + 1) * table_acc.nCycleCov + key.qual * (table_acc.nContextCov + 1) * table_acc.nCycleCov +
                                        (key.contextCovariate + 1) * table_acc.nCycleCov + (-1 * key.cycleCovariate + bqsr::defaults::MAX_CYCLE - 1);
                            } else {
                                buffer_idx = key.readGroupCovariate * table_acc.nQualScore * (table_acc.nContextCov + 1) * table_acc.nCycleCov + key.qual * (table_acc.nContextCov + 1) * table_acc.nCycleCov +
                                        (key.contextCovariate + 1) * table_acc.nCycleCov + key.cycleCovariate;
                            }
                            atom_cached_value.store(recalibratedQualityScore);
                            //lld_recalibrated_quality_scores[buffer_idx] = static_cast<int8_t>(recalibratedQualityScore);
                            d_qualities[base_i] = static_cast<char>(recalibratedQualityScore);
                        }
                    }, tid, lbs.merge_range.a_count());
            });
        });
    }
    batch->_q_process.wait_and_throw();
    Logger::getLogger().log(LOG_DEBUG, "    [UpdateBaseQualityScore kernel] Time: " + std::to_string(t.elapsed()) + " sec");
     }