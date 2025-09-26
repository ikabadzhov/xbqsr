#include "ComputeCovariateVectorsFunctions.h"

#include <boost/iostreams/categories.hpp>

#include "../../sycl_parallel_library/sorting/mergesort.h"
#include "../../timer/Timer.h"
#include "../../sycl_parallel_library/segmented_reduce/segreduce.h"
#include "../../sycl_parallel_library/intervalsmove/intervalmove.h"
#include "ComputeCovariateVectorsData.h"
#include "../../logger/Logger.h"
#include "../../timer/Profiler.h"

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
    auto d_cigar_start = batch->_d_data.d_cigar_start;
    auto d_cigar_len = batch->_d_data.d_cigar_len;
    auto d_pos = batch->_d_data.d_pos;
    auto d_rname = batch->_d_data.d_rname;
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

    auto d_sz_qual_score = batch->_d_data.d_sz_qual_score_cov;
    auto d_sz_cycle_cov = batch->_d_data.d_sz_cycle_cov;
    auto d_sz_context_cov = batch->_d_data.d_sz_context_cov;

    auto d_quality_table_indexes = batch->_d_data.d_quality_table_indexes;
    auto d_cycle_table_indexes = batch->_d_data.d_cycle_table_indexes;
    auto d_context_table_indexes = batch->_d_data.d_context_table_indexes;
    auto d_quality_table_vals = batch->_d_data.d_quality_table_vals;
    auto d_cycle_table_vals = batch->_d_data.d_cycle_table_vals;
    auto d_context_table_vals = batch->_d_data.d_context_table_vals;

    auto d_segments_quality = batch->_d_data.d_segments_quality;
    auto d_segments_cycle = batch->_d_data.d_segments_cycle;
    auto d_segments_context = batch->_d_data.d_segments_context;

    auto d_quality_table = batch->_d_data.d_quality_table;
    auto d_cycle_table = batch->_d_data.d_cycle_table;
    auto d_context_table = batch->_d_data.d_context_table;


    auto d_temporary_keys = batch->_d_data.d_temporary_keys;
    auto d_temporary_vals = batch->_d_data.d_temporary_vals;
    auto d_temporary_segments_sequences_offsets = batch->_d_data.d_temporary_segments_sequences_offsets;

    auto d_read_index = batch->_d_data.d_read_index;
    auto d_i_offset = batch->_d_data.d_i_offset;
    auto d_j_offset = batch->_d_data.d_j_offset;

    size_t offset_sequences = batch->_d_data.d_sequence_span_offset;
    size_t offset_cigar = batch->_d_data.d_cigar_span_offset;


    t.reset();
    batch->_q_process.parallel_for<class SetOffset>(sycl::range<1>(batch->_n), [=](sycl::id<1> idx) {
        d_cigar_start[idx] -= offset_cigar;

        d_sequences_ptr[idx].start -= offset_sequences;
        d_stranded_seq_ptr[idx].start -= offset_sequences;
        d_base_context_ptr[idx].start -= offset_sequences;
    });
    batch->_q_process.wait_and_throw();
    Logger::getLogger().log(LOG_DEBUG, "    [SetOffset kernel] Time: " + std::to_string(t.elapsed()) + " sec");

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
            auto cigar = d_cigar + d_cigar_start[idx];
            auto start = d_pos[idx];
            for (size_t i = 0; i < d_cigar_len[idx]; i++) {
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
        for (size_t i = 0; i < d_cigar_len[idx]; i++) {
            CigarOperation op = d_cigar[d_cigar_start[idx] + i];
            // We can comment out this
            if (CigarOperatorConsumesReferenceBases[static_cast<uint8_t>(op.Operation)] > 0) {
                length += static_cast<int32_t>(CigarOperatorConsumesReferenceBases[static_cast<uint8_t>(op.Operation)] * op.Length);
            }
        }

        int64_t end = d_pos[idx] + length - 1;

        int64_t sftEnd = end;
        {
            bool saw_nonHS = false;
            auto cigar = d_cigar + d_cigar_start[idx];

            for (int64_t i = static_cast<int64_t>(d_cigar_len[idx]) - 1; i >= 0; i--) {
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
            SpanField cigar_ptr(d_cigar_start[idx], d_cigar_len[idx]);
            calculateSkipSlice(d_sequences, d_sequences_ptr[idx], d_qualities, d_sequences_ptr[idx]/*d_qualities_ptr[idx]*/,
                               d_cigar, cigar_ptr, d_pos[idx], sftStart, sftEnd,
                               d_known_site, d_known_site_ptr[d_rname[idx]],
                               d_skip_slice + d_sequences_ptr[idx].start/*d_skip_slice_ptr[idx].start*/);
        }
    });
    batch->_q_process.wait_and_throw();
    Logger::getLogger().log(LOG_DEBUG, "        [ComputeSkipSliceVector kernel] Time: " + std::to_string(t2.elapsed()) + " sec");
    profile.stop(StageId::ComputeSkipVector);
    //  cigar :                 || 14I  13I  9I   2D   10I  3I  ||  6I   9I   1M   8I ||  12D  13M  3D   3I   8I ||  15M  13M  15M ||
    //  cigars_op_segments :    || 0    14   27   36   38   48      51   57   66   67     75 ...
    //  i_offset :              || 14   27   36   36   46   49      55   64 ....

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
                             d_cigar, d_cigar_start[idx], d_cigar_len[idx], d_reference, d_reference_ptr[d_rname[idx]],
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

    t.reset();
    {
        enum { nt = 256, vt = 1, vt0 = 5, nv = nt*vt };
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

                        //for (size_t i = 0; i < l; i++) {
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

                        if (qual < minInterestingQual) {
                            //    continue;
                            return;
                        }

                        int32_t errVal = d_snps[output_i];

                        bqsrTableKey key1(qual, 0, readGrpCov);
                        auto idx1 = static_cast<int32_t>(BaseRecalibratorTables::getQualityScoresIdx(key1, d_sz_qual_score));

                        d_quality_table_indexes[output_i] = idx1;
                        d_quality_table_vals[output_i] = errVal;

                        auto baseCycleCovariate = computeBaseCycleCovariate(cycleFactor, cycleIncrement, i, MAX_CYCLE);

                        bqsrTableKey key2(qual, baseCycleCovariate, readGrpCov);
                        auto idx2 = static_cast<int32_t>(BaseRecalibratorTables::getCycleCovIdx(key2, d_sz_qual_score, d_sz_cycle_cov));

                        d_cycle_table_indexes[output_i] = idx2;
                        d_cycle_table_vals[output_i] = errVal;


                        if (d_base_context_ptr[idx].len > 0 && d_base_context[output_i] >= 0) {
                            bqsrTableKey key3(qual, d_base_context[output_i], readGrpCov);
                            auto idx3 = static_cast<int32_t>(BaseRecalibratorTables::getContextCovIdx(key3, d_sz_qual_score, d_sz_cycle_cov, d_sz_context_cov));

                            d_context_table_indexes[output_i] = idx3;
                            d_context_table_vals[output_i] = errVal;
                        }
                    }, tid, lbs.merge_range.a_count());
            });
        });
    }
    batch->_q_process.wait_and_throw();
    Logger::getLogger().log(LOG_DEBUG, "    [ComputeCovariateTableIndexes kernel] Time: " + std::to_string(t.elapsed()) + " sec");

    t.reset();
    const size_t batch_size = 2000000000;
    size_t total = batch->_d_data.table_sz;
    for (size_t offset = 0; offset < total; offset += batch_size) {
        //size_t current_size = min(batch_size, total - offset);

        /**
         * Quality
         */
        mergesort_developer_version(d_quality_table_indexes, d_quality_table_vals, batch->_d_data.table_sz, less_t<key_t>(),
                  batch->_q_process, d_temporary_keys, d_temporary_vals);

        Logger::getLogger().log(LOG_DEBUG, "    [Mergesort] Time: " + std::to_string(t.elapsed()) + " sec");

        int32_t num_segments_quality = 0;
        int32_t num_segments_cycle = 0;
        int32_t num_segments_context = 0;
        t.reset();

        reduce_by_key_preprocessing(d_quality_table_indexes, batch->_d_data.table_sz, d_segments_quality,
                                    d_temporary_keys, num_segments_quality, plus_t<key_t>(), equal_to_t<key_t>(),
                                    val_t(0), batch->_q_process);

        Logger::getLogger().log(LOG_DEBUG,
                                "    [Reduce by Key Preprocessing] Time: " + std::to_string(t.elapsed()) + " sec");

        t.reset();

        segreduce(d_quality_table_vals, batch->_d_data.table_sz, d_segments_quality, num_segments_quality,
                  d_quality_table_vals, plus_t<val_t>(), 0, batch->_q_process);

        Logger::getLogger().log(LOG_DEBUG, "    [Segreduce] Time: " + std::to_string(t.elapsed()) + " sec");

        if (num_segments_quality > 1) {

            batch->_q_process.parallel_for(sycl::range<1>(num_segments_quality - 1), [=](sycl::id<1> segment_idx) {

                int32_t covariate_idx = d_temporary_keys[segment_idx];
                int32_t err_val = d_quality_table_vals[segment_idx];
                int32_t segment_sz = d_segments_quality[segment_idx + 1] - d_segments_quality[segment_idx];

                if (covariate_idx == INVALID_INDEX) {
                    return;
                }
                if (d_quality_table[covariate_idx].Observations == -1) {
                    d_quality_table[covariate_idx].Observations = 0;
                }
                d_quality_table[covariate_idx].Mismatches += err_val;
                d_quality_table[covariate_idx].Observations += segment_sz;
            });
        }
        batch->_q_process.wait_and_throw();

        /**
         * Cycle
         */
        t.reset();

        mergesort_developer_version(d_cycle_table_indexes, d_cycle_table_vals, batch->_d_data.table_sz, less_t<key_t>(),
                                    batch->_q_process, d_temporary_keys, d_temporary_vals);

        Logger::getLogger().log(LOG_DEBUG, "    [Mergesort] Time: " + std::to_string(t.elapsed()) + " sec");

        t.reset();

        reduce_by_key_preprocessing(d_cycle_table_indexes, batch->_d_data.table_sz, d_segments_cycle,
                                    d_temporary_keys, num_segments_cycle, plus_t<key_t>(), equal_to_t<key_t>(),
                                    val_t(0), batch->_q_process);

        Logger::getLogger().log(LOG_DEBUG,
                                "    [Reduce by Key Preprocessing] Time: " + std::to_string(t.elapsed()) + " sec");

        t.reset();

        segreduce(d_cycle_table_vals, batch->_d_data.table_sz, d_segments_cycle, num_segments_cycle, d_cycle_table_vals,
                  plus_t<val_t>(), 0, batch->_q_process);

        Logger::getLogger().log(LOG_DEBUG, "    [Segreduce] Time: " + std::to_string(t.elapsed()) + " sec");

        if (num_segments_cycle > 1) {
            batch->_q_process.parallel_for(sycl::range<1>(num_segments_cycle - 1), [=](sycl::id<1> segment_idx) {

                int32_t covariate_idx = d_temporary_keys[segment_idx];
                int32_t err_val = d_cycle_table_vals[segment_idx];
                int32_t segment_sz = d_segments_cycle[segment_idx + 1] - d_segments_cycle[segment_idx];

                if (covariate_idx == INVALID_INDEX) {
                    return;
                }
                if (d_cycle_table[covariate_idx].Observations == -1) {
                    d_cycle_table[covariate_idx].Observations = 0;
                }
                d_cycle_table[covariate_idx].Mismatches += err_val;
                d_cycle_table[covariate_idx].Observations += segment_sz;
            });
        }
        batch->_q_process.wait_and_throw();

        /**
        * Context
        */
        t.reset();

        mergesort_developer_version(d_context_table_indexes, d_context_table_vals, batch->_d_data.table_sz, less_t<key_t>(),
                                    batch->_q_process, d_temporary_keys, d_temporary_vals);

        Logger::getLogger().log(LOG_DEBUG, "    [Mergesort] Time: " + std::to_string(t.elapsed()) + " sec");

        t.reset();
        reduce_by_key_preprocessing(d_context_table_indexes, batch->_d_data.table_sz, d_segments_context,
                                    d_temporary_keys, num_segments_context,
                                    plus_t<key_t>(), equal_to_t<key_t>(), val_t(0), batch->_q_process);

        Logger::getLogger().log(LOG_DEBUG,
                                "    [Reduce by Key Preprocessing] Time: " + std::to_string(t.elapsed()) + " sec");

        t.reset();
        segreduce(d_context_table_vals, batch->_d_data.table_sz, d_segments_context, num_segments_context,
                  d_context_table_vals, plus_t<val_t>(), 0, batch->_q_process);

        batch->_q_process.wait_and_throw();

        Logger::getLogger().log(LOG_DEBUG, "    [Segreduce] Time: " + std::to_string(t.elapsed()) + " sec");

        t.reset();

        if (num_segments_context > 1) {

            batch->_q_process.parallel_for(sycl::range<1>(num_segments_context - 1), [=](sycl::id<1> segment_idx) {

                int32_t covariate_idx = d_temporary_keys[segment_idx];
                int32_t err_val = d_context_table_vals[segment_idx];
                int32_t segment_sz = d_segments_context[segment_idx + 1] - d_segments_context[segment_idx];

                if (covariate_idx == INVALID_INDEX) {
                    return;
                }
                if (d_context_table[covariate_idx].Observations == -1) {
                    d_context_table[covariate_idx].Observations = 0;
                }
                d_context_table[covariate_idx].Mismatches += err_val;
                d_context_table[covariate_idx].Observations += segment_sz;
            });
        }
        batch->_q_process.wait_and_throw();
    }
    Logger::getLogger().log(LOG_DEBUG, "    [Tables updates] Time: " + std::to_string(t.elapsed()) + " sec");
    profile.stop(StageId::UpdateIndexes);
}
