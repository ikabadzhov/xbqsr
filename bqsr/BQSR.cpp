#include "BQSR.h"
#include "../filters/UtilsFilters.h"
#include "TablesBQSR.h"
#include "../sam_data_structures/sam_fields/Sequence.h"
#include "../timer/Timer.h"
#include <iterator>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <iomanip>
#include <execution>
#include <fstream>
#include <pstl/glue_execution_defs.h>
#include <tbb/parallel_for.h>

#include "BaseRecalibratorTables.h"
#include "UtilsBQSR.h"
#include "../Config.h"
#include "compute_vectors_for_covariates/ComputeCovariateVectorsData.h"
#include "compute_vectors_for_covariates/ComputeCovariateVectorsFunctions.h"
#include "../device/batch.h"
#include "../sycl_parallel_library/segmented_reduce/segreduce.h"
#include "../sycl_parallel_library/sorting/mergesort.h"
#include "../sycl_parallel_library/sorting/external_sort.h"
#include "../logger/Logger.h"
#include "../timer/Profiler.h"
#include "apply_bqsr/ApplyBqsrFunctions.h"
#include "apply_bqsr/ApplyBqsrKernelData.h"

const uint32_t lowQualityTail = 2;

//const int32_t globalQualityScorePrior = -1;


int32_t fastRound(double d) {
    return (d > 0.0) ? static_cast<int>(d + 0.5) : static_cast<int>(d - 0.5);
}

SYCL_EXTERNAL int32_t keyFromContext(char* dna, SpanField dnaPtrs, int64_t start, int64_t end) {
    auto key = static_cast<int32_t>(end - start);

    uint32_t bitOffset = LENGTH_BITS;

    for (size_t i = start; i < end; i++) {
        int32_t baseIndex = simpleBaseToBaseIndexTable[static_cast<uint8_t>(dna[dnaPtrs.start + i])];
        if( baseIndex == -1 ){
            return -1;
        }
        key |= baseIndex << bitOffset;
        bitOffset += 2;
    }
    return key;
}

void contextWith(char* bases, SpanField basesPtrs, int32_t contextSize, int32_t* keys, uint16_t& keysLength) {
    int32_t mask = ((1u << (contextSize * 2)) - 1) << LENGTH_BITS;
    size_t readLength = basesPtrs.len;

    int32_t fillCount = std::min(contextSize - 1, static_cast<int32_t>(readLength));
    //std::fill_n(keys, fillCount, -1);
    std::fill(keys, keys + fillCount, -1);
    size_t keyVectorIdx = fillCount;

    if (readLength < contextSize) {
        keysLength = 0;
        return;
    }

    // keys.resize(bases.size());
    keysLength = basesPtrs.len;

    auto newBaseOffset = static_cast<uint32_t>(2 * (contextSize - 1) + LENGTH_BITS);
    int32_t currentKey = keyFromContext(bases, basesPtrs, 0, contextSize);

    keys[keyVectorIdx] = currentKey;
    keyVectorIdx++;
    auto tmp_bases = bases;
    bases += basesPtrs.start;
    int32_t currentNPenalty = 0;
    if (currentKey == -1) {
        currentKey = 0;
        currentNPenalty = contextSize - 1;
        uint32_t offset = newBaseOffset;
        while (bases[currentNPenalty] != 'N') {
            int32_t baseIndex = simpleBaseToBaseIndexTable[static_cast<uint8_t>(bases[currentNPenalty])];
            currentKey |= baseIndex << offset;
            offset -= 2;
            currentNPenalty--;
        }
    }
    auto contextSize32 = static_cast<int32_t>(contextSize);
    for (size_t currentIndex = contextSize; currentIndex < readLength; currentIndex++) {
        int32_t baseIndex = simpleBaseToBaseIndexTable[static_cast<uint8_t>(bases[currentIndex])];
        if (baseIndex == -1) {
            currentNPenalty = contextSize;
            currentKey = 0;
        } else {
            currentKey = (currentKey >> 2) & mask;
            currentKey |= baseIndex << newBaseOffset;
            currentKey |= contextSize32;
        }
        if (currentNPenalty == 0) {
            keys[keyVectorIdx] = currentKey;
            keyVectorIdx++;
        } else {
            currentNPenalty--;
            keys[keyVectorIdx]=-1;
            keyVectorIdx++;
        }
    }
    keysLength = keyVectorIdx;
    bases = tmp_bases;
    tmp_bases = nullptr;
}
//TODO: FIX MISMATCHED_COTEXT_SIZE
void computeBaseContextCovariate(bool isReversed, char* strandedClippedSeq, SpanField& strandedClippedSeqPtrs, int32_t* keys, uint16_t& keyLength) {
    contextWith(strandedClippedSeq, strandedClippedSeqPtrs, MISMATCHES_CONTEXT_SIZE, keys, keyLength);
    //assert(keyLength <= strandedClippedSeqPtrs.len);
    if (isReversed){
        std::reverse(keys, keys + keyLength);
    }
}



bool recalibrateAln(Header& hdr, size_t idx, uint8_t MAPQ, uint16_t FLAG, ssize_t RNAME, int32_t PNEXT, int64_t POS, size_t len_sequence, size_t len_qual, CigarOperation* cigar, SpanField cigarPtrs, int32_t RG) {
    /*if (!reads.GetSR(idx).empty()) {
        return false;
    }*/
    auto isEmptyRG = CompressedStringCollection::isEmpty(RG);

    if (MAPQ > 0 && MAPQ < 255 &&
        SamAlignments::FlagNotAny( Secondary | Duplicate | QCFailed, FLAG) &&
        !isStrictUnmapped(SamAlignments::IsUnmapped(FLAG), RNAME, PNEXT) &&
        POS > 0 && len_sequence > 0 &&
        len_sequence == len_qual &&
        !isEmptyRG && alignmentAgreesWithHeader(hdr, RNAME, POS)) {

        return !cigarContainsN(cigar, cigarPtrs) && ReferenceLengthFromCigar(cigar, cigarPtrs) >= 0 &&
               len_sequence == static_cast<uint32_t>(ReadLengthFromCigar(cigar, cigarPtrs));
    }
    return false;
}
//[[clang::optnone]]
void computeSnpEvents(char* sequence, SpanField sequencePtrs, int32_t POS, CigarOperation* cigar, size_t cigar_start, uint16_t cigar_len, uint8_t* ref, SpanFieldLong& refPtrs, int8_t* snps) {
    //size_t readLength = sequencePtrs.len;

    size_t snps_i=0;

    size_t i = sequencePtrs.start;
    auto j = POS - 1;

    for (size_t c = 0; c < cigar_len; c++){
        CigarOperation cigarOp = cigar[cigar_start + c];
        auto ln = static_cast<int32_t>(cigarOp.Length);
        switch (cigarOp.Operation) {
            case 'M':
            case '=':
            case 'X':
                for (size_t k = 0; k < ln; k++) {
                    if (baseToIntMap[static_cast<uint8_t>(sequence[i])] != baseToIntMap[ref[refPtrs.start + j]]) {
                        snps[snps_i] = 1;
                    }
                    i++;
                    j++;
                    snps_i++;
                }
                break;
            case 'D':
            case 'N':
                j += ln;
                break;
            case 'I':
            case 'S':
                i += ln;
                snps_i+=ln;
                break;
        }
    }
}


std::pair<int32_t, bool> readStartsWithInsertion( std::vector<CigarOperation>& cigarVec ) {

    for( auto& element : cigarVec ){
        switch( element.Operation ){
            case 'I':
                return {element.Length, true};
            case 'H':
            case 'S':
                continue;
            default:
                return {-1, false};
        }
    }
    return {-1, false};
}

void computeStrandedClippedSeq(char* SEQ, SpanField seqPtrs, char* QUAL, SpanField& qualPtrs, bool isReversed, char* newSeq, SpanField& newSeqPtrs) {
    ssize_t seqLength = seqPtrs.len;
    ssize_t leftPos = seqLength;
    for (ssize_t i = 0; i < leftPos; i++) {
        if (QUAL[qualPtrs.start + i] > lowQualityTail) {
            leftPos = i;
            break;
        }
    }
    ssize_t rightPos = leftPos - 1;
    for (ssize_t i = seqLength - 1; i >= leftPos; i--) {
        if (QUAL[qualPtrs.start + i] > lowQualityTail) {
            rightPos = i;
            break;
        }
    }
    if (leftPos > rightPos) {
        newSeqPtrs.len = 0;
        return;
    }
    assert(seqLength <= seqPtrs.len);
    newSeqPtrs.len = seqLength;
    //newSeq.resize(seqLength);

    if (isReversed) {
        std::fill(newSeq + newSeqPtrs.start, newSeq + newSeqPtrs.start + (seqLength - 1 - rightPos), 'N');
        std::fill(newSeq + newSeqPtrs.start + seqLength - leftPos, newSeq + newSeqPtrs.start + seqLength, 'N');
        std::transform(
            std::make_reverse_iterator(SEQ + seqPtrs.start + rightPos + 1),
            std::make_reverse_iterator(SEQ + seqPtrs.start + leftPos),
            newSeq + newSeqPtrs.start + (seqLength - 1 - rightPos),
            [](uint8_t base) {
                return baseComplementTable[base];
            }
        );
    } else {
        std::fill(newSeq + newSeqPtrs.start, newSeq + newSeqPtrs.start + leftPos, 'N');
        std::copy(SEQ + seqPtrs.start + leftPos, SEQ + seqPtrs.start + rightPos + 1, newSeq + newSeqPtrs.start + leftPos);
        std::fill(newSeq + newSeqPtrs.start + rightPos + 1, newSeq + newSeqPtrs.start + seqLength, 'N');
    }
}

int32_t checkCycleCovariate(int32_t cycle, int32_t maxCycle) {
    int32_t result = std::abs(cycle);
    if( result > maxCycle){
        std::cout<<"cycle value exceeds maximum cycle value"<<std::endl;
    }
    return (result << 1) | (cycle < 0);
}

const int32_t reversedShift = 4;
const int32_t lastShift = 7;



CycleCovariatesResult prepareCycleCovariates(uint16_t sequence_len, uint16_t FLAG) {
    auto reversed = static_cast<int32_t>((FLAG & Reversed) >> reversedShift);
    auto last = static_cast<int32_t>((FLAG & Last) >> lastShift);
    int32_t readOrderFactor = 1 - 2*last;
    int32_t cycleFactor = readOrderFactor + reversed * (sequence_len - 1) * readOrderFactor;
    int32_t increment = (1 - 2 * reversed) * readOrderFactor;
    return {cycleFactor, increment};
}

int32_t computeBaseCycleCovariate( int32_t cycleFactor, int32_t increment, uint32_t index, int32_t maxCycle ) {
    int32_t cycle = cycleFactor+index*increment;
    int32_t result = std::abs(cycle);
    assert(result < maxCycle);
    result <<= 1;
    result = result | static_cast<int32_t>(cycle<0);
    return result;
    //return checkCycleCovariate(cycleFactor+index*increment, maxCycle);
}

void calculateSkipSlice(char* seq, SpanField& sequencePtrs, char* qual, SpanField qualPtrs, CigarOperation* cigar, SpanField cigar_ptr, int32_t POS, int32_t sftStart, int32_t sftEnd, Interval* knownSites, SpanFieldLong& knownSitesPtrs, uint8_t* skipSlice) {
    size_t seqLength = sequencePtrs.len;
    if (knownSitesPtrs.len == 0) return;
    // [[0, 5] [7, 10] [12,25]]
    auto intervals = Intersect(knownSites, knownSitesPtrs, sftStart, sftEnd);

    for (size_t interval_idx = intervals.first; interval_idx < intervals.second; interval_idx++) {
        auto site = knownSites[interval_idx];
        auto t = getReadCoordinateForReferenceCoordinate(cigar, cigar_ptr, sftStart, static_cast<int32_t>(site.Start));
        int32_t featureStartOnRead = (t.first >= 0 ? t.first : 0);
        CigarOperation op = t.second;
        if (op.Operation == 'D') {
            featureStartOnRead--;
        }

        t = getReadCoordinateForReferenceCoordinate(cigar, cigar_ptr, sftStart, int32_t(site.End));
        auto featureEndOnRead = (t.first>=0 ? t.first : seqLength);

        if( featureStartOnRead > seqLength ) {
            featureStartOnRead = featureEndOnRead = seqLength;
        }
        for (size_t i = sycl::max(0, featureStartOnRead); i < sycl::min(seqLength, featureEndOnRead+1); i++) {
            skipSlice[i] = 1;
        }
    }
}



void BaseRecalibrator::close() {
}



// Recalibrate implements the first step of base recalibration.
//[[clang::optnone]]
void BaseRecalibrator::Recalibrate(SamAlignments& reads) {
    Timer timerAll;
    if (reads.GetAlnsNum() == 0) return;
    auto& profiler = ProfilerBQSR::instance();
    timerAll.reset();
    Timer t;

    t.reset();

    int32_t maxIdx = reads.GetAlnsNum();
    auto& hdr = reads.GetHdr();

    profiler.start(StageId::AddReadGroup);

    //reads.TMP.reads_group_covariate_.resize(maxIdx);
    std::vector<int32_t> reads_group_covariate_(maxIdx);
    tbb::parallel_for(0, maxIdx, [&](size_t idx) {
        reads_group_covariate_[idx] = readGroupCovariateId(hdr, reads, idx);
    });

    profiler.stop(StageId::AddReadGroup);

    auto *clippedSequences = reads.SEQ.Data();
    auto *clippedSequencesPtrs = reads.SEQ.Ptrs().data();
    auto *clippedQualities = reads.QUAL.Data();
    auto *clippedQualitiesPtrs = reads.QUAL.Ptrs().data();

    std::vector<uint8_t> recalibrate_aln_flag(reads.GetAlnsNum());

    Logger::getLogger().log(LOG_INFO, "Time to declare and initialize variables: " + std::to_string(t.elapsed()));

    t.reset();

    tbb::parallel_for(0UL, reads.GetAlnsNum(), [&] (size_t idx) {
        recalibrate_aln_flag[idx] = recalibrateAln(hdr, idx, reads.MAPQ[idx], reads.FLAG[idx], reads.RNAME[idx], reads.PNEXT[idx], reads.POS[idx], reads.GetSEQ(idx).size(), reads.GetQUAL(idx).size(), reads.CIGAR.Data(), reads.CIGAR.Ptrs()[idx], reads.TAGS.RG[idx]);
    });

    Logger::getLogger().log(LOG_INFO, "Time to set recalibrate flag: " + std::to_string(t.elapsed()));

    t.reset();

    size_t batch_size = reads.GetAlnsNum();

    profiler.start(StageId::RecalibrateBatches);

    t.reset();
    auto* cigar = reads.CIGAR.Data();
    auto* cigar_ptr = reads.CIGAR.Ptrs().data();
    BaseRecalibratorTablesAccessor recalibration_tables_acc = recalibration_tables.make_accessor();
    SequencesInput in(local_known_intervals, local_known_intervals_ptrs,
                      fasta, fasta_ptrs,
                      tables_quality.data(), tables_cycle.data(), tables_context.data(), tables_quality.size(), tables_cycle.size(), tables_context.size(),
                      recalibration_tables_acc,
                      clippedSequences,
                      clippedSequencesPtrs,
                      clippedQualities,
                      clippedQualitiesPtrs,
                      cigar,
                      cigar_ptr,
                      reads.POS.data(),
                      reads.RNEXT.data(), reads.PNEXT.data(), reads.TLEN.data(),
                      reads.RNAME.data(),
                      recalibrate_aln_flag.data(),
                      reads.FLAG.data(),
                      reads_group_covariate_.data(),
                      batch_size);


    TablesOutput out; // Not actullay used, but here for compatibility with my API


    // Set up your configuration using available device memory
    /*ComputeSkipSliceConfiguration config;
    //size_t available_mem = ConfigDeviceMemoryUsage::max_memory;//32*std::pow(2,30); //q.get_device().get_info<sycl::info::device::global_mem_size>() - 2 * pow(2, 30);
    size_t available_mem = q_.get_device().get_info<sycl::info::device::global_mem_size>() - constant_allocation_;
    config.max_memory = available_mem / 2;
*/
    // Process the current batch with your kernel
    batch_processing_variable_length<ComputeSkipSliceKernelData, ComputeSkipSliceConfiguration, SequencesInput, TablesOutput>(q_process_, q_transfer_, in, out, compute_skip_slice_kernel, config, batch_process, batch_transfer);
    q_.wait_and_throw();
    Logger::getLogger().log(LOG_INFO, "Time to batch process alignments: " + std::to_string(t.elapsed()));

    profiler.stop(StageId::RecalibrateBatches);

    Logger::getLogger().log(LOG_INFO, "Total time recalibrate: " + std::to_string(timerAll.elapsed()) + "sec");
    profiler.report();
}

BaseRecalibrator::BaseRecalibrator(MappedFasta &referenceGenome, KnownSites& knownSites,
                                   sycl::queue &queue, Header &hdr, BaseRecalibratorTables& recalibration_tables) : recalibration_tables(recalibration_tables), q_(queue), q_process_(queue.get_context(), queue.get_device(), sycl::async_handler{}), q_transfer_(queue.get_context(), queue.get_device(), sycl::async_handler{}) {
    Timer t;
    t.reset();

    this->fasta = mem_t<uint8_t>(referenceGenome.getFullFasta().size(), q_, memory_space_t::memory_space_device);
    this->fasta_ptrs = mem_t<SpanFieldLong>(referenceGenome.getFullFastaPtrs().size(), q_, memory_space_t::memory_space_device);

    this->local_known_intervals = mem_t<Interval>(knownSites.getIntervals().size(), q_, memory_space_t::memory_space_device);
    this->local_known_intervals_ptrs = mem_t<SpanFieldLong>(knownSites.getIntervalsPtrs().size(), q_, memory_space_t::memory_space_device);

    // Copy FASTA bases
    q_.memcpy(this->fasta.data(),
             referenceGenome.getFullFasta().data(),
             referenceGenome.getFullFasta().size() * sizeof(uint8_t));

    // Copy FASTA pointers
    q_.memcpy(this->fasta_ptrs.data(),
             referenceGenome.getFullFastaPtrs().data(),
             referenceGenome.getFullFastaPtrs().size() * sizeof(SpanFieldLong));

    // Copy known-site intervals
    q_.memcpy(this->local_known_intervals.data(),
             knownSites.getIntervals().data(),
             knownSites.getIntervals().size() * sizeof(Interval));

    // Copy interval pointers
    q_.memcpy(this->local_known_intervals_ptrs.data(),
             knownSites.getIntervalsPtrs().data(),
             knownSites.getIntervalsPtrs().size() * sizeof(SpanFieldLong));

    Logger::getLogger().log(LOG_INFO, "Time to flatten fasta and intervals: " + std::to_string(t.elapsed()));
    this->tables_quality = mem_t<bqsrEntry>(recalibration_tables.qualityScores().size(), q_, memory_space_t::memory_space_device);
    this->tables_context = mem_t<bqsrEntry>(recalibration_tables.contexts().size(), q_, memory_space_t::memory_space_device);
    this->tables_cycle = mem_t<bqsrEntry>(recalibration_tables.cycles().size(), q_, memory_space_t::memory_space_device);


    bqsrEntry defaultEntry;
    defaultEntry.Observations = 0;
    defaultEntry.Mismatches = 0;

    q_.fill(tables_quality.data(), defaultEntry, recalibration_tables.qualityScores().size());
    q_.fill(tables_cycle.data(), defaultEntry, recalibration_tables.cycles().size());
    q_.fill(tables_context.data(), defaultEntry, recalibration_tables.contexts().size());
    q_.wait();
    constant_allocation_ =
        referenceGenome.getFullFasta().size()                       * sizeof(uint8_t)       +
        referenceGenome.getFullFastaPtrs().size()                   * sizeof(SpanFieldLong) +
        knownSites.getIntervals().size()                            * sizeof(Interval)      +
        knownSites.getIntervalsPtrs().size()                        * sizeof(SpanFieldLong) +
        recalibration_tables.qualityScores().size()                 * sizeof(bqsrEntry)     +
        recalibration_tables.contexts().size()                      * sizeof(bqsrEntry)     +
        recalibration_tables.cycles().size()                        * sizeof(bqsrEntry)     +
            515 * 1024 * 1024; // 512 MB as margin

    Logger::getLogger().log(LOG_INFO,
        "Total device memory allocated: " + std::to_string(constant_allocation_) +
        " bytes (" + std::to_string(constant_allocation_ / (1024.0 * 1024.0 * 1024.0)) + " GB)");

    Logger::getLogger().log(LOG_INFO, "Time to allocate and initialize device memory: " + std::to_string(t.elapsed()));
}

void BaseRecalibrator::initializeDeviceMemory() {
    size_t available_mem = q_.get_device().get_info<sycl::info::device::global_mem_size>() - constant_allocation_;
    config.max_memory = available_mem / 3;

    batch_process = std::make_shared<Batch<ComputeSkipSliceKernelData, ComputeSkipSliceConfiguration>>(q_transfer_, q_process_, config);
    batch_transfer = std::make_shared<Batch<ComputeSkipSliceKernelData, ComputeSkipSliceConfiguration>>(q_transfer_, q_process_, config);
}

//[[clang::optnone]]
void BaseRecalibrator::exportTables(BaseRecalibratorTables &recalibration_tables) {
    tables_quality.context().memcpy(recalibration_tables.qualityScores().table.data(), tables_quality.data(), sizeof(bqsrEntry)*recalibration_tables.qualityScores().size());
    tables_cycle.context().memcpy(recalibration_tables.cycles().table.data(), tables_cycle.data(), sizeof(bqsrEntry)*recalibration_tables.cycles().size());
    tables_context.context().memcpy(recalibration_tables.contexts().table.data(), tables_context.data(), sizeof(bqsrEntry)*recalibration_tables.contexts().size());
}

void BaseRecalibrator::mergeTables(BaseRecalibratorTables &recalibration_tables) {
    auto merge_one_table = [&](auto &dstTableVec,
                               auto &srcMem /* mem_t<bqsrEntry>& */)
    {
        const size_t n = dstTableVec.size();
        if (n == 0) return;
        // Safety: assume srcMem.size() == n. If paranoid, you can assert.
        // assert(srcMem.size() == n && "BQSR table size mismatch across devices");

        const bqsrEntry* srcHostPtr = nullptr;
        std::vector<bqsrEntry> staging; // only used if we need to copy from device

        if (srcMem.space() == memory_space_t::memory_space_shared ||
            srcMem.space() == memory_space_t::memory_space_host) {
            // Host can directly read srcMem.data()
            srcHostPtr = srcMem.data();
            } else {
                // Device-only: copy to staging buffer on host
                staging.resize(n);
                // enqueue memcpy and wait
                srcMem.context().memcpy(staging.data(),
                                      srcMem.data(),
                                      n * sizeof(bqsrEntry)).wait();
                srcHostPtr = staging.data();
            }

        // Now merge into dst
        bqsrEntry* dstPtr = dstTableVec.data();
        for (size_t i = 0; i < n; ++i) {
            dstPtr[i].Observations += srcHostPtr[i].Observations;
            dstPtr[i].Mismatches   += srcHostPtr[i].Mismatches;
        }
    };

    // qualityScores()
    merge_one_table(recalibration_tables.qualityScores().table, tables_quality);

    // cycles()
    merge_one_table(recalibration_tables.cycles().table, tables_cycle);

    // contexts()
    merge_one_table(recalibration_tables.contexts().table, tables_context);
}

double estimateHierarchicalBayesianQuality(double epsilon, combinedBqsrEntry& empiricalReadGroupEntry, bqsrEntry& empiricalBqEntry, bqsrEntry& empiricalContextEntry, bqsrEntry& empiricalCycleEntry) {
    double deltaGlobalQuality = 0.0;
    double deltaReportedQuality = 0.0;

    deltaGlobalQuality = static_cast<double>(empiricalReadGroupEntry.entry.calculateEmpiricalQuality(epsilon)) - epsilon;

    deltaReportedQuality = static_cast<double>(empiricalBqEntry.calculateEmpiricalQuality(deltaGlobalQuality+epsilon)) - deltaGlobalQuality - epsilon;

    double deltaCovariatesQuality = 0.0;
    double conditionalPrior = deltaReportedQuality + deltaGlobalQuality + epsilon;

    deltaCovariatesQuality += static_cast<double>(empiricalCycleEntry.calculateEmpiricalQuality(conditionalPrior)) - conditionalPrior;

    if (empiricalContextEntry.Observations > 0) {
        deltaCovariatesQuality += static_cast<double>(empiricalContextEntry.calculateEmpiricalQuality(conditionalPrior)) - conditionalPrior;
    }
    return conditionalPrior + deltaCovariatesQuality;
}
