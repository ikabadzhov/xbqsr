#ifndef GENOMICPIPELINE_BQSR_H
#define GENOMICPIPELINE_BQSR_H

#include <string>
#include <memory>

#include "BaseRecalibratorTables.h"
#include "../sam_data_structures/SamHeader.h"
#include "../file_format/MappedFasta.h"
#include "../intervals/Interval.h"
#include "TablesBQSR.h"
#include "../sam_data_structures/SamAlignments.h"
#include "../sycl_parallel_library/memory/context.h"
#include "compute_vectors_for_covariates/ComputeCovariateVectorsData.h"

/*
 * BaseRecalibrator implements the first step of base recalibration.
 */
class BaseRecalibrator {
public:
    BaseRecalibratorTables& recalibration_tables;
    sycl::queue q_;
    sycl::queue q_process_;
    sycl::queue q_transfer_;
    mem_t<Interval> local_known_intervals;
    mem_t<SpanFieldLong> local_known_intervals_ptrs;
    mem_t<SpanFieldLong> fasta_ptrs;
    mem_t<uint8_t> fasta;

    ComputeSkipSliceConfiguration config;
    std::shared_ptr<Batch<ComputeSkipSliceKernelData, ComputeSkipSliceConfiguration>> batch_process;
    std::shared_ptr<Batch<ComputeSkipSliceKernelData, ComputeSkipSliceConfiguration>> batch_transfer;

    mem_t<bqsrEntry> tables_quality;
    mem_t<bqsrEntry> tables_cycle;
    mem_t<bqsrEntry> tables_context;

    size_t constant_allocation_;

    explicit BaseRecalibrator(MappedFasta &referenceGenome, KnownSites &knownSites,
                              sycl::queue &queue, Header &hdr, BaseRecalibratorTables& recalibration_tables);// : fastaMap(fastaMap) {};

    void initializeDeviceMemory();

    void exportTables(BaseRecalibratorTables& recalibration_tables);
    void mergeTables(BaseRecalibratorTables& recalibration_tables);
    void Recalibrate(SamAlignments& reads);
    void close();
};

struct DeviceBqsrEntry {
    int64_t Observations;
    int64_t Mismatches;
};

struct applyKey {
    int32_t readGroupCovariate;
    uint8_t qual;
    int32_t cycleCovariate;
    int32_t contextCovariate;

    bool operator==(const applyKey &rhs) const {
        return readGroupCovariate == rhs.readGroupCovariate &&
               qual == rhs.qual &&
               cycleCovariate == rhs.cycleCovariate &&
               contextCovariate == rhs.contextCovariate;
    };
};

struct applyBuffers {
    std::string strandedClippedSeq;
    std::vector<int32_t> baseContextCovariate;
    /* recalibratedQualityScores vector must allow to have context covariate as -1. In order to support that, every value is shifted by 1*/
    std::vector<int8_t> recalibratedQualityScores; // Max possible value 93. So, we can use int8_t as type and use -1 as null value
    size_t nReadGroups;
    size_t szQual;
    size_t szCycle;
    size_t szContext;
    size_t MAX_CYCLE;

    applyBuffers(size_t sizeReadGroup, size_t sizeQuality, size_t maxCycle, size_t sizeContext) : nReadGroups(sizeReadGroup), szQual(sizeQuality), szCycle(maxCycle*2), szContext(sizeContext+1), MAX_CYCLE(maxCycle) {
        recalibratedQualityScores.resize(nReadGroups * szQual * szCycle * szContext, -1);
    }

    size_t getIdxRecalQualScore(applyKey& key) {
        if(key.cycleCovariate < 0) {
            return key.readGroupCovariate * szQual * szContext * szCycle + key.qual * szContext * szCycle +
                    (key.contextCovariate + 1) * szCycle + (-1 * key.cycleCovariate + MAX_CYCLE - 1);
        } else {
            return key.readGroupCovariate * szQual * szContext * szCycle + key.qual * szContext * szCycle +
                    (key.contextCovariate + 1) * szCycle + key.cycleCovariate;
        }
    }

};

struct CycleCovariatesResult{
    int32_t cycleFactor;
    int32_t increment;
};

SYCL_EXTERNAL void calculateSkipSlice(char* seq, SpanField& sequencePtrs, char* qual, SpanField qualPtrs, CigarOperation* cigar, SpanField cigar_ptr, int32_t POS, int32_t sftStart, int32_t sftEnd, Interval* knownSites, SpanFieldLong& knownSitesPtrs, uint8_t* skipSlice);
SYCL_EXTERNAL void computeStrandedClippedSeq(char* SEQ, SpanField seqPtrs, char* QUAL, SpanField& qualPtrs, bool isReversed, char* newSeq, SpanField& newSeqPtrs);
SYCL_EXTERNAL void computeBaseContextCovariate(bool isReversed, char* strandedClippedSeq, SpanField& strandedClippedSeqPtrs, int32_t* keys, uint16_t& keyLength);
SYCL_EXTERNAL void contextWith(char* bases, SpanField basesPtrs, int32_t contextSize, int32_t* keys, uint16_t& keysLength);
SYCL_EXTERNAL CycleCovariatesResult prepareCycleCovariates(uint16_t sequence_len, uint16_t FLAG);
SYCL_EXTERNAL int32_t computeBaseCycleCovariate( int32_t cycleFactor, int32_t increment, uint32_t index, int32_t maxCycle );
SYCL_EXTERNAL void computeSnpEvents(char* sequence, SpanField sequencePtrs, int32_t POS, CigarOperation* cigar, size_t cigar_start, uint16_t cigar_len, uint8_t* ref, SpanFieldLong& refPtrs, int8_t* snps);

// --- Inlined Functions for GPU visibility ---

SYCL_EXTERNAL inline int32_t fastRound(float d) {
    return (d > 0.0f) ? static_cast<int>(d + 0.5f) : static_cast<int>(d - 0.5f);
}

SYCL_EXTERNAL inline float estimateHierarchicalBayesianQuality(float epsilon, combinedBqsrEntry& empiricalReadGroupEntry, bqsrEntry& empiricalBqEntry, bqsrEntry& empiricalContextEntry, bqsrEntry& empiricalCycleEntry) {
    float deltaGlobalQuality = 0.0f;
    float deltaReportedQuality = 0.0f;

    deltaGlobalQuality = static_cast<float>(empiricalReadGroupEntry.entry.calculateEmpiricalQuality(epsilon)) - epsilon;

    deltaReportedQuality = static_cast<float>(empiricalBqEntry.calculateEmpiricalQuality(deltaGlobalQuality+epsilon)) - deltaGlobalQuality - epsilon;

    float deltaCovariatesQuality = 0.0f;
    float conditionalPrior = deltaReportedQuality + deltaGlobalQuality + epsilon;

    deltaCovariatesQuality += static_cast<float>(empiricalCycleEntry.calculateEmpiricalQuality(conditionalPrior)) - conditionalPrior;

    if (empiricalContextEntry.Observations > 0) {
        deltaCovariatesQuality += static_cast<float>(empiricalContextEntry.calculateEmpiricalQuality(conditionalPrior)) - conditionalPrior;
    }
    return conditionalPrior + deltaCovariatesQuality;
}

#endif //GENOMICPIPELINE_BQSR_H
