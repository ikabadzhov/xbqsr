#ifndef GENOMICPIPELINE_BQSR_H
#define GENOMICPIPELINE_BQSR_H

#include <string>
#include <memory>
#include "../sam_data_structures/SamHeader.h"
#include "../file_format/MappedFasta.h"
#include "../intervals/Interval.h"
#include "TablesBQSR.h"
#include "../sam_data_structures/SamAlignments.h"
#include "../sycl_parallel_library/memory/context.h"


const int32_t INVALID_INDEX = INT32_MAX - 1;

const uint8_t minInterestingQual = 6;

const int quantizationLevel = 16;
const std::string countString = "Count";
const std::string covariateNameString = "CovariateName";
const std::string covariateValueString = "CovariateValue";
const std::string empiricalQualityString = "EmpiricalQuality";
const std::string errorsString = "Errors";
const std::string estimatedQReportedString = "EstimatedQReported";
const std::string eventTypeString = "EventType";
const std::string observationsString = "Observations";
const std::string qualityScoreString = "QualityScore";
const std::string quantizedScoreString = "QuantizedScore";
const std::string readGroupString = "ReadGroup";

const int32_t MAX_CYCLE = 500;
const std::string BQSRTablenamePrefix = "GATK"; // To be defined the actual prefix


struct covariateBqsrTableKey {
    bool cycleNotContext;
    std::string text;
    bqsrTableKey key;
};

/* BaseRecalibratorTables is the result of the base recalibration.
*  All subsequent steps, including ApplyBQSR, are based on these tables.
*/
class BaseRecalibratorTables {
public:
    bqsrTable QualityScores;
    bqsrTable Cycles;
    bqsrTable Contexts;
    size_t szQualScore;
    size_t szCycleCov;
    size_t szContextCov;
    size_t nReadGroups;
    size_t MAX_CYCLE;

    BaseRecalibratorTables(size_t nReadGroups, size_t szQual, size_t maxCycle);
    SYCL_EXTERNAL static ssize_t getQualityScoresIdx(bqsrTableKey &key, size_t szQualScore);
    SYCL_EXTERNAL static ssize_t getCycleCovIdx(bqsrTableKey &key, size_t szQualScore, size_t szCycleCov);
    SYCL_EXTERNAL static ssize_t getContextCovIdx(bqsrTableKey &key, size_t szQualScore, size_t szCycleCov, size_t szContextCov);

    bqsrTableKey idxToQualityScoresKey(uint8_t idx);
    bqsrTableKey idxToCycleCovKey(size_t idx);
    bqsrTableKey idxToContextCovKey(size_t idx);

    void PrintTables(Header& samHeader);

    void initializeCombinedBQSRTable(bqsrTable& bqsrTable, std::vector<combinedBqsrEntry>& table);
    void FinalizeBQSRTables();
    void ApplyBQSR(sycl::device& device, SamAlignments &reads, int32_t quantizeLevels, std::vector<uint8_t> &sqqList,
                   int32_t maxCycle);

    void
    PrintBQSRTables(const std::string &filename, CompressedStringCollection &RG, Header &hdr);

private:
    void printOtherCovariateTable(std::ostream& file, Header &hdr);
    void printBQSRTable(std::ofstream& filename, Header &hdr);
    void printQuantizationTable(std::ofstream &file);
    void printCombinedBQSRTable(std::ostream& file, CompressedStringCollection& RG, Header &hdr);
};

uint32_t MaxContextCovariateKey();


/*
 * BaseRecalibrator implements the first step of base recalibration.
 */
class BaseRecalibrator {
public:
    BaseRecalibratorTables& recalibration_tables;
    MappedFasta& fastaMap;
    std::unordered_map<std::string, std::vector<Interval>> knownIntervals;
    sycl::queue q;
    mem_t<Interval> local_known_intervals;
    mem_t<SpanFieldLong> local_known_intervals_ptrs;
    mem_t<SpanFieldLong> fasta_ptrs;
    mem_t<uint8_t> fasta;

    mem_t<bqsrEntry> tables_quality;//(result.QualityScores.size(), q, memory_space_t::memory_space_shared);
    mem_t<bqsrEntry> tables_cycle;//(result.Cycles.size(), q, memory_space_t::memory_space_shared);
    mem_t<bqsrEntry> tables_context;//(result.Contexts.size(), q, memory_space_t::memory_space_shared);

    explicit BaseRecalibrator(MappedFasta &fastaMap, std::vector<std::string> &knownSites,
                              sycl::device &device, Header &hdr, BaseRecalibratorTables& recalibration_tables);// : fastaMap(fastaMap) {};
    void getTables(BaseRecalibratorTables& recalibration_tables);
    void Recalibrate(sycl::device& device, SamAlignments& reads, int32_t maxCycle);
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

    applyBuffers(size_t sizeReadGroup, size_t sizeQuality, size_t maxCycle, size_t sizeContext) : szQual(sizeQuality), szCycle(maxCycle*2), szContext(sizeContext+1), nReadGroups(sizeReadGroup), MAX_CYCLE(maxCycle) {
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

std::unique_ptr<BaseRecalibrator> NewBaseRecalibrator(std::vector<std::string>& knownSites, MappedFasta& referenceFasta);

SYCL_EXTERNAL void calculateSkipSlice(char* seq, SpanField& sequencePtrs, char* qual, SpanField qualPtrs, CigarOperation* cigar, SpanField cigar_ptr, int32_t POS, int32_t sftStart, int32_t sftEnd, Interval* knownSites, SpanFieldLong& knownSitesPtrs, uint8_t* skipSlice);
SYCL_EXTERNAL void computeStrandedClippedSeq(char* SEQ, SpanField seqPtrs, char* QUAL, SpanField& qualPtrs, bool isReversed, char* newSeq, SpanField& newSeqPtrs);
SYCL_EXTERNAL void computeBaseContextCovariate(bool isReversed, char* strandedClippedSeq, SpanField& strandedClippedSeqPtrs, int32_t* keys, uint16_t& keyLength);
SYCL_EXTERNAL void contextWith(char* bases, SpanField basesPtrs, int32_t contextSize, int32_t* keys, uint16_t& keysLength);
SYCL_EXTERNAL CycleCovariatesResult prepareCycleCovariates(uint16_t sequence_len, uint16_t FLAG);
SYCL_EXTERNAL int32_t computeBaseCycleCovariate( int32_t cycleFactor, int32_t increment, uint32_t index, int32_t maxCycle );
SYCL_EXTERNAL void computeSnpEvents(char* sequence, SpanField sequencePtrs, int32_t POS, CigarOperation* cigar, size_t cigar_start, uint16_t cigar_len, uint8_t* ref, SpanFieldLong& refPtrs, int8_t* snps);
SYCL_EXTERNAL double estimateHierarchicalBayesianQuality(double epsilon, combinedBqsrEntry& empiricalReadGroupEntry, bqsrEntry& empiricalBqEntry, bqsrEntry& empiricalContextEntry, bqsrEntry& empiricalCycleEntry);
SYCL_EXTERNAL int32_t fastRound(double d);

#endif //GENOMICPIPELINE_BQSR_H
