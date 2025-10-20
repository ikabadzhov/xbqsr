#pragma once
#include <sycl/device.hpp>
#include <sys/types.h>

#include "BaseRecalibratorTablesAccessor.h"
#include "TablesBQSR.h"
#include "../sam_data_structures/SamAlignments.h"
#include "../sam_data_structures/SamHeader.h"
#include "../sycl_parallel_library/memory/memory.h"

constexpr int32_t INVALID_INDEX = INT32_MAX - 1;

constexpr uint8_t MIN_INTERESTING_QUAL = 6;
constexpr uint32_t MAX_STATIC_QUANTIZE_QUALITY = 254;
constexpr uint32_t LENGTH_BITS = 4;
constexpr uint32_t MISMATCHES_CONTEXT_SIZE = 2;
constexpr uint32_t QUANTIZATION_LEVEL = 16;

namespace bqsr::defaults {
    inline constexpr size_t MAX_QUALITY_SCORE = 94;      // 0–93
    inline constexpr size_t MAX_CYCLE        = 500;
    inline constexpr size_t CONTEXT_LEN = 2;   // 4^4 = 256
}

const std::string BQSRTablenamePrefix = "GATK"; // To be defined the actual prefix
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


struct covariateBqsrTableKey {
    bool cycleNotContext;
    std::string text;
    bqsrTableKey key;
};

struct BQSRTableConfig {
    size_t n_read_groups;
    size_t n_qual_score;
    size_t max_cycle;
    size_t context_len;
};

class BaseRecalibratorTables {
public:
    explicit BaseRecalibratorTables(const BQSRTableConfig& cfg);

    [[nodiscard]] BaseRecalibratorTablesAccessor make_accessor() const noexcept {
        return { nReadGroups, nQualScore, maxCycle, contextLen, nCycleCov, nContextCov };
    }
    // --- Getters (const + mutable) ---
    [[nodiscard]] inline const bqsrTable& qualityScores() const noexcept { return QualityScores_; }
    [[nodiscard]] inline bqsrTable&       qualityScores()       noexcept { return QualityScores_; }

    [[nodiscard]] inline const bqsrTable& cycles() const noexcept { return Cycles_; }
    [[nodiscard]] inline bqsrTable&       cycles()       noexcept { return Cycles_; }

    [[nodiscard]] inline const bqsrTable& contexts() const noexcept { return Contexts_; }
    [[nodiscard]] inline bqsrTable&       contexts()       noexcept { return Contexts_; }

    // --- helper: total element counts ---
    [[nodiscard]] inline size_t qualityScoresCount() const noexcept { return nReadGroups * nQualScore; }
    [[nodiscard]] inline size_t cyclesCount() const noexcept        { return nReadGroups * nQualScore * nCycleCov; }
    [[nodiscard]] inline size_t contextsCount() const noexcept      { return nReadGroups * nQualScore * nContextCov; }

    static constexpr size_t MaxContextCovariateKey(size_t contextLen) noexcept {
        // GATK-style encoding:
        // lowest LENGTH_BITS bits store length; each base adds 2 bits.
        uint32_t length = contextLen;  // the length of the context
        uint32_t key = contextLen;
        uint32_t bitOffset = LENGTH_BITS;
        for (int i = 0; i <length ; i++) {
            key |= (3 << bitOffset);
            bitOffset += 2;
        }
        return key;
    }

    void initializeCombinedBQSRTable(bqsrTable& bqsrTable, std::vector<combinedBqsrEntry>& table);
    void FinalizeBQSRTables();
    void ApplyBQSR(sycl::device& device, SamAlignments &reads, int32_t quantizeLevels, std::vector<uint8_t> &sqqList);

    void PrintBQSRTables(const std::string &filename, std::vector<int32_t>& RG, Header &hdr);
    void LoadFromReport(const std::string &filename, Header &hdr);

private:
    void printOtherCovariateTable(std::ostream& file, Header &hdr);


    void printBQSRTable(std::ofstream& filename, Header &hdr);
    void printQuantizationTable(std::ofstream &file);
    void printCombinedBQSRTable(std::ostream &file, Header &hdr);

private:
    // --- configuration (runtime-const) ---
    const size_t nReadGroups;   // RG
    const size_t nQualScore;    // Q
    const size_t maxCycle;      // M
    const size_t contextLen;    // L

    // --- derived sizes ---
    const size_t nCycleCov;     // 2 * maxCycle (forward/backward)
    const size_t nContextCov;   // 4^contextLen

    // tables (unchanged)
    bqsrTable QualityScores_;    // size = RG * Q
    bqsrTable Cycles_;           // size = RG * Q * nCycleCov
    bqsrTable Contexts_;         // size = RG * Q * nContextCov

    std::vector<combinedBqsrEntry> combinedBQSRTable;
    std::vector<uint8_t> staticQuantizedScores;
    std::vector<int64_t> quantizationMap;
    std::vector<uint8_t> quantizedQualityScores;
    mem_t<int> recalibrated_quality_scores;
    bool is_apply_bqsr_initialized = false;
};


void initializeQuantizedQualityScores(bqsrTable& bqsrTable, int32_t quantLevel, std::vector<int64_t>& quantizationMap, std::vector<uint8_t>& quantizedScores);
std::vector<uint8_t> initializeStaticQuantizedScores(std::vector<uint8_t>& quals);

inline char simpleBaseIndexToBase(int32_t baseIndex) {
    switch( baseIndex ) {
        case -1:
            return 'N';
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
        default:
            std::cout<<"invalid baseIndex"<<std::endl;
            exit(-1);
    }
}
