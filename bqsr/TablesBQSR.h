#ifndef GENOMICPIPELINE_TABLESBQSR_H
#define GENOMICPIPELINE_TABLESBQSR_H

#include <cstdint>
#include <string>
#include <cmath>
#include <tbb/concurrent_unordered_map.h>
#include <memory>
#include <algorithm>
#include <sycl/sycl.hpp>

const float Ln10 = 2.30258509299404568401799145468436420760110148862877297603332790f;
const float Log10E = 1.0f / Ln10;

const int32_t smoothingConstant = 1;
const int32_t maxRecalibratedQualScore = 93;
const int32_t maxNumberOfObservations = INT32_MAX - 1;

const int32_t maxQualityScore = 93;
const int32_t maxReasonableQualScore = 60;

const int32_t LOG10_QUAL_EMPIRICAL_PRIOR_CACHE_SIZE = 21;
static const float LOG10_QUAL_EMPIRICAL_PRIOR_CACHE[LOG10_QUAL_EMPIRICAL_PRIOR_CACHE_SIZE] = {
            -0.045757490560675115f,
            -0.9143464543671788f,
            -3.5201133457866898f,
            -7.863058164819208f,
            -13.943180911464733f,
            -21.760481585723266f,
            -31.314960187594806f,
            -42.606616717079355f,
            -55.63545117417691f,
            -70.40146355888747f,
            -86.90465387121104f,
            -105.14502211114761f,
            -125.1225682786972f,
            -146.83729237385978f,
            -170.2891943966354f,
            -195.47827434702398f,
            -222.4045322250256f,
            -251.06796803064023f,
            -281.46858176386786f,
            -313.60637342472336f,
            -3.402823466e+38f // Roughly -FLT_MAX
};

struct bqsrTableKey {
    uint8_t Qual;
    int32_t Covariate;
    int32_t ReadGroup;

    bqsrTableKey(uint8_t qual, int32_t covariate, int32_t readGroup) : Qual(qual), Covariate(covariate), ReadGroup(readGroup) {};
    bqsrTableKey() : Qual(0), Covariate(0), ReadGroup(0) {};

    bool operator==(const bqsrTableKey &rhs) const;
    bool operator!=(const bqsrTableKey &rhs) const;
};

struct bqsrTableKeyHasher
{
    std::size_t operator()(const bqsrTableKey& key) const
    {
        return std::hash<uint8_t>()(key.Qual) ^ std::hash<int32_t>()(key.Covariate) ^ std::hash<int32_t>()(key.ReadGroup);
    }
};

struct bqsrEntry {
    int64_t Observations{};
    int64_t Mismatches{};

    bqsrEntry() = default;
    bqsrEntry(const bqsrEntry&) = default;
    bqsrEntry& operator=(const bqsrEntry&) = default;

    explicit bqsrEntry(int64_t observ, int64_t mis) : Observations(observ), Mismatches(mis) {}

    SYCL_EXTERNAL uint8_t calculateEmpiricalQuality(float priorQual);
};

static_assert(std::is_trivially_copyable_v<bqsrEntry>,
              "bqsrEntry must be trivially copyable for SYCL device use");


class bqsrTable {
public:
    explicit bqsrTable(size_t sz) : table(sz), empirical_qualities(sz) {};
    bqsrTable() = default;
    SYCL_EXTERNAL static void update(bqsrEntry* table, size_t idx, int64_t errorVal);
    void merge(bqsrTable &other);

    std::vector<bqsrEntry> table;
    std::vector<uint8_t> empirical_qualities;
    [[nodiscard]] size_t size() const noexcept { return table.size(); }
};

struct combinedBqsrEntry {
    float reportedQuality{};
    bqsrEntry entry;
    uint8_t empiricalQuality;
    bool isValid = false;
};

// --- Math Functions (Inlined for GPU visibility) ---

SYCL_EXTERNAL inline float qualityToErrorProbability(float phred) {
    return sycl::pow(10.0f, phred / -10.0f);
}

SYCL_EXTERNAL inline float qualityToProbability(float phred) {
    return 1.0f - sycl::pow(10.0f, phred / -10.0f);
}

SYCL_EXTERNAL inline float calculateExpectedErrors(int64_t observations, float reportedQuality) {
    return static_cast<float>(observations) * qualityToErrorProbability(reportedQuality);
}

SYCL_EXTERNAL inline float log10QualEmpiricalPrior(float empiricalQual, float reportedQual) {
    auto diff = std::min(std::abs(static_cast<int32_t>(empiricalQual - reportedQual)), 
                         static_cast<int32_t>(LOG10_QUAL_EMPIRICAL_PRIOR_CACHE_SIZE - 1));
    return LOG10_QUAL_EMPIRICAL_PRIOR_CACHE[diff];
}

SYCL_EXTERNAL inline float log10Gamma(int64_t n) {
#if defined(__ADAPTIVECPP__)
    // Workaround for AdaptiveCpp builtin_interface.hpp where __acpp_lgamma errantly
    // dispatches to __acpp_acos, producing NaN for arguments > 1.
    int sign = 0;
    return sycl::lgamma_r(static_cast<float>(n), sycl::private_ptr<int>(&sign)) * Log10E;
#else
    return sycl::lgamma(static_cast<float>(n)) * Log10E;
#endif
}

SYCL_EXTERNAL inline float log10BinomialCoefficient(int64_t n, int64_t k) {
    return log10Gamma(n + 1) - log10Gamma(k + 1) - log10Gamma(n - k + 1);
}

SYCL_EXTERNAL inline float log10BinomialProbability(int64_t n, int64_t k, float log10p) {
    if (log10p == 0.0f) {
        return -3.402823466e+38f; // -FLT_MAX
    }
    float log10MinP = sycl::log10(1.0f - sycl::pow(10.0f, log10p));
    return log10BinomialCoefficient(n, k) + log10p * static_cast<float>(k) + log10MinP * static_cast<float>(n - k);
}

SYCL_EXTERNAL inline float log10QualEmpiricalLikelihood(float empiricalQual, float observations, int64_t mismatches) {
    if (observations == 0.0f) {
        return 0.0f;
    }
    float qualToErrorProbLog10 = empiricalQual / -10.0f;
    return log10BinomialProbability(static_cast<int64_t>(observations), mismatches, qualToErrorProbLog10);
}

SYCL_EXTERNAL inline uint8_t calculateBayesianEstimateOfEmpiricalQuality(int64_t observations, int64_t mismatches, float priorQuality) {
    if (observations > maxNumberOfObservations) {
        mismatches = static_cast<int64_t>(sycl::round(static_cast<float>(mismatches) * (static_cast<float>(maxNumberOfObservations) / static_cast<float>(observations))));
        observations = maxNumberOfObservations;
    }
    float max_val = -3.402823466e+38f;
    uint8_t maxI = 0;
    for (uint8_t i = 0; i <= maxReasonableQualScore; i++) {
        float floati = static_cast<float>(i);
        float p1 = log10QualEmpiricalPrior(floati, priorQuality);
        float p2 = log10QualEmpiricalLikelihood(floati, static_cast<float>(observations), mismatches);
        float log10Posterior = p1 + p2;
        if (max_val < log10Posterior) {
            max_val = log10Posterior;
            maxI = i;
        }
    }
    return maxI;
}

SYCL_EXTERNAL inline uint8_t bqsrEntry::calculateEmpiricalQuality(float priorQual) {
    int64_t smoothMismatches = this->Mismatches + smoothingConstant;
    int64_t smoothObservations = this->Observations + smoothingConstant + smoothingConstant;
    uint8_t empiricalQual = calculateBayesianEstimateOfEmpiricalQuality(smoothObservations, smoothMismatches, priorQual);
    return std::min(empiricalQual, static_cast<uint8_t>(maxRecalibratedQualScore));
}

#endif //GENOMICPIPELINE_TABLESBQSR_H
