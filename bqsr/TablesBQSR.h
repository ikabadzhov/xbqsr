#ifndef GENOMICPIPELINE_TABLESBQSR_H
#define GENOMICPIPELINE_TABLESBQSR_H

#include <cstdint>
#include <string>
#include <cmath>
#include <tbb/concurrent_unordered_map.h>
#include <memory>


const double Ln10 = 2.30258509299404568401799145468436420760110148862877297603332790;
const double Log10E = 1 / Ln10;

const int32_t smoothingConstant = 1;
const int32_t maxRecalibratedQualScore = 93;
const int32_t maxNumberOfObservations = INT32_MAX - 1;

const int32_t maxQualityScore = 93;
const int32_t maxReasonableQualScore = 60;


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
    uint8_t EmpiricalQuality{};
    int64_t Observations{};
    int64_t Mismatches{};

    explicit bqsrEntry(uint8_t emp_qual, int64_t observ, int64_t mis) : EmpiricalQuality(emp_qual), Observations(observ), Mismatches(mis) {}
    bqsrEntry() = default;

    SYCL_EXTERNAL uint8_t calculateEmpiricalQuality(double priorQual);
};

class bqsrTable {
public:
    explicit bqsrTable(size_t sz) : table(sz) {};
    bqsrTable() = default;
    SYCL_EXTERNAL static void update(bqsrEntry* table, size_t idx, int64_t errorVal);
    void merge(bqsrTable &other);

    std::vector<bqsrEntry> table;
    size_t size();
};

struct combinedBqsrEntry {
    double reportedQuality{};
    bqsrEntry entry;
    bool isValid = false;
};


SYCL_EXTERNAL double qualityToErrorProbability(double phred);
SYCL_EXTERNAL double qualityToProbability(double phred);
SYCL_EXTERNAL double calculateExpectedErrors(int64_t observations, double reportedQuality);
SYCL_EXTERNAL double log10QualEmpiricalLikelihood(double empiricalQual, double observations, int64_t mismatches);
SYCL_EXTERNAL double log10BinomialProbability(int64_t n, int64_t k, double log10p);
SYCL_EXTERNAL double log10BinomialCoefficient(int64_t n, int64_t k);
SYCL_EXTERNAL double log10Gamma(int64_t n);
SYCL_EXTERNAL double log10QualEmpiricalPrior(double empiricalQual, double reportedQual);

#endif //GENOMICPIPELINE_TABLESBQSR_H
