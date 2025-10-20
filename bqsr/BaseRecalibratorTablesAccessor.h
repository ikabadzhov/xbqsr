#pragma once
#include "TablesBQSR.h"

class BaseRecalibratorTablesAccessor {
public:
    // --- Immutable metadata ---
    size_t nReadGroups;   // RG
    size_t nQualScore;    // Q
    size_t maxCycle;      // M
    size_t contextLen;    // L
    size_t nCycleCov;     // = 2 * M
    size_t nContextCov;   // = 4^L

    // --- Constructors ---
    BaseRecalibratorTablesAccessor() = default;

    BaseRecalibratorTablesAccessor(size_t rg, size_t q, size_t mc, size_t ctx,
                                  size_t cycle_cov, size_t context_cov)
       : nReadGroups(rg), nQualScore(q), maxCycle(mc), contextLen(ctx),
         nCycleCov(cycle_cov), nContextCov(context_cov) {}

    // --- Indexing helpers ---
    [[nodiscard]] SYCL_EXTERNAL
    inline ssize_t getQualityScoresIdx(const bqsrTableKey& key) const noexcept {
        return static_cast<ssize_t>(
            key.ReadGroup * nQualScore + key.Qual
        );
    }

    [[nodiscard]] SYCL_EXTERNAL
    inline ssize_t getCycleCovIdx(const bqsrTableKey& key) const noexcept {
        return static_cast<ssize_t>(
            key.ReadGroup * nQualScore * nCycleCov +
            key.Qual      * nCycleCov +
            key.Covariate
        );
    }

    [[nodiscard]] SYCL_EXTERNAL
    inline ssize_t getContextCovIdx(const bqsrTableKey& key) const noexcept {
        return static_cast<ssize_t>(
            key.ReadGroup * nQualScore * nContextCov +
            key.Qual      * nContextCov +
            key.Covariate
        );
    }

    // --- Reverse mapping ---
    [[nodiscard]] SYCL_EXTERNAL
    inline bqsrTableKey idxToQualityScoresKey(size_t idx) const noexcept {
        bqsrTableKey key;
        key.ReadGroup = idx / nQualScore;
        key.Qual      = idx % nQualScore;
        key.Covariate = 0;
        return key;
    }

    [[nodiscard]] SYCL_EXTERNAL
    inline bqsrTableKey idxToCycleCovKey(size_t idx) const noexcept {
        bqsrTableKey key;
        key.ReadGroup = idx / (nQualScore * nCycleCov);
        const auto tmp = idx % (nQualScore * nCycleCov);
        key.Qual      = tmp / nCycleCov;
        key.Covariate = tmp % nCycleCov;
        return key;
    }

    [[nodiscard]] SYCL_EXTERNAL
    inline bqsrTableKey idxToContextCovKey(size_t idx) const noexcept {
        bqsrTableKey key;
        key.ReadGroup = idx / (nQualScore * nContextCov);
        const auto tmp = idx % (nQualScore * nContextCov);
        key.Qual      = tmp / nContextCov;
        key.Covariate = tmp % nContextCov;
        return key;
    }
};
