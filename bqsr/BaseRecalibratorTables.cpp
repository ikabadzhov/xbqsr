#include "BQSR.h"
#include <sycl/queue.hpp>
#include "../timer/Timer.h"
#include "BaseRecalibratorTables.h"

#include "UtilsBQSR.h"
#include "../Config.h"

#include "apply_bqsr/ApplyBqsrFunctions.h"
#include "apply_bqsr/ApplyBqsrKernelData.h"


BaseRecalibratorTables::BaseRecalibratorTables(const BQSRTableConfig& cfg, const sycl::queue& q)
    : nReadGroups(cfg.n_read_groups),
      nQualScore(cfg.n_qual_score),
      maxCycle  (cfg.max_cycle),
      contextLen(cfg.context_len),
      nCycleCov (cfg.max_cycle * 2),
      nContextCov(MaxContextCovariateKey(cfg.context_len) + 1),
      q_(q), ctx_(q_.get_context()), dev_(q_.get_device()), q_process_(q.get_context(), q.get_device(), sycl::async_handler{}), q_transfer_(q.get_context(), q.get_device(), sycl::async_handler{}) {

    bqsrEntry defaultEntry(0, 0);
    QualityScores_.table.resize(nReadGroups * nQualScore, defaultEntry);
    QualityScores_.empirical_qualities.resize(nReadGroups * nQualScore, 0);
    Cycles_.table.resize(nReadGroups * nQualScore * nCycleCov, defaultEntry);
    Cycles_.empirical_qualities.resize(nReadGroups * nQualScore * nCycleCov, 0);
    Contexts_.table.resize(nReadGroups * nQualScore * nContextCov, defaultEntry);
    Contexts_.empirical_qualities.resize(nReadGroups * nQualScore * nContextCov, 0);



   /* batch_process = std::make_shared<Batch<AppyBqsrKernelData, AppyBqsrKernelConfiguration>>(q_transfer_, q_process_, config);
    batch_transfer = std::make_shared<Batch<AppyBqsrKernelData, AppyBqsrKernelConfiguration>>(q_transfer_, q_process_, config);*/
}

void BaseRecalibratorTables::initializeDeviceMemory() {
    recalibrated_quality_scores = mem_t<int>(nReadGroups * nQualScore * nCycleCov * (nContextCov+1), q_, memory_space_device);
    q_.fill(recalibrated_quality_scores.data(), -1, recalibrated_quality_scores.size()).wait();

    size_t constant_memory = nReadGroups * nQualScore * nCycleCov * (nContextCov+1) * sizeof(int);
    size_t available_mem = q_.get_device().get_info<sycl::info::device::global_mem_size>() - constant_memory - 515 * 1024 * 1024;
    size_t max_alloc = q_.get_device().get_info<sycl::info::device::max_mem_alloc_size>();
    config.max_memory = std::min(available_mem / 3, max_alloc - 1024 * 1024); // Cap it with some margin

    batch_process = std::make_shared<Batch<AppyBqsrKernelData, AppyBqsrKernelConfiguration>>(q_transfer_, q_process_, config);
    batch_transfer = std::make_shared<Batch<AppyBqsrKernelData, AppyBqsrKernelConfiguration>>(q_transfer_, q_process_, config);
}


void BaseRecalibratorTables::merge(const BaseRecalibratorTables& other) {
    // 1. sanity: sizes must match
    const auto qs_size = qualityScores().size();
    const auto oc_qs_size = other.qualityScores().size();
    if (qs_size != oc_qs_size) {
        throw std::runtime_error("merge(): qualityScores size mismatch");
    }

    const auto cy_size = cycles().size();
    const auto oc_cy_size = other.cycles().size();
    if (cy_size != oc_cy_size) {
        throw std::runtime_error("merge(): cycles size mismatch");
    }

    const auto cx_size = contexts().size();
    const auto oc_cx_size = other.contexts().size();
    if (cx_size != oc_cx_size) {
        throw std::runtime_error("merge(): contexts size mismatch");
    }

    // 2. merge qualityScores
    {
        auto *dst = qualityScores().table.data();        // mutable bqsrEntry*
        const auto *src = other.qualityScores().table.data(); // const bqsrEntry*
        tbb::parallel_for(
            std::size_t(0), qs_size, [&](std::size_t i) {
                dst[i].Observations += src[i].Observations;
                dst[i].Mismatches   += src[i].Mismatches;
                // if you later track insertions/deletions/etc, add them here too
            }
        );
    }

    // 3. merge cycles
    {
        auto *dst = cycles().table.data();
        const auto *src = other.cycles().table.data();
        tbb::parallel_for(
            std::size_t(0), cy_size, [&](std::size_t i) {
                dst[i].Observations += src[i].Observations;
                dst[i].Mismatches   += src[i].Mismatches;
            }
        );
    }

    // 4. merge contexts
    {
        auto *dst = contexts().table.data();
        const auto *src = other.contexts().table.data();
        tbb::parallel_for(
            std::size_t(0), cx_size, [&](std::size_t i) {
                dst[i].Observations += src[i].Observations;
                dst[i].Mismatches   += src[i].Mismatches;
            }
        );
    }
}

void BaseRecalibratorTables::initializeCombinedBQSRTable(bqsrTable& bqsrTable, std::vector<combinedBqsrEntry>& table) {

    auto table_acc = make_accessor();

    for (size_t idx = 0; idx < bqsrTable.table.size(); idx++) {
        if (bqsrTable.table[idx].Observations > 0) {
            auto key = table_acc.idxToQualityScoresKey(idx);
            auto entry = bqsrTable.table[idx];

            if (table[key.ReadGroup].isValid) {
                auto& tableEntry = table[key.ReadGroup];
                auto sumErrors = calculateExpectedErrors(tableEntry.entry.Observations, tableEntry.reportedQuality) +
                                 calculateExpectedErrors(entry.Observations, float(key.Qual));
                tableEntry.entry.Observations += entry.Observations;
                tableEntry.entry.Mismatches += entry.Mismatches;
                tableEntry.reportedQuality = -10.0f * log10f(sumErrors / float(tableEntry.entry.Observations));
            } else {
                table[key.ReadGroup] = combinedBqsrEntry{static_cast<float>(key.Qual), entry};
                table[key.ReadGroup].isValid = true;
            }
        }
    }
    for (auto& tableEntry : table) {
        tableEntry.empiricalQuality = tableEntry.entry.calculateEmpiricalQuality(tableEntry.reportedQuality);
    }
}


void BaseRecalibratorTables::FinalizeBQSRTables() {
    auto table_acc = make_accessor();

    for (size_t idx = 0; idx < QualityScores_.size(); idx++) {
        if (QualityScores_.table[idx].Observations > 0) {
            auto key = table_acc.idxToQualityScoresKey(idx);
            QualityScores_.empirical_qualities[idx] = QualityScores_.table[idx].calculateEmpiricalQuality(
                    float(key.Qual));
        }
    }
    for (size_t idx = 0; idx < Contexts_.size(); idx++) {
        if (Contexts_.table[idx].Observations > 0) {
            auto key = table_acc.idxToContextCovKey(idx);
            Contexts_.empirical_qualities[idx] = Contexts_.table[idx].calculateEmpiricalQuality(
                    float(key.Qual));
        }
    }
    for (size_t idx = 0; idx < Cycles_.size(); idx++) {
        if (Cycles_.table[idx].Observations > 0) {
            auto key = table_acc.idxToCycleCovKey(idx);
            Cycles_.empirical_qualities[idx] = Cycles_.table[idx].calculateEmpiricalQuality(
                    float(key.Qual));
        }
    }
}

void BaseRecalibratorTables::ApplyBQSR(SamAlignments &reads, int32_t quantizeLevels, std::vector<uint8_t> &sqqList) {
    Timer t;

    std::vector<int32_t> reads_group_covariate(reads.GetAlnsNum());
    tbb::parallel_for(0UL, reads.GetAlnsNum(), [&](size_t idx) {
    //for (size_t idx = 0; idx < reads.GetAlnsNum(); idx++) {
        reads_group_covariate[idx] = readGroupCovariateId(reads.GetHdr(), reads, idx);
    });/*
    std::vector<int32_t> RG(reads.GetAlnsNum());
    tbb::parallel_for(0UL, reads.GetAlnsNum(), [&](size_t idx) {
        RG[idx] = reads.TAGS.RG[idx];
    });*/

    //if (!is_apply_bqsr_initialized) {
    std::call_once(apply_bqsr_once_, [&]{
        std::vector<uint8_t> sqqCopy = sqqList;
        combinedBQSRTable.resize(reads.GetHdr().getReadGroupCovariates().size());
        initializeCombinedBQSRTable(this->QualityScores_, combinedBQSRTable);

        if (!sqqCopy.empty()) {
            staticQuantizedScores = initializeStaticQuantizedScores(sqqCopy);
        }

        initializeQuantizedQualityScores(this->QualityScores_, quantizeLevels, quantizationMap, quantizedQualityScores);

        is_apply_bqsr_initialized = true;
    });

    t.reset();

    size_t flat_offset_sequences = 0;
    size_t flat_offset_qualities = 0;
    size_t offset = 0;
    size_t current_batch_size = reads.GetAlnsNum();

    auto table_acc = make_accessor();

    // Update the SequencesInput instance to process only the current batch.
    ApplyBqsrSequencesInput in(this->QualityScores_.table.data(), this->Cycles_.table.data(), this->Contexts_.table.data(),
                        this->QualityScores_.table.size(), this->Cycles_.table.size(), this->Contexts_.table.size(),
                          table_acc,
                          reads.SEQ.Data() + flat_offset_sequences,    // correct flat offset
                          reads.SEQ.Ptrs().data() + offset,                 // pointer array remains sliced at batch start
                          reinterpret_cast<uint8_t*>(reads.QUAL.Data() + flat_offset_qualities),      // adjust using pointer array
                          reads.FLAG.data() + offset,
                          reads_group_covariate.data(),
                          combinedBQSRTable.data(), combinedBQSRTable.size(), quantizedQualityScores.data(), quantizedQualityScores.size(),
                          staticQuantizedScores.data(), staticQuantizedScores.size(), recalibrated_quality_scores.data(), current_batch_size);


    ApplyBqsrOutput out(reinterpret_cast<uint8_t*>(reads.QUAL.Data())); // Not actullay used, but here for compatibility with my API


    // // Set up your configuration using available device memory
    // AppyBqsrKernelConfiguration config;
    // size_t available_mem = q_.get_device().get_info<sycl::info::device::global_mem_size>();
    // config.max_memory = available_mem / 3;

    // Process the current batch with your kernel
    batch_processing_variable_length<AppyBqsrKernelData, AppyBqsrKernelConfiguration, ApplyBqsrSequencesInput, ApplyBqsrOutput>(q_process_, q_transfer_, in, out, apply_bqsr_kernel, config, batch_process, batch_transfer);
    q_.wait();
    Logger::getLogger().log(LOG_INFO, "Time to batch process alignments: " + std::to_string(t.elapsed()));

}


/*
ssize_t BaseRecalibratorTables::getQualityScoresIdx(bqsrTableKey &key, size_t szQualScore) {
    return key.ReadGroup * szQualScore + key.Qual;
}

ssize_t BaseRecalibratorTables::getCycleCovIdx(bqsrTableKey &key, size_t szQualScore, size_t szCycleCov) {
    return key.ReadGroup * szQualScore * szCycleCov + key.Qual * szCycleCov + key.Covariate;//(-1 * key.Covariate + MAX_CYCLE);
}

ssize_t BaseRecalibratorTables::getContextCovIdx(bqsrTableKey &key, size_t szQualScore, size_t szCycleCov, size_t szContextCov) {
    return key.ReadGroup * szQualScore * szContextCov + key.Qual * szContextCov + key.Covariate;
}

bqsrTableKey BaseRecalibratorTables::idxToQualityScoresKey(uint8_t idx) {
    bqsrTableKey key;
    key.ReadGroup = idx / nQualScore;
    key.Qual = idx % nQualScore;
    key.Covariate = 0;
    return key;
}

bqsrTableKey BaseRecalibratorTables::idxToCycleCovKey(size_t idx) {
    bqsrTableKey key;
    key.ReadGroup = idx / (nQualScore * nCycleCov);
    auto tmpVal = idx % (nQualScore * nCycleCov);
    key.Qual = tmpVal / nCycleCov;
    key.Covariate = tmpVal % nCycleCov;
//    if (key.Covariate >= MAX_CYCLE) {
//        key.Covariate -= MAX_CYCLE;
//        key.Covariate *= -1;
//    }
    return key;
}

bqsrTableKey BaseRecalibratorTables::idxToContextCovKey(size_t idx) {
    bqsrTableKey key;
    key.ReadGroup = idx / (nQualScore * nContextCov);
    auto tmpVal = idx % (nQualScore * nContextCov);
    key.Qual = tmpVal / nContextCov;
    key.Covariate = tmpVal % nContextCov;
    return key;
}
*/

int32_t errorProbabilityToQuality(float prob) {
    if (prob <= 0.0f) {
        return maxQualityScore;
    }
    return fastRound(-10.0f * log10f(prob));
}

struct quantizationInterval  {
    int32_t next;
    double errorRate;
    int64_t nobs;
    int64_t leafNobs;
    int64_t nerrors;
};

bool leafInterval(int32_t index, quantizationInterval& interval) {
    if (interval.next < 0) {
        return index == maxQualityScore;
    }
    return interval.next == index+1;
}

std::vector<quantizationInterval> initializeQuantizationIntervals(std::vector<int64_t>& quantizationMap) {
    std::vector<quantizationInterval> intervals(quantizationMap.size());
    for (size_t i=0; i<quantizationMap.size(); i++) {
        auto nobs = quantizationMap[i];
        auto errorRate = qualityToErrorProbability(float(i));
        auto nerrors = int32_t(float(nobs) * errorRate);
        auto j = i + 1;
        if (j == quantizationMap.size()) {
            j = -1;
        }
        intervals[i].next = j;
        intervals[i].errorRate = errorRate;
        intervals[i].nobs = nobs;
        intervals[i].leafNobs = nobs;
        intervals[i].nerrors = nerrors;
    }
    return intervals;
}

float leafPenalty(int32_t k, std::vector<quantizationInterval>& quantizationIntervals, float globalErrorRate) {
    auto interval = &quantizationIntervals[k];
    if (k <= MIN_INTERESTING_QUAL) {
        return 0.0f;
    }
    return std::abs(log10f(interval->errorRate)-log10f(globalErrorRate)) * static_cast<float>(interval->leafNobs);
}

float calculateErrorRate(int64_t nobs, int64_t nerrors) {
    if (nobs == 0) {
        return 0.0f;
    }
    return static_cast<float>(nerrors+1) / static_cast<float>(nobs+1);
}

float computeMergePenalty(int32_t i, int32_t j, std::vector<quantizationInterval>& quantizationIntervals) {
    auto intervalI = &quantizationIntervals[i];
    auto intervalJ = &quantizationIntervals[j];

    int64_t mergedNobs = intervalI->nobs + intervalJ->nobs;
    int64_t mergedNerrors = intervalI->nerrors + intervalJ->nerrors;
    float mergedErrorRate = calculateErrorRate(mergedNobs, mergedNerrors);
    if (mergedErrorRate == 0.0f) {
        return 0.0f;
    }
    float sumI = 0.0f, sumJ = 0.0f;
    for (int32_t k = i; k < j; k++) {
        sumI += leafPenalty(k, quantizationIntervals, mergedErrorRate);
    }

    ssize_t kend;
    if (intervalJ->next >= 0) {
        kend = intervalJ->next;
    } else {
        kend = static_cast<ssize_t>(quantizationIntervals.size());
    }
    for (size_t k = j; k < kend; k++) {
        sumJ += leafPenalty(k, quantizationIntervals, mergedErrorRate);
    }

    return sumI + sumJ;
}

bool mergeMinimalPenaltyQuantizationIntervals(std::vector<quantizationInterval>& quantizationIntervals) {
    auto i = 0;
    auto intervalI = &quantizationIntervals[0];
    auto j = intervalI->next;
    if (j < 0) {
        return false;
    }
    auto minI = i;
    auto mergePenalty = computeMergePenalty(i, j, quantizationIntervals);
    while(true) {
        i = j;
        intervalI = &quantizationIntervals[i];
        j = intervalI->next;
        if (j < 0) {
            break;
        }
        auto penalty = computeMergePenalty(i, j, quantizationIntervals);
        if (penalty < mergePenalty) {
            minI = i;
            mergePenalty = penalty;
        }
    }
    auto intervalMinI = &quantizationIntervals[minI];
    auto intervalMinJ = &quantizationIntervals[intervalMinI->next];
    auto mergedNobs = intervalMinI->nobs + intervalMinJ->nobs;
    auto mergedNerrors = intervalMinI->nerrors + intervalMinJ->nerrors;
    intervalMinI->next = intervalMinJ->next;
    intervalMinI->nobs = mergedNobs;
    intervalMinI->nerrors = mergedNerrors;
    return true;
}

void mergeQuantizationIntervals(std::vector<quantizationInterval>& quantizationIntervals, int32_t levels) {
    auto n = quantizationIntervals.size();
    while (n > levels) {
        if (mergeMinimalPenaltyQuantizationIntervals(quantizationIntervals)) {
            n--;
        } else {
            break;
        }
    }
}

void initializeQuantizedQualityScores(bqsrTable& bqsrTable, int32_t quantLevel, std::vector<int64_t>& quantizationMap, std::vector<uint8_t>& quantizedScores) {
    quantizationMap.resize(maxQualityScore+1);
    quantizedScores.resize(maxQualityScore+1);
    if (quantLevel == 0) {
        for (size_t i = 0; i < maxQualityScore+1; i++) {
            quantizedScores[i] = uint8_t(i);
        }
        return;
    }

    for (size_t i = 0; i < bqsrTable.table.size(); i++) {
        if (bqsrTable.table[i].Observations >=0) {
            quantizationMap[bqsrTable.empirical_qualities[i]] += bqsrTable.table[i].Observations;
        }
    }
    auto quantizationIntervals = initializeQuantizationIntervals(quantizationMap);
    mergeQuantizationIntervals(quantizationIntervals, quantLevel);
    for (ssize_t i = 0; i >= 0;) {
        auto& intervalI = quantizationIntervals[i];
        uint8_t quantizedScore;
        if (leafInterval(i, intervalI)) {
            quantizedScore = uint8_t(i);
        } else {
            auto errorRate = calculateErrorRate(intervalI.nobs, intervalI.nerrors);
            quantizedScore = uint8_t(errorProbabilityToQuality(errorRate));
        }
        auto j = intervalI.next;
        int32_t kend;
        if (j >= 0) {
            kend = j;
        } else {
            kend = quantizationIntervals.size();
        }
        for (size_t k = i; k < kend; k++) {
            quantizedScores[k] = quantizedScore;
        }
        i = intervalI.next;
    }
}

std::vector<uint8_t> initializeStaticQuantizedScores(std::vector<uint8_t>& quals) {
    std::vector<uint8_t> staticScores(MAX_STATIC_QUANTIZE_QUALITY);
    for (uint8_t i = 0; i < MIN_INTERESTING_QUAL; i++) {
        staticScores[i] = uint8_t(i);
    }
    if (quals.size() == 1) {
        auto qual = quals[0];
        for (auto i = MIN_INTERESTING_QUAL; i < MAX_STATIC_QUANTIZE_QUALITY; i++) {
            staticScores[i] = qual;
        }
        return staticScores;
    }
    sort(quals.begin(), quals.end(), [](uint8_t a, uint8_t b){
        return a < b;
    });

    auto firstQual = uint8_t(MIN_INTERESTING_QUAL);
    auto prevQual = firstQual;
    auto prevProb = qualityToProbability(float(prevQual));

    for (auto nextQual :quals) {
        for (auto i = prevQual; i < nextQual; i++) {
            auto nextProb = qualityToProbability(float(nextQual));
            auto iProb = qualityToProbability(float(i));
            if (iProb-prevProb > nextProb-iProb) {
                staticScores[i] = nextQual;
            } else {
                staticScores[i] = prevQual;
            }
            prevProb = nextProb;
            prevQual = nextQual;
        }
    }
    for (auto i = prevQual; i < MAX_STATIC_QUANTIZE_QUALITY; i++) {
        staticScores[i] = prevQual;
    }
    return staticScores;
}

std::string keyToString(int32_t key) {
    auto length = key & 0xF;
    int32_t realKey = key >> 4;
    std::string result;
    for( size_t i = 0; i < length; i++ ){
        result += simpleBaseIndexToBase(realKey&0x3);
        realKey >>= 2;
    }

    if( realKey != 0 ){
        Logger::getLogger().log(LOG_FATAL, "keyToString(): invalid key");
        abort();
    }
    return result;
}
//alns.GetHdr().getReadGroupRecord().getString(alns.TAGS.RG[idx])
void BaseRecalibratorTables::PrintBQSRTables(const std::string &filename, Header &hdr) {
    std::ofstream file(filename);

    file << "#:" << BQSRTablenamePrefix << "Report.v1.1:5\n";
    file << "#:" << BQSRTablenamePrefix << "Table:2:17:%s:%s:;\n";
    file << "#:" << BQSRTablenamePrefix << "Table:Arguments:Recalibration argument collection values used in this run\n";
    file << "Argument                    Value                                                                   \n";
    file << "binary_tag_name             null                                                                    \n";
    file << "covariate                   ReadGroupCovariate,QualityScoreCovariate,ContextCovariate,CycleCovariate\n";
    file << "default_platform            null                                                                    \n";
    file << "deletions_default_quality   45                                                                      \n";
    file << "force_platform              null                                                                    \n";
    file << "indels_context_size         3                                                                       \n";
    file << "insertions_default_quality  45                                                                      \n";
    file << "low_quality_tail            2                                                                       \n";
    file << "maximum_cycle_value         500                                                                     \n";
    file << "mismatches_context_size     2                                                                       \n";
    file << "mismatches_default_quality  -1                                                                      \n";
    file << "no_standard_covs            false                                                                   \n";
    file << "quantizing_levels           16                                                                      \n";
    file << "recalibration_report        null                                                                    \n";
    file << "run_without_dbsnp           false                                                                   \n";
    file << "solid_nocall_strategy       THROW_EXCEPTION                                                         \n";
    file << "solid_recal_mode            SET_Q_ZERO                                                              \n";
    file << "\n";

    printQuantizationTable(file);
    printCombinedBQSRTable(file, hdr);
    printBQSRTable(file, hdr);
    printOtherCovariateTable(file, hdr);
}


void BaseRecalibratorTables::printQuantizationTable(std::ofstream &file) {

    std::vector<int64_t> observations; // placeholder type
    std::vector<uint8_t> scores;       // placeholder type

    initializeQuantizedQualityScores(this->QualityScores_, QUANTIZATION_LEVEL, observations, scores);

    file << "#:" << BQSRTablenamePrefix << "Table:3:" << observations.size() << ":%d:%d:%d:;\n";
    file << "#:" << BQSRTablenamePrefix << "Table:Quantized:Quality quantization map\n";

    size_t maxLenQualityScore = qualityScoreString.size();
    size_t maxLenCount = countString.size();
    size_t maxLenQuantizedScore = quantizedScoreString.size();

    for (size_t i = 0; i < observations.size(); i++) {
        maxLenQualityScore = std::max(maxLenQualityScore, std::to_string(i).size());
        maxLenCount = std::max(maxLenCount, std::to_string(observations[i]).size());
        maxLenQuantizedScore = std::max(maxLenQuantizedScore, std::to_string(scores[i]).size());
    }

    file << std::left << std::setw(static_cast<int32_t>(maxLenQualityScore)) << qualityScoreString;
    file << "  " << std::left << std::setw(static_cast<int32_t>(maxLenCount)) << countString;
    file << "  " << std::left << std::setw(static_cast<int32_t>(maxLenQuantizedScore)) << quantizedScoreString << "\n";

    for (size_t i = 0; i < observations.size(); i++) {
        file << std::right << std::setw(static_cast<int32_t>(maxLenQualityScore)) << i;
        file << "  " << std::right << std::setw(static_cast<int32_t>(maxLenCount)) << observations[i];
        file << "  " << std::right << std::setw(static_cast<int32_t>(maxLenQuantizedScore)) << static_cast<int>(scores[i]) << "\n";
    }

    file << "\n";
}

void BaseRecalibratorTables::printCombinedBQSRTable(std::ostream &file, Header &hdr) {
    std::vector<combinedBqsrEntry> table(hdr.getReadGroupCovariates().size()); // Table mapping from read groups to entries
    auto table_acc = make_accessor();
    initializeCombinedBQSRTable(QualityScores_, table);

    file << "#:" << BQSRTablenamePrefix << "Table:6:" << table.size() << ":%s:%s:%.4f:%.4f:%d:%.2f:;\n";
    file << "#:" << BQSRTablenamePrefix << "Table:RecalTable0:\n";

    size_t maxLenReadGroup = readGroupString.size();
    size_t maxLenEventType = eventTypeString.size();
    size_t maxLenEmpiricalQuality = empiricalQualityString.size();
    size_t maxLenEstimatedQReported = estimatedQReportedString.size();
    size_t maxLenObservations = observationsString.size();
    size_t maxLenErrors = errorsString.size();

    std::vector<int32_t> readGroups;

    for (size_t i = 0; i < table.size(); i++) {
        auto& entry = table[i].entry;
        auto& empiricalQuality = table[i].empiricalQuality;
        auto key = table_acc.idxToQualityScoresKey(i);
        readGroups.push_back(key.ReadGroup);
        maxLenReadGroup = std::max(maxLenReadGroup, hdr.covLabelForKey(key.ReadGroup).size());
        maxLenEmpiricalQuality = std::max(maxLenEmpiricalQuality, std::to_string(empiricalQuality).size() + 5);
        maxLenEstimatedQReported = std::max(maxLenEstimatedQReported, std::to_string(table[i].reportedQuality).size());
        maxLenObservations = std::max(maxLenObservations, std::to_string(entry.Observations).size());
        maxLenErrors = std::max(maxLenErrors, std::to_string(entry.Mismatches).size() + 3);
    }

    file << std::left << std::setw(static_cast<int32_t>(maxLenReadGroup)) << readGroupString
         << "  " << std::left << std::setw(static_cast<int32_t>(maxLenEventType)) << eventTypeString
         << "  " << std::left << std::setw(static_cast<int32_t>(maxLenEmpiricalQuality)) << empiricalQualityString
         << "  " << std::left << std::setw(static_cast<int32_t>(maxLenEstimatedQReported)) << estimatedQReportedString
         << "  " << std::left << std::setw(static_cast<int32_t>(maxLenObservations)) << observationsString
         << "  " << std::left << std::setw(static_cast<int32_t>(maxLenErrors)) << errorsString << "\n";

    std::sort(readGroups.begin(), readGroups.end(), [&](auto& e1, auto& e2){
        return hdr.covLabelForKey(e1) < hdr.covLabelForKey(e2);
    });

    for ( auto &rg : readGroups) {
        auto& combinedEntry = table[rg];
        auto& empiricalQuality = table[rg].empiricalQuality;
        auto& entry = table[rg].entry;
        file << std::left << std::setw(static_cast<int32_t>(maxLenReadGroup)) << hdr.getReadGroupCovariates().getString(rg)
             << "  " << std::left << std::setw(static_cast<int32_t>(maxLenEventType)) << "M"
             << "  " << std::right << std::setw(static_cast<int32_t>(maxLenEmpiricalQuality) - 5) << (int)empiricalQuality << ".0000"
             << "  " << std::fixed << std::setprecision(4) << std::right << std::setw(static_cast<int32_t>(maxLenEstimatedQReported)) << combinedEntry.reportedQuality
             << "  " << std::right << std::setw(static_cast<int32_t>(maxLenObservations)) << entry.Observations
             << "  " << std::right << std::setw(static_cast<int32_t>(maxLenErrors) - 3) << entry.Mismatches << ".00\n";
    }

    file << "\n";
}


void BaseRecalibratorTables::printBQSRTable(std::ofstream& file, Header &hdr) {
    size_t szQualityTable = 0;
    size_t maxLenReadGroup = readGroupString.size();
    size_t maxLenQualityScore = qualityScoreString.size();
    size_t maxLenEventType = eventTypeString.size();
    size_t maxLenEmpiricalQuality = empiricalQualityString.size();
    size_t maxLenObservations = observationsString.size();
    size_t maxLenErrors = errorsString.size();
    auto table_acc = make_accessor();
    std::vector<bqsrTableKey> keys;

    for (size_t i = 0; i < QualityScores_.table.size(); i++) {
        if (QualityScores_.table[i].Observations > 0) {
            auto key = table_acc.idxToQualityScoresKey(i);
            auto &entry = QualityScores_.table[i];
            auto &empiricalQuality = QualityScores_.empirical_qualities[i];
            keys.push_back(key);
            maxLenReadGroup = std::max(maxLenReadGroup, hdr.covLabelForKey(key.ReadGroup).size());
            maxLenQualityScore = std::max(maxLenQualityScore, std::to_string(key.Qual).size());
            maxLenEmpiricalQuality = std::max(maxLenEmpiricalQuality,
                                              std::to_string(empiricalQuality).size() + 5);
            maxLenObservations = std::max(maxLenObservations, std::to_string(entry.Observations).size());
            maxLenErrors = std::max(maxLenErrors, std::to_string(entry.Mismatches).size() + 3);
            szQualityTable++;
        }
    }

    file << "#:" << BQSRTablenamePrefix << "Table:6:" << szQualityTable << ":%s:%d:%s:%.4f:%d:%.2f:;\n";
    file << "#:" << BQSRTablenamePrefix << "Table:RecalTable1:\n";

    file << std::left << std::setw(maxLenReadGroup) << readGroupString
         << "  " << std::left << std::setw(maxLenQualityScore) << qualityScoreString
         << "  " << std::left << std::setw(maxLenEventType) << eventTypeString
         << "  " << std::left << std::setw(maxLenEmpiricalQuality) << empiricalQualityString
         << "  " << std::left << std::setw(maxLenObservations) << observationsString
         << "  " << std::left << std::setw(maxLenErrors) << errorsString << "\n";

    std::sort(keys.begin(), keys.end(), [](const bqsrTableKey &e1, const bqsrTableKey &e2) {
        if (e1.ReadGroup < e2.ReadGroup) return true;
        if (e1.ReadGroup == e2.ReadGroup) return e1.Qual < e2.Qual;
        return false;
    });

    for (auto &key : keys) {
        auto& entry = QualityScores_.table[table_acc.getQualityScoresIdx(key)];
        auto& empiricalQuality = QualityScores_.empirical_qualities[table_acc.getQualityScoresIdx(key)];
        file << std::left << std::setw(static_cast<int32_t>(maxLenReadGroup)) << hdr.covLabelForKey(key.ReadGroup)
             << "  " << std::right << std::setw(static_cast<int32_t>(maxLenQualityScore)) << (int)key.Qual
             << "  " << std::left << std::setw(static_cast<int32_t>(maxLenEventType)) << "M"
             << "  " << std::right << std::setw(static_cast<int32_t>(maxLenEmpiricalQuality) - 5) << (int)empiricalQuality << ".0000"
             << "  " << std::right << std::setw(static_cast<int32_t>(maxLenObservations)) << entry.Observations
             << "  " << std::right << std::setw(static_cast<int32_t>(maxLenErrors) - 3) << entry.Mismatches << ".00\n";
    }
    file << "\n";
}

void BaseRecalibratorTables::printOtherCovariateTable(std::ostream& file, Header &hdr) {
    size_t szCycleTable = 0;
    size_t szContextTable = 0;

    size_t maxLenReadGroup = readGroupString.size();
    size_t maxLenQualityScore = qualityScoreString.size();
    size_t maxLenCovariateValue = covariateValueString.size();
    size_t maxLenCovariateName = std::max(std::max(covariateNameString.size(), std::string("Cycle").size()), std::string("Context").size());
    size_t maxLenEventType = eventTypeString.size();
    size_t maxLenEmpiricalQuality = empiricalQualityString.size();
    size_t maxLenObservations = observationsString.size();
    size_t maxLenErrors = errorsString.size();

    auto table_acc = make_accessor();
    std::vector<covariateBqsrTableKey> keys;

    for (size_t i = 0; i < Cycles_.size(); i++) {
        if (Cycles_.table[i].Observations > 0) {
            auto key = table_acc.idxToCycleCovKey(i);
            auto entry = Cycles_.table[i];
            auto empiricalQuality = Cycles_.empirical_qualities[i];
            int32_t sign = (key.Covariate & 0x1);
            int32_t realCycle = (sign==0) ? (key.Covariate>>1) : (-1)*(key.Covariate>>1);
            std::string cycleText = std::to_string(realCycle);
            keys.push_back({true, cycleText, key});
            maxLenReadGroup = std::max(maxLenReadGroup, hdr.covLabelForKey(key.ReadGroup).size());
            maxLenQualityScore = std::max(maxLenQualityScore, std::to_string(key.Qual).size());
            maxLenCovariateValue = std::max(maxLenCovariateValue, cycleText.size());
            maxLenEmpiricalQuality = std::max(maxLenEmpiricalQuality,
                                              std::to_string(empiricalQuality).size() + 5);
            maxLenObservations = std::max(maxLenObservations, std::to_string(entry.Observations).size());
            maxLenErrors = std::max(maxLenErrors, std::to_string(entry.Mismatches).size() + 3);
            szCycleTable++;
        }
    }

    for (size_t i = 0; i < Contexts_.size(); i++) {
        if (Contexts_.table[i].Observations > 0) {
            auto key = table_acc.idxToContextCovKey(i);
            auto entry = Contexts_.table[i];
            auto empiricalQuality = Contexts_.empirical_qualities[i];

            std::string text = keyToString(key.Covariate);
            keys.push_back({false, text, key});
            maxLenReadGroup = std::max(maxLenReadGroup, hdr.covLabelForKey(key.ReadGroup).size());
            maxLenQualityScore = std::max(maxLenQualityScore, std::to_string(key.Qual).size());
            maxLenCovariateValue = std::max(maxLenCovariateValue, text.size());
            maxLenEmpiricalQuality = std::max(maxLenEmpiricalQuality,
                                              std::to_string(empiricalQuality).size() + 5);
            maxLenObservations = std::max(maxLenObservations, std::to_string(entry.Observations).size());
            maxLenErrors = std::max(maxLenErrors, std::to_string(entry.Mismatches).size() + 3);
            szContextTable++;
        }
    }

    file << "#:" << BQSRTablenamePrefix << "Table:8:" << szCycleTable + szContextTable << ":%s:%d:%s:%s:%s:%.4f:%d:%.2f:;\n";
    file << "#:" << BQSRTablenamePrefix << "Table:RecalTable2:\n";

    file << std::left << std::setw(static_cast<int32_t>(maxLenReadGroup)) << readGroupString
         << "  " << std::left << std::setw(static_cast<int32_t>(maxLenQualityScore)) << qualityScoreString
         << "  " << std::left << std::setw(static_cast<int32_t>(maxLenCovariateValue)) << covariateValueString
         << "  " << std::left << std::setw(static_cast<int32_t>(maxLenCovariateName)) << covariateNameString
         << "  " << std::left << std::setw(static_cast<int32_t>(maxLenEventType)) << eventTypeString
         << "  " << std::left << std::setw(static_cast<int32_t>(maxLenEmpiricalQuality)) << empiricalQualityString
         << "  " << std::left << std::setw(static_cast<int32_t>(maxLenObservations)) << observationsString
         << "  " << std::left << std::setw(static_cast<int32_t>(maxLenErrors)) << errorsString << "\n";

    std::sort(keys.begin(), keys.end(), [&](const covariateBqsrTableKey& e1, const covariateBqsrTableKey& e2) {
        if (hdr.covLabelForKey(e1.key.ReadGroup) < hdr.covLabelForKey(e2.key.ReadGroup)) return true;
        if (hdr.covLabelForKey(e1.key.ReadGroup) == hdr.covLabelForKey(e2.key.ReadGroup)) {
            if (e1.key.Qual < e2.key.Qual) return true;
            if (e1.key.Qual == e2.key.Qual) {
                return e1.text < e2.text;
            }
        }
        return false;
    });

    for (auto& key : keys) {
        bqsrEntry entry;
        uint8_t empiricalQuality;
        std::string name;
        if (key.cycleNotContext) {
            entry = Cycles_.table[table_acc.getCycleCovIdx(key.key)];
            empiricalQuality = Cycles_.empirical_qualities[table_acc.getCycleCovIdx(key.key)];
            name = "Cycle";
        } else {
            entry = Contexts_.table[table_acc.getContextCovIdx(key.key)];
            empiricalQuality = Contexts_.empirical_qualities[table_acc.getContextCovIdx(key.key)];
            name = "Context";
        }

        file << std::left << std::setw(static_cast<int32_t>(maxLenReadGroup)) << hdr.covLabelForKey(key.key.ReadGroup)
             << "  " << std::right << std::setw(static_cast<int32_t>(maxLenQualityScore)) << static_cast<int>(key.key.Qual)
             << "  " << std::left << std::setw(static_cast<int32_t>(maxLenCovariateValue)) << key.text
             << "  " << std::left << std::setw(static_cast<int32_t>(maxLenCovariateName)) << name
             << "  " << std::left << std::setw(static_cast<int32_t>(maxLenEventType)) << "M"
             << "  " << std::right << std::setw(static_cast<int32_t>(maxLenEmpiricalQuality) - 5) << static_cast<int>(empiricalQuality) << ".0000"
             << "  " << std::right << std::setw(static_cast<int32_t>(maxLenObservations)) << entry.Observations
             << "  " << std::right << std::setw(static_cast<int32_t>(maxLenErrors) - 3) << entry.Mismatches << ".00" << "\n";
    }
    file << "\n";
}

// --- small helpers ----------------------------------------------------------
static inline uint32_t baseToIndex(char b) {
    switch (b) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:
            Logger::getLogger().log(LOG_FATAL, std::string("LoadFromReport: invalid base '") + b + "'");
            std::abort();
    }
}

// pack a context string to the same integer key used by keyToString()
static inline int32_t contextStringToKey(const std::string& ctx) {
    const size_t L = ctx.size();
    if (L == 0 || L > 15) {
        Logger::getLogger().log(LOG_FATAL, "LoadFromReport: context length out of range (1..15)");
        std::abort();
    }
    int32_t realKey = 0;
    // NOTE: keyToString() encodes least-significant 2 bits as the first character (i==0).
    // When printing, it iterates i=0..L-1 doing realKey&0x3 then >>2.
    // To invert, we place ctx[i] at bits [2*i +: 2].
    for (size_t i = 0; i < L; ++i) {
        realKey |= static_cast<int32_t>(baseToIndex(ctx[i])) << (2 * i);
    }
    return (realKey << 4) | static_cast<int32_t>(L & 0xF);
}

static inline int32_t cycleToEncoded(int32_t cycle) {
    // same encoding you used when printing:
    // sign bit in LSB: 0 => non-negative, 1 => negative; value in the higher bits, shifted by 1
    if (cycle >= 0) return (cycle << 1) | 0;
    return ((-cycle) << 1) | 1;
}

static inline std::string trim(const std::string& s) {
    size_t b = 0, e = s.size();
    while (b < e && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(s[e-1]))) --e;
    return s.substr(b, e - b);
}

// Read a possibly “12.00” number into int64_t safely
static inline int64_t parseIntLike(const std::string& tok) {
    // try integer fast path
    bool all_digits = !tok.empty() && std::all_of(tok.begin(), tok.end(), [](char c){ return std::isdigit((unsigned char)c); });
    if (all_digits) return static_cast<int64_t>(std::stoll(tok));
    // fallback: double -> llround
    return static_cast<int64_t>(std::llround(std::stod(tok)));
}

void BaseRecalibratorTables::LoadFromReport(const std::string& filename, Header& hdr) {
    std::ifstream in(filename);
    if (!in) {
        Logger::getLogger().log(LOG_FATAL, "LoadFromReport: cannot open file: " + filename);
        std::abort();
    }

    // map RG string -> RG index (as in header)
    std::unordered_map<std::string, int32_t> rgIndex;
    const auto rgCount = hdr.getReadGroupCovariates().size();
    rgIndex.reserve(rgCount);
    for (int32_t rg = 0; rg < static_cast<int32_t>(rgCount); ++rg) {
        rgIndex.emplace(hdr.covLabelForKey(rg), rg);
    }

    enum class Section { None, RecalTable1, RecalTable2 };
    Section section = Section::None;

    auto acc = make_accessor();

    std::string line;
    while (std::getline(in, line)) {
        const std::string s = trim(line);
        if (s.empty()) { section = Section::None; continue; }

        // detect section starts
        if (s.rfind("#:", 0) == 0) {
            // table sentinels that matter for ApplyBQSR
            if (s.find("Table:RecalTable1:") != std::string::npos) {
                section = Section::RecalTable1;
            } else if (s.find("Table:RecalTable2:") != std::string::npos) {
                section = Section::RecalTable2;
            } else {
                // any other meta line (quantization / arguments / RecalTable0, etc.)
                // skip but reset to None so we don't parse their headers as data
                if (s.find("Table:") != std::string::npos) section = Section::None;
            }
            continue;
        }

        // skip column headers
        if (s.find("ReadGroup") != std::string::npos && s.find("QualityScore") != std::string::npos) {
            continue;
        }
        if (s.find("ReadGroup") != std::string::npos && s.find("CovariateValue") != std::string::npos) {
            continue;
        }

        // parse rows
        if (section == Section::RecalTable1) {
            // Columns:
            // RG  Q  EventType(M)  EmpQ  Observations  Errors
            std::istringstream iss(s);
            std::string rgStr, eventType;
            int q = 0;
            float empQ = 0.0f;
            std::string obsTok, errTok;

            if (!(iss >> rgStr >> q >> eventType >> empQ >> obsTok >> errTok)) {
                // be tolerant: some writers pad differently, try again reading all tokens
                continue;
            }
            if (eventType != "M") continue; // only mismatches used by your code

            auto rgIt = rgIndex.find(rgStr);
            if (rgIt == rgIndex.end()) continue; // unknown RG; skip

            bqsrTableKey key{};
            key.ReadGroup = static_cast<uint16_t>(rgIt->second);
            key.Qual      = static_cast<uint8_t>(q);
            key.Covariate = 0;

            const auto idx = acc.getQualityScoresIdx(key);
            if (idx < 0 || static_cast<size_t>(idx) >= QualityScores_.table.size()) continue;

            bqsrEntry& e = QualityScores_.table[static_cast<size_t>(idx)];
            e.Observations = parseIntLike(obsTok);
            e.Mismatches   = parseIntLike(errTok);
            // empirical_qualities will be recomputed in FinalizeBQSRTables()
        }
        else if (section == Section::RecalTable2) {
            // Columns:
            // RG  Q  CovValue  CovName(Context|Cycle)  M  EmpQ  Obs  Err
            std::istringstream iss(s);
            std::string rgStr, covValue, covName, eventType;
            int q = 0;
            float empQ = 0.0f;
            std::string obsTok, errTok;

            if (!(iss >> rgStr >> q >> covValue >> covName >> eventType >> empQ >> obsTok >> errTok)) {
                continue;
            }
            if (eventType != "M") continue;

            auto rgIt = rgIndex.find(rgStr);
            if (rgIt == rgIndex.end()) continue;

            bqsrTableKey key{};
            key.ReadGroup = static_cast<uint16_t>(rgIt->second);
            key.Qual      = static_cast<uint8_t>(q);

            if (covName == "Cycle") {
                // covValue is a signed integer string (e.g., "-12", "34")
                const int cyc = std::stoi(covValue);
                key.Covariate = static_cast<uint32_t>(cycleToEncoded(cyc));

                const auto idx = acc.getCycleCovIdx(key);
                if (idx < 0 || static_cast<size_t>(idx) >= Cycles_.table.size()) continue;

                bqsrEntry& e = Cycles_.table[static_cast<size_t>(idx)];
                e.Observations = parseIntLike(obsTok);
                e.Mismatches   = parseIntLike(errTok);
            } else if (covName == "Context") {
                // covValue is context string like "AC", "TT", etc.
                key.Covariate = static_cast<uint32_t>(contextStringToKey(covValue));

                const auto idx = acc.getContextCovIdx(key);
                if (idx < 0 || static_cast<size_t>(idx) >= Contexts_.table.size()) continue;

                bqsrEntry& e = Contexts_.table[static_cast<size_t>(idx)];
                e.Observations = parseIntLike(obsTok);
                e.Mismatches   = parseIntLike(errTok);
            } else {
                // Unknown covariate — ignore
                continue;
            }
        }
    }

    // Compute empirical qualities from what we loaded
    FinalizeBQSRTables();
    Logger::getLogger().log(LOG_INFO, "Loaded BQSR tables from file: " + filename);
}