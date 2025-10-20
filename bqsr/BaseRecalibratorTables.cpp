#include <sycl/queue.hpp>
#include "../timer/Timer.h"
#include "BaseRecalibratorTables.h"
#include "../Config.h"

#include "apply_bqsr/ApplyBqsrFunctions.h"
#include "apply_bqsr/ApplyBqsrKernelData.h"



BaseRecalibratorTables::BaseRecalibratorTables(const BQSRTableConfig& cfg)
    : nReadGroups(cfg.n_read_groups),
      nQualScore(cfg.n_qual_score),
      maxCycle  (cfg.max_cycle),
      contextLen(cfg.context_len),
      nCycleCov (cfg.max_cycle * 2),
      nContextCov(MaxContextCovariateKey(cfg.context_len) + 1) {

    bqsrEntry defaultEntry(-1, 0);
    QualityScores_.table.resize(nReadGroups * nQualScore, defaultEntry);
    QualityScores_.empirical_qualities.resize(nReadGroups * nQualScore, 0);
    Cycles_.table.resize(nReadGroups * nQualScore * nCycleCov, defaultEntry);
    Cycles_.empirical_qualities.resize(nReadGroups * nQualScore * nCycleCov, 0);
    Contexts_.table.resize(nReadGroups * nQualScore * nContextCov, defaultEntry);
    Contexts_.empirical_qualities.resize(nReadGroups * nQualScore * nContextCov, 0);
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
                                 calculateExpectedErrors(entry.Observations, double(key.Qual));
                tableEntry.entry.Observations += entry.Observations;
                tableEntry.entry.Mismatches += entry.Mismatches;
                tableEntry.reportedQuality = -10 * log10(sumErrors / double(tableEntry.entry.Observations));
            } else {
                table[key.ReadGroup] = combinedBqsrEntry{static_cast<double>(key.Qual), entry};
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
                    double(key.Qual));
        }
    }
    for (size_t idx = 0; idx < Contexts_.size(); idx++) {
        if (Contexts_.table[idx].Observations > 0) {
            auto key = table_acc.idxToContextCovKey(idx);
            Contexts_.empirical_qualities[idx] = Contexts_.table[idx].calculateEmpiricalQuality(
                    double(key.Qual));
        }
    }
    for (size_t idx = 0; idx < Cycles_.size(); idx++) {
        if (Cycles_.table[idx].Observations > 0) {
            auto key = table_acc.idxToCycleCovKey(idx);
            Cycles_.empirical_qualities[idx] = Cycles_.table[idx].calculateEmpiricalQuality(
                    double(key.Qual));
        }
    }
}

void BaseRecalibratorTables::ApplyBQSR(sycl::device& device, SamAlignments &reads, int32_t quantizeLevels, std::vector<uint8_t> &sqqList,
                                       int32_t maxCycle) {
    Timer t;

    sycl::queue q(device);

    std::vector<combinedBqsrEntry> combinedBQSRTable(reads.GetHdr().getReadGroupRecord().size());
    initializeCombinedBQSRTable(this->QualityScores_, combinedBQSRTable);

    std::vector<uint8_t> staticQuantizedScores;
    if (!sqqList.empty()) {
        staticQuantizedScores = initializeStaticQuantizedScores(sqqList);
    }

    std::vector<int64_t> quantizationMap;
    std::vector<uint8_t> quantizedQualityScores;
    initializeQuantizedQualityScores(this->QualityScores_, quantizeLevels, quantizationMap, quantizedQualityScores);


    mem_t<int> recalibrated_quality_scores(nReadGroups * nQualScore * nCycleCov * (nContextCov+1), q, memory_space_device);
    q.fill(recalibrated_quality_scores.data(), -1, recalibrated_quality_scores.size()).wait();

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
                          reads.TMP.reads_group_covariate_.toVec().data() + offset,
                          combinedBQSRTable.data(), combinedBQSRTable.size(), quantizedQualityScores.data(), quantizedQualityScores.size(),
                          staticQuantizedScores.data(), staticQuantizedScores.size(), recalibrated_quality_scores.data(), current_batch_size);


    ApplyBqsrOutput out(reinterpret_cast<uint8_t*>(reads.QUAL.Data())); // Not actullay used, but here for compatibility with my API


    // Set up your configuration using available device memory
    AppyBqsrKernelConfiguration config;
    size_t available_mem = ConfigDeviceMemoryUsage::max_memory;//32*std::pow(2,30); //q.get_device().get_info<sycl::info::device::global_mem_size>() - 2 * pow(2, 30);
    config.max_memory = available_mem / 2;

    // Create temporary kernel data for this batch to compute the max elements
    AppyBqsrKernelData tmp_kernel_data(q, config, memory_space_t::none);
    size_t max_n_elements = apply_bqsr::update_batch_num_elements(in, tmp_kernel_data, 0);
    config.max_n_elements_first_batch = max_n_elements;

    // Process the current batch with your kernel
    batch_processing_variable_length<AppyBqsrKernelData, AppyBqsrKernelConfiguration, ApplyBqsrSequencesInput, ApplyBqsrOutput>(device, in, out, apply_bqsr_kernel, config);
    q.wait();
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

int32_t errorProbabilityToQuality(double prob) {
    if (prob == 0.0) {
        return maxQualityScore;
    }
    return std::max(std::min(int32_t(round(-10*log10(prob))), maxQualityScore), 1);
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
        auto errorRate = qualityToErrorProbability(double(i));
        auto nerrors = int32_t(double(nobs) * errorRate);
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

double leafPenalty(int32_t k, std::vector<quantizationInterval>& quantizationIntervals, double globalErrorRate) {
    auto interval = &quantizationIntervals[k];
    if (k <= MIN_INTERESTING_QUAL) {
        return 0.0;
    }
    return abs(log10(interval->errorRate)-log10(globalErrorRate)) * double(interval->leafNobs);
}

double calculateErrorRate(int64_t nobs, int64_t nerrors) {
    if (nobs == 0) {
        return 0.0;
    }
    return double(nerrors+1) / double(nobs+1);
}

double computeMergePenalty(int32_t i, int32_t j, std::vector<quantizationInterval>& quantizationIntervals) {
    auto intervalI = &quantizationIntervals[i];
    auto intervalJ = &quantizationIntervals[j];

    int64_t mergedNobs = intervalI->nobs + intervalJ->nobs;
    int64_t mergedNerrors = intervalI->nerrors + intervalJ->nerrors;
    double mergedErrorRate = calculateErrorRate(mergedNobs, mergedNerrors);
    if (mergedErrorRate == 0) {
        return 0.0;
    }
    double sumI = 0, sumJ = 0;
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
    sort(quals.begin(), quals.end(), [](auto& a, auto& b){
        return a < b;
    });

    auto firstQual = uint8_t(MIN_INTERESTING_QUAL);
    auto prevQual = firstQual;
    auto prevProb = qualityToProbability(double(prevQual));

    for (auto nextQual :quals) {
        for (auto i = prevQual; i < nextQual; i++) {
            auto nextProb = qualityToProbability(double(nextQual));
            auto iProb = qualityToProbability(double(i));
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
void BaseRecalibratorTables::PrintBQSRTables(const std::string &filename, std::vector<int32_t>& RG, Header &hdr) {
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
    printCombinedBQSRTable(file, RG, hdr);
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

void BaseRecalibratorTables::printCombinedBQSRTable(std::ostream &file, std::vector<int32_t>& RG, Header &hdr) {
    std::vector<combinedBqsrEntry> table(hdr.getReadGroupRecord().size()); // Table mapping from read groups to entries
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
        maxLenReadGroup = std::max(maxLenReadGroup, hdr.getReadGroupRecord().getString(key.ReadGroup).size());
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
        return hdr.getReadGroupRecord().getString(RG[e1]) < hdr.getReadGroupRecord().getString(RG[e2]);
    });

    for ( auto &rg : readGroups) {
        auto& combinedEntry = table[rg];
        auto& empiricalQuality = table[rg].empiricalQuality;
        auto& entry = table[rg].entry;
        file << std::left << std::setw(static_cast<int32_t>(maxLenReadGroup)) << hdr.getReadGroupRecord().getString(rg)
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
            maxLenReadGroup = std::max(maxLenReadGroup, hdr.getReadGroupRecord().getString(key.ReadGroup).size());
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
        file << std::left << std::setw(static_cast<int32_t>(maxLenReadGroup)) << hdr.getReadGroupRecord().getString(key.ReadGroup)
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
            maxLenReadGroup = std::max(maxLenReadGroup, hdr.getReadGroupRecord().getString(key.ReadGroup).size());
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
            maxLenReadGroup = std::max(maxLenReadGroup, hdr.getReadGroupRecord().getString(key.ReadGroup).size());
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
        if (hdr.getReadGroupRecord().getString(e1.key.ReadGroup) < hdr.getReadGroupRecord().getString(e2.key.ReadGroup)) return true;
        if (hdr.getReadGroupRecord().getString(e1.key.ReadGroup) == hdr.getReadGroupRecord().getString(e2.key.ReadGroup)) {
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

        file << std::left << std::setw(static_cast<int32_t>(maxLenReadGroup)) << hdr.getReadGroupRecord().getString(key.key.ReadGroup)
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

