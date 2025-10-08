#include "BQSR.h"
#include "../filters/UtilsFilters.h"
#include "TablesBQSR.h"
#include "../sam_data_structures/sam_fields/Sequence.h"
#include "../timer/Timer.h"
#include <iterator>
#include <tbb/parallel_reduce.h>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <iomanip>
#include <fstream>
#include <tbb/parallel_for.h>

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

using namespace std;

const uint32_t lengthBits = 4;
const uint32_t mismatchesContextSize = 2;
const uint32_t lowQualityTail = 2;

//const int32_t globalQualityScorePrior = -1;

const int maxStaticQuantizedQuality = 254;

int32_t fastRound(double d) {
    return (d > 0.0) ? static_cast<int>(d + 0.5) : static_cast<int>(d - 0.5);
}

int32_t keyFromContext(char* dna, SpanField dnaPtrs, int64_t start, int64_t end) {
    auto key = static_cast<int32_t>(end - start);

    uint32_t bitOffset = lengthBits;

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
    int32_t mask = ((1u << (contextSize * 2)) - 1) << lengthBits;
    size_t readLength = basesPtrs.len;

    int32_t fillCount = std::min(contextSize - 1, static_cast<int32_t>(readLength));
    std::fill_n(keys, fillCount, -1);
    size_t keyVectorIdx = fillCount;

    if (readLength < contextSize) {
        keysLength = 0;
        return;
    }

    //keys.resize(bases.size());
    keysLength = basesPtrs.len;

    auto newBaseOffset = static_cast<uint32_t>(2 * (contextSize - 1) + lengthBits);
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

void computeBaseContextCovariate(bool isReversed, char* strandedClippedSeq, SpanField& strandedClippedSeqPtrs, int32_t* keys, uint16_t& keyLength) {
    contextWith(strandedClippedSeq, strandedClippedSeqPtrs, mismatchesContextSize, keys, keyLength);
    assert(keyLength <= strandedClippedSeqPtrs.len);
    if (isReversed){
        std::reverse(keys, keys + keyLength);
    }
}

char simpleBaseIndexToBase( int32_t baseIndex ) {
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
            cout<<"invalid baseIndex"<<std::endl;
            exit(-1);
    }
}

std::string keyToString(int32_t key) {
    auto length = static_cast<int32_t>(key & 0xF);
    int32_t realKey = key >> 4;
    string result;
    for( size_t i = 0; i < length; i++ ){
        result += simpleBaseIndexToBase(realKey&0x3);
        realKey >>= 2;
    }

    if( realKey != 0 ){
        cout<<"invalid key"<<std::endl;
        exit(-1);
    }
    return result;
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
        !isEmptyRG && alignmentAgreesWithHeader(hdr, RNAME, POS) ) {

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
        cout<<"cycle value exceeds maximum cycle value"<<std::endl;
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
    // This condition is moved later
    /*for (size_t skipIdx = 0; skipIdx < skipSlice.size(); skipIdx++) {
        skipSlice[skipIdx] = !simpleBaseToBaseIndex(SEQ.Base(skipIdx)) || ...
    }*/
}

// NewBaseRecalibrator returns a struct for the first step of base recalibration.
std::unique_ptr<BaseRecalibrator> NewBaseRecalibrator(std::vector<std::string>& knownSites, MappedFasta& referenceFasta) {

    //std::unique_ptr<BaseRecalibrator> recal = std::make_unique<BaseRecalibrator>(referenceFasta, knownSites);

    return {};
}

// NewBaseRecalibratorTables returns a struct for storing the result
// of the base recalibration.
//BaseRecalibratorTables NewBaseRecalibratorTables()  {}

void BaseRecalibrator::close() {
    this->knownIntervals.clear();
}

std::string readGroupCovariate(Header &hdr, SamAlignments& reads, size_t idx) {

    std::string rg = hdr.getReadGroupRecord().getString(reads.TAGS.RG.Get(idx));
    if(!rg.empty()){

        for( auto& record : hdr.getReadGroupRecord().getAll() ){
            std::string id = record.getReadGroupId();

            if(!id.empty()) {
                id = record.getReadGroupId();
                if (rg == id) {
                    bool ok = record.getAllAttributes().containsKey("PU");
                    std::string unit;
                    if(ok) {
                        unit = record.getAttribute("PU");
                        return unit;
                    }
                    break;
                }
            }
        }
        return rg;
    } else {
        std::cout<<"Error: bqsr requires input with read groups. Alignment at pos " << idx <<" has no read group."/*Please fix input, e.g. rerun with the --replace-read-group option.*/<<std::endl;
        exit(-1);
        std::cout<<"unreachable code" << std::endl;
    }
}

// Recalibrate implements the first step of base recalibration.
//[[clang::optnone]]
void BaseRecalibrator::Recalibrate(sycl::device& device, SamAlignments& reads, int32_t maxCycle) {
    Timer timerAll;
    if (reads.GetAlnsNum()==0) return;
    auto& profiler = ProfilerBQSR::instance();
    timerAll.reset();
    Timer t;

    t.reset();

    int32_t maxIdx = reads.GetAlnsNum();
    auto& hdr = reads.GetHdr();

    profiler.start(StageId::AddReadGroup);

    reads.TMP.reads_group_covariate_.Resize(maxIdx);
    tbb::parallel_for(0, maxIdx, [&](size_t idx) {
        reads.TMP.reads_group_covariate_.Set(idx, readGroupCovariate(hdr, reads, idx));
    });

    profiler.stop(StageId::AddReadGroup);

    //std::vector<char> clippedSequences(reads.SEQ.Data(), reads.SEQ.Data() + reads.SEQ.size());
    //std::vector<SpanField> clippedSequencesPtrs(reads.SEQ.Ptrs().begin(), reads.SEQ.Ptrs().end());

    //std::string clippedQualities(reads.QUAL.Data(), reads.QUAL.Data() + reads.QUAL.size());
    //std::vector<SpanField> clippedQualitiesPtrs(reads.QUAL.Ptrs().begin(), reads.QUAL.Ptrs().end());

    profiler.start(StageId::Transform);

    std::vector<int64_t> pos(reads.POS.begin(), reads.POS.end());

    //auto *pos = reads.POS.data();
    auto *clippedSequences = reads.SEQ.Data();
    auto *clippedSequencesPtrs = reads.SEQ.Ptrs().data();
    auto *clippedQualities = reads.QUAL.Data();
    auto *clippedQualitiesPtrs = reads.QUAL.Ptrs().data();


    std::vector<uint8_t> recalibrate_aln_flag(reads.GetAlnsNum());

    Logger::getLogger().log(LOG_INFO, "Time to declare and initialize variables: " + std::to_string(t.elapsed()));

    t.reset();

    tbb::parallel_for(0UL, reads.GetAlnsNum(), [&] (size_t idx) {
        recalibrate_aln_flag[idx] = recalibrateAln(hdr, idx, reads.MAPQ[idx], reads.FLAG[idx], reads.RNAME[idx], reads.PNEXT[idx], reads.POS[idx], reads.GetSEQ(idx).size(), reads.GetQUAL(idx).size(), reads.CIGAR.Data(), reads.CIGAR.Ptrs()[idx], reads.TAGS.RG.Get(idx));
    });

    Logger::getLogger().log(LOG_INFO, "Time to set recalibrate flag: " + std::to_string(t.elapsed()));

    t.reset();

    size_t max_thread = std::thread::hardware_concurrency();
    size_t work_sz = (reads.GetAlnsNum() + max_thread - 1) / max_thread;

    t.reset();

    std::vector<std::vector<std::vector<CigarOperation>>> local_clipped_cigar(max_thread, std::vector<std::vector<CigarOperation>>(work_sz));
    std::vector<std::vector<SpanField>> local_clipped_cigar_ptrs(max_thread);

    tbb::parallel_for(0UL, max_thread, [&](size_t threadIdx) {

        for (size_t idx = threadIdx * work_sz, local_idx = 0; idx < (threadIdx * work_sz + work_sz); idx++, local_idx++) {
            if (idx >= maxIdx) {
                return;
            }
            if (!recalibrate_aln_flag[idx]) {
                continue;
            }

            std::vector<CigarOperation> CIGAR(reads.GetCIGAR(idx).begin(), reads.GetCIGAR(idx).end());


            hardClipAdaptorSequence(CIGAR, pos[idx], clippedSequences, clippedSequencesPtrs[idx],
                                    clippedQualities, clippedQualitiesPtrs[idx], idx, reads.GetRNAME(idx),
                                    reads.GetFLAG(idx),
                                    reads.GetRNEXT(idx), reads.GetTLEN(idx),
                                    reads.GetPNEXT(idx)); // It call a function that modifies the aln

            if (clippedSequencesPtrs[idx].len == 0) {
                continue;
            }
            hardClipSoftClippedBases(CIGAR, pos[idx], clippedSequences, clippedSequencesPtrs[idx],
                                     clippedQualities, clippedQualitiesPtrs[idx], reads.GetRNAME(idx),
                                     reads.GetFLAG(idx)); // It call a function that modifies the aln

            if (clippedSequencesPtrs[idx].len == 0) {
                continue;
            }
            local_clipped_cigar[threadIdx][local_idx] = CIGAR; // TODO: Change this when doing multithreading because idx will not start from 0 when number of threads is greater than 1
        }
    });

    Logger::getLogger().log(LOG_INFO, "Time to clip sequences: " + std::to_string(t.elapsed()));

    t.reset();
    size_t totSz = 0;
    for (const auto& cigars_thread : local_clipped_cigar) {
        for (const auto& cigar : cigars_thread) {
            totSz += cigar.size();
        }
    }

    std::vector<CigarOperation> clipped_cigar_flat(totSz);
    std::vector<size_t> clipped_cigar_flat_start(reads.GetAlnsNum());
    std::vector<uint16_t> clipped_cigar_flat_len(reads.GetAlnsNum());
    size_t i = 0;
    size_t ptr = 0;
    for (const auto& cigars_thread : local_clipped_cigar) {
        for (const auto& cigar : cigars_thread) {
            if (i >= reads.GetAlnsNum()) {
                break;
            }
            if (cigar.size() > 0) {
                std::memcpy(clipped_cigar_flat.data() + ptr, cigar.data(), cigar.size() * sizeof(CigarOperation));
            }
            clipped_cigar_flat_start[i] = ptr;
            clipped_cigar_flat_len[i] = cigar.size();
            ptr += cigar.size();
            i++;
        }
    }

    Logger::getLogger().log(LOG_INFO, "Time to flatten clipped cigar vectors: " + std::to_string(t.elapsed()));

    t.reset();

    size_t total_sequences = reads.GetAlnsNum();
    size_t batch_size = reads.GetAlnsNum();

    profiler.stop(StageId::Transform);

    profiler.start(StageId::RecalibrateBatches);

    t.reset();
    // Update the SequencesInput instance to process only the current batch.
    SequencesInput in(local_known_intervals, local_known_intervals_ptrs,
                      fasta, fasta_ptrs,
                      tables_quality.data(), tables_cycle.data(), tables_context.data(), tables_quality.size(), tables_cycle.size(), tables_context.size(),
                      recalibration_tables.szQualScore, recalibration_tables.szCycleCov, recalibration_tables.szContextCov,
                      clippedSequences,
                      clippedSequencesPtrs,
                      clippedQualities,
                      clippedQualitiesPtrs,
                      clipped_cigar_flat.data(),
                      clipped_cigar_flat_start.data(),
                      clipped_cigar_flat_len.data(),
                      pos.data(),
                      reads.RNAME.data(),
                      recalibrate_aln_flag.data(),
                      reads.FLAG.data(),
                      reads.TMP.reads_group_covariate_.toVec().data(),
                      batch_size);


    TablesOutput out; // Not actullay used, but here for compatibility with my API


    // Set up your configuration using available device memory
    ComputeSkipSliceConfiguration config;
    size_t available_mem = ConfigDeviceMemoryUsage::max_memory;//32*std::pow(2,30); //q.get_device().get_info<sycl::info::device::global_mem_size>() - 2 * pow(2, 30);
    config.max_memory = available_mem / 2;

    // Create temporary kernel data for this batch to compute the max elements
    ComputeSkipSliceKernelData tmp_kernel_data(q, config, memory_space_t::none);
    size_t max_n_elements = bqsr_covariate_vector_kernel::update_batch_num_elements(in, tmp_kernel_data, 0);
    config.max_n_elements_first_batch = max_n_elements;

    // Process the current batch with your kernel
    batch_processing_variable_length<ComputeSkipSliceKernelData, ComputeSkipSliceConfiguration, SequencesInput, TablesOutput>(device, in, out, compute_skip_slice_kernel, config);
    q.wait();
    Logger::getLogger().log(LOG_INFO, "Time to batch process alignments: " + std::to_string(t.elapsed()));

    profiler.stop(StageId::RecalibrateBatches);

    Logger::getLogger().log(LOG_INFO, "Total time recalibrate: " + std::to_string(timerAll.elapsed()) + "sec");
    profiler.report();
}

BaseRecalibrator::BaseRecalibrator(MappedFasta &fastaMap, vector<std::string> &knownSites,
                                   sycl::device &device, Header &hdr, BaseRecalibratorTables& recalibration_tables) : recalibration_tables(recalibration_tables), fastaMap(fastaMap), q(device) {
    Timer t;
    t.reset();
    for (std::string& knownSitesFile : knownSites){
        std::unordered_map<std::string, std::vector<Interval>> knownSitesIntervals;
        Interval::ImportIntervals(knownSitesFile, knownSitesIntervals);
        for (auto& p : knownSitesIntervals) {
            auto chrom = p.first;
            auto ivals = p.second;
            auto s1 = knownIntervals[chrom].size();
            auto s2 = ivals.size();
            knownIntervals[chrom].reserve(s1+s2);
            std::move(ivals.begin(), ivals.end(), std::back_inserter(knownIntervals[chrom]));
        }
    }
    for (auto& [chrom, ivals] : knownIntervals) {
        Interval::ParallelSortByStart(ivals);
        knownIntervals[chrom] = Interval::ParallelFlatten(ivals);
    }

    size_t sz_local_intervals_ptr = hdr.getSequenceDictionary().getAll().size();

    local_known_intervals_ptrs = mem_t<SpanFieldLong>(sz_local_intervals_ptr, q, memory_space_shared);

    size_t totSz = 0;
    for (const auto& chrIntervals : this->knownIntervals) {
        if (hdr.getSequenceDictionary().containsSequence(chrIntervals.first)) {
            auto chrIdx = hdr.getSequenceDictionary().getIdx(chrIntervals.first);
            totSz += chrIntervals.second.size();
            local_known_intervals_ptrs.data()[chrIdx] = SpanFieldLong(0, chrIntervals.second.size());
        }
    }

    local_known_intervals = mem_t<Interval>(totSz, q, memory_space_shared);

    for (size_t idx = 1; idx < local_known_intervals_ptrs.size(); idx++) {
        local_known_intervals_ptrs.data()[idx].start += local_known_intervals_ptrs.data()[idx-1].start + local_known_intervals_ptrs.data()[idx - 1].len;
    }

    for (const auto& chrIntervals : knownIntervals) {
        if (hdr.getSequenceDictionary().containsSequence(chrIntervals.first)) {
            auto chrIdx = hdr.getSequenceDictionary().getIdx(chrIntervals.first);
            assert(local_known_intervals_ptrs.data()[chrIdx].len == chrIntervals.second.size());
            for (size_t intervIdx = 0; intervIdx < chrIntervals.second.size(); intervIdx++) {
                local_known_intervals.data()[local_known_intervals_ptrs.data()[chrIdx].start + intervIdx] = chrIntervals.second[intervIdx];
            }
        }
    }

    fasta_ptrs = mem_t<SpanFieldLong>(hdr.getSequenceDictionary().getAll().size(), q, memory_space_shared);

    totSz = 0;
    for (auto& references : fastaMap.getFullFasta()) {
        if (hdr.getSequenceDictionary().containsSequence(references.first)) {
            auto chrIdx = hdr.getSequenceDictionary().getIdx(references.first);
            totSz += references.second.size();
            fasta_ptrs.data()[chrIdx] = SpanFieldLong(0, references.second.size());
        }
    }

    fasta = mem_t<uint8_t>(totSz, q, memory_space_shared);

    for (size_t idx = 1; idx < fasta_ptrs.size(); idx++) {
        fasta_ptrs.data()[idx].start += fasta_ptrs.data()[idx-1].start + fasta_ptrs.data()[idx - 1].len;
    }

    for (auto& references : this->fastaMap.getFullFasta()) {
        if (hdr.getSequenceDictionary().containsSequence(references.first)) {
            auto chrIdx = hdr.getSequenceDictionary().getIdx(references.first);
            assert(fasta_ptrs.data()[chrIdx].len == references.second.size());
            for (size_t intervIdx = 0; intervIdx < references.second.size(); intervIdx++) {
                fasta.data()[fasta_ptrs.data()[chrIdx].start + intervIdx] = references.second[intervIdx];
            }
        }
    }

    Logger::getLogger().log(LOG_INFO, "Time to flatten fasta and intervals: " + std::to_string(t.elapsed()));
    this->tables_quality = mem_t<bqsrEntry>(recalibration_tables.QualityScores.size(), q, memory_space_t::memory_space_shared);
    this->tables_context = mem_t<bqsrEntry>(recalibration_tables.Contexts.size(), q, memory_space_t::memory_space_shared);
    this->tables_cycle = mem_t<bqsrEntry>(recalibration_tables.Cycles.size(), q, memory_space_t::memory_space_shared);


    bqsrEntry defaultEntry;
    defaultEntry.Observations = 0;
    defaultEntry.Mismatches = 0;

    q.fill(tables_quality.data(), defaultEntry, recalibration_tables.QualityScores.size());
    q.fill(tables_cycle.data(), defaultEntry, recalibration_tables.Cycles.size());
    q.fill(tables_context.data(), defaultEntry, recalibration_tables.Contexts.size());
    q.wait();
    Logger::getLogger().log(LOG_INFO, "Time to allocate and initialize device tables in shared memory: " + std::to_string(t.elapsed()));
}

//[[clang::optnone]]
void BaseRecalibrator::getTables(BaseRecalibratorTables &recalibration_tables) {
    memcpy(recalibration_tables.QualityScores.table.data(), tables_quality.data(), sizeof(bqsrEntry)*recalibration_tables.QualityScores.size());
    memcpy(recalibration_tables.Cycles.table.data(), tables_cycle.data(), sizeof(bqsrEntry)*recalibration_tables.Cycles.size());
    memcpy(recalibration_tables.Contexts.table.data(), tables_context.data(), sizeof(bqsrEntry)*recalibration_tables.Contexts.size());
}

void BaseRecalibratorTables::initializeCombinedBQSRTable(bqsrTable& bqsrTable, std::vector<combinedBqsrEntry>& table) {

    for (size_t idx = 0; idx < bqsrTable.table.size(); idx++) {
        if (bqsrTable.table[idx].Observations > 0) {
            auto key = idxToQualityScoresKey(idx);
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

// FinalizeBQSRTables finalizes the first step of base recalibration.
//[[clang::optnone]]
void BaseRecalibratorTables::FinalizeBQSRTables() {
    for (size_t idx = 0; idx < QualityScores.size(); idx++) {
        if (QualityScores.table[idx].Observations > 0) {
            auto key = idxToQualityScoresKey(idx);
            QualityScores.empirical_qualities[idx] = QualityScores.table[idx].calculateEmpiricalQuality(
                    double(key.Qual));
        }
    }
    for (size_t idx = 0; idx < Contexts.size(); idx++) {
        if (Contexts.table[idx].Observations > 0) {
            auto key = idxToContextCovKey(idx);
            Contexts.empirical_qualities[idx] = Contexts.table[idx].calculateEmpiricalQuality(
                    double(key.Qual));
        }
    }
    for (size_t idx = 0; idx < Cycles.size(); idx++) {
        if (Cycles.table[idx].Observations > 0) {
            auto key = idxToCycleCovKey(idx);
            Cycles.empirical_qualities[idx] = Cycles.table[idx].calculateEmpiricalQuality(
                    double(key.Qual));
        }
    }
}

int32_t errorProbabilityToQuality(double prob) {
    if (prob == 0.0) {
        return maxQualityScore;
    }
    return std::max(std::min(int32_t(round(-10*log10(prob))), maxQualityScore), 1);
}

std::vector<uint8_t> initializeStaticQuantizedScores(std::vector<uint8_t>& quals) {
    std::vector<uint8_t> staticScores(maxStaticQuantizedQuality);
    for (uint8_t i = 0; i < minInterestingQual; i++) {
        staticScores[i] = uint8_t(i);
    }
    if (quals.size() == 1) {
        auto qual = quals[0];
        for (auto i = minInterestingQual; i < maxStaticQuantizedQuality; i++) {
            staticScores[i] = qual;
        }
        return staticScores;
    }
    sort(quals.begin(), quals.end(), [](auto& a, auto& b){
        return a < b;
    });

    auto firstQual = uint8_t(minInterestingQual);
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
    for (auto i = prevQual; i < maxStaticQuantizedQuality; i++) {
        staticScores[i] = prevQual;
    }
    return staticScores;
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
    if (k <= minInterestingQual) {
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


// ApplyBQSR applies the base recalibration result to the QUAL strings of the given reads.
//[[clang::optnone]]
void BaseRecalibratorTables::ApplyBQSR(sycl::device& device, SamAlignments &reads, int32_t quantizeLevels, std::vector<uint8_t> &sqqList,
                                       int32_t maxCycle) {
    Timer t;

    sycl::queue q(device);

    std::vector<combinedBqsrEntry> combinedBQSRTable(reads.TMP.reads_group_covariate_.NumEntriesTable());
    initializeCombinedBQSRTable(this->QualityScores, combinedBQSRTable);

    std::vector<uint8_t> staticQuantizedScores;
    if (!sqqList.empty()) {
        staticQuantizedScores = initializeStaticQuantizedScores(sqqList);
    }

    std::vector<int64_t> quantizationMap;
    std::vector<uint8_t> quantizedQualityScores;
    initializeQuantizedQualityScores(this->QualityScores, quantizeLevels, quantizationMap, quantizedQualityScores);


    mem_t<int> recalibrated_quality_scores(nReadGroups * szQualScore * szCycleCov * (szContextCov+1), q, memory_space_device);
    q.fill(recalibrated_quality_scores.data(), -1, recalibrated_quality_scores.size()).wait();

    t.reset();

    size_t flat_offset_sequences = 0;
    size_t flat_offset_qualities = 0;
    size_t offset = 0;
    size_t current_batch_size = reads.GetAlnsNum();

    // Update the SequencesInput instance to process only the current batch.
    ApplyBqsrSequencesInput in(this->QualityScores.table.data(), this->Cycles.table.data(), this->Contexts.table.data(),
                        this->QualityScores.table.size(), this->Cycles.table.size(), this->Contexts.table.size(),
                          szQualScore, szCycleCov, szContextCov,
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


/*
      size_t workload = reads.GetAlnsNum() / 12;

      applyBuffers buffers(nReadGroups, szQualScore, maxCycle, szContextCov);

      size_t count = 0;
      size_t countValid = 0;
      size_t nReads = reads.GetAlnsNum();
      for (size_t idx = 0; idx < nReads; idx++, count++) {
          if (count % 100000==0) {
              std::cout << count << "\r"<<std::flush;
          }

          auto result = true;

          int32_t readGrpCov = reads.TMP.reads_group_covariate_.Get(idx);

          if (!combinedBQSRTable[readGrpCov].isValid) {
              return;
          }
          countValid++;

          auto empiricalReadGroupEntry = combinedBQSRTable[readGrpCov];

          auto [cycleFactor, cycleIncrement] = prepareCycleCovariates((uint16_t)reads.GetSEQ(idx).Len(), reads.GetFLAG(idx));

          computeStrandedClippedSeq(reads.GetSEQ(idx), reads.GetQUAL(idx), SamAlignments::IsReversed(reads.GetFLAG(idx)), buffers.strandedClippedSeq);

          computeBaseContextCovariate(SamAlignments::IsReversed(reads.GetFLAG(idx)), buffers.strandedClippedSeq, buffers.baseContextCovariate);

          double epsilon;
          if (globalQualityScorePrior > 0.0) {
              epsilon = globalQualityScorePrior;
          } else {
              epsilon = empiricalReadGroupEntry.reportedQuality;
          }
          // Adapt each base's quality score
          auto alnLength = static_cast<uint32_t>(reads.GetSEQ(idx).Len());
          applyKey key{};
          key.readGroupCovariate = readGrpCov;

          for (uint32_t i = 0; i < alnLength; i++) {
              uint8_t qual = reads.GetQUAL(idx)[i];
              if (qual >= minInterestingQual) {
                  key.qual = qual;
                  key.cycleCovariate = computeBaseCycleCovariate(cycleFactor, cycleIncrement, i, maxCycle);
                  key.contextCovariate = buffers.baseContextCovariate[i];


                  if (buffers.recalibratedQualityScores[buffers.getIdxRecalQualScore(key)] > 0) {
                          auto recalibratedQualityScore = buffers.recalibratedQualityScores[buffers.getIdxRecalQualScore(
                                  key)];
                          reads.GetQUAL(idx)[i] = static_cast<char>(recalibratedQualityScore);
                          continue;
                  }

                  bqsrTableKey bqKey;
                  bqKey.ReadGroup = readGrpCov;
                  bqKey.Qual = qual;

                  auto empiricalBqEntry = this->QualityScores.table[getQualityScoresIdx(bqKey)];

                  bqsrTableKey cycleKey;
                  cycleKey.Covariate = key.cycleCovariate;
                  cycleKey.ReadGroup = readGrpCov;
                  cycleKey.Qual = qual;

                  auto empiricalCycleEntry = this->Cycles.table[getCycleCovIdx(cycleKey)];

                  bqsrTableKey contextKey;
                  contextKey.Covariate = key.contextCovariate;
                  contextKey.ReadGroup = readGrpCov;
                  contextKey.Qual = qual;

                  bqsrEntry empiricalContextEntry(0,-1,0);

                  if (contextKey.Covariate >= 0) {
                      empiricalContextEntry = this->Contexts.table[getContextCovIdx(contextKey)];
                  }

                  double recalibratedQual = estimateHierarchicalBayesianQuality(epsilon, empiricalReadGroupEntry,
                                                                              empiricalBqEntry,
                                                                              empiricalContextEntry,
                                                                              empiricalCycleEntry);

                  uint8_t recalibratedQualityScore = quantizedQualityScores[std::max(1,
                                                                                std::min(fastRound(recalibratedQual),
                                                                                       maxRecalibratedQualScore))];
                  if (!staticQuantizedScores.empty()) {
                      recalibratedQualityScore = staticQuantizedScores[recalibratedQualityScore];
                  }

                  size_t bufferIdx = buffers.getIdxRecalQualScore(key);
                  buffers.recalibratedQualityScores[bufferIdx] = static_cast<int8_t>(recalibratedQualityScore);

                  reads.GetQUAL(idx)[i] = static_cast<char>(recalibratedQualityScore);
              }
          }
      }*/
}

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
    key.ReadGroup = idx / szQualScore;
    key.Qual = idx % szQualScore;
    key.Covariate = 0;
    return key;
}

bqsrTableKey BaseRecalibratorTables::idxToCycleCovKey(size_t idx) {
    bqsrTableKey key;
    key.ReadGroup = idx / (szQualScore * szCycleCov);
    auto tmpVal = idx % (szQualScore * szCycleCov);
    key.Qual = tmpVal / szCycleCov;
    key.Covariate = tmpVal % szCycleCov;
//    if (key.Covariate >= MAX_CYCLE) {
//        key.Covariate -= MAX_CYCLE;
//        key.Covariate *= -1;
//    }
    return key;
}

bqsrTableKey BaseRecalibratorTables::idxToContextCovKey(size_t idx) {
    bqsrTableKey key;
    key.ReadGroup = idx / (szQualScore * szContextCov);
    auto tmpVal = idx % (szQualScore * szContextCov);
    key.Qual = tmpVal / szContextCov;
    key.Covariate = tmpVal % szContextCov;
    return key;
}

BaseRecalibratorTables::BaseRecalibratorTables(size_t nReadGroups, size_t szQual, size_t maxCycle) : szQualScore(szQual), szCycleCov(maxCycle * 2), nReadGroups(nReadGroups), MAX_CYCLE(maxCycle) {
    bqsrEntry defaultEntry;
    defaultEntry.Observations = -1;
    defaultEntry.Mismatches = 0;
    szContextCov = MaxContextCovariateKey() + 1;

    QualityScores.table.resize(nReadGroups * szQualScore, defaultEntry);
    QualityScores.empirical_qualities.resize(nReadGroups * szQualScore, 0);
    Cycles.table.resize(nReadGroups * szQualScore * szCycleCov, defaultEntry);
    Cycles.empirical_qualities.resize(nReadGroups * szQualScore * szCycleCov, 0);
    Contexts.table.resize(nReadGroups * szQualScore * szContextCov, defaultEntry);
    Contexts.empirical_qualities.resize(nReadGroups * szQualScore * szContextCov, 0);
}

void BaseRecalibratorTables::PrintTables(Header& samHeader) {
    //std::vector<std::string> reverseReadGroupDict(samHeader.getReadGroupRecord().getAll().size());
    std::string readGroupStr;
//    for(auto& readGroup : samHeader.getReadGroupRecord().getAll()) {
//        readGroup.getReadGroupId();
//        reverseReadGroupDict[v] = k;
//    }
    {
        std::ofstream qualityScoreTable("quality-score-table");
        for (size_t qIdx = 0; qIdx < szQualScore; qIdx++) {
            if(this->QualityScores.table[qIdx].Observations!=-1) {
                bqsrTableKey key = idxToQualityScoresKey(qIdx);
                qualityScoreTable << "[" << unsigned(qIdx) << ", " << 0 << ", " << samHeader.getReadGroupRecord().getString(key.ReadGroup) << "] -> ";
                qualityScoreTable << unsigned(qIdx) << " "
                                  << this->QualityScores.table[qIdx].Observations << " "
                                  << this->QualityScores.table[qIdx].Mismatches
                                  << std::endl;
            }
        }
    }
    {
        std::ofstream cyclesTable("cycles-table");
        for (size_t idx = 0; idx < this->Cycles.table.size(); idx++) {
            if(this->Cycles.table[idx].Observations!=-1) {
                bqsrTableKey key = idxToCycleCovKey(idx);
                cyclesTable << "[" << unsigned(key.Qual) << ", " << key.Covariate << ", " << samHeader.getReadGroupRecord().getString(key.ReadGroup) << "] -> ";
                cyclesTable << unsigned(this->Cycles.empirical_qualities[idx]) << " "
                            << this->Cycles.table[idx].Observations << " " << this->Cycles.table[idx].Mismatches
                            << std::endl;
            }
        }
    }

    {
        std::ofstream contextTable("context-table");
        for (size_t idx = 0; idx < this->Contexts.table.size(); idx++) {
            if(this->Contexts.table[idx].Observations!=-1) {
                bqsrTableKey key = idxToContextCovKey(idx);
                contextTable << "[" << unsigned(key.Qual) << ", " << key.Covariate << ", " << samHeader.getReadGroupRecord().getString(key.ReadGroup) << "] -> ";
                contextTable << unsigned(this->Contexts.empirical_qualities[idx]) << " "
                             << this->Contexts.table[idx].Observations << " " << this->Contexts.table[idx].Mismatches
                             << std::endl;
            }
        }
    }
}


uint32_t MaxContextCovariateKey() {
    // the maximum value is T (11 in binary) for each base in the context
    uint32_t length = mismatchesContextSize;  // the length of the context
    uint32_t key = length;
    uint32_t bitOffset = lengthBits;
    for (int i = 0; i <length ; i++) {
        key |= (3 << bitOffset);
        bitOffset += 2;
    }
    return key;
}

void BaseRecalibratorTables::PrintBQSRTables(const std::string &filename, CompressedStringCollection &RG, Header &hdr) {
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

    initializeQuantizedQualityScores(this->QualityScores, quantizationLevel, observations, scores);

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

void BaseRecalibratorTables::printCombinedBQSRTable(std::ostream &file, CompressedStringCollection& RG, Header &hdr) {
    std::vector<combinedBqsrEntry> table(hdr.getReadGroupRecord().size()); // Table mapping from read groups to entries

    initializeCombinedBQSRTable(QualityScores, table);

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
        auto key = this->idxToQualityScoresKey(i);
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
        return RG.IntToString(e1) < RG.IntToString(e2);
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

    std::vector<bqsrTableKey> keys;

    for (size_t i = 0; i < QualityScores.table.size(); i++) {
        if (QualityScores.table[i].Observations > 0) {
            auto key = this->idxToQualityScoresKey(i);
            auto &entry = QualityScores.table[i];
            auto &empiricalQuality = QualityScores.empirical_qualities[i];
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
        auto& entry = QualityScores.table[getQualityScoresIdx(key, szQualScore)];
        auto& empiricalQuality = QualityScores.empirical_qualities[getQualityScoresIdx(key, szQualScore)];
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

    std::vector<covariateBqsrTableKey> keys;

    for (size_t i = 0; i < Cycles.size(); i++) {
        if (Cycles.table[i].Observations > 0) {
            auto key = this->idxToCycleCovKey(i);
            auto entry = Cycles.table[i];
            auto empiricalQuality = Cycles.empirical_qualities[i];
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

    for (size_t i = 0; i < Contexts.size(); i++) {
        if (Contexts.table[i].Observations > 0) {
            auto key = this->idxToContextCovKey(i);
            auto entry = Contexts.table[i];
            auto empiricalQuality = Contexts.empirical_qualities[i];

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
            entry = Cycles.table[getCycleCovIdx(key.key, szQualScore, szCycleCov)];
            empiricalQuality = Cycles.empirical_qualities[getCycleCovIdx(key.key, szQualScore, szCycleCov)];
            name = "Cycle";
        } else {
            entry = Contexts.table[getContextCovIdx(key.key, szQualScore, szCycleCov, szContextCov)];
            empiricalQuality = Contexts.empirical_qualities[getContextCovIdx(key.key, szQualScore, szCycleCov, szContextCov)];
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

