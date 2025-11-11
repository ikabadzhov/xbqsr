#ifndef GENOMICPIPELINE_SEQUENCEFIELD_H
#define GENOMICPIPELINE_SEQUENCEFIELD_H

#include <string>
#include <utility>
#include <unordered_set>
#include <regex>
#include "SamHeaderRecord.h"
#include "../../data_types/LinkedHashMap.h"
#include <sycl/sycl.hpp>

const int RESERVED_SEQUENCE_NAMES_SIZE = 64;

class SequenceRecord : public SamHeaderRecord {
    std::string sequenceName;
    int32_t length{};
    std::unordered_set<std::string> alternativeSequenceName;
public:
    static const int64_t serialVersionUID = 1L;
    static const int32_t UNAVAILABLE_SEQUENCE_INDEX = -3;
    static const int32_t UNKNOWN_SEQUENCE_LENGTH = 0;
   // static const std::unordered_map<std::string, ssize_t> RESERVED_SEQUENCE_NAMES;
    constexpr static ssize_t RESERVED_SEQUENCE_NAMES[RESERVED_SEQUENCE_NAMES_SIZE] = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 0-15
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 16-31
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, // 32-47, '*' at index 42
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, // 48-63, '=' at index 61
    };
    static const std::unordered_map<ssize_t, std::string> RESERVED_SEQUENCE_NAMES_VALS;

    static const std::regex LEGAL_RNAME_PATTERN;

    explicit SequenceRecord(const std::string &name, int sequenceLength = UNKNOWN_SEQUENCE_LENGTH)
            : SamHeaderRecord("SQ"), sequenceName(name), length(sequenceLength) {
        validateSequenceName(name);
    }

    const std::string &getSequenceName() const;
    void setSequenceName(const std::string &name);

    int32_t getSequenceLength() const;
    void setSequenceLength(int value);

    std::string getAssembly() const;
    void setAssembly(const std::string &value);

    std::string getSpecies() const;
    void setSpecies(const std::string &valu);

    std::string getMd5() const;
    void setMd5(const std::string &valu);

    std::string getDescription() const;
    void setDescription(const std::string &valu);

    std::unordered_set<std::string> getAlternativeSequenceName() const;
    void addAlternativeSequenceName(const std::string &name);
    void setAlternativeSequenceName(std::unordered_set<std::string>& alternativeSequences);
    bool hasAlternativeSequenceNames();

    SYCL_EXTERNAL static bool isEmpty(int32_t val);

private:
    static void validateSequenceName(const std::string &name);
    static void validateAltRegExp(const std::string &basicString);
    bool isSameSequence(const SequenceRecord &that) const;
    static void truncateSequenceName(std::string& sequenceName);

};

class SequenceDictionary {
    std::vector<SequenceRecord> sequences;
    std::unordered_map<std::string, ssize_t> seqToIdx;
    std::vector<std::string> idxToSeq;
public:
    void add(SequenceRecord& sequence);

    std::vector<SequenceRecord>& getAll();
    ssize_t getIdx(const std::string &sequence_name);

    bool hasIndex(ssize_t idx) const {
        return idx >= 0 && static_cast<size_t>(idx) < sequences.size();
    }

    const SequenceRecord& byIndex(ssize_t idx) const {
        return sequences[static_cast<size_t>(idx)];
    }

    const std::string& getSequenceName(ssize_t id);

    bool containsSequence(const std::string &sequence_name);
};

#endif //GENOMICPIPELINE_SEQUENCEFIELD_H
