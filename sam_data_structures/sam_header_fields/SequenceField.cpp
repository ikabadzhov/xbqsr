#include "SequenceField.h"

//const ssize_t SequenceRecord::RESERVED_SEQUENCE_NAMES[RESERVED_SEQUENCE_NAMES_SIZE] = {
//        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 0-15
//        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // 16-31
//        0, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 32-47, '*' at index 42
//        0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, // 48-63, '=' at index 61
//};
const std::unordered_map<ssize_t, std::string> SequenceRecord::RESERVED_SEQUENCE_NAMES_VALS = {{-1, "="}, {-2, "*"}};

const std::regex SequenceRecord::LEGAL_RNAME_PATTERN("[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*");

/**
 * SEQUENCE DICTIONARY METHODS
 * */

void SequenceDictionary::add(SequenceRecord &sequence) {
    this->sequences.emplace_back(sequence);
    this->seqToIdx[sequence.getSequenceName()] = idxToSeq.size();
    this->idxToSeq.emplace_back(sequence.getSequenceName());
}

std::vector<SequenceRecord> &SequenceDictionary::getAll() {
    return sequences;
}

ssize_t SequenceDictionary::getIdx(const std::string& sequence_name){
    for (uint8_t i = 0; i < RESERVED_SEQUENCE_NAMES_SIZE; i++) {
        if (static_cast<char>(i) == sequence_name[0]) {
            return SequenceRecord::RESERVED_SEQUENCE_NAMES[i];
        }
    }
    return this->seqToIdx.at(sequence_name);
}

const std::string& SequenceDictionary::getSequenceName(ssize_t id) {
    for (auto& vals : SequenceRecord::RESERVED_SEQUENCE_NAMES_VALS) {
        if (vals.first == id) {
            return SequenceRecord::RESERVED_SEQUENCE_NAMES_VALS.at(vals.first);
        }
    }
    return sequences[id].getSequenceName();
}

bool SequenceDictionary::containsSequence(const std::string &sequence_name) {
    return this->seqToIdx.contains(sequence_name);
}


/**
 * SEQUENCE RECORD METHODS
 * */
void SequenceRecord::truncateSequenceName(std::string &sequenceName) {
    size_t truncateAt = sequenceName.length();
    const char WHITESPACE_CHARS[] = {' ', '\t', '\n', '\v', '\f', '\r'};

    for (char c : WHITESPACE_CHARS) {
        auto index = sequenceName.find(c);
        if (index != std::string::npos && index < truncateAt) {
            truncateAt = index;
        }
    }
    sequenceName = sequenceName.substr(0, truncateAt);
}

void SequenceRecord::setSequenceName(const std::string &name) {
    sequenceName = name;
}

void SequenceRecord::setSequenceLength(int value) {
    length = value;
}

int32_t SequenceRecord::getSequenceLength() const {
    return length;
}

const std::string &SequenceRecord::getSequenceName() const  {
    return sequenceName;
}

void SequenceRecord::validateSequenceName(const std::string &name) {
    if (!std::regex_match(name, LEGAL_RNAME_PATTERN)) {
        throw std::invalid_argument("Sequence name doesn't match regex");
    }
}

std::string SequenceRecord::getAssembly() const {
    return getAttribute("AS");
}

void SequenceRecord::setAssembly(const std::string& value) {
    putAttribute("AS", value);
}

std::string SequenceRecord::getSpecies() const {
    return getAttribute("SP");
}

void SequenceRecord::setSpecies(const std::string& value) {
    putAttribute("SP", value);
}

std::string SequenceRecord::getMd5() const {
    return getAttribute("M5");
}

void SequenceRecord::setMd5(const std::string& value) {
    putAttribute("M5", value);
}

std::string SequenceRecord::getDescription() const {
    return getAttribute("DS");
}

void SequenceRecord::setDescription(const std::string& value) {
    return putAttribute("DS", value);
}

void SequenceRecord::addAlternativeSequenceName(const std::string &name) {
    if (!alternativeSequenceName.contains(name)) {
        validateAltRegExp(name);
        alternativeSequenceName.emplace(name);
    }
}


std::unordered_set<std::string> SequenceRecord::getAlternativeSequenceName() const {
    return alternativeSequenceName;
}


void SequenceRecord::validateAltRegExp(const std::string& name) {
    if (!std::regex_match(name, LEGAL_RNAME_PATTERN)) {
        throw std::invalid_argument("Invalid alternative sequence name '" + name + "': does not match the pattern");
    }
}

bool SequenceRecord::hasAlternativeSequenceNames() {
    return !alternativeSequenceName.empty();
}

bool SequenceRecord::isSameSequence(const SequenceRecord& that) const {
    if (length != 0 && that.length != 0 && length != that.length) {
        return false;
    }

    std::string thisMd5 = getAttribute("M5");
    std::string thatMd5 = that.getAttribute("M5");

    if (!thisMd5.empty() && !thatMd5.empty()) {
        // Assuming you have a way to compare MD5 strings or convert them to comparable values
        if (thisMd5 != thatMd5) {
            return false;
        }
    } else if (sequenceName != that.sequenceName) {
        if (!this->getAlternativeSequenceName().contains(that.sequenceName) && !that.getAlternativeSequenceName().contains(sequenceName)) {
            return false;
        }
        return true;
    }
    return true;
}

void SequenceRecord::setAlternativeSequenceName(std::unordered_set<std::string> &alternativeSequences) {
    this->alternativeSequenceName = alternativeSequences;
}

bool SequenceRecord::isEmpty(int32_t val) {
    return val==UNAVAILABLE_SEQUENCE_INDEX;
}
