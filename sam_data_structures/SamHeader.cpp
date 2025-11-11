#include <iostream>
#include "SamHeader.h"


std::string SortingOrderToString(SortingOrder so) {
    if (so == SortingOrder::Keep){
        return "Keep";
    }
    else if (so == SortingOrder::Unknown) {
        return "Unknown";
    }
    else if (so == SortingOrder::Unsorted) {
        return "Unsorted";
    }
    else if (so == SortingOrder::Queryname) {
        return "Queryname";
    }
    else if (so == SortingOrder::Coordinate) {
        return "Coordinate";
    }
    else{
        std::cerr << "Error in conversion for string to SortingOrder" << std::endl;
        exit(-1);
    }
};

inline bool iequals(const std::string& a, const std::string& b)
{
    size_t sz = a.size();
    if (b.size() != sz)
        return false;
    for (size_t i = 0; i < sz; ++i)
        if (tolower(a[i]) != tolower(b[i]))
            return false;
    return true;
}

SortingOrder StringToSortingOrder(const std::string& so){
    if(iequals(so, "Keep")){
        return SortingOrder::Keep;
    }
    else if(iequals(so,"Unknown")){
        return SortingOrder::Unknown;
    }
    else if(iequals(so,"Unsorted")){
        return SortingOrder::Unsorted;
    }
    else if(iequals(so, "Queryname")){
        return SortingOrder::Queryname;
    }
    else if(iequals(so,"Coordinate")){
        return SortingOrder::Coordinate;
    }else{
        std::cerr<<"Error in conversion for string to SortingOrder"<<std::endl;
        exit(-1);
    };
}

GroupingOrder StringToGroupingOrder(const std::string& go){
    if (go == "None") {
        return GroupingOrder::None;
    }
    else if (go == "Query") {
        return GroupingOrder::Query;
    }
    else if (go == "Reference") {
        return GroupingOrder::Reference;
    }
    else {
        std::cerr<<"Error in conversion for string to GroupingOrder"<<std::endl;
        exit(-1);
    };
}

std::string GroupingOrderToString(GroupingOrder go){
    if(go==GroupingOrder::None){
        return "None";
    }
    else if(go==GroupingOrder::Query){
        return "Query";
    }
    else if(go==GroupingOrder::Reference){
        return "Reference";
    }
    else{
        std::cerr<<"Error in conversion for string to GroupingOrder"<<std::endl;
        exit(-1);
    };
}


void Header::addProgramRecord(ProgramRecord &program) {
    this->programGroups.add(program);
}

void Header::addComment(const std::string& line) {
    this->comments.add(line);
}

void Header::addSequenceRecord(SequenceRecord &sequence) {
    this->sequenceDictionary.add(sequence);
}

void Header::addReadGroupRecord(const ReadGroupRecord &readGroup) {
    this->readGroups.add(readGroup);
}

int32_t Header::covKeyForRGId(int32_t rg_id) const {
    if (rg_id < 0 || static_cast<size_t>(rg_id) >= rgIdToCovKey_.size()) return -1;
    return rgIdToCovKey_[static_cast<size_t>(rg_id)];
}

const std::string & Header::covLabelForKey(int32_t key) const {
    return rgCov_.getString(static_cast<size_t>(key));
}

std::string Header::covLabelForRGId(int32_t rg_id) const {
    const int32_t k = covKeyForRGId(rg_id);
    return (k < 0) ? std::string{} : rgCov_.getString(static_cast<size_t>(k));
}

QnameStringMapper &Header::QnameMapper() {
    return qnameMapper;
}

ProgramGroups& Header::getProgramGroup() {
    return programGroups;
}

ReadGroupsDictionary & Header::getReadGroupRecord() noexcept { return readGroups; }

const ReadGroupsDictionary & Header::getReadGroupRecord() const noexcept { return readGroups; }

Comments& Header::getComments() {
    return comments;
}

SequenceDictionary &Header::getSequenceDictionary() {
    return sequenceDictionary;
}

const SequenceDictionary & Header::getSequenceDictionary() const noexcept {
    return sequenceDictionary;
}

void Header::buildReadGroupCovariateIndex() {
    std::unique_lock lk(mtx);

    // Rebuild safely (idempotent)
    rgCov_.clear();
    rgIdToCovKey_.clear();
    rgIdToCovKey_.resize(static_cast<size_t>(readGroups.size()), -1);

    const auto& all = readGroups.getAll();
    for (size_t rg_id = 0; rg_id < all.size(); ++rg_id) {
        const auto& rg = all[rg_id];

        // Prefer PU; fallback to RG ID
        std::string ident;
        if (rg.containsAttribute("PU")) {
            ident = rg.getAttribute("PU");
            if (ident.empty()) ident = rg.getReadGroupId();   // guard empty PU
        } else {
            ident = rg.getReadGroupId();
        }

        // Insert (StringDictionary::add should dedupe & return stable key)
        const int32_t key = rgCov_.add(ident);
        rgIdToCovKey_[rg_id] = key;
    }
}

FileMetadata &Header::getHD() {
    return attributesHD;
}

void FileMetadata::put(const std::string &key, const std::string &value) {
    attributes.put(key, value);
}

std::string FileMetadata::get(const std::string& key) {
    return attributes.get(key);
}

void Comments::add(const std::string &line) {
    this->comments.emplace_back(line);
}

std::vector<std::string> &Comments::getAll() {
    return comments;
}
