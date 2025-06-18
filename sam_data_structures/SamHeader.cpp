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

void Header::addReadGroupRecord(ReadGroupRecord& readGroup) {
    this->readGroups.add(readGroup);
}

QnameStringMapper &Header::QnameMapper() {
    return qnameMapper;
}

ProgramGroups& Header::getProgramGroup() {
    return programGroups;
}

ReadGroupsDictionary &Header::getReadGroupRecord() {
    return readGroups;
}

Comments& Header::getComments() {
    return comments;
}

SequenceDictionary &Header::getSequenceDictionary() {
    return sequenceDictionary;
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
