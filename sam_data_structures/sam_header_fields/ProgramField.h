#ifndef GENOMICPIPELINE_PROGRAMFIELD_H
#define GENOMICPIPELINE_PROGRAMFIELD_H


#include <string>
#include "../../data_types/LinkedHashMap.h"
#include "SamHeaderRecord.h"

class ProgramRecord : public SamHeaderRecord{
    std::string programGroupId;
public:
    explicit ProgramRecord(std::string progId) : SamHeaderRecord("PG"), programGroupId(std::move(progId)) {}
    const std::string & getProgramId() const;
};

class ProgramGroups {
    std::vector<ProgramRecord> programs;
    std::unordered_map<std::string, int32_t> progToIdx;
    std::vector<std::string> idxToProg;
public:
    void add(ProgramRecord &record);
    std::vector<ProgramRecord>& getAll();
};

#endif //GENOMICPIPELINE_PROGRAMFIELD_H
