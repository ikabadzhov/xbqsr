
#include "ProgramField.h"

void ProgramGroups::add(ProgramRecord &record) {
    programs.emplace_back(record);
    idxToProg.emplace_back(record.getProgramId());
    progToIdx[record.getProgramId()] = static_cast<int32_t>(idxToProg.size());
}

std::vector<ProgramRecord>& ProgramGroups::getAll() {
    return programs;
}

const std::string & ProgramRecord::getProgramId() const {
    return programGroupId;
}
