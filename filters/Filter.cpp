#include <unordered_map>
#include "Filter.h"
#include "../sam_data_structures/sam_header_fields/SequenceField.h"
#include "../file_format/Global.h"

using namespace std;

std::function<bool (SamAlignments&, size_t)> Filter::AddREFID(SamAlignments& reads) {

    reads.TMP.REFID.resize(reads.GetAlnsNum());
    reads.TMP.nextREFID.resize(reads.GetAlnsNum());

    return [=](SamAlignments& sam, size_t idx){


        //if (dictTable.count(std::string{sam.RNAME[idx]})<=0){
          //  value=-1;
        //}
        //else {
          //  value = dictTable.at(std::string{sam.RNAME[idx]});
        //}
        int32_t value = sam.RNAME[idx];
        sam.TMP.REFID[idx] = value;

        if( sam.RNEXT[idx] == SequenceRecord::RESERVED_SEQUENCE_NAMES['=']){
            sam.TMP.nextREFID[idx] = value;
        } else {
            sam.TMP.nextREFID[idx] = SequenceRecord::UNAVAILABLE_SEQUENCE_INDEX;
        }
        return true;
    };
}

std::function<bool(SamAlignments &, size_t)>
Filter::AddNewTAG(TagKeyAndType tagKeyAndType, const std::string &tagValue) {
    return [=](SamAlignments& sam, size_t idx){

        OptionalFields::addTag(sam.TAGS, idx, tagKeyAndType, tagValue);

        return true;
    };
}
