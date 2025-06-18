#include "SamHeaderFormatter.h"
#include "../streams/BufferWriter.h"

void SamHeaderFormatter::FormatSamHeader(Header &hdr, std::string& line) {
    writeHDLine(hdr, line);
    line.append(1,'\n');
    for (auto& sequenceRecord : hdr.getSequenceDictionary().getAll()) {
        getSQLine(sequenceRecord, line);
        line.append(1,'\n');
    }
    for (int32_t id = 0; id < hdr.getReadGroupRecord().size(); id++) {
        getRGLine(hdr.getReadGroupRecord().getReadGroupRecord(id), line);
        line.append(1,'\n');
    }
    for (auto& programRecord : hdr.getProgramGroup().getAll()) {
        getPGLine(programRecord, line);
        line.append(1,'\n');
    }
    for (auto& comment : hdr.getComments().getAll()) {
        line.append(comment);
        line.append(1,'\n');
    }
}

void SamHeaderFormatter::writeHDLine(Header &hdr, std::string& line) {
    std::vector<std::string> fields(hdr.getHD().getAllAttributes().size() + 1);
    fields[0].append("@HD");
    encodeTags(hdr.getHD(), fields, 1);

    join('\t', fields, line);
}


void SamHeaderFormatter::getSQLine(const SequenceRecord& sequenceRecord, std::string& line) {
    size_t numAttributes = sequenceRecord.getAllAttributes().size();
    std::vector<std::string> fields(3 + numAttributes);

    fields[0] = "@SQ";
    fields[1] = "SN:" + sequenceRecord.getSequenceName();
    fields[2] = "LN:" + std::to_string(sequenceRecord.getSequenceLength());

    encodeTags(sequenceRecord, fields, 3);
    join('\t', fields, line);
}

void SamHeaderFormatter::getPGLine(const ProgramRecord& program, std::string& line) {
    size_t numAttributes = program.getAllAttributes().size();
    std::vector<std::string> fields(2 + numAttributes);

    fields[0] = "@PG";
    fields[1] = "ID:" + program.getProgramId();

    encodeTags(program, fields, 2);
    join('\t', fields, line);
}

void SamHeaderFormatter::getRGLine(const ReadGroupRecord& readGroupRecord, std::string& line) {
    size_t numAttributes = readGroupRecord.getAllAttributes().size();
    std::vector<std::string> fields(2 + numAttributes);

    fields[0] = "@RG";
    fields[1] = "ID:" + readGroupRecord.getReadGroupId();

    encodeTags(readGroupRecord, fields, 2);
    join('\t', fields, line);
}

void SamHeaderFormatter::encodeTags(const SamHeaderRecord& record, std::vector<std::string>& fields, size_t offset) {
    record.getAllAttributes().for_each([&](const std::string& key, const std::string& value){
        fields[offset].append(key);
        fields[offset].append(":");
        fields[offset].append(value);
        offset++;
    });
}

void SamHeaderFormatter::join(char sep, std::vector<std::string>& fields, std::string& line) {
    bool first_field = true;
    for (auto& s : fields) {
        if (!first_field) {
            line.append(1, sep);
        }
        line.append(s);
        first_field = false;
    }
}
