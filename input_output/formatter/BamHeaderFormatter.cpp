#include "BamHeaderFormatter.h"
#include "HeaderFormatter.h"
#include "../utils/ByteFormatter.h"
#include "../utils/Constant.h"

void BamHeaderFormatter::FormatBamHeader(Header &hdr, std::string &line) {
    std::string textHeaderLine;
    std::string referenceSequencesLine;
    ByteFormatter bf(line);

    //Field: magic
    bf.appendString("BAM\001");
    //line.append(BAM_MAGIC);

    HeaderFormatter::FormatSamHeader(hdr, textHeaderLine);

    //Field: l_text
    bf.appendInt32(textHeaderLine.size());

    //Field: text
    if (!textHeaderLine.empty())
        bf.appendString(textHeaderLine);

    //Field: n_ref
    bf.appendInt32(hdr.getSequenceDictionary().getAll().size());

    //Field: for each reference information:
    //l_name, name, l_ref
    for (auto& sequences : hdr.getSequenceDictionary().getAll()) {
        bf.appendInt32(sequences.getSequenceName().size() + 1);
        bf.appendString(sequences.getSequenceName() + '\0'); //Null terminated string
        bf.appendInt32(sequences.getSequenceLength());
    }
}
