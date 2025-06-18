#ifndef GENOMICPIPELINE_SAMHEADER_H
#define GENOMICPIPELINE_SAMHEADER_H

#include <vector>
#include <map>
#include "../data_types/RGStringMapper.h"
#include "../data_types/QnameStringMapper.h"
#include "sam_header_fields/ReadGroupField.h"
#include "sam_header_fields/SequenceField.h"
#include "sam_header_fields/ProgramField.h"
#include "sam_header_fields/FileMetadata.h"

// The SAM file format version and date strings supported by this
// library. This is entered by default in an @HD line in the header
// section of a SAM file, unless user code explicitly asks for a
// different version number. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3.
extern const std::string FileFormatVersion;// = "1.6";
extern const std::string FileFormatDate;//    = "22 May 2018";



enum SortingOrder {
    Keep,
    Unknown,
    Unsorted,
    Queryname,
    Coordinate
};

enum GroupingOrder {
    None,
    Query,
    Reference
};

class Comments {
    std::vector<std::string> comments;
public:
    void add(const std::string &line);
    std::vector<std::string>& getAll();
};

/*
 * The class Header stores the information read in a sam/bam header file.
 * We arrange those information as follows.
 * Every line in a header file has a header record type (HD, SQ, RG, PG, CO)
 * and some attributes (mandatory and/or optional).
 *
 * The line whose record type is HD (there can be only one HD line) is stored in the class FileMetadata.
 *
 * All the lines whose record type is SQ are stored using the class SequenceRecord and are grouped in the class
 * variable sequenceDictionary;
 *
 * All the lines whose record type is RG are stored using the class ReadGroupRecord and are grouped in the class
 * ReadGroupsDictionary.
 *
 * All the lines whose record type is PG are stored using the class ProgramRecord and are grouped in the class
 * ProgramRecordDictionary.
 *
 * To sum up:
 *
 * HD -> FileMetadata class
 * SQ -> SequenceRecord class
 * RG -> ReadGroupRecord
 * PG -> ProgramRecord
 * CO -> simple string
 *
 * All the class we defined above hinerit from a common, generic class called SamHeaderRecord
 *
 * Many SQ -> SequenceDictionary class
 * Many RG -> ReadGroupsDictionary class
 * Many PG -> ProgramGroups class
 * Many CO -> Comments class
 *
 * As every line has severa attributes, those are stored as class variables
 * of the corresponding record type class. For example, the attributes (tags)
 * of the line @SQ, which are SN, LN, AN, etc are stored as class variables.
 *
 * Moreover, since these attributes are more frequently accessed than other,
 * those are explicitly saved as class variables of the specialized class
 * (e.g. SequenceRecord will have SN as sequenceName variable, etc), while
 * the rest is pushed up in the parent class and save in a hash map.
 *
 * */
class Header {

    FileMetadata attributesHD;

    // Save the lines @SQ of a sam/bam header file
    SequenceDictionary sequenceDictionary;

    // Save the lines @RG of a sam/bam header file
    ReadGroupsDictionary readGroups;

    // Save the lines @PG of a sam/bam header file
    ProgramGroups programGroups;

    // The @CO lines in the order they occur in the header.
    Comments comments;

public:
    Header() = default;

    void addProgramRecord(ProgramRecord& program);
    void addReadGroupRecord(ReadGroupRecord& program);
    void addComment(const std::string& line);
    void addSequenceRecord(SequenceRecord& sequence);

    FileMetadata& getHD();
    ProgramGroups& getProgramGroup();
    ReadGroupsDictionary& getReadGroupRecord();
    Comments& getComments();
    SequenceDictionary& getSequenceDictionary();

private:
    QnameStringMapper qnameMapper;
public:
    QnameStringMapper& QnameMapper();
};


/*int32_t SQLN(HeaderFieldMap& record);

std::string SortingOrderToString(SortingOrder so);

SortingOrder StringToSortingOrder(const std::string& so);

GroupingOrder StringToGroupingOrder(const std::string& go);

std::string GroupingOrderToString(GroupingOrder go);
*/
#endif //GENOMICPIPELINE_SAMHEADER_H
