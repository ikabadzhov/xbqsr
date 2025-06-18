#include "BamReader.h"
#include "../parser/SamHeaderParser.h"
#include "../parser/BamAlignmentParser.h"

/*
 * This function is designed to parse the header utilizing the nested BAMUtilsIO
 * function. Subsequently, upon extracting the information, a validation function
 * is called to verify the accuracy of the quantity of information
 * (not the content itself), ensuring that all header lines have been
 * successfully read.
*/
void BamReader::ParseHeader(Header &hdr) {

    std::string headerData(sizeof(uint32_t), 0);
    uint32_t lengthToRead = 0;
    size_t numBytesRead = in->read(headerData, 0, sizeof(uint32_t));
    assert(numBytesRead == sizeof(uint32_t));
    std::memcpy(&lengthToRead, headerData.data(), sizeof(uint32_t));

    headerData.resize(lengthToRead, 0);
    numBytesRead = in->read(headerData, 0, lengthToRead);
    assert(numBytesRead == lengthToRead);
    size_t pos = 0;

    while (pos < headerData.size()) {
        auto endPos = std::find(headerData.begin() + pos, headerData.end(), '\n');

        // Check if a newline character was found
        if (endPos != headerData.end()) {
            // Extract the line and add it to the lines vector
            std::string line(headerData.begin() + pos, endPos);
            SamHeaderParser::ParseSamHeaderFromSlice(hdr, line);
            // Move pos to the character after the found newline
            pos = std::distance(headerData.begin(), endPos) + 1;
        } else {
            break; // Exit the loop
        }
    }

    headerData.clear();
    //headerData.resize(pos);
    headerData.resize(sizeof(uint32_t), 0);

    uint32_t checkNREF = 0;
    numBytesRead = in->read(headerData, 0, sizeof(uint32_t));
    assert(numBytesRead == sizeof(uint32_t));
    std::memcpy(&checkNREF, headerData.data(), sizeof(uint32_t));
    /**
     *
     * Save these in the header
     * */
    for(int i=0; i<checkNREF; i++) {
        uint32_t tempLengthToRead = 0;
        numBytesRead = in->read(headerData, 0, sizeof(uint32_t));
        assert(numBytesRead == sizeof(uint32_t));
        std::memcpy(&tempLengthToRead, headerData.data(), sizeof(uint32_t));

        std::string tempBuffer(tempLengthToRead, 0);
        numBytesRead = in->read(tempBuffer, 0, tempLengthToRead);

        std::string_view name(tempBuffer.data(), tempBuffer.data() + tempLengthToRead);
        //auto null_terminator_pos = std::find(name.begin(), name.end(), '\0');

        uint32_t tempLengthValue = 0;
        numBytesRead = in->read(headerData, 0, sizeof(uint32_t));
        assert(numBytesRead == sizeof(uint32_t));
        std::memcpy(&tempLengthValue, headerData.data(), sizeof(uint32_t));

        if (!hdr.getSequenceDictionary().getAll().empty()) {
            std::string_view hdrSN = hdr.getSequenceDictionary().getSequenceName(i);
            assert(hdrSN.compare(name.data()) == 0);
        } else {
            SequenceRecord sr(name.data(), tempLengthValue);
            hdr.getSequenceDictionary().add(sr);
        }
    }
}

/*
 * Once the file analysis is complete, there's no requirement to maintain the stream
 * and the gzFile in an open state. Closing them aids in conserving resources.
*/
void BamReader::Close() {
    this->in->close();
}

ssize_t BamReader::ParseAlignment(SamAlignments &reads) {
    if (buffer.is_empty()) {
        buffer.fill(in);
    }

    std::string remaining;
    size_t alnsRead = 0;
    size_t numBytesRead = 0;
    //size_t bytesToRead = 0;
    uint32_t recordLength = 0;

    while (buffer.has_data()) {

        if (alnsRead == reads.GetAlnsNum()) {
            reads.Resize(reads.GetAlnsNum() + BATCH_SIZE);
        }

        // Check if we have enough data to read the record length
        if (remaining.size() + buffer.available_data() < sizeof(int32_t)) {
            std::string_view slice = buffer.getslice(static_cast<size_t>(sizeof(int32_t)), numBytesRead);
            remaining.append(slice);
            buffer.fill(in);
            continue;
        }

        /*
         * If the remaining string (from previous reading) is not enough yet to read a uint32_t size
         * */
        if (remaining.size() < sizeof(int32_t)) {
            std::string_view slice = buffer.getslice(sizeof(int32_t) - remaining.size(), numBytesRead);
            remaining.append(slice);
           // recordLength = reinterpret_cast<uint32_t*>(remaining.data());
            std::memcpy(&recordLength, remaining.c_str(), sizeof(recordLength));
            remaining.clear();
        } else {
            std::memcpy(&recordLength, remaining.c_str(), sizeof(recordLength));
        }

        while (buffer.available_data() < recordLength) {
            remaining.append(buffer.getslice(recordLength, numBytesRead));
            recordLength -= numBytesRead;
            buffer.fill(in);
            //if the buffer is empty even after filling it up, that means there is nothing else to read
            if(buffer.is_empty())
                break;
        }

        std::string_view slice = buffer.getslice(recordLength, numBytesRead);

        try {
            if (!remaining.empty()) {
                remaining.append(slice);
                // Process the record
                BamAlignmentParser::ParseBamAlignmentFromSlice(remaining, reads, alnsRead);
                alnsRead++;
                remaining.clear();
            } else {
                // Process the record
                BamAlignmentParser::ParseBamAlignmentFromSlice(slice, reads, alnsRead);
                alnsRead++;
            }
        } catch (const std::runtime_error& e) {
            std::cerr << "EOF condition encountered!" << std::endl;
        }

        if (!buffer.has_data()) {
            buffer.fill(in);
        }
    }

    reads.Resize(alnsRead);
    reads.GetHdr().QnameMapper().finalize(reads.QNAME);
    return static_cast<ssize_t>(alnsRead);
}

std::string_view BamReader::getslice(size_t numBytesToRead, size_t &numBytesRead) {
    return buffer.getslice(numBytesToRead, numBytesRead);
}


