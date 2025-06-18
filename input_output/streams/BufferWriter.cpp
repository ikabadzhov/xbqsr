#include "BufferWriter.h"

bool BufferWriter::addSlice(std::string content){
    if(content.size()<=(size-pos)) {
        memcpy(&buffer[pos], content.data(), content.size());
        this->pos += content.size();
        //if(pos>count) //otherwise it means we flushed and we are overwriting
        this->count += content.size();
        return true;
    }
    return false;
}

bool BufferWriter::addNewLine(){
    if(remaining_space()>sizeof(std::to_string('\n'))) {
        buffer[pos] = '\n';
        this->pos += 1;
        this->count += 1;
        return true;
    }
    return false;
}

std::string_view BufferWriter::getslice(char delim, bool& isDelimFound) {
    if (pos >= count) {
        isDelimFound = false;
        return {}; // Return an empty string_view if no more data
    }
    // Find the delim in the line or the end of the buffer
    auto lineEnd = static_cast<ssize_t>(buffer.find(delim, pos));
    if (lineEnd != std::string::npos) { // Delim found
        ssize_t oldPos = pos;
        pos = lineEnd + 1; // Move past the newline for the next call
        isDelimFound = true;
        return std::string_view{buffer.data() + oldPos, static_cast<std::string_view::size_type>(lineEnd - oldPos)};
    } else if (pos < count) {
        // Delim not found but some chars are available
        // Return this first part of the line
        ssize_t oldPos = pos;
        pos = count;
        isDelimFound = false;
        return std::string_view{buffer.data() + oldPos, static_cast<std::string_view::size_type>(count - oldPos)};
    } else {
        // Delim not found and no chars avaialble in the buffer
        // I think it should never arrive here
        isDelimFound = false;
        return {};
    }
}

size_t BufferWriter::getLastCount(){
    return count;
}

size_t BufferWriter::getLastPos(){
    return pos;
}

size_t BufferWriter::getSize(){
    return buffer.size();
}

std::string_view BufferWriter::getslice(size_t numBytesToRead, size_t& numBytesRead) {
    if (pos >= count) {
        numBytesRead = 0;
        return {}; // Return an empty string_view if no more data
    }

    // Calculate the number of bytes that can be read
    numBytesRead = std::min(numBytesToRead, count - pos);

    // Update the position
    size_t oldPos = pos;
    pos += numBytesRead;

    // Create a string_view for the slice
    return {buffer.data() + oldPos, numBytesRead};
}

void BufferWriter::getContent(size_t startPos, size_t length, std::string &value) {
    if(startPos+length<=buffer.size()) { //startPos>=0 useless, size_t it's >=0
        //value.append(buffer.substr(startPos, length));
        if (value.capacity() < startPos + length){
            //buffer.resize(getLastPos());
            value.resize(startPos + length);
        }
        //value.append(buffer.data() + startPos, length);
        value.replace(0, length, buffer.data() + startPos, length);
    }
}

void BufferWriter::flushContent(){
    //buffer = std::string(); No need, just overwriting existing values
    pos=0;
    //count=0;
}

/*std::string_view BufferReader::getslice(size_t numBytesToRead, size_t& numBytesRead) {
    if (pos >= count) {
        numBytesRead = 0;
        return {}; // Return an empty string_view if no more data
    }

    // Calculate the number of bytes that can be read
    numBytesRead = std::min(numBytesToRead, count - pos);

    // Update the position
    size_t oldPos = pos;
    pos += numBytesRead;

    // Create a string_view for the slice
    return {buffer.data() + oldPos, numBytesRead};
}*/

char BufferWriter::peek() {
    if (pos >= count) {
        return -1;  // Return EOF if no more data
    }
    return static_cast<unsigned char>(buffer[pos]); // Return the next character without advancing pos
}

[[maybe_unused]] bool BufferWriter::is_empty() const {
    if (pos >= count) {
        return true;
    }
    return false;
}

[[maybe_unused]] bool BufferWriter::has_data() const {
    if (count > 0 && pos < count && pos != 0) { //TODO: To verify
        return true;
    }
    return false;
}

size_t BufferWriter::available_data() {
    return size - pos;
}

size_t BufferWriter::remaining_space() {
    return this->size-count;
}
