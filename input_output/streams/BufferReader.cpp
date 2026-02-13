#include "BufferReader.h"

void BufferReader::fill(std::unique_ptr<InputStream> &in) {
    pos = 0;
    count_ = in->read(buffer, pos, buffer.size() - pos);
}

std::string_view BufferReader::getslice(char delim, bool& isDelimFound) {
    if (pos >= count_) {
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
    } else if (pos < count_) {
        // Delim not found but some chars are available
        // Return this first part of the line
        ssize_t oldPos = pos;
        pos = count_;
        isDelimFound = false;
        return std::string_view{buffer.data() + oldPos, static_cast<std::string_view::size_type>(count_ - oldPos)};
    } else {
        // Delim not found and no chars avaialble in the buffer
        // I think it should never arrive here
        isDelimFound = false;
        return {};
    }
}

std::string_view BufferReader::getslice(size_t numBytesToRead, size_t& numBytesRead) {
    if (pos >= count_) {
        numBytesRead = 0;
        return {}; // Return an empty string_view if no more data
    }

    // Calculate the number of bytes that can be read
    numBytesRead = std::min(numBytesToRead, count_ - pos);

    // Update the position
    size_t oldPos = pos;
    pos += numBytesRead;

    // Create a string_view for the slice
    return {buffer.data() + oldPos, numBytesRead};
}

std::string_view BufferReader::getslicerange(size_t begin, size_t end) {
    if (begin > end) {
        throw std::out_of_range("BufferReader::getslice: begin > end");
    }
    if (end > count_) {
        throw std::out_of_range("BufferReader::getslice: end beyond valid data");
    }
    // Create a string_view for the slice
    size_t sz = end - begin;
    return {buffer.data() + begin, sz};
}

/*std::string_view BufferReader::getslice(size_t numBytesToRead, size_t& numBytesRead) {
    if (pos >= count_) {
        numBytesRead = 0;
        return {}; // Return an empty string_view if no more data
    }

    // Calculate the number of bytes that can be read
    numBytesRead = std::min(numBytesToRead, count_ - pos);

    // Update the position
    size_t oldPos = pos;
    pos += numBytesRead;

    // Create a string_view for the slice
    return {buffer.data() + oldPos, numBytesRead};
}*/

char BufferReader::peek() {
    if (pos >= count_) {
        return -1;  // Return EOF if no more data
    }
    return static_cast<unsigned char>(buffer[pos]); // Return the next character without advancing pos
}

[[maybe_unused]] bool BufferReader::is_empty() const {
    if (pos >= count_) {
        return true;
    }
    return false;
}

[[maybe_unused]] bool BufferReader::has_data() const {
    if (count_ > 0 && pos < count_) {
        return true;
    }
    return false;
}

size_t BufferReader::available_data() {
    return count_ - pos;
}
