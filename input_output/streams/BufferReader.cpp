#include "BufferReader.h"

void BufferReader::fill(std::unique_ptr<InputStream> &in) {
    pos = 0;
    count = in->read(buffer, pos, buffer.size() - pos);
}

std::string_view BufferReader::getslice(char delim, bool& isDelimFound) {
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

std::string_view BufferReader::getslice(size_t numBytesToRead, size_t& numBytesRead) {
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

char BufferReader::peek() {
    if (pos >= count) {
        return -1;  // Return EOF if no more data
    }
    return static_cast<unsigned char>(buffer[pos]); // Return the next character without advancing pos
}

[[maybe_unused]] bool BufferReader::is_empty() const {
    if (pos >= count) {
        return true;
    }
    return false;
}

[[maybe_unused]] bool BufferReader::has_data() const {
    if (count > 0 && pos < count) {
        return true;
    }
    return false;
}

size_t BufferReader::available_data() {
    return count - pos;
}
