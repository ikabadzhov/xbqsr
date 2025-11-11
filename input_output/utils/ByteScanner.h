#pragma once
#include <string_view>
#include <cstdint>
#include <cstddef>
#include <cstring>
#include <string>
#include <stdexcept>

class ByteScanner {
    std::string_view bytes_;
    std::size_t index_ = 0;

    inline void checkAvailable(std::size_t size) const {
        if (index_ + size > bytes_.size()) {
            throw std::out_of_range("ByteScanner: not enough data");
        }
    }

public:
    ByteScanner() = default;
    explicit ByteScanner(std::string_view s) noexcept { reset(s); }

    // Caller must ensure 's' stays alive as long as ByteScanner is used.
    inline void reset(std::string_view s) noexcept {
        bytes_ = s;
        index_ = 0;
    }

    inline bool canRead(std::size_t size) const noexcept {
        return index_ + size <= bytes_.size();
    }

    inline uint8_t readByte() noexcept {
        // assume checked at higher level
        uint8_t v = static_cast<uint8_t>(bytes_[index_]);
        ++index_;
        return v;
    }

    inline uint32_t readUint32() noexcept {
        const uint8_t* p = reinterpret_cast<const uint8_t*>(bytes_.data() + index_);
        uint32_t v = (uint32_t)p[0]
                   | (uint32_t)p[1] << 8
                   | (uint32_t)p[2] << 16
                   | (uint32_t)p[3] << 24;
        index_ += 4;
        return v;
    }

    inline int32_t readInt32() noexcept {
        const uint8_t* p = reinterpret_cast<const uint8_t*>(bytes_.data() + index_);
        int32_t v = (int32_t)p[0]
                  | (int32_t)p[1] << 8
                  | (int32_t)p[2] << 16
                  | (int32_t)p[3] << 24;
        index_ += 4;
        return v;
    }

    inline int16_t readInt16() noexcept {
        const uint8_t* p = reinterpret_cast<const uint8_t*>(bytes_.data() + index_);
        int16_t v = (int16_t)p[0]
                  | (int16_t)p[1] << 8;
        index_ += 2;
        return v;
    }

    inline uint16_t readUint16() noexcept {
        const uint8_t* p = reinterpret_cast<const uint8_t*>(bytes_.data() + index_);
        uint16_t v = (uint16_t)p[0]
                   | (uint16_t)p[1] << 8;
        index_ += 2;
        return v;
    }

    inline float readFloat() noexcept {
        float v;
        std::memcpy(&v, bytes_.data() + index_, sizeof(float));
        index_ += sizeof(float);
        return v;
    }

    // This one can throw → do NOT mark noexcept.
    inline std::string readString(std::size_t len) {
        if (index_ + len > bytes_.size()) {
            throw std::out_of_range("ByteScanner::readString out of range");
        }
        std::string s(bytes_.substr(index_, len));
        index_ += len;
        return s;
    }

    inline std::string readUntil(char delim) {
        std::size_t start = index_;
        while (index_ < bytes_.size() && bytes_[index_] != delim) {
            ++index_;
        }
        std::size_t end = index_;
        if (index_ < bytes_.size() && bytes_[index_] == delim) {
            ++index_; // consume delimiter
        }
        return std::string(bytes_.substr(start, end - start));
    }

    inline const char* data() const noexcept {
        return bytes_.data() + index_;
    }

    inline std::size_t getIndex() const noexcept {
        return index_;
    }

    inline void incrementIndex(std::size_t n) noexcept {
        index_ += n;
    }

    inline void setIndex(std::size_t pos) noexcept {
        index_ = pos;
    }
};
