#include "tileInfo.h"

tileInfo::tileInfo(int16_t tile, int16_t xCoord, int16_t yCoord) : t(tile), x(xCoord), y(yCoord){}

bool tileInfo::operator<(const tileInfo &rhs) const {
    if (t < rhs.t)
        return true;
    if (rhs.t < t)
        return false;
    if (x < rhs.x)
        return true;
    if (rhs.x < x)
        return false;
    return y < rhs.y;
}

bool tileInfo::operator>(const tileInfo &rhs) const {
    return rhs < *this;
}

bool tileInfo::operator<=(const tileInfo &rhs) const {
    return !(rhs < *this);
}

bool tileInfo::operator>=(const tileInfo &rhs) const {
    return !(*this < rhs);
}

bool tileInfo::operator==(const tileInfo &rhs) const {
    return t == rhs.t &&
           x == rhs.x &&
           y == rhs.y;
}

bool tileInfo::operator!=(const tileInfo &rhs) const {
    return !(rhs == *this);
}

inline bool parse_int_fast(std::string_view s, int& out) {
    if (s.empty()) return false;
    int sign = 1;
    size_t i = 0;
    if (s[0] == '-') {
        sign = -1;
        i = 1;
        if (i == s.size()) return false;
    }

    int val = 0;
    for (; i < s.size(); ++i) {
        unsigned char c = static_cast<unsigned char>(s[i]);
        if (c < '0' || c > '9') return false;
        int d = c - '0';
        val = val * 10 + d;
        // you can clamp / early-fail if you care about overflow
    }
    out = sign * val;
    return true;
}

tileInfo computeTileInfo(std::string_view qname) {
    // We only ever need up to 7 fields
    std::string_view fields[7];
    int field_count = 0;

    size_t start = 0;
    for (size_t i = 0; i <= qname.size(); ++i) {
        if (i == qname.size() || qname[i] == ':') {
            if (field_count < 7) {
                fields[field_count++] = std::string_view(qname.data() + start, i - start);
            } else {
                // More than 7 columns: we don't support this format
                return tileInfo{-1, -1, -1};
            }
            start = i + 1;
        }
    }

    int t = 0, x = 0, y = 0;
    bool ok = false;

    switch (field_count) {
        case 7: {
            // instrument:run:flowcell:lane:tile:x:y
            int it, ix, iy;
            ok =
                parse_int_fast(fields[4], it) &&
                parse_int_fast(fields[5], ix) &&
                parse_int_fast(fields[6], iy);
            if (!ok) return tileInfo{-1, -1, -1};
            t = it; x = ix; y = iy;
            break;
        }
        case 5: {
            // (legacy style) instrument:run:tile:x:y or similar
            int it, ix, iy;
            ok =
                parse_int_fast(fields[2], it) &&
                parse_int_fast(fields[3], ix) &&
                parse_int_fast(fields[4], iy);
            if (!ok) return tileInfo{-1, -1, -1};
            t = it; x = ix; y = iy;
            break;
        }
        default:
            return tileInfo{-1, -1, -1};
    }

    // clamp / validate for int16_t if needed
    if (t < std::numeric_limits<int16_t>::min() || t > std::numeric_limits<int16_t>::max() ||
        x < std::numeric_limits<int16_t>::min() || x > std::numeric_limits<int16_t>::max() ||
        y < std::numeric_limits<int16_t>::min() || y > std::numeric_limits<int16_t>::max()) {
        return tileInfo{-1, -1, -1};
    }

    return tileInfo{static_cast<int16_t>(t),
                    static_cast<int16_t>(x),
                    static_cast<int16_t>(y)};
}