#ifndef TILEINFO_H
#define TILEINFO_H

#include <cstdint>
#include <string_view>
#include <sycl/sycl.hpp>

class tileInfo {
public:
    int16_t t, x, y;
    // TODO: This is for compatibility with GATK.. They treat these variables as uint16_t (to check), causing overflow in many dataset
    tileInfo() = default;
    tileInfo(int16_t tile, int16_t xCoord, int16_t yCoord);

    SYCL_EXTERNAL bool operator<(const tileInfo &rhs) const;

    SYCL_EXTERNAL bool operator>(const tileInfo &rhs) const;

    bool operator<=(const tileInfo &rhs) const;

    bool operator>=(const tileInfo &rhs) const;

    bool operator==(const tileInfo &rhs) const;

    bool operator!=(const tileInfo &rhs) const;
};

class TileRGKey {
public:
    int32_t rg;
    int t;
    friend class TileRGHasher;
    TileRGKey(int32_t RG, int t) : rg(RG), t(t) {}

    bool operator==(const TileRGKey &rhs) const {
        return rg == rhs.rg &&
               t == rhs.t;
    }

    bool operator!=(const TileRGKey &rhs) const {
        return !(rhs == *this);
    }
};

class TileRGHasher {
public:
    std::size_t operator()(const TileRGKey& tile) const
    {
        return static_cast<int32_t>(tile.rg) ^ static_cast<uint64_t>(tile.t);
    }
};

tileInfo computeTileInfo(std::string_view QNAME);


#endif //TILEINFO_H
