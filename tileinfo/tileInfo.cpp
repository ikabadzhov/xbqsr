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

tileInfo computeTileInfo(const std::string &QNAME) {
    std::vector<std::string> qnameInfo;
    qnameInfo.reserve(7);

    ssize_t pos = 0;
    ssize_t endPos = 0;
    while((endPos = QNAME.find_first_of(':', pos))!=std::string::npos) {
        qnameInfo.push_back(QNAME.substr(pos, endPos-pos));
        pos = endPos + 1;
    }
    qnameInfo.push_back(QNAME.substr(pos, endPos-pos));

    int32_t nColumns = static_cast<int32_t>(qnameInfo.size());

    int32_t t = 0, x = 0, y = 0;

    switch (nColumns) {
        case 7:
            t = stoi(qnameInfo[4]);
        x = stoi(qnameInfo[5]);
        y = stoi(qnameInfo[6]);
        break;
        case 5:
            t = stoi(qnameInfo[2]);
        x = stoi(qnameInfo[3]);
        y = stoi(qnameInfo[4]);
        break;
        default:
            // std::cout << "Warning: Unsupported qname format for extracting tile info. Cannot properly perform optical duplicate marking. "
            //           << nColumns << " columns." << std::endl;
        return tileInfo{-1, -1, -1};
    }
    return tileInfo{static_cast<int16_t>(t), static_cast<int16_t>(x), static_cast<int16_t>(y)};
}