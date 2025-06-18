#include <fstream>
#include <boost/random.hpp>
#include <boost/crc.hpp>
#include "BitVector.h"

/**
 * BitVector class allows to create a vector of "_nBits" bits.
 * Bits are stored as vector of bytes, and the addressable information unit is the bit.
 * */
BitVector::BitVector(size_t nBits, bool initVal) {
    this->_nBits = nBits;
    if (initVal) {
        this->_buffer.resize((nBits + BITS_PER_BYTE - 1) / BITS_PER_BYTE, 0xFF);
    } else {
        this->_buffer.resize((nBits + BITS_PER_BYTE - 1) / BITS_PER_BYTE, 0);
    }
}


bool BitVector::GetBit(uint64_t idx) const {
    if (idx >= this->_nBits) {
        throw out_of_bound();
    }
    size_t bytePos = idx / BITS_PER_BYTE;
    uint8_t bitPos = BITS_PER_BYTE - 1 - idx % BITS_PER_BYTE;
    auto b = static_cast<byte_t>(this->_buffer[bytePos]);
    return static_cast<bool>( (b >> bitPos) & static_cast<byte_t>(1) );
}

void BitVector::SetBit(uint64_t idx, bool value) {
    size_t bytePos = idx / BITS_PER_BYTE;
    uint8_t bitPos = BITS_PER_BYTE - 1 - idx % BITS_PER_BYTE;
    if (idx >= _nBits) {
        throw out_of_bound();
    }
    if (value) {
        _buffer[bytePos] |= (static_cast<byte_t>(1) << bitPos);
    } else {
        _buffer[bytePos] &= ~(static_cast<byte_t>(1) << bitPos);
    }
}

size_t BitVector::SizeInBits() const {
    return _nBits;
}

size_t BitVector::SizeInBytes() const {
    return this->_buffer.size();
}

void BitVector::Clear() {
    _buffer.clear();
    _nBits = 0;
}

void BitVector::ResetValues(){
    for (auto& b : _buffer){
        b=0;
    }
}

uint32_t BitVector::ComputeCRC32(){
    boost::crc_32_type  crc32code;
    crc32code.process_bytes(_buffer.data(), _nBits );
    return crc32code.checksum();
}

uint32_t BitVector::GetRangeBits(uint64_t startIdx, uint64_t nBits) const {
    uint32_t rangeValue = 0;
    if( nBits >= BITS_PER_BYTE * sizeof(rangeValue) || (startIdx + nBits) > _nBits ){
        throw out_of_bound();
    }

    size_t end = std::min(startIdx + nBits, _nBits);
    for (size_t idx = startIdx; idx < end; idx++) {
        rangeValue <<= 1;
        rangeValue += static_cast<byte_t>(this->GetBit(idx));
    }
    return rangeValue;
}

void BitVector::FillBuffer(const std::vector<byte_t>& data, uint64_t startIdx, std::uint64_t nBitsToCopy) {
    BitsView dataView(data);
    uint64_t locPtr = 0; //TODO: Optimize this bits insertion
    for (uint64_t bIdx = startIdx; bIdx < startIdx + nBitsToCopy; bIdx++, locPtr++) {
        SetBit(locPtr, dataView.GetBit(bIdx));
    }
}

void BitVector::Resize(size_t nBits) {
    this->_nBits = nBits;
    this->_buffer.resize((nBits + BITS_PER_BYTE - 1) / BITS_PER_BYTE, 0);
}

std::istream &operator>>(std::istream &is, BitVector &vector) {
    uint64_t idx = 0;
    std::string line;
    std::getline(is, line);
    for (auto bit : line) {
        if (bit == '1' || bit == '0') {
            bool b = (bit == '1');
            vector.SetBit(idx, b);
            idx++;
            if (idx >= vector.SizeInBits()) {
                break;
            }
        }
    }
    return is;
}

BitVector BitVector::ExtractBits(uint64_t startIdx, uint64_t nBits) {
    BitVector copiedVec(nBits);
    copiedVec.FillBuffer(_buffer, startIdx, nBits);
    return copiedVec;
}

std::string BitVector::ToString() const {
    std::string out(_nBits, '0');
    for (size_t bIdx = 0; bIdx < _nBits; bIdx++) {
        out[bIdx] = GetBit(bIdx) ? '1' : '0';
    }
    return out;
}

const std::vector<byte_t> &BitVector::GetBuffer() {
    return _buffer;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * BitsView class allows to create a "view" over the bits in a vector of bytes.
 * This means that it allowas to access the single bit rather than the byte, without allocating a BitVector
 */
BitsView::BitsView(const std::vector<byte_t>& data) : _ptr(data.data()), _len(data.size() * BITS_PER_BYTE) {}

BitsView::BitsView(const byte_t* ptr, size_t sz) : _ptr(ptr), _len(sz) {}

bool BitsView::GetBit(uint64_t idx){
    if (idx >= _len) {
        throw out_of_bound();
    }
    size_t bytePos = idx / BITS_PER_BYTE;
    uint8_t bitPos = BITS_PER_BYTE - 1 - idx % BITS_PER_BYTE;
    auto b = static_cast<byte_t>(_ptr[bytePos]);
    return static_cast<bool>( (b >> bitPos) & static_cast<byte_t>(1) );
}
