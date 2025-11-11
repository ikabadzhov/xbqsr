#ifndef GENOMICPIPELINE_COMPRESSEDINPUTSTREAM_H
#define GENOMICPIPELINE_COMPRESSEDINPUTSTREAM_H

#include <cstring>
#include <iostream>
#include <vector>
#include <array>
#include <fcntl.h>
#include <cinttypes>
#include "PlainInputStream.h"
#include "../../zlib.h"

#define USE_LIBDEFLATE 1
#ifdef USE_LIBDEFLATE
#include <libdeflate.h>
#endif


struct BlockMeta {
    uint64_t comp_off;
    uint64_t decomp_off;
    size_t comp_size;
    size_t decomp_size;
};

struct BGZFBlockInfo {
    std::size_t comp_size;  // compressed size (bytes on disk)
    std::size_t decomp_size; // uncompressed size (bytes after inflate)
};


const size_t BGZF_BLOCK_SIZE = 65536;

class CompressedInputStream : public InputStream {
    struct Options {
        // Max bytes we’re willing to keep in the decompressed stash at once.
        // (Protects memory use on huge read() requests; can be large if you want.)
        std::size_t stash_cap_bytes = 1ull << 30; // 1 GiB default
        // Optional: preallocate stash to this many bytes (avoid first growth cost)
        std::size_t stash_reserve_bytes = 1ull << 30;
    };

    gzFile fileDescriptor;
    std::string path;
    //z_stream strm{};

    static bool isValidBlockHeader(std::string& buffer);

public:
    static bool isValidFile(PlainInputStream& stream);

    explicit CompressedInputStream(const std::string &input);


    ~CompressedInputStream() override;

    ssize_t read(std::string& buffer, ssize_t offset, ssize_t length) override;
    void seekg(ssize_t pos) override;
    const std::string & getFilePath() override;
    int peek() override;
    bool good() override;
    bool bad() override;
    bool eof() override;
    void close() override;

    private:
    // ------------ Helpers & state ------------

    static inline uint16_t le16_(const uint8_t* p) {
        return static_cast<uint16_t>(p[0] | (p[1] << 8));
    }
    static inline uint32_t le32_(const uint8_t* p) {
        return static_cast<uint32_t>(p[0]) |
               (static_cast<uint32_t>(p[1]) << 8) |
               (static_cast<uint32_t>(p[2]) << 16) |
               (static_cast<uint32_t>(p[3]) << 24);
    }

    // Parse GZIP header extra field at 'off' to get BGZF BSIZE.
    // Returns compressed block size in bytes (= BSIZE + 1).
    // Throws on invalid header/short read.
    BGZFBlockInfo probe_bgzf_block_size_(uint64_t off);

    // Inflate one BGZF block at comp_off_.
    // Returns true if block produced >0 decompressed bytes;
    // Returns false if this was the BGZF EOF block (0 bytes decompressed).
    bool inflate_one_bgzf_block_();
    bool inflate_bgzf_block_(BlockMeta& block);
    bool get_bgzf_blocks_indexes(std::vector<BlockMeta>& blocks);
    // Move unconsumed bytes to the front so we can append efficiently.
    void compact_stash_();

private:
    // File + position (compressed)
    int       fd_;
    uint64_t  comp_off_;
    uint64_t  decomp_off_;
    bool      seen_eof_block_;

    // Decompressed stash
    std::vector<uint8_t> stash_;
    std::size_t          stash_pos_;
    std::size_t          stash_cap_;

    // Temp buffers per block
    std::vector<uint8_t> tmp_comp_;
    std::vector<uint8_t> tmp_decomp_;
    std::vector<std::array<uint8_t, BGZF_BLOCK_SIZE>> comp_pool;

    // Inflaters
    //z_stream zstrm_;


    libdeflate_decompressor* decomp_ = nullptr;

};


#endif //GENOMICPIPELINE_COMPRESSEDINPUTSTREAM_H
