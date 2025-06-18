#include "BamHeader.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "../zlib.h"
#include "../sam_data_structures/sam_fields/Sequence.h"
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

/* CONSTANTS */
// int32_t|uint32_t == 4, uint8_t == 1, uint16_t == 2
const int refIDIndex = 0;
const int posIndex = 4;
const int lReadNameIndex = posIndex + 4; //8
const int mapqIndex = lReadNameIndex + 1; //9
const int binIndex = mapqIndex + 1; //10
const int nCigarOpIndex = binIndex + 2; //12
const int flagIndex = nCigarOpIndex + 2; //14
const int lSeqIndex = flagIndex + 2; //16
const int nextRefIDIndex = lSeqIndex + 4; //20
const int nextPosIndex = nextRefIDIndex + 4; //24
const int tlenIndex = nextPosIndex + 4; //28
const int readNameIndex = tlenIndex + 4; //32


struct AlignmentRecord {
    int32_t block_size;
    int32_t refID; /* Reference sequence ID, −1 ≤ refID < n ref; -1 for a read without a mapping position*/
    int32_t pos; /* 0-based leftmost coordinate (= POS − 1) */
    uint8_t l_read_name; /* Length of read name below (= length(QNAME) + 1) */
    uint8_t mapq; /* Mapping quality */
    uint16_t bin; /* BAI index bin, TODO Section 4.2.1 */
    uint16_t n_cigar_op; /* Number of operations in CIGAR */
    uint16_t flag; /* Bitwise flag */
    uint32_t l_seq; /* Length of SEQ */
    int32_t next_refID; /* Ref-ID of the next segment (−1 ≤ next refID < n ref) */
    int32_t next_pos; /* 0-based leftmost pos of the next segment (= PNEXT − 1) */
    int32_t tlen; /* Template length */
    //std::string read_name; /* It is the same thing as QNAME - Read name, NB: NUL-terminated (QNAME with trailing ‘\0’) */
    std::string QNAME; /* read_name == QNAME - The Reference sequence NAME*/
    std::vector<CigarOperation> cigar; /* CIGAR: op len<<4|op. ‘MIDNSHP=X’→‘012345678’ */
    std::vector<uint8_t> seq; /* 4-bit encoded read: ‘=ACMGRSVTWYHKDBN’→ [0, 15].*/
    std::string seq_clear;
    std::string qual; /* Phred-scaled base qualities. TODO Section 4.2.3 */
};

/* How to understand machine-indianess - My machine is Little-endian friendly
Endianness_old checkEndiannessOfHostMachine() {
    union {
        uint32_t value;
        char bytes[sizeof(uint32_t)];
    } u;

    u.value = 0x12345678;

    if (u.bytes[0] == 0x78)
        return Endianness_old::LittleEndian;
    else if (u.bytes[sizeof(uint32_t) - 1] == 0x78)
        return Endianness_old::BigEndian;
    else
        return Endianness_old::Unknown;
}*/


/*
 * The function utilizes zlib's compression functions to convert a file into the BGZF format.
 * It establishes the output file for writing the compressed data and the input file for reading, setting the window size, chunk size, and buffer size.
 * In a loop, it compresses segments of the input data using the deflate algorithm and writes these compressed segments to the output file.
 *
 * Default values of 1024B for Zlib buffer size was used.
*/
void compressToBGZF(const std::string &inputFilename, const std::string &outputFilename) {
    // 16 for BGZF + MAX_WBITS for default gzip [From doc]
    const int GZIP_WINDOW_BITS = 16 + MAX_WBITS;
    const int GZIP_CHUNK_SIZE = 16384;
    const int ZLIB_BUFFER_SIZE = 1024;

    gzFile file_out = gzopen(outputFilename.c_str(), "wb");
    if (!file_out) {
        std::cerr << "Error opening output file for writing compressed data." << std::endl;
        return;
    }

    std::ifstream input_file(inputFilename, std::ios::binary);
    if (!input_file) {
        std::cerr << "Error opening input file." << std::endl;
        gzclose(file_out);
        return;
    }

    char in_buffer[ZLIB_BUFFER_SIZE];
    int ret;
    z_stream zs;
    memset(&zs, 0, sizeof(zs));

    //Initialize zlib for default compression
    if (deflateInit2(&zs, Z_DEFAULT_COMPRESSION, Z_DEFLATED, GZIP_WINDOW_BITS, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
        std::cerr << "Error initializing zlib for compression." << std::endl;
        input_file.close();
        gzclose(file_out);
        return;
    }

    zs.avail_in = 0;
    zs.next_in = nullptr;

    do {
        //Read data from the input file into the buffer to avoid too many opening of the same file
        input_file.read(in_buffer, ZLIB_BUFFER_SIZE);
        zs.avail_in = static_cast<uInt>(input_file.gcount());
        zs.next_in = reinterpret_cast<Bytef*>(in_buffer);

        do {
            //compression of the file until Z_FINISH [end-of-file for zlib]
            zs.avail_out = GZIP_CHUNK_SIZE;
            zs.next_out = reinterpret_cast<Bytef*>(in_buffer);
            ret = deflate(&zs, Z_FINISH);
            if (ret < 0) {
                std::cerr << "Error during compression." << std::endl;
                deflateEnd(&zs);
                input_file.close();
                gzclose(file_out);
                return;
            }
            int compressed_size = GZIP_CHUNK_SIZE - zs.avail_out;
            gzwrite(file_out, in_buffer, compressed_size);
        } while (zs.avail_out == 0);
    } while (!input_file.eof());

    // Clean up: Close the input and output files, and free zlib resources
    deflateEnd(&zs);
    input_file.close();
    gzclose(file_out);
}

/*
 * The function utilizes zlib's decompression functions to convert a file from a BGZF format to a human readable one.
 * It establishes the output file for writing the compressed data and the input file for reading, setting the window size, chunk size, and buffer size.
 * Sequentially, it reads 512-byte blocks of compressed data from the input file, considering the initial 10 bytes as the gzip header information.
 * Utilizing zlib's inflate function, it decompresses these blocks and writes the decompressed data into the output file in 64-byte segments.
 *
 * Decompressed chunks of 64Kb are read for each loop, getting closer to the definition of the biggest side of the BAM file's chunks.
*/
void decompressBGZF(const std::string &inputFilename, const std::string &outputFilename) {
    gzFile file_in = gzopen(inputFilename.c_str(), "rb");
    if (!file_in) {
        std::cerr << "Error opening input file for reading compressed data." << std::endl;
        return;
    }

    std::ofstream output_file(outputFilename, std::ios::binary | std::ios::out);
    if (!output_file) {
        std::cerr << "Error opening output file for writing decompressed data." << std::endl;
        gzclose(file_in);
        return;
    }

    z_stream zs;
    memset(&zs, 0, sizeof(zs));

    // Use 32 for GZIP decoding [Read doc]
    if (inflateInit2(&zs, 32 + MAX_WBITS) != Z_OK) {
        std::cerr << "Error initializing zlib for decompression." << std::endl;
        gzclose(file_in);
        output_file.close();
        return;
    }

    unsigned char in_buffer[512]; //This is the buffer that decides the size of each block (in our case, 64)

    int ret;
    do {
        // The gzip header is located in the first 10 bytes, the one containing @SQ for instance
        // TODO: To confirm why is 10 bytes, it works but not specified in the original doc
        gzread(file_in, in_buffer, 10);
        // Read the remaining block (>>10 --> buff+10) after the header
        int block_size = gzread(file_in, in_buffer + 10, sizeof(in_buffer) - 10);
        if (block_size <= 0) {
            std::cerr << "Error reading compressed data." << std::endl;
            inflateEnd(&zs);
            gzclose(file_in);
            output_file.close();
            return;
        }

        zs.avail_in = block_size + 10;
        zs.next_in = reinterpret_cast<Bytef*>(in_buffer);

        do {
            char out_buffer[64]; // This is the buffer that will contain the decompressed file data
            zs.avail_out = sizeof(out_buffer);
            zs.next_out = reinterpret_cast<Bytef*>(out_buffer);

            ret = inflate(&zs, Z_NO_FLUSH);
            if (ret < 0 && ret != Z_STREAM_END) {
                std::cerr << "Error during decompression." << std::endl;
                inflateEnd(&zs);
                gzclose(file_in);
                output_file.close();
                return;
            }

            int decompressed_size = sizeof(out_buffer) - zs.avail_out;
            output_file.write(out_buffer, decompressed_size);
        } while (zs.avail_out == 0);

    } while (ret != Z_STREAM_END);

    // Clean up: Close the input and output files, and free zlib resources
    inflateEnd(&zs);
    gzclose(file_in);
    output_file.close();
}

void printField(const std::string& fieldName, const std::string& value) {
    std::cout << fieldName << " --> " << value << std::endl;
}

/*
 * This function is designed to organize information from the header data in the BAM file, represented in plain text.
 * Its goal is to structure the details stored in the headerData variable into a vector of the ReferenceSequence structure.
 * This approach allows for storing the data in a manner resembling a dictionary, facilitating the establishment of correspondence between
 * the names and lengths of each reference.
 *
 * Note: For a line to be considered part of the sequence dictionary, it must commence with @SQ.
*/
std::vector<ReferenceSequence_old> parseHeader(const std::string& headerData) {
    std::vector<ReferenceSequence_old> sequences;
    std::istringstream iss(headerData);
    std::string line;

    while (std::getline(iss, line)) {
        if (line.compare(0, 4, "@SQ\t") == 0) {
            ReferenceSequence_old sequence;
            std::istringstream lineStream(line.substr(4)); // Remove "@SQ\t" from the line

            std::string token;
            while (lineStream >> token) {
                size_t colonPos = token.find(':');
                if (colonPos != std::string::npos) {
                    std::string key = token.substr(0, colonPos);
                    std::string value = token.substr(colonPos + 1);

                    if (key == "SN") {
                        sequence.name = value;
                    } else if (key == "LN") {
                        sequence.length = std::stoull(value);
                    }
                }
            }

            sequences.push_back(sequence);
        }
    }

    return sequences;
}

void decompressBGZF_working(const std::string &inputFilename) {
    gzFile file_in = gzopen(inputFilename.c_str(), "rb");
    if (!file_in) {
        std::cerr << "Error opening input file for reading compressed data." << std::endl;
        return;
    }

    z_stream zs;
    memset(&zs, 0, sizeof(zs));
    // Use 32 for GZIP decoding, for more details go here: https://www.zlib.net/manual.html#inflateInit2
    if (inflateInit2(&zs, 32 + MAX_WBITS) != Z_OK) {
        std::cerr << "Error initializing zlib for decompression." << std::endl;
        gzclose(file_in);
        return;
    }

    HeaderRecord_old headerRecord;
    unsigned char in_buffer_header[4];
    gzread(file_in, in_buffer_header, sizeof(in_buffer_header));

    /* if printing is needed
     * std::cout << "magicNumber[0] --> " << in_buffer[0] << std::endl;
     * std::cout << "magicNumber[1] --> " << in_buffer[1] << std::endl;
     * std::cout << "magicNumber[2] --> " << in_buffer[2] << std::endl;
     * std::cout << "magicNumber[3] --> " << in_buffer[3] << std::endl; */

    // Check if it's a valid BAM file
    if (in_buffer_header[0] != 'B' || in_buffer_header[1] != 'A' || in_buffer_header[2] != 'M' || in_buffer_header[3] != '\1') {
        std::cerr << "Invalid BAM file format." << std::endl;
        gzclose(file_in);
        return;
    }

    std::cout << " ---------- START  HEADER PARCING ---------- " << std::endl;
    // Header parcing
    gzread(file_in, reinterpret_cast<char *>(&headerRecord.header_length), sizeof(uint32_t));
    std::cout << "Header length: " << headerRecord.header_length << std::endl;

    char *headerData = new char[headerRecord.header_length];
    gzread(file_in, headerData, headerRecord.header_length); //length to be modified!
    headerRecord.references = parseHeader(headerData);
    delete[] headerData;
    std::cout << "Parsed header:\n";
    for (const auto& sequence : headerRecord.references) {
        std::cout << "Name: " << sequence.name << ", Length: " << sequence.length << std::endl;
    }
    headerRecord.n_references = headerRecord.references.size();
    uint32_t checkNREF;
    gzread(file_in, reinterpret_cast<char *>(&checkNREF), sizeof(uint32_t));
    if(std::to_string(checkNREF) != std::to_string(headerRecord.n_references)){
        std::cout << "Missmatch between #references calculated and in the file" << std::endl;
        printField("CheckNREF", std::to_string(checkNREF));
        printField("HeaderRecord", std::to_string(headerRecord.n_references));
        return;
    }
    std::cout << " ---------- END HEADER PARCING ---------- " << std::endl;
    std::cout << " ---------- START LIST OF REFERENCE INFORMATION PARCING ---------- " << std::endl;

    //list of reference information parcing
    for(int i=0; i<headerRecord.n_references; i++) {
        ReferenceInformation_old referenceInformation;
        gzread(file_in, reinterpret_cast<char *>(&headerRecord.references[i].referenceInformation.length_name), sizeof(uint32_t));
        printField("Length name", std::to_string(headerRecord.references[i].referenceInformation.length_name));

        std::vector<char> tempBuffer(headerRecord.references[i].referenceInformation.length_name);
        gzread(file_in, tempBuffer.data(), headerRecord.references[i].referenceInformation.length_name);
        headerRecord.references[i].referenceInformation.name = std::string(tempBuffer.data(), tempBuffer.data() + headerRecord.references[i].referenceInformation.length_name);
        auto null_terminator_pos = std::find(headerRecord.references[i].referenceInformation.name.begin(), headerRecord.references[i].referenceInformation.name.end(), '\0');
        headerRecord.references[i].referenceInformation.name = std::string(headerRecord.references[i].referenceInformation.name.data(), null_terminator_pos - headerRecord.references[i].referenceInformation.name.begin());
        printField("Name", headerRecord.references[i].referenceInformation.name);

        gzread(file_in, reinterpret_cast<char *>(&headerRecord.references[i].referenceInformation.length_sequence), sizeof(uint32_t));
        printField("Reference sequence length", std::to_string(headerRecord.references[i].referenceInformation.length_sequence));
        std::cout << "---------- reference information["<< i << "] ---------- " << std::endl;
    }

    std::cout << " ---------- END LIST OF REFERENCE INFORMATION PARCING ---------- " << std::endl;
    std::cout << " ---------- START LIST OF ALIGNMENTS PARCING ---------- " << std::endl;

    //we shouldn't do.
    // 1f 8b 08 04 00 00 00 00 00 ff 06 00 42 43 02 00 1b 00 03 00 00 00 00 00 00 00 00 00
    int cont = 1;
    while (true) {
        AlignmentRecord record;
        gzread(file_in, &record.block_size, sizeof(uint32_t));
        printField("Block Size", std::to_string(record.block_size));

        /*We should read 64Kb blocks, but we trust block_size which is more specific*/
        //unsigned char in_buffer[record.block_size];
        unsigned char *in_buffer = new unsigned char[record.block_size];
        gzread(file_in, in_buffer, sizeof(in_buffer));

        std::memcpy(&record.refID, in_buffer + refIDIndex, sizeof(int32_t));
        if (record.refID < 0)
            record.refID = 42; //ASCII code for STAR --> "*"
        printField("refID", std::to_string(record.refID));
        //The RefID must be either STAR (unknown) or −1 ≤ refID < n ref
        if (record.refID != 42 && (record.refID < -1 && record.refID > headerRecord.n_references)) {
            std::cout << "refID not valid" << std::endl;
            return;
        }

        // 0-based leftmost coordinate, the real value is POS-1
        std::memcpy(&record.pos, in_buffer + posIndex, sizeof(int32_t));
        record.pos -= 1;
        if (record.pos == -1) //equal to --> 1f 8b 08 04 00 00 00 00 00 ff 06 00 42 43 02 00 1b 00 03 00 00 00 00 00 00 00 00 00
            break;
        else  printField("Pos", std::to_string(record.pos));

        std::memcpy(&record.l_read_name, in_buffer + lReadNameIndex, sizeof(uint8_t));
        printField("Length read Name", std::to_string(record.l_read_name));

        std::memcpy(&record.mapq, in_buffer + mapqIndex, sizeof(uint8_t));
        printField("MAPQ", std::to_string(record.mapq));

        std::memcpy(&record.bin, in_buffer + binIndex, sizeof(uint16_t));
        // Note unmapped reads with POS 0 (which becomes −1 in BAM) therefore use reg2bin(−1, 0) which is computed as 4680
        record.bin -= 1;
        printField("BIN", std::to_string(record.bin)); //TODO: Read Section 4.2.1

        //for instance here: 43S13M6872N45M there are 4 different operations.
        // S (soft clipping), M (match), and N (gap or intron).
        std::memcpy(&record.n_cigar_op, in_buffer + nCigarOpIndex, sizeof(uint16_t));
        printField("# cigar operations", std::to_string(record.n_cigar_op));

        //TODO: As noted in Section 1.4, reserved FLAG bits should be written as zero and ignored on reading by current software.
        std::memcpy(&record.flag, in_buffer + flagIndex, sizeof(uint16_t));
        printField("Flag", std::to_string(record.flag));

        std::memcpy(&record.l_seq, in_buffer + lSeqIndex, sizeof(uint32_t));
        printField("Length sequence", std::to_string(record.l_seq));


        std::memcpy(&record.next_refID, in_buffer + nextRefIDIndex, sizeof(int32_t));
        if (record.next_refID < 0)
            record.next_refID = 42;
        //The RefID must be either STAR (unknown) or −1 ≤ refID < n ref
        if (record.next_refID != 42 && (record.next_refID < -1 && record.next_refID > headerRecord.n_references)) {
            std::cout << "Next refID not valid" << std::endl;
            return;
        }
        printField("Next refID", std::to_string(record.next_refID));

        // 0-based leftmost coordinate, the real value is NEXT_POS-1
        std::memcpy(&record.next_pos, in_buffer + nextPosIndex, sizeof(int32_t));
        record.next_pos -= 1;
        printField("Next pos", std::to_string(record.next_pos));

        std::memcpy(&record.tlen, in_buffer + tlenIndex, sizeof(int32_t));
        printField("Template length", std::to_string(record.tlen));

        //the last character is unreadable
        if (record.l_read_name-1 > 0){
            record.QNAME = std::string(in_buffer + readNameIndex, in_buffer + readNameIndex + (record.l_read_name-1));
            printField("QNAME", record.QNAME);
        }else std::cout << "Name length == 0, not printable!" << std::endl;


        int index = readNameIndex + record.l_read_name;

        record.cigar.resize(record.n_cigar_op);

        //std::cout << "Index before: " << index << std::endl;
        for (uint32_t i = 0; i < record.n_cigar_op; ++i, index += 4) {
            uint32_t cigar = 0;
            std::memcpy(&cigar, in_buffer + index, sizeof(uint32_t));

            record.cigar[i] = CigarOperation{
                    static_cast<uint32_t>(cigar >> 4),
                    cigarOps[static_cast<int>(0xF & cigar)]
            };
        }

        std::cout << "CIGAR Operations: " << std::endl;

        for (const auto& op : record.cigar) {
            std::cout << op.Operation << "[" << op.Length << "] - ";
        }
        std::cout << std::endl;

        std::cout << "Length seq: " << record.l_seq << std::endl;
        record.seq.resize(record.l_seq);
        std::string decoded_sequence;

        for (size_t i = 0; i < record.l_seq; i += 2) {
            uint8_t nibbles = 0;
            // Copies 1 byte (8 bits) from in_buffer starting at the index into the nibbles variable
            // which represents the encoded nucleotide sequence in 4-bit chunks.
            std::memcpy(&nibbles, in_buffer + index, sizeof(uint8_t));
            index++;

            uint8_t first_nibble = nibbles >> 4;        //first 4 bits
            uint8_t second_nibble = nibbles & 0x0F;     //last 4 bits

            decoded_sequence += "=ACMGRSVTWYHKDBN"[first_nibble]; // first nucleotide

            if (i < record.l_seq - 1) {
                decoded_sequence += "=ACMGRSVTWYHKDBN"[second_nibble]; // second nucleotide
            }
        }
        record.seq_clear = decoded_sequence;

        /*  Could be written easily like this but the for below avoid printing the dash at the last element
         * for (const auto& op : record.cigar) {
            std::cout << op.Operation << "[" << op.Length << "] - ";
        }
        std::cout << std::endl; */

        for (size_t i = 0; i < record.cigar.size(); ++i) {
            const auto& op = record.cigar[i];
            std::cout << op.Operation << "/" << std::to_string(op.Operation) <<"[" << op.Length << "]";

            if (i != record.cigar.size() - 1) {
                std::cout << " - ";
            }
        }
        std::cout << std::endl;
        printField("SEQ", record.seq_clear);
        //if (cont>2)
        //    break;
        cont += 1;
        std::cout << "This is the " << cont << "th alignment read" << std::endl;
        std::cout << "---------------------------------------" << std::endl;
        delete[] in_buffer;
    }

    std::cout << " ---------- END LIST OF ALIGNMENTS PARCING ---------- " << std::endl;

    inflateEnd(&zs);
    gzclose(file_in);
}

void BAMHeaderReader::readHeader(const std::string& bamFilePath) {

    /*
     * compressToBGZF(inputFile, compressedFile);
     * decompressBGZF(compressedFile, decompressedFile);
    */

    decompressBGZF_working(bamFilePath);
    std::cout << "Done parcing!" << std::endl;
}

