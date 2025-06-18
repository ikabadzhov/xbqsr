#include "AlignmentIO.h"


using namespace std;


bool AlignmentReader::Good() {
    return in->good();
}

bool AlignmentReader::Bad() {
    return in->bad();
}

bool AlignmentReader::isEOF() {
    return in->eof();
}

bool AlignmentWriter::Good() {
    return this->out->good();
}

bool AlignmentWriter::Bad() {
    return this->out->bad();
}

bool AlignmentWriter::isEOF() {
    return this->out->eof();
}

AlignmentWriter::~AlignmentWriter() = default;
