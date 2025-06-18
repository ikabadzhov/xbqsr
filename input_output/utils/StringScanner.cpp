#include <vector>
#include <string>
#include <iostream>
#include "StringScanner.h"
#include "../../sam_data_structures/sam_fields/Sequence.h"
#include "../../file_format/Global.h"


using namespace std;


// reset resets the scanner, and initializes it with the given string.
void StringScanner::reset(std::string_view s) {
    this->index = 0;
    this->data = s;
}

// len returns the number of ASCII characters that still need to be
// scanned/parsed. Returns 0 if Err() would return a non-nil value.
size_t StringScanner::len() {
    return this->data.size() - this->index;
}

bool StringScanner::readByteUntil(char c, char &val) {

    size_t start = this->index;
    size_t next = start + 1;

    if ( next >= this->data.size() ) {
        this->index = this->data.size();
        val=this->data[start];
        return false;
    } else if ( this->data[next] != c ){
        cerr<<"unexpected character " << this->data[next] <<" in stringScanner.ReadByteUntil";
        exit(-1);
    } else {
        this->index = next + 1;
        val=this->data[start];
        return true;
    }
}

bool StringScanner::readUntil(char c, string& s) {
    s.clear();
    size_t n=this->data.find(c, this->index);

    if(n!=std::string::npos) {
        size_t offset = n - this->index;
        s.append(data.substr(this->index, offset));
        this->index = n + 1;
        return true;
    }
    else{
        //s.resize(this->data.size());
        s.append(data.substr(this->index, this->data.size()));
        this->index = this->data.size();
        return false;
    }
}

bool StringScanner::readBytes(std::string& bytes) {
    //this->readUntil('\t', bytes);
    if (!readUntil('\t', bytes)) {
        return false;
    }else   return true;
}

void StringScanner::readUntil2(char c1, char c2, std::string& s, char &b) {// (s string, b byte) {
    size_t start = this->index;
    for( size_t end = this->index; end < this->data.size(); end++ ) {
        b = this->data[end];
        if ( (b == c1) || (b == c2) ) {
            this->index = end + 1;
            s.resize(end-start);
            std::copy(this->data.begin()+start, this->data.begin()+end, s.begin());
            return;
        }
    }
    this->index = this->data.size();
    s.resize(this->data.size()-start);
    std::copy(this->data.begin()+start, this->data.end(), s.begin());
}

void StringScanner::doString( std::string& value) {
    value.clear();
    this->readUntil('\t', value);
}

bool StringScanner::doSeq(std::string& seq) {
    if (this->index >= this->data.size()) {
        cerr << "Panic: index out of range" << std::endl;
        return false;
    }

    size_t count = 0;
    bool found = false;
    size_t end = 0;

    for (end = this->index; end < this->data.size(); ++end) {
        if (data[end] == '\t') {
            found = true;
            //this->index = end + 1; // Update index for the next call
            break;
        }
        count++;
    }

    if (!found) {
        cerr << "Panic: missing tabulator in SAM alignment line" << std::endl;
        return false;
    }

    seq.resize(count);

    std::copy_n(data.data()+index, count, seq.begin());
    //seq.SetSequenceFromString(sequence);
    /*for (size_t i = 0; i < count; i++, index++) {
        seq.SetBase(static_cast<int>(i), data[this->index]);
    }*/
    this->index = end + 1;
    return true;
}

int32_t StringScanner::doInt32(){
    string value;
    this->doString(value);
    return int32_t(stoll(value));
}

uint64_t StringScanner::doUint() {
    string value;
    this->doString(value);
    return uint64_t(stoull(value));
}
