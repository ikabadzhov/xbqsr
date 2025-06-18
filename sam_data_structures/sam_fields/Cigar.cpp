#include "Cigar.h"
#include "../../file_format/Global.h"
#include "../../data_types/VectorCollection.h"

using namespace std;

std::unordered_map<std::string, std::vector<CigarOperation>> cigarCache = { {"*" , vector<CigarOperation>(1, CigarOperation{0,0}) } };


bool CigarOperation::operator==(const CigarOperation &rhs) const {
    return Length == rhs.Length &&
           Operation == rhs.Operation;
}

bool CigarOperation::operator!=(const CigarOperation &rhs) const {
    return !(rhs == *this);
};

CigarOperationsTable::CigarOperationsTable() {
    this->cigar_operations_table.resize(256);
    for(const char c:cigarOperations){
        this->cigar_operations_table[static_cast<uint8_t>(c)]=std::toupper(c);
    }
}

char &CigarOperationsTable::operator[](size_t i) {
    return this->cigar_operations_table[i];
}


CigarOperation new_cigar_operation(std::string& cigar, size_t i, size_t& j){

    for(j=i; ;j++){
        char ch = cigar[j];

        if(!isdigit(ch)){
            uint32_t length = static_cast<uint32_t>(std::stoull(cigar.substr(i,(j-i))));
            char operation = cigarOperationsTable[ch];
            if(operation != 0){
                //op = CigarOperation(length, operation);
                j++;
            }
            else{
                std::cerr<<"invalid CIGAR operation " << operation<<std::endl;
            }
            return CigarOperation(length, operation);
        }
    }
}

bool slow_scan_cigar_string(std::string& cigar, std::vector<CigarOperation>& out) {

    if( cigar.empty() ) {
        return false;
    }

    size_t i;
    CigarOperation cigar_operation = new_cigar_operation(cigar, 0, i);

    out.emplace_back(cigar_operation);

    while ( i < cigar.size() ) {
        size_t j;
        CigarOperation nextCigarOperation = new_cigar_operation(cigar, i, j);
        if( nextCigarOperation.Operation == cigar_operation.Operation ) {
            out[out.size()-1].Length += nextCigarOperation.Length;
        } else {
            out.emplace_back(nextCigarOperation);
        }
        cigar_operation = nextCigarOperation;
        i = j;
    }

    cigarCache[cigar]=out;

    /*
     *

    if value, found := cigarSliceCache[cigar]; found {
            slice = value
    } else {
        cigarSliceCache[cigar] = slice
    }
    cigarSliceCacheMutex.Unlock()*/
    return true;
}

/*var (
        cigarSliceCache      = map[string][]CigarOperation{"*": []CigarOperation{}}
cigarSliceCacheMutex = sync.RWMutex{}
)*/

// ScanCigarString converts a CIGAR string to a slice of
// CigarOperation. See http://samtools.github.io/hts-specs/SAMv1.pdf -
// Section 1.4.6.
//
// Uses an internal cache to reduce memory overhead. It is safe for
// multiple goroutines to call ScanCigarString concurrently.
void  ScanCigarString( std::string& cigar_string, std::vector<CigarOperation>& out) {

    if(cigar_string=="*"){
        return;
    }

    if( cigarCache.count(cigar_string) > 0 ){
        out=cigarCache[cigar_string];
        return;
    }
    /* value, found = cigarSliceCache[cigar]
    cigarSliceCacheMutex.RUnlock()
    if found {
        return value;
    }*/
    slow_scan_cigar_string(cigar_string, out);
}


// ReadLengthFromCigar sums the lengths of all CIGAR operations that consume read bases.
uint32_t ReadLengthFromCigar( const std::vector<CigarOperation>& cigar ){
    int32_t length=0;
    for(const CigarOperation& op:cigar){
        if(CigarOperatorConsumesReadBases[static_cast<uint8_t>(op.Operation)] > 0) {
            length += CigarOperatorConsumesReadBases[static_cast<uint8_t>(op.Operation)] * op.Length;
        }
    }
    return length;
}
uint32_t ReadLengthFromCigar(CigarOperation* cigar, SpanField cigarPtrs){
    int32_t length=0;
    for (size_t i = cigarPtrs.start; i < cigarPtrs.start + cigarPtrs.len; i++) {
        CigarOperation op = cigar[i];
        if(CigarOperatorConsumesReadBases[static_cast<uint8_t>(op.Operation)]>0) {
            length += CigarOperatorConsumesReadBases[static_cast<uint8_t>(op.Operation)] * op.Length;
        }
    }
    return length;
}

// ReferenceLengthFromCigar sums the lengths of all CIGAR operations that consume reference bases
uint32_t ReferenceLengthFromCigar(const vector<CigarOperation>& cigar) {
    uint32_t length = 0;
    for(const CigarOperation& op:cigar){
        if(CigarOperatorConsumesReferenceBases[static_cast<uint8_t>(op.Operation)] > 0) {
            length += CigarOperatorConsumesReferenceBases[static_cast<uint8_t>(op.Operation)] * op.Length;
        }
    }
    return length;
}
uint32_t ReferenceLengthFromCigar(CigarOperation* cigar, SpanField cigarPtrs) {
    uint32_t length = 0;
    for (size_t i = cigarPtrs.start; i < cigarPtrs.start + cigarPtrs.len; i++) {
        CigarOperation op = cigar[i];
        if(CigarOperatorConsumesReferenceBases[static_cast<uint8_t>(op.Operation)] > 0) {
            length += CigarOperatorConsumesReferenceBases[static_cast<uint8_t>(op.Operation)] * op.Length;
        }
    }
    return length;
}

uint32_t cigarStringToValue(const std::string& cigarString) {
    uint32_t length = std::stoi(cigarString.substr(0, cigarString.size() - 1));
    char operation = cigarString[cigarString.size() - 1];
    int operationValue = 0;

    for (int i = 0; i < sizeof(cigarOps); ++i) {
        if (cigarOps[i] == operation) {
            operationValue = i;
            break;
        }
    }

    return (length << 4) | operationValue;
}

void cigarToString(std::span<CigarOperation> cigar, string& res_string) {
   if( cigar.empty() ){
        res_string+="*";
        return;
    }
    for(auto& op:cigar){
       // if(op.Length==0){
         //   res_string+='*';
        //}else {
        res_string += (to_string(op.Length) + op.Operation);
        //}//res_string += op.Operation;
    }
}

//const std::string cigarOperations = "MmIiDdNnSsHhPpXx=";





//const int32_t CigarOperatorConsumesReadBases[CIGAR_CONSUME_BASES_SZ] = {
//        /* 0-31 */
//        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//        /* 32-63*/
//        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/* = */, 0, 0,
//        /* 64-79*/
//        0, 0, 0, 0, 0, 0, 0, 0, 0, 1/* I */, 0, 0, 0, 1 /* M */, 0, 0,
//        /* 80-95 */
//        0, 0, 0, 1 /* S */, 0, 0, 0, 0, 1 /* X */, 0, 0, 0, 0, 0, 0, 0,
//        /* 96-127*/
//        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
//};

//const int32_t CigarOperatorConsumesReferenceBases[CIGAR_CONSUME_BASES_SZ] = {
//        /* 0-31 */
//        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//        /* 32-63 */
//        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/* = */, 0, 0,
//        /* 64-79 */
//        0, 0, 0, 0, 1/* D */, 0, 0, 0, 0, 0, 0, 0, 0, 1 /* M */, 1/* N */, 0,
//        /* 80-95 */
//        0, 0, 0, 0, 0, 0, 0, 0, 1 /* X */, 0, 0, 0, 0, 0, 0, 0,
//        /* 96-127 */
//        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
//};

/*
const std::unordered_map<char,int32_t> CigarOperatorConsumesReadBases={
        {'M',1},
        {'I',1},
        {'S',1},
        {'=',1},
        {'X',1}
};

const std::unordered_map<char,int32_t> CigarOperatorConsumesReferenceBases={
        {'M',1},
        {'D',1},
        {'N',1},
        {'=',1},
        {'X',1}
};

*/
