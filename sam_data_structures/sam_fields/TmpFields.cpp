#include "TmpFields.h"
#include "../../Utils.h"

void TmpFields::Permute(std::vector<size_t> &indices) {
    permute(indices, this->adaptPos);
    permute(indices, this->adaptScore);
    permute(indices, this->nextREFID);
    permute(indices, this->REFID);
    permute(indices, this->LIBID);
}



