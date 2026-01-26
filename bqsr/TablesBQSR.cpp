
#include "TablesBQSR.h"

// Log10QualEmpiricalPriorCache is now in the header for visibility

bool bqsrTableKey::operator==(const bqsrTableKey &rhs) const {
    return Qual == rhs.Qual &&
           Covariate == rhs.Covariate &&
           ReadGroup == rhs.ReadGroup;
}

bool bqsrTableKey::operator!=(const bqsrTableKey &rhs) const {
    return !(rhs == *this);
}

/**
 * bqsrTable functions
 * */

void bqsrTable::update(bqsrEntry* table, size_t idx, int64_t errorVal){
    if (table[idx].Observations == -1) {
        table[idx].Observations++;
    }
    table[idx].Observations++;
    table[idx].Mismatches += errorVal;
}

void bqsrTable::merge(bqsrTable &other) {
    if (table.size() != other.table.size()) {
       // Should not happen
       return;
    }
    for (size_t i = 0; i < table.size(); ++i) {
        table[i].Observations += other.table[i].Observations;
        table[i].Mismatches += other.table[i].Mismatches;
    }
}
