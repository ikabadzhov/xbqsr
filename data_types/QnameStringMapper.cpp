#include "QnameStringMapper.h"
#include "../Utils.h"
#include "../sycl_parallel_library/sorting/external_sort.h"
#include <tbb/parallel_sort.h>
#include <numeric>
#include <iostream>


ssize_t QnameStringMapper::addString(const std::string &str) {
    if (nextPosition==sortedStrings.size()) {
        sortedStrings.resize(sortedStrings.size() + REALLOCATION);
        firstStringsIdxs.resize(firstStringsIdxs.size() + REALLOCATION);
    }
    ssize_t oldIdx = nextPosition;
    sortedStrings[nextPosition] = str;
    firstStringsIdxs[nextPosition] = nextPosition;
    nextPosition++;
    return oldIdx;
}

void QnameStringMapper::resize(size_t sz) {
    sortedStrings.resize(sz);
    firstStringsIdxs.resize(sz);
}

void QnameStringMapper::finalize(std::vector<int32_t> &qname) { // TODO: To optimize
    sortedStrings.resize(qname.size());
    std::vector<size_t> indexes;
    std::vector<size_t> uniqIndexes;

    indexes.resize(sortedStrings.size());

    std::iota(indexes.begin(), indexes.end(), 0);

    tbb::parallel_sort(indexes.begin(), indexes.end(), [&](auto& p1, auto& p2){
        return sortedStrings[p1] < sortedStrings[p2];
    });

    permute(indexes, sortedStrings);
    permute(indexes, firstStringsIdxs);

    std::iota(indexes.begin(), indexes.end(), 0);

    uniqIndexes.resize(sortedStrings.size());
    std::vector<uint8_t> duplicates_mask(sortedStrings.size(), 0);

    duplicates_mask[0] = 1;
    tbb::parallel_for (1UL, duplicates_mask.size(), [&](size_t i) {
        duplicates_mask[i] = (sortedStrings[i] == sortedStrings[i - 1]) ? 0 : 1;
    });

//    std::inclusive_scan(duplicates_mask.begin(), duplicates_mask.end(), uniqIndexes.begin(),
//                        [](size_t a, uint8_t b) { return a + b; });

    uniqIndexes[0] = static_cast<size_t>(duplicates_mask[0]);  // Start with the first element

    for (size_t i = 1; i < duplicates_mask.size(); ++i) {
        uniqIndexes[i] = uniqIndexes[i - 1] + static_cast<size_t>(duplicates_mask[i]);
    }


    tbb::parallel_sort(indexes.begin(), indexes.end(), [&](auto& idx1, auto& idx2){
        return firstStringsIdxs[idx1] < firstStringsIdxs[idx2];
    });

    tbb::parallel_for(0UL, qname.size(), [&](size_t i) {
        size_t oldStringID = qname[i];
        size_t ptrToNewStringID = indexes[oldStringID];
        qname[i] = uniqIndexes[ptrToNewStringID] - 1;
    });

    auto lastItr = std::unique(sortedStrings.begin(), sortedStrings.end(), [&](auto& p1, auto& p2) {
        return p1 == p2;
    });

    sortedStrings.resize(std::distance(sortedStrings.begin(), lastItr));
    nextPosition = sortedStrings.size();
    firstStringsIdxs.clear();
}

const std::string& QnameStringMapper::GetAt(int32_t id) {
    return this->sortedStrings[id];
}

