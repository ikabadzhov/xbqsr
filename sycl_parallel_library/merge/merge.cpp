#include "merge.h"

merge_range_t compute_merge_range(int a_count, int b_count, int partition, int spacing, int mp0, int mp1) {

    int diag0 = spacing * partition;
    int diag1 = min(a_count + b_count, diag0 + spacing);

    return merge_range_t {mp0, mp1, diag0 - mp0, diag1 - mp1};
}


merge_range_64_t compute_merge_range_long(int64_t a_count, int64_t b_count, int64_t partition, int64_t spacing, int64_t mp0, int64_t mp1) {
    int64_t diag0 = spacing * partition;
    int64_t diag1 = min(a_count + b_count, diag0 + spacing);

    return merge_range_64_t { mp0, mp1, diag0 - mp0, diag1 - mp1 };
}
