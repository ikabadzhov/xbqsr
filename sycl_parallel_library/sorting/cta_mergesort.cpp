#include "cta_mergesort.h"

int out_of_range_flags(int first, int vt, int count) {
    int out_of_range = min(vt, first + vt - count);
    int head_flags = 0;
    if(out_of_range > 0) {
        const int mask = (1<< vt) - 1;
        head_flags = mask & (~mask>> out_of_range);
    }
    return head_flags;
}

merge_range_t compute_mergesort_frame(int partition, int coop, int spacing){
    int size = spacing * (coop / 2);
    int start = ~(coop - 1) & partition;
    int a_begin = spacing * start;
    int b_begin = spacing * start + size;

    return merge_range_t {
            a_begin,
            a_begin + size,
            b_begin,
            b_begin + size
    };
}

merge_range_t compute_mergesort_range(int count, int partition, int coop, int spacing) {

    merge_range_t frame = compute_mergesort_frame(partition, coop, spacing);

    return merge_range_t {
            frame.a_begin,
            min(count, frame.a_end),
            min(count, frame.b_begin),
            min(count, frame.b_end)
    };
}

merge_range_t compute_mergesort_range(int count, int partition, int coop, int spacing, int mp0, int mp1) {

    merge_range_t range = compute_mergesort_range(count, partition, coop, spacing);

    // Locate the diagonal from the start of the A sublist.
    int diag = spacing * partition - range.a_begin;

    // The end partition of the last cta for each merge operation is computed
    // and stored as the begin partition for the subsequent merge. i.e. it is
    // the same partition but in the wrong coordinate system, so its 0 when it
    // should be listSize. Correct that by checking if this is the last cta
    // in this merge operation.
    if(coop - 1 != ((coop - 1) & partition)) {
        range.a_end = range.a_begin + mp1;
        range.b_end = min(count, range.b_begin + diag + spacing - mp1);
    }

    range.a_begin = range.a_begin + mp0;
    range.b_begin = min(count, range.b_begin + diag - mp0);

    return range;
}
