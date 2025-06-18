#include "types.h"
#include "meta.h"

range_t get_tile(int cta, int nv, int count) {
    return range_t { nv * cta, min(count, nv * (cta + 1)) };
}
