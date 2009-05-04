/* This source file is our test case for sse-2 support. If this can be
 * compiled, then we're happy because sse-2 will be usable in cantor/ ;
 * it also means that SSE_NORM_INIT will be usable in sieve/las, although
 * the latter is currently disabled because it's slower than the
 * mainstream code anyway (!).
 */

#include <stdint.h>
#include <emmintrin.h>

__v2di x;
__v2df g;

void main() {
    x = (__v2di) { (uint64_t) 42, (uint64_t) 17 };
    g = (__v2df) { 42.0, 17.0 };
    __v2di shift = { (uint64_t) 42, (uint64_t) 17 };
    x = _mm_srl_epi64(x, shift);
}

