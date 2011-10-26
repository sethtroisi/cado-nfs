/* This source file is our test case for sse-2 support. If this can be
 * compiled, then we're happy because sse-2 will be usable in cantor/ ;
 * it also means that SSE_NORM_INIT will be usable in sieve/las.
 */

#include <stdint.h>
#include <emmintrin.h>

volatile __v2di x;
volatile __v2df g;

int main() {
    x = (__v2di) { (uint64_t) 42, (uint64_t) 17 };
    g = (__v2df) { 42.0, 17.0 };
    __v2di shift = { (uint64_t) 42, (uint64_t) 17 };
    x = _mm_srl_epi64(x, shift);
    return 0;
}

