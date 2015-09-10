/* This source file is our test case for sse-2 support. If this can be
 * compiled, then we're happy because sse-2 will be usable in cantor/ ;
 * it also means that SSE_NORM_INIT will be usable in sieve/las.
 */

#include <stdint.h>
#include <emmintrin.h>

int main() {
    __m128i x = _mm_setr_epi64((__m64) (uint64_t) 42, (__m64) (uint64_t) 17);
    __m128d g = _mm_set_pd(42.0, 17.0);
    __m128i shift = _mm_setr_epi64((__m64) (uint64_t) 42, (__m64) (uint64_t) 17);
    x = _mm_srl_epi64(x, shift);
    return 0;
}

