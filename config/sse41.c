#include <stdint.h>
#include <stdlib.h>
#include <smmintrin.h>

int main() {
    __m128i x = _mm_setr_epi32(42, 0, 17, 0);
    __m128i y = _mm_setr_epi32(41, 0, 17, 0);
    x = _mm_cmpeq_epi64(x, y);
    // x = 0....0 1....1
    y =  _mm_insert_epi64(y, _mm_extract_epi64(x, 1), 0); 
    // y = 1....1 0x11
    y = _mm_andnot_si128(y, y);
    // y = 0....0 0x11
    x = _mm_cmpeq_epi64(x, y);
    // x = 1....1 1....1
    /* the following test is for emulated 32-bit on physical 64-bit */
    if (sizeof(unsigned long) != 8)
      abort ();
    unsigned int z = _mm_extract_epi32(x, 2);
    return
        (
       (_mm_extract_epi32(x, 0) != 0)
    && (_mm_extract_epi32(x, 1) != 0)
    && (_mm_extract_epi32(x, 2) == 0)
    && (_mm_extract_epi32(x, 3) == 0)) ? EXIT_SUCCESS : EXIT_FAILURE;
}
