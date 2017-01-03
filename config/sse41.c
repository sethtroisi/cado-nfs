#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <smmintrin.h>

int main(int argc, char * argv[])
{
    /* the following test is for emulated 32-bit on physical 64-bit */
    if (sizeof(unsigned long) != 8)
      abort ();
    volatile int a0 = 17;
    volatile int a1 = 42;
    __m128i x = _mm_setr_epi32(a1, 0, a0, 0);
    // x = 0 0x2a 0 0x11
    __m128i y = _mm_setr_epi32(42, 0, 17, 0);
    // y = 0 0x2a 0 0x11
    __m128i ma = _mm_max_epi32(x, y);
    __m128i mi = _mm_min_epi32(x, y);
    __m128i z = _mm_xor_si128(mi, ma);
    int ok0 = _mm_testz_si128(z, z);
    __m128i c = _mm_cmpeq_epi64(x, y);
    int ok1 = _mm_extract_epi32(c, 0) && _mm_extract_epi32(c, 1);
    return (ok0 && ok1) ? EXIT_SUCCESS : EXIT_FAILURE;
}
