#include <stdint.h>
#include <smmintrin.h>

int main() {
    __m128i x = _mm_setr_epi32(42, 0, 17, 0);
    __m128i y = _mm_setr_epi32(41, 0, 17, 0);
    x = _mm_cmpeq_epi64(x, y);
    /* the following test is for emulated 32-bit on physical 64-bit */
    if (sizeof(unsigned long) != 8)
      abort ();
    return 0;
}
