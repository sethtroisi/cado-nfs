/* This source file is our test case for ssse3 support. */
#include <stdint.h>
#include <string.h>
#include <tmmintrin.h>

int main() {
    __m128i x = _mm_setr_epi64(
            (__m64) UINT64_C(0x0706050403020100),
            (__m64) UINT64_C(0x0F0E0D0C0B0A0908));
    __m128i y = _mm_setr_epi64(
            (__m64) UINT64_C(0x1716151413121110),
            (__m64) UINT64_C(0x1F1E1D1C1B1A1918));
    uint64_t a[2], b[2] = { 0x0C0B0A0908070605, 0x14131211100F0E0D };
    y = _mm_alignr_epi8(y, x, 0x5);
    memcpy (a, &y, 16);
    return(a[0] != b[0] || a[1] != b[1]);
}

