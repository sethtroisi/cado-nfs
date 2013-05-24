/* This source file is our test case for ssse3 support. */
#include <stdint.h>
#include <string.h>
#include <tmmintrin.h>

__v2di x = { 0x0706050403020100, 0x0F0E0D0C0B0A0908 },
  y = { 0x1716151413121110, 0x1F1E1D1C1B1A1918 };
uint64_t a[2], b[2] = { 0x0C0B0A0908070605, 0x14131211100F0E0D };

int main() {
  y = _mm_alignr_epi8(y, x, 0x5);
  memcpy (a, &y, 16);
  return(a[0] != b[0] || a[1] != b[1]);
}

