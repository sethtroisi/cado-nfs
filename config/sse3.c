/* This source file is our test case for sse-3 support. */
#include <stdint.h>
#include <string.h>
#include <pmmintrin.h>

__v2df x = (__v2df) { 12.34, 34.12 }, y =(__v2df) { 78.56, 56.78 };
double a[2], b[2] = { 78.56+56.78, 12.34+34.12 };

int main() {
  y = _mm_hadd_pd(y, x);
  memcpy(a, &y, 16);
  return(a[0] != b[0] || a[1] != b[1]);
}

