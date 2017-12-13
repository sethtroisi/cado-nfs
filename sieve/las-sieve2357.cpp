#include "cado.h"
#include <cstddef>
#include <cstdint>

#ifdef HAVE_SSE2
#include "emmintrin.h"
#endif
#ifdef HAVE_AVX2
#include "immintrin.h"
#endif

#include "ularith.h"
#include "las-sieve2357.hpp"

static inline unsigned long ATTRIBUTE((__always_inline__, __artificial__))
modsub(const unsigned long a, const unsigned long b, const unsigned long m)
{
  unsigned long r;
  ularith_submod_ul_ul(&r, a, b, m);
  return r;
}

template <typename SIMDTYPE, typename ELEMTYPE>
static SIMDTYPE set0();

template <typename SIMDTYPE, typename ELEMTYPE>
static SIMDTYPE adds(SIMDTYPE, SIMDTYPE);

template<typename SIMDTYPE, typename ELEMTYPE>
static SIMDTYPE bcaststride_inl(ELEMTYPE, unsigned int, unsigned int);

template<>
inline unsigned long ATTRIBUTE((__always_inline__, __artificial__))
adds<unsigned long, unsigned char>(const unsigned long a, const unsigned long b)
{
  return a + b; /* WARNING, we assume no carry between elements here! */
}

template<>
inline unsigned long ATTRIBUTE((__always_inline__, __artificial__))
set0<unsigned long, unsigned char>()
{
  return 0UL;
}

template<>
inline unsigned long
bcaststride_inl<unsigned long, unsigned char>(const unsigned char v,
  const unsigned int offset, const unsigned int stride)
{
  unsigned long r = 0;
  for (unsigned int i = 0; i < sizeof(unsigned long); i += stride) {
    r += (unsigned long) v << (i * 8);
  }
  /* offset is usually not a compile-time constant; handle it separately so
     the above loop can be compile-time expanded */
  r <<= (offset * 8);
  return r;
}

/* A non-static wrapper for testing */
template<>
unsigned long bcaststride<unsigned long, unsigned char>(const unsigned char v,
  const unsigned int offset, const unsigned int stride)
{
  return bcaststride_inl<unsigned long, unsigned char>(v, offset, stride);
}

#ifdef HAVE_SSE2

template<>
inline __m128i ATTRIBUTE((__always_inline__, __artificial__))
set0<__m128i, unsigned char>()
{
  return _mm_setzero_si128 ();
}

template<>
inline __m128i ATTRIBUTE((__always_inline__, __artificial__))
adds<__m128i, unsigned char>(const __m128i a, const __m128i b)
{
  return _mm_adds_epu8(a, b);
}

static const __m128i shuffle_constants[17] {
  _mm_set1_epi8(128), /* Stride = 0: undefined, result gets set to 0 everywhere */
  _mm_set_epi8 ( 15,  14,  13,  12,  11,  10,   9,   8,   7,   6,   5,   4,   3,   2,   1, 0), /* 1 */
  _mm_set_epi8 (128,   7, 128,   6, 128,   5, 128,   4, 128,   3, 128,   2, 128,   1, 128, 0), /* 2 */
  _mm_set_epi8 (  5, 128, 128,   4, 128, 128,   3, 128, 128,   2, 128, 128,   1, 128, 128, 0), /* 3 */
  _mm_set_epi8 (128, 128, 128,   3, 128, 128, 128,   2, 128, 128, 128,   1, 128, 128, 128, 0), /* 4 */
  _mm_set_epi8 (  3, 128, 128, 128, 128,   2, 128, 128, 128, 128,   1, 128, 128, 128, 128, 0), /* 5 */
  _mm_set_epi8 (128, 128, 128,   2, 128, 128, 128, 128, 128,   1, 128, 128, 128, 128, 128, 0), /* 6 */
  _mm_set_epi8 (128,   2, 128, 128, 128, 128, 128, 128,   1, 128, 128, 128, 128, 128, 128, 0), /* 7 */
  _mm_set_epi8 (128, 128, 128, 128, 128, 128, 128,   1, 128, 128, 128, 128, 128, 128, 128, 0), /* 8 */
  _mm_set_epi8 (128, 128, 128, 128, 128, 128,   1, 128, 128, 128, 128, 128, 128, 128, 128, 0), /* 9 */
  _mm_set_epi8 (128, 128, 128, 128, 128,   1, 128, 128, 128, 128, 128, 128, 128, 128, 128, 0), /* 10 */
  _mm_set_epi8 (128, 128, 128, 128,   1, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 0), /* 11 */
  _mm_set_epi8 (128, 128, 128,   1, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 0), /* 12 */
  _mm_set_epi8 (128, 128,   1, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 0), /* 13 */
  _mm_set_epi8 (128,   1, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 0), /* 14 */
  _mm_set_epi8 (  1, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 0), /* 15 */
  _mm_set_epi8 (128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 0)  /* 16 */
};

static const __m128i shift_constants[2] = {
  _mm_set1_epi8(128),
  _mm_set_epi8(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0),
};

template<>
inline __m128i bcaststride_inl<__m128i, unsigned char>(const unsigned char v,
  const unsigned int offset, const unsigned int stride)
{
  __m128i r = _mm_set1_epi8(v);

  r = _mm_shuffle_epi8(r, shuffle_constants[stride]);

  const char *sc = (const char *)&shift_constants + 16 - offset;
  __m128i s2 = _mm_loadu_si128((__m128i *)sc);
  r = _mm_shuffle_epi8(r, s2);

  return r;
}

template<>
__m128i bcaststride<__m128i, unsigned char>(const unsigned char v,
  const unsigned int offset, const unsigned int stride)
{
  return bcaststride_inl<__m128i, unsigned char>(v, offset, stride);
}
#endif

#ifdef HAVE_AVX2

template<>
inline __m256i ATTRIBUTE((__always_inline__, __artificial__))
set0<__m256i, unsigned char>()
{
  return _mm256_setzero_si256 ();
}

template<>
inline __m256i ATTRIBUTE((__always_inline__, __artificial__))
adds<__m256i, unsigned char>(const __m256i a, const __m256i b)
{
  return _mm256_adds_epu8(a, b);
}

template<>
inline __m256i bcaststride_inl<__m256i, unsigned char>(const unsigned char v,
  const unsigned int offset, const unsigned int stride)
{
  /* In the low 16-byte word, we hit those locations x where x = k*stride + offset.
     In the high 16-byte word, we hit in those locations x where x = k*stride + offset - 16
     x == offset - 16 (mod stride) */

  __m128i r0, r1;
  
  r0 = bcaststride_inl<__m128i, unsigned char>(v, offset, stride);
  r1 = bcaststride_inl<__m128i, unsigned char>(v, modsub(offset, 16 % stride, stride), stride);

  return _mm256_insertf128_si256(_mm256_castsi128_si256(r0), r1, 1);
}

template<>
__m256i bcaststride<__m256i, unsigned char>(const unsigned char v,
  const unsigned int offset, const unsigned int stride)
{
  return bcaststride_inl<__m256i, unsigned char>(v, offset, stride);
}
#endif

template <typename SIMDTYPE, typename ELEMTYPE>
static inline void
sieve_one_prime(SIMDTYPE * const result, const size_t len, const ELEMTYPE logp, const unsigned long idx)
{
  for (size_t i = 0; i < len; i++) {
    unsigned int offset = modsub(idx, sizeof(SIMDTYPE) * i % len, len);
    result[i] = adds<SIMDTYPE, ELEMTYPE>(result[i], bcaststride_inl<SIMDTYPE, ELEMTYPE>(logp, offset, len));
  }
}

template <typename SIMDTYPE, typename ELEMTYPE>
void
sieve2357(SIMDTYPE * const sievearray, const size_t arraylen, const sieve2357_prime_t *primes)
{
  const SIMDTYPE zero = set0<SIMDTYPE, ELEMTYPE>();
  const size_t N = sizeof(SIMDTYPE) / sizeof(ELEMTYPE);
  SIMDTYPE pattern2 = zero;

  /* Sieve projective primes with q=1 */
  for ( ; primes->q == 1 ; primes++) {
    pattern2 = adds<SIMDTYPE, ELEMTYPE>(pattern2, bcaststride_inl<SIMDTYPE, ELEMTYPE>(primes->logp, 0, 1));
  }

  /* Initialise the first word of the pattern with the sieve entries that are
     powers of 2 */
  for ( ; primes->p == 2 ; primes++) {
    pattern2 = adds<SIMDTYPE, ELEMTYPE>(pattern2, bcaststride_inl<SIMDTYPE, ELEMTYPE>(primes->logp, primes->idx, primes->q));
  }

  SIMDTYPE pattern23[3] = {pattern2, pattern2, pattern2};
  for ( ; primes->p == 3 ; primes++) {
    sieve_one_prime(pattern23, 3, primes->logp, primes->idx);
  }

  SIMDTYPE pattern235[15 + 6] = {
    pattern23[0], pattern23[1], pattern23[2],
    pattern23[0], pattern23[1], pattern23[2],
    pattern23[0], pattern23[1], pattern23[2],
    pattern23[0], pattern23[1], pattern23[2],
    pattern23[0], pattern23[1], pattern23[2],
    pattern23[0], pattern23[1], pattern23[2],
    pattern23[0], pattern23[1], pattern23[2],
  };
  SIMDTYPE pattern5[5] = {zero, zero, zero, zero, zero};
  for ( ; primes->p == 5 ; primes++) {
    sieve_one_prime(pattern5, 5, primes->logp, primes->idx);
  }
 
  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 5; j++) {
      pattern235[i*5 + j] = adds<SIMDTYPE, ELEMTYPE>(pattern235[i*5 + j], pattern5[j]);
    }
  }

  /* Pad the list at the end so that we don't have to check for mod 15
     wrap-around after each increment. */
  pattern235[15 + 0] = pattern235[0];
  pattern235[15 + 1] = pattern235[1];
  pattern235[15 + 2] = pattern235[2];
  pattern235[15 + 3] = pattern235[3];
  pattern235[15 + 4] = pattern235[4];
  pattern235[15 + 5] = pattern235[5];

  SIMDTYPE pattern7[7] = {zero, zero, zero, zero, zero, zero, zero};
  for ( ; primes->p == 7 ; primes++) {
    sieve_one_prime(pattern7, 7, primes->logp, primes->idx);
  }

  const size_t l7 = (arraylen / 7 / N) * 7;
  size_t im15 = 0;

  for (size_t i = 0; i < l7; i += 7) {
    sievearray[i]     = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15],     pattern7[0]);
    sievearray[i + 1] = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15 + 1], pattern7[1]);
    sievearray[i + 2] = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15 + 2], pattern7[2]);
    sievearray[i + 3] = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15 + 3], pattern7[3]);
    sievearray[i + 4] = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15 + 4], pattern7[4]);
    sievearray[i + 5] = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15 + 5], pattern7[5]);
    sievearray[i + 6] = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15 + 6], pattern7[6]);
    im15 += 7;
    if (im15 >= 15) im15 -= 15;
  }
  for (size_t i = l7; i < arraylen / N; i++) {
    sievearray[i] = adds<SIMDTYPE, ELEMTYPE>(pattern235[i % 15], pattern7[i % 7]);
  }
}

template
void sieve2357<unsigned long, unsigned char>(unsigned long * const sievearray, const size_t arraylen, const sieve2357_prime_t *primes);
#ifdef HAVE_SSE2
template
void sieve2357<__m128i, unsigned char>(__m128i * const sievearray, const size_t arraylen, const sieve2357_prime_t *primes);
#endif
#ifdef HAVE_AVX2
template
void sieve2357<__m256i, unsigned char>(__m256i * const sievearray, const size_t arraylen, const sieve2357_prime_t *primes);
#endif
