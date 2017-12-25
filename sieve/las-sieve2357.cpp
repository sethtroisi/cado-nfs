#include "cado.h"
#include <cstddef>
#include <cstdint>

#ifdef HAVE_SSSE3
#include "tmmintrin.h"
#endif
#ifdef HAVE_AVX2
#include "immintrin.h"
#endif
#ifdef HAVE_ARM_NEON
#include <arm_neon.h>
#endif

#include "ularith.h"
#include "cado-endian.h"
#include "las-sieve2357.hpp"

/* Shift "b" so that each byte that does not get shifted out moves to a
   memory address that is "idx" bytes higher */
template<typename T>
static inline T shift_in_memory(T b, int idx)
{
#if defined(CADO_LITTLE_ENDIAN)
    return b << (idx * 8);
#elif defined(CADO_BIG_ENDIAN)
    return b >> (idx * 8);
#else
#error PDP endianness not implemented
#endif
}

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
static SIMDTYPE set1(ELEMTYPE);

template <typename SIMDTYPE, typename ELEMTYPE>
static SIMDTYPE adds(SIMDTYPE, SIMDTYPE);

template <typename SIMDTYPE, typename ELEMTYPE>
static SIMDTYPE _and(SIMDTYPE, SIMDTYPE);

template <typename SIMDTYPE>
static SIMDTYPE loadu(const SIMDTYPE *);

template<typename SIMDTYPE, typename ELEMTYPE>
inline SIMDTYPE ATTRIBUTE((__always_inline__, __artificial__))
adds(const SIMDTYPE a, const SIMDTYPE b)
{
  return a + b; /* WARNING, we assume no carry between elements here! */
}

template<typename SIMDTYPE, typename ELEMTYPE>
inline SIMDTYPE ATTRIBUTE((__always_inline__, __artificial__))
set0()
{
  return 0;
}

template<typename SIMDTYPE, typename ELEMTYPE>
static SIMDTYPE set1(const ELEMTYPE);

template<>
inline uint32_t ATTRIBUTE((__always_inline__, __artificial__))
set1<uint32_t, uint8_t>(const uint8_t c)
{
  uint32_t r = (uint32_t) c;
  r += r << 8;
  r += r << 16;
  return r;
}

template<>
inline uint64_t ATTRIBUTE((__always_inline__, __artificial__))
set1<uint64_t, uint8_t>(const uint8_t c)
{
  uint64_t r = (uint32_t) c;
  r += r << 8;
  r += r << 16;
  r += r << 32;
  return r;
}

template<typename SIMDTYPE, typename ELEMTYPE>
inline SIMDTYPE ATTRIBUTE((__always_inline__, __artificial__))
_and(const SIMDTYPE a, const SIMDTYPE b)
{
  return a & b;
}

template<typename T>
inline T ATTRIBUTE((__always_inline__, __artificial__))
loadu(const T *p)
{
  return *p;
}

static const uint8_t ff = ~(uint8_t)0;
static const uint8_t mask2[64] =  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                   ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0};
static const uint8_t mask4[64] =  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                   ff, 0, 0, 0, ff, 0, 0, 0, ff, 0, 0, 0, ff, 0, 0, 0, ff, 0, 0, 0, ff, 0, 0, 0, ff, 0, 0, 0, ff, 0, 0, 0};
static const uint8_t mask8[64] =  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                   ff, 0, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0, 0, 0, 0, 0};
static const uint8_t mask16[64] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                   ff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
static const uint8_t mask32[64] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                   ff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
static const uint8_t mask3[64] =  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                   ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0};
static const uint8_t mask5[64] =  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                   ff, 0, 0, 0, 0, ff, 0, 0, 0, 0, ff, 0, 0, 0, 0, ff, 0, 0, 0, 0, ff, 0, 0, 0, 0, ff, 0, 0, 0, 0, ff, 0};
static const uint8_t mask7[64] =  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                   ff, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0};

#ifdef HAVE_SSSE3

template<>
inline __m128i ATTRIBUTE((__always_inline__, __artificial__))
set0<__m128i, uint8_t>()
{
  return _mm_setzero_si128 ();
}

template<>
inline __m128i ATTRIBUTE((__always_inline__, __artificial__))
set1<__m128i, uint8_t>(const uint8_t c)
{
  return _mm_set1_epi8(c);
}

template<>
inline __m128i ATTRIBUTE((__always_inline__, __artificial__))
adds<__m128i, uint8_t>(const __m128i a, const __m128i b)
{
  return _mm_adds_epu8(a, b);
}

template<>
inline __m128i ATTRIBUTE((__always_inline__, __artificial__))
_and<__m128i, uint8_t>(const __m128i a, const __m128i b)
{
  return _mm_and_si128(a, b);
}

template<>
inline __m128i ATTRIBUTE((__always_inline__, __artificial__))
loadu<__m128i>(const __m128i *p)
{
  return _mm_loadu_si128(p);
}
#endif

#ifdef HAVE_AVX2

template<>
inline __m256i ATTRIBUTE((__always_inline__, __artificial__))
set0<__m256i, uint8_t>()
{
  return _mm256_setzero_si256 ();
}

template<>
inline __m256i ATTRIBUTE((__always_inline__, __artificial__))
set1<__m256i, uint8_t>(const uint8_t c)
{
  __m128i t = _mm_cvtsi64x_si128(c);
  return _mm256_broadcastb_epi8(t);
}

template<>
inline __m256i ATTRIBUTE((__always_inline__, __artificial__))
adds<__m256i, uint8_t>(const __m256i a, const __m256i b)
{
  return _mm256_adds_epu8(a, b);
}

template<>
inline __m256i ATTRIBUTE((__always_inline__, __artificial__))
_and<__m256i, uint8_t>(const __m256i a, const __m256i b)
{
  return _mm256_and_si256(a, b);
}

template<>
inline __m256i ATTRIBUTE((__always_inline__, __artificial__))
loadu<__m256i>(const __m256i *p)
{
  return _mm256_loadu_si256(p);
}
#endif

#ifdef HAVE_ARM_NEON

template<>
inline uint8x16_t ATTRIBUTE((__always_inline__, __artificial__))
set0<uint8x16_t, uint8_t>()
{
  return vdupq_n_u8(0);
}

template<>
inline uint8x16_t ATTRIBUTE((__always_inline__, __artificial__))
set1<uint8x16_t, uint8_t>(const uint8_t c)
{
   return vdupq_n_u8(c);
}

template<>
inline uint8x16_t ATTRIBUTE((__always_inline__, __artificial__))
adds<uint8x16_t, uint8_t>(const uint8x16_t a, const uint8x16_t b)
{
  return vqaddq_u8(a, b);
}

template<>
inline uint8x16_t ATTRIBUTE((__always_inline__, __artificial__))
_and<uint8x16_t, uint8_t>(const uint8x16_t a, const uint8x16_t b)
{
  return vandq_u8(a, b);
}

template<>
inline uint8x16_t ATTRIBUTE((__always_inline__, __artificial__))
loadu<uint8x16_t>(const uint8x16_t *p)
{
  return vld1q_u8(p);
}
#endif

/* A demultiplexer that returns the correct mask array for a given "STRIDE"
   value */
template <typename SIMDTYPE, typename ELEMTYPE, unsigned int STRIDE>
static inline const ELEMTYPE *get_mask() {
    switch (STRIDE) {
        case 2: return mask2;
        case 4: return mask4;
        case 8: return mask8;
        case 16: return mask16;
        case 32: return mask32;
        case 3: return mask3;
        case 5: return mask5;
        case 7: return mask7;
        default: abort();
    }
}

/* Return a mask of stride "STRIDE", shifted by "shift", via an unaligned
   memory read */
template <typename SIMDTYPE, typename ELEMTYPE, unsigned int STRIDE>
static inline SIMDTYPE
get_shifted_mask(const fbprime_t shift)
{
    const ELEMTYPE * mask = get_mask<SIMDTYPE, ELEMTYPE, STRIDE>();
    const ELEMTYPE * p = &mask[32 - shift];
    return loadu<SIMDTYPE>((const SIMDTYPE *) p);
}

/* Return a sieving pattern of stride "STRIDE", shifted by "shift", with
   value "elem" in locations where it hits */
template <typename SIMDTYPE, typename ELEMTYPE, unsigned int STRIDE>
static inline SIMDTYPE
get_pattern(const fbprime_t offset, ELEMTYPE elem)
{
    const SIMDTYPE shifted_mask = get_shifted_mask<SIMDTYPE, ELEMTYPE, STRIDE>(offset);
    return _and<SIMDTYPE, ELEMTYPE>(shifted_mask, set1<SIMDTYPE, ELEMTYPE>(elem));
}

template <typename SIMDTYPE, typename ELEMTYPE, unsigned int STRIDE>
static inline void
sieve_odd_prime(SIMDTYPE * const result, const ELEMTYPE logp, const fbprime_t idx)
{
    fbprime_t offset = idx;
    for (size_t i = 0; i < STRIDE; i++) {
        const SIMDTYPE pattern = get_pattern<SIMDTYPE, ELEMTYPE, STRIDE>(offset, logp);
        result[i] = adds<SIMDTYPE, ELEMTYPE>(result[i], pattern);
        offset = modsub(offset, sizeof(SIMDTYPE) % STRIDE, STRIDE);
    }
}

template <typename SIMDTYPE, typename ELEMTYPE>
static inline void
big_loop(SIMDTYPE * __restrict__ sievearray, const SIMDTYPE * const __restrict__ pattern235,
         const SIMDTYPE * const __restrict__ pattern7, size_t l7)
{
  size_t im15 = 0;
  for (size_t i = l7; i > 0; i--) {
    *sievearray++ = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15], pattern7[0]);
    *sievearray++ = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15 + 1], pattern7[1]);
    *sievearray++ = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15 + 2], pattern7[2]);
    *sievearray++ = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15 + 3], pattern7[3]);
    *sievearray++ = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15 + 4], pattern7[4]);
    *sievearray++ = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15 + 5], pattern7[5]);
    *sievearray++ = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15 + 6], pattern7[6]);
    if (im15 >= 8) {im15 -= 8;} else {im15 += 7;}
  }
}

#if 0 && defined(HAVE_AVX2) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
template <>
inline void
big_loop<__m256i, uint8_t>(__m256i * __restrict__ sievearray, const __m256i * __restrict__ pattern235,
         const __m256i * const __restrict__ pattern7, const size_t l7)
{
  __m256i r, pattern235_0 = pattern235[0],
             pattern235_1 = pattern235[1],
             pattern235_2 = pattern235[2],
             pattern235_3 = pattern235[3],
             pattern235_4 = pattern235[4],
             pattern235_5 = pattern235[5],
             pattern235_6 = pattern235[6];

  const __m256i * pattern235_plus_7 = pattern235 + 7;
  const __m256i * const pattern235_plus_8 = pattern235 + 8;

  if (l7 == 0)
      return;

#define LOOP_BLOCK(i) \
               "vpaddusb %[p235_"CADO_STRINGIZE(i)"], %[p7_"CADO_STRINGIZE(i)"], %[r] \n\t"	/* p015 */ \
               "vmovdqa "CADO_STRINGIZE(i)"*32(%[pattern235]), %[p235_"CADO_STRINGIZE(i)"] \n\t"/* p23  */ \
               "vmovdqa %[r], "CADO_STRINGIZE(i)"*32(%[sievearray]) \n\t"			/* p237 p4 */

  for (size_t i = l7 - 1; i > 0; i--) {
      size_t t;
      __asm__ volatile (
               LOOP_BLOCK(0)
               "cmp %[pp8], %[pattern235] \n\t"				// p0156 carry iff pp8 > pattern235
               LOOP_BLOCK(1)
               LOOP_BLOCK(2)
               LOOP_BLOCK(3)
               LOOP_BLOCK(4)
               LOOP_BLOCK(5)
               LOOP_BLOCK(6)
               "lea 7*32(%[pattern235]), %[t] \n\t"			// p15
               "lea -8*32(%[pattern235]), %[pattern235] \n\t"		// p15
               "cmovc %[t], %[pattern235] \n\t"				// p06
               "add $7*32, %[sievearray] \n\t"				// p0156
                                                                        // p6 for dec/jnz
               : [sievearray] "+r" (sievearray), [r] "=&x" (r), [t] "=&r" (t),
                 [pattern235] "+r" (pattern235_plus_7), [p235_0] "+x" (pattern235_0),
                 [p235_1] "+x" (pattern235_1), [p235_2] "+x" (pattern235_2),
                 [p235_3] "+x" (pattern235_3), [p235_4] "+x" (pattern235_4),
                 [p235_5] "+x" (pattern235_5), [p235_6] "+x" (pattern235_6)
               : [pp8] "r" (pattern235_plus_8), [p7_0] "x" (pattern7[0]),
                 [p7_1] "x" (pattern7[1]), [p7_2] "x" (pattern7[2]),
                 [p7_3] "x" (pattern7[3]), [p7_4] "x" (pattern7[4]),
                 [p7_5] "x" (pattern7[5]), [p7_6] "x" (pattern7[6])
               : "cc" //, "memory"
      );
  }
  *sievearray++ = _mm256_adds_epu8(pattern235_0, pattern7[0]);
  *sievearray++ = _mm256_adds_epu8(pattern235_1, pattern7[1]);
  *sievearray++ = _mm256_adds_epu8(pattern235_2, pattern7[2]);
  *sievearray++ = _mm256_adds_epu8(pattern235_3, pattern7[3]);
  *sievearray++ = _mm256_adds_epu8(pattern235_4, pattern7[4]);
  *sievearray++ = _mm256_adds_epu8(pattern235_5, pattern7[5]);
  *sievearray++ = _mm256_adds_epu8(pattern235_6, pattern7[6]);
}
#endif

template <typename SIMDTYPE, typename ELEMTYPE>
inline
SIMDTYPE sieve2(const fbprime_t q, const fbprime_t idx, const uint8_t logp)
{
    const size_t N MAYBE_UNUSED = sizeof(SIMDTYPE) / sizeof(ELEMTYPE);
    ASSERT(q <= N);
    switch (q) {
        case 2: return get_pattern<SIMDTYPE, ELEMTYPE, 2>(idx, logp);
        case 4: return get_pattern<SIMDTYPE, ELEMTYPE, 4>(idx, logp);
        case 8: return get_pattern<SIMDTYPE, ELEMTYPE, 8>(idx, logp);
        case 16: return get_pattern<SIMDTYPE, ELEMTYPE, 16>(idx, logp);
        case 32: return get_pattern<SIMDTYPE, ELEMTYPE, 32>(idx, logp);
        default: abort();
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
    pattern2 = adds<SIMDTYPE, ELEMTYPE>(pattern2, set1<SIMDTYPE, ELEMTYPE>(primes->logp));
  }

  /* Sieve powers of 2 */
  for ( ; primes->p == 2 ; primes++) {
    pattern2 = adds<SIMDTYPE, ELEMTYPE>(pattern2, sieve2<SIMDTYPE, ELEMTYPE>(primes->q, primes->idx, primes->logp));
  }

  alignas(sizeof(SIMDTYPE)) 
  SIMDTYPE pattern23[3] = {pattern2, pattern2, pattern2};
  for ( ; primes->p == 3 ; primes++) {
    sieve_odd_prime<SIMDTYPE, ELEMTYPE, 3>(pattern23, primes->logp, primes->idx);
  }

  alignas(sizeof(SIMDTYPE)) 
  SIMDTYPE pattern235[15 + 6] = {
    pattern23[0], pattern23[1], pattern23[2],
    pattern23[0], pattern23[1], pattern23[2],
    pattern23[0], pattern23[1], pattern23[2],
    pattern23[0], pattern23[1], pattern23[2],
    pattern23[0], pattern23[1], pattern23[2],
    pattern23[0], pattern23[1], pattern23[2],
    pattern23[0], pattern23[1], pattern23[2],
  };
  alignas(sizeof(SIMDTYPE)) 
  SIMDTYPE pattern5[5] = {zero, zero, zero, zero, zero};
  for ( ; primes->p == 5 ; primes++) {
    sieve_odd_prime<SIMDTYPE, ELEMTYPE, 5>(pattern5, primes->logp, primes->idx);
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

  alignas(sizeof(SIMDTYPE))
  SIMDTYPE pattern7[7] = {zero, zero, zero, zero, zero, zero, zero};
  for ( ; primes->p == 7 ; primes++) {
    sieve_odd_prime<SIMDTYPE, ELEMTYPE, 7>(pattern7, primes->logp, primes->idx);
  }

  const size_t l7 = arraylen / 7 / N;
  
  big_loop<SIMDTYPE, ELEMTYPE>(sievearray, pattern235, pattern7, l7);

  for (size_t i = l7 * 7; i < arraylen / N; i++) {
    sievearray[i] = adds<SIMDTYPE, ELEMTYPE>(pattern235[i % 15], pattern7[i % 7]);
  }
}

template
void sieve2357<uint32_t, uint8_t>(uint32_t * const sievearray, const size_t arraylen, const sieve2357_prime_t *primes);
template
void sieve2357<uint64_t, uint8_t>(uint64_t * const sievearray, const size_t arraylen, const sieve2357_prime_t *primes);
#ifdef HAVE_SSSE3
template
void sieve2357<__m128i, uint8_t>(__m128i * const sievearray, const size_t arraylen, const sieve2357_prime_t *primes);
#endif
#ifdef HAVE_AVX2
template
void sieve2357<__m256i, uint8_t>(__m256i * const sievearray, const size_t arraylen, const sieve2357_prime_t *primes);
#endif
#ifdef HAVE_ARM_NEON
template
void sieve2357<uint8x16_t, uint8_t>(uint8x16_t * const sievearray, const size_t arraylen, const sieve2357_prime_t *primes);
#endif
