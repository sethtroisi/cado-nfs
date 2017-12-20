#include "cado.h"
#include <cstddef>
#include <cstdint>

#ifdef HAVE_SSSE3
#include "tmmintrin.h"
#endif
#ifdef HAVE_AVX2
#include "immintrin.h"
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

/* This class defines constants of type SIMDTYPE, used as a mask for sieving
   patterns, depending on template parameters. These are prefixed by zero-
   valued entries so that a shifted mask can be read by an (unaligned)
   memory access. E.g., for a SIMD type of 8 one-byte elements with STRIDE=3,
   the constants will contain, in order of increasing memory address, the
   values 0, 0, 0, 0, 0, 0, 0, 0,
   0xff, 0, 0, 0xff, 0, 0, 0xff, 0, 0 */
template<typename SIMDTYPE, typename ELEMTYPE, unsigned int STRIDE>
class patterns {
public:
    static const SIMDTYPE pattern[2];

    static SIMDTYPE
    get_shifted_mask(const fbprime_t shift) {
        const ELEMTYPE * p = (const ELEMTYPE *)(&pattern[1]) - shift;
        return loadu<SIMDTYPE>((const SIMDTYPE *) p);
   }
};

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

template<typename SIMDTYPE, typename ELEMTYPE>
inline SIMDTYPE ATTRIBUTE((__always_inline__, __artificial__))
set1(const ELEMTYPE c)
{
  SIMDTYPE r = (SIMDTYPE) c;
  if (sizeof(SIMDTYPE) > sizeof(ELEMTYPE))
      r += r << (8*sizeof(ELEMTYPE));
  if (sizeof(SIMDTYPE) > 2*sizeof(ELEMTYPE))
      r += r << (2*8*sizeof(ELEMTYPE));
  if (sizeof(SIMDTYPE) > 4*sizeof(ELEMTYPE))
      r += r << (4*8*sizeof(ELEMTYPE));
  if (sizeof(SIMDTYPE) > 8*sizeof(ELEMTYPE))
      r += r << (8*8*sizeof(ELEMTYPE));
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

/* Here come a crapload of magic constants to set up the patterns.
   Maybe there is a way to have the compiler generate them rather than
   defining each one individually. */
#if defined(CADO_LITTLE_ENDIAN)
template <>
const uint32_t patterns<uint32_t, unsigned char, 2>::pattern[2] = {
    0, UINT32_C(0x00ff00ff)
};
template <>
const uint32_t patterns<uint32_t, unsigned char, 4>::pattern[2] = {
    0, UINT32_C(0x000000ff)
};
template <>
const uint32_t patterns<uint32_t, unsigned char, 3>::pattern[2] = {
    0, UINT32_C(0xff0000ff)
};
template <>
const uint32_t patterns<uint32_t, unsigned char, 5>::pattern[2] = {
    0, UINT32_C(0x000000ff)
};
template <>
const uint32_t patterns<uint32_t, unsigned char, 7>::pattern[2] = {
    0, UINT32_C(0x000000ff)
};
template <>
const uint64_t patterns<uint64_t, unsigned char, 2>::pattern[2] = {
    0, UINT64_C(0x00ff00ff00ff00ff)
};
template <>
const uint64_t patterns<uint64_t, unsigned char, 4>::pattern[2] = {
    0, UINT64_C(0x000000ff000000ff)
};
template <>
const uint64_t patterns<uint64_t, unsigned char, 8>::pattern[2] = {
    0, UINT64_C(0x00000000000000ff)
};
template <>
const uint64_t patterns<uint64_t, unsigned char, 3>::pattern[2] = {
    0, UINT64_C(0x00ff0000ff0000ff)
};
template <>
const uint64_t patterns<uint64_t, unsigned char, 5>::pattern[2] = {
    0, UINT64_C(0x0000ff00000000ff)
};
template <>
const uint64_t patterns<uint64_t, unsigned char, 7>::pattern[2] = {
    0, UINT64_C(0xff000000000000ff)
};
#else
template <>
const uint32_t patterns<uint32_t, unsigned char, 2>::pattern[2] = {
    0, UINT32_C(0xff00ff00)
};
template <>
const uint32_t patterns<uint32_t, unsigned char, 4>::pattern[2] = {
    0, UINT32_C(0xff000000)
};
template <>
const uint32_t patterns<uint32_t, unsigned char, 3>::pattern[2] = {
    0, UINT32_C(0xff0000ff)
};
template <>
const uint32_t patterns<uint32_t, unsigned char, 5>::pattern[2] = {
    0, UINT32_C(0xff000000)
};
template <>
const uint32_t patterns<uint32_t, unsigned char, 7>::pattern[2] = {
    0, UINT32_C(0xff000000)
};
template <>
const uint64_t patterns<uint64_t, unsigned char, 2>::pattern[2] = {
    0, UINT64_C(0xff00ff00ff00ff00)
};
template <>
const uint64_t patterns<uint64_t, unsigned char, 4>::pattern[2] = {
    0, UINT64_C(0xff000000ff000000)
};
template <>
const uint64_t patterns<uint64_t, unsigned char, 8>::pattern[2] = {
    0, UINT64_C(0xff00000000000000)
};
template <>
const uint64_t patterns<uint64_t, unsigned char, 3>::pattern[2] = {
    0, UINT64_C(0xff0000ff0000ff00)
};
template <>
const uint64_t patterns<uint64_t, unsigned char, 5>::pattern[2] = {
    0, UINT64_C(0xff00000000ff0000)
};
template <>
const uint64_t patterns<uint64_t, unsigned char, 7>::pattern[2] = {
    0, UINT64_C(0xff000000000000ff)
};
#endif

#ifdef HAVE_SSSE3

template<>
inline __m128i ATTRIBUTE((__always_inline__, __artificial__))
set0<__m128i, unsigned char>()
{
  return _mm_setzero_si128 ();
}

template<>
inline __m128i ATTRIBUTE((__always_inline__, __artificial__))
set1<__m128i, unsigned char>(const unsigned char c)
{
  return _mm_set1_epi8(c);
}

template<>
inline __m128i ATTRIBUTE((__always_inline__, __artificial__))
adds<__m128i, unsigned char>(const __m128i a, const __m128i b)
{
  return _mm_adds_epu8(a, b);
}

template<>
inline __m128i ATTRIBUTE((__always_inline__, __artificial__))
_and<__m128i, unsigned char>(const __m128i a, const __m128i b)
{
  return _mm_and_si128(a, b);
}

template<>
inline __m128i ATTRIBUTE((__always_inline__, __artificial__))
loadu<__m128i>(const __m128i *p)
{
  return _mm_loadu_si128(p);
}

static const unsigned char ff = ~(unsigned char)0;
template <>
const __m128i patterns<__m128i, unsigned char, 2>::pattern[2] = {
        _mm_set1_epi8(0),
        _mm_setr_epi8(ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0)
};
template <>
const __m128i patterns<__m128i, unsigned char, 4>::pattern[2] = {
        _mm_set1_epi8(0),
        _mm_setr_epi8(ff, 0, 0, 0, ff, 0, 0, 0, ff, 0, 0, 0, ff, 0, 0, 0)
};
template <>
const __m128i patterns<__m128i, unsigned char, 8>::pattern[2] = {
        _mm_set1_epi8(0),
        _mm_setr_epi8(ff, 0, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0, 0, 0, 0, 0)
};
template <>
const __m128i patterns<__m128i, unsigned char, 16>::pattern[2] = {
        _mm_set1_epi8(0),
        _mm_setr_epi8(ff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
};
template <>
const __m128i patterns<__m128i, unsigned char, 3>::pattern[2] = {
        _mm_set1_epi8(0),
        _mm_setr_epi8(ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0, 0, ff)
};
template <>
const __m128i patterns<__m128i, unsigned char, 5>::pattern[2] = {
        _mm_set1_epi8(0),
        _mm_setr_epi8(ff, 0, 0, 0, 0, ff, 0, 0, 0, 0, ff, 0, 0, 0, 0, ff)
};
template <>
const __m128i patterns<__m128i, unsigned char, 7>::pattern[2] = {
        _mm_set1_epi8(0),
        _mm_setr_epi8(ff, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0, 0, 0, 0, ff, 0)
};
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
set1<__m256i, unsigned char>(const unsigned char c)
{
  __m128i t = _mm_cvtsi64x_si128(c);
  return _mm256_broadcastb_epi8(t);
}

template<>
inline __m256i ATTRIBUTE((__always_inline__, __artificial__))
adds<__m256i, unsigned char>(const __m256i a, const __m256i b)
{
  return _mm256_adds_epu8(a, b);
}

template<>
inline __m256i ATTRIBUTE((__always_inline__, __artificial__))
_and<__m256i, unsigned char>(const __m256i a, const __m256i b)
{
  return _mm256_and_si256(a, b);
}

template<>
inline __m256i ATTRIBUTE((__always_inline__, __artificial__))
loadu<__m256i>(const __m256i *p)
{
  return _mm256_loadu_si256(p);
}

template <>
const __m256i patterns<__m256i, unsigned char, 2>::pattern[2] = {
        _mm256_set1_epi8(0),
        _mm256_setr_epi8(ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0)
};
template <>
const __m256i patterns<__m256i, unsigned char, 4>::pattern[2] = {
        _mm256_set1_epi8(0),
        _mm256_setr_epi8(ff, 0, 0, 0, ff, 0, 0, 0, ff, 0, 0, 0, ff, 0, 0, 0, ff, 0, 0, 0, ff, 0, 0, 0, ff, 0, 0, 0, ff, 0, 0, 0)
};
template <>
const __m256i patterns<__m256i, unsigned char, 8>::pattern[2] = {
        _mm256_set1_epi8(0),
        _mm256_setr_epi8(ff, 0, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0, 0, 0, 0, 0)
};
template <>
const __m256i patterns<__m256i, unsigned char, 16>::pattern[2] = {
        _mm256_set1_epi8(0),
        _mm256_setr_epi8(ff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
};
template <>
const __m256i patterns<__m256i, unsigned char, 32>::pattern[2] = {
        _mm256_set1_epi8(0),
        _mm256_setr_epi8(ff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
};
template <>
const __m256i patterns<__m256i, unsigned char, 3>::pattern[2] = {
    _mm256_set1_epi8(0),
    _mm256_setr_epi8(ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0)
};
template <>
const __m256i patterns<__m256i, unsigned char, 5>::pattern[2] = {
    _mm256_set1_epi8(0),
    _mm256_setr_epi8(ff, 0, 0, 0, 0, ff, 0, 0, 0, 0, ff, 0, 0, 0, 0, ff, 0, 0, 0, 0, ff, 0, 0, 0, 0, ff, 0, 0, 0, 0, ff, 0)
};
template <>
const __m256i patterns<__m256i, unsigned char, 7>::pattern[2] = {
    _mm256_set1_epi8(0),
    _mm256_setr_epi8(ff, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0)
};
#endif

template <typename SIMDTYPE, typename ELEMTYPE, unsigned int STRIDE>
static inline SIMDTYPE
get_shifted_pattern(const fbprime_t offset, ELEMTYPE elem)
{
    const SIMDTYPE mask = patterns<SIMDTYPE, ELEMTYPE, STRIDE>::get_shifted_mask(offset);
    return _and<SIMDTYPE, ELEMTYPE>(mask, set1<SIMDTYPE, ELEMTYPE>(elem));
}

template <typename SIMDTYPE, typename ELEMTYPE, unsigned int STRIDE>
static inline void
sieve_one_prime(SIMDTYPE * const result, const ELEMTYPE logp, const fbprime_t idx)
{
    fbprime_t offset = idx;
    for (size_t i = 0; i < STRIDE; i++) {
        const SIMDTYPE pattern = get_shifted_pattern<SIMDTYPE, ELEMTYPE, STRIDE>(offset, logp);
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
static
SIMDTYPE sieve2(fbprime_t, fbprime_t, ELEMTYPE);

template <>
inline
uint32_t sieve2<uint32_t, unsigned char>(const fbprime_t q, const fbprime_t idx, const unsigned char logp)
{
    switch (q) {
        case 2: return get_shifted_pattern<uint32_t, unsigned char, 2>(idx, logp);
        case 4: return get_shifted_pattern<uint32_t, unsigned char, 4>(idx, logp);
        default: abort();
    }
}

template <>
inline
uint64_t sieve2<uint64_t, unsigned char>(const fbprime_t q, const fbprime_t idx, const unsigned char logp)
{
    switch (q) {
        case 2: return get_shifted_pattern<uint64_t, unsigned char, 2>(idx, logp);
        case 4: return get_shifted_pattern<uint64_t, unsigned char, 4>(idx, logp);
        case 8: return get_shifted_pattern<uint64_t, unsigned char, 8>(idx, logp);
        default: abort();
    }
}

template <>
inline
__m128i sieve2<__m128i, unsigned char>(const fbprime_t q, const fbprime_t idx, const unsigned char logp)
{
    switch (q) {
        case 2: return get_shifted_pattern<__m128i, unsigned char, 2>(idx, logp);
        case 4: return get_shifted_pattern<__m128i, unsigned char, 4>(idx, logp);
        case 8: return get_shifted_pattern<__m128i, unsigned char, 8>(idx, logp);
        case 16: return get_shifted_pattern<__m128i, unsigned char, 16>(idx, logp);
        default: abort();
    }
}

template <>
inline
__m256i sieve2<__m256i, unsigned char>(const fbprime_t q, const fbprime_t idx, const unsigned char logp)
{
    switch (q) {
        case 2: return get_shifted_pattern<__m256i, unsigned char, 2>(idx, logp);
        case 4: return get_shifted_pattern<__m256i, unsigned char, 4>(idx, logp);
        case 8: return get_shifted_pattern<__m256i, unsigned char, 8>(idx, logp);
        case 16: return get_shifted_pattern<__m256i, unsigned char, 16>(idx, logp);
        case 32: return get_shifted_pattern<__m256i, unsigned char, 32>(idx, logp);
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
    sieve_one_prime<SIMDTYPE, ELEMTYPE, 3>(pattern23, primes->logp, primes->idx);
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
    sieve_one_prime<SIMDTYPE, ELEMTYPE, 5>(pattern5, primes->logp, primes->idx);
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
    sieve_one_prime<SIMDTYPE, ELEMTYPE, 7>(pattern7, primes->logp, primes->idx);
  }

  const size_t l7 = arraylen / 7 / N;
  
  big_loop<SIMDTYPE, ELEMTYPE>(sievearray, pattern235, pattern7, l7);

  for (size_t i = l7 * 7; i < arraylen / N; i++) {
    sievearray[i] = adds<SIMDTYPE, ELEMTYPE>(pattern235[i % 15], pattern7[i % 7]);
  }
}

template
void sieve2357<unsigned long, unsigned char>(unsigned long * const sievearray, const size_t arraylen, const sieve2357_prime_t *primes);
#ifdef HAVE_SSSE3
template
void sieve2357<__m128i, unsigned char>(__m128i * const sievearray, const size_t arraylen, const sieve2357_prime_t *primes);
#endif
#ifdef HAVE_AVX2
template
void sieve2357<__m256i, unsigned char>(__m256i * const sievearray, const size_t arraylen, const sieve2357_prime_t *primes);
#endif
