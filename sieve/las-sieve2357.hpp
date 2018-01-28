#ifndef LAS_SIEVE2357_HPP
#define LAS_SIEVE2357_HPP

#include "las-smallsieve-types.hpp"
#include "las-debug.hpp"

namespace sieve2357 {

#if defined(HAVE_AVX2)
typedef __m256i preferred_simd_type;
#elif defined(HAVE_SSSE3)
typedef __m128i preferred_simd_type;
#elif defined(HAVE_ARM_NEON)
typedef uint8x16_t preferred_simd_type;
#elif LONG_BIT == 64 /* FIXME: this is false on, e.g., MinGW. What's a better condition here? */
typedef uint64_t preferred_simd_type;
#else
typedef uint32_t preferred_simd_type;
#endif

/* We hit sievearray[x] for x = idx + i*q with 0 <= x < arraylen */
struct prime_t {
  fbprime_t q, idx;
  unsigned char logp;
  /* Make prime_t sortable in the order in which sieve2357 expects them:
     first q=1, then powers of 2, then odd primes in increasing order */
  bool operator<(const sieve2357::prime_t &other) const {
    if (q == 1 && other.q != 1)
      return true;
    if (q != 1 && other.q == 1)
      return false;
    if (q % 2 == 0 && other.q % 2 == 1)
      return true;
    if (q % 2 == 1 && other.q % 2 == 0)
      return false;
    return q < other.q;
  }
};

/* A predicate that tells whether a prime power q = p^k can be sieved by
   sieve2357 with a given SIMD and element data type */
template<typename SIMDTYPE, typename ELEMTYPE>
static inline bool can_sieve(const fbprime_t q)
{
  const size_t N = sizeof(SIMDTYPE) / sizeof(ELEMTYPE);
  return q == 1 || (q % 2 == 0 && N % q  == 0) || q == 3 || q == 5
      || q == 7;
}

enum {
    update_set,
    update_add
};

template <typename SIMDTYPE, typename ELEMTYPE>
void sieve(SIMDTYPE *sievearray, size_t arraylen, const prime_t *primes,
    bool only_odd, int update_operation, where_am_I & w);

}
#endif
