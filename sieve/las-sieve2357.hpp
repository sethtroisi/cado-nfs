#ifndef LAS_SIEVE2357_HPP
#define LAS_SIEVE2357_HPP

#include "las-smallsieve-types.hpp"

typedef struct {
  fbprime_t p, q, idx;
  unsigned char logp;
} sieve2357_prime_t;

/* A predicate that tells whether a prime power q = p^k can be sieved by
   sieve2357 with a given SIMD and element data type */
template<typename SIMDTYPE, typename ELEMTYPE>
static bool sieve2357_can_sieve(const fbprime_t p, const fbprime_t q)
{
  return (p == 2 && q <= sizeof(SIMDTYPE) / sizeof(ELEMTYPE) && q <= 16)
    // currently, bcaststride() can't do stride >16 */
         || (p == 3 && q <= 3)
         || (p == 5 && q <= 5)
         || (p == 7 && q <= 7);
}

template<typename SIMDTYPE, typename ELEMTYPE>
SIMDTYPE bcaststride(const ELEMTYPE v, unsigned int offset,
    unsigned int stride);

template <typename SIMDTYPE, typename ELEMTYPE>
void sieve2357(SIMDTYPE *, size_t, const sieve2357_prime_t *);
#endif
