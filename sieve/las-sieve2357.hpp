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
static inline bool sieve2357_can_sieve(const fbprime_t p, const fbprime_t q)
{
  const size_t N = sizeof(SIMDTYPE) / sizeof(ELEMTYPE);
  return (p == 2 && q <= N)
         || (p == 3 && q <= 3)
         || (p == 5 && q <= 5)
         || (p == 7 && q <= 7);
}

template <typename SIMDTYPE, typename ELEMTYPE>
void sieve2357(SIMDTYPE *, size_t, const sieve2357_prime_t *);
#endif
