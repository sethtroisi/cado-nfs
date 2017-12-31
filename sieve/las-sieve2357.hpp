#ifndef LAS_SIEVE2357_HPP
#define LAS_SIEVE2357_HPP

#include "las-smallsieve-types.hpp"

namespace sieve2357 {

/* We hit sievearray[x] for x = idx + i*q with 0 <= x < arraylen */
typedef struct {
  fbprime_t q, idx;
  unsigned char logp;
} prime_t;

/* A predicate that tells whether a prime power q = p^k can be sieved by
   sieve2357 with a given SIMD and element data type */
template<typename SIMDTYPE, typename ELEMTYPE>
static inline bool can_sieve(const fbprime_t q)
{
  const size_t N = sizeof(SIMDTYPE) / sizeof(ELEMTYPE);
  return q == 1 || (q % 2 == 0 && N % q  == 0) || q == 3 || q == 5
      || q == 7;
}

template <typename SIMDTYPE, typename ELEMTYPE>
void sieve(SIMDTYPE *sievearray, size_t arraylen, const prime_t *primes,
    bool only_odd);

}
#endif
