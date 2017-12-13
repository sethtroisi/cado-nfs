#ifndef LAS_SIEVE2357_HPP
#define LAS_SIEVE2357_HPP

#include "las-smallsieve-types.hpp"

typedef struct {
  fbprime_t p, q, idx;
  unsigned char logp;
} sieve2357_prime_t;


template<typename SIMDTYPE, typename ELEMTYPE>
SIMDTYPE bcaststride(const ELEMTYPE v, unsigned int offset,
    unsigned int stride);

template <typename SIMDTYPE, typename ELEMTYPE>
void sieve2357(SIMDTYPE *, size_t, const sieve2357_prime_t *);
#endif
