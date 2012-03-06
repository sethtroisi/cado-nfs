#ifndef __FB_H__
#define __FB_H__

#include "types.h"



// Element of the factor base as an ideal (p, r).
// Given the two vectors (a0, b0) and (a1, b1) of the reduced q-lattice,
// lambda is precomputed as lambda = - (a1 - r*b1) / (a0 - r*b0) mod p.
// In the case of a projective root (i.e., a0 - r*b0 = 0 mod p),
// proj is set to 1.
// TODO: with alignements, we loose a lot of space
// See the types in cado-nfs and try to imitate ?
typedef struct {
  fbprime_t p;
  fbprime_t r;
  fbprime_t lambda;
  uint8_t   degp;
  _Bool     proj;
} fbideal_t;

typedef struct {
  unsigned   n;  // nb of entries in the factor base
  fbideal_t *elts;
} __factor_base_struct;

typedef       __factor_base_struct  factor_base_t[1];
typedef       __factor_base_struct *factor_base_ptr;
typedef const __factor_base_struct *factor_base_srcptr;

// Initialize a factor base, reading the ideals from a file.
// Return 1 if successful.
int factor_base_init(factor_base_ptr FB, const char *filename);

// Precompute lambda for each element of the factor base.
void factor_base_precomp_lambda(factor_base_ptr FB, qlat_srcptr qlat);

// Clean up memory.
void factor_base_clear(factor_base_ptr FB);

#endif   /* __FB_H__ */
