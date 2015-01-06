#ifndef MAKEFB_H
#define MAKEFB_H

#include <stdint.h>
#include "factor_base.h"

/*
  Set an ideal_1 at an index. Do not forget to define LINESIEVE if you
   want to set an ideal mod r.

  fb: the factor base.
  index: index in the factor base.
  r: the r of the ideal (r, h).
  h: the h of the ideal (r, h).
  fbb: factor base bound for this side.
*/
void add_ideal_1(factor_base_ptr fb, uint64_t * index, uint64_t r,
                 mpz_poly_srcptr h, uint64_t fbb, unsigned int t);

/*
  Set an ideal_1 at an index. Do not forget to define LINESIEVE if you
   want to set an ideal mod r.

  fb: the factor base.
  index: index in the factor base.
  r: the r of the ideal (r, h).
  h: the h of the ideal (r, h).
  fbb: factor base bound for this side.
  lpb: large prime bound.
*/
void add_ideal_u(factor_base_ptr fb, uint64_t * index, uint64_t r,
                 mpz_poly_srcptr h, uint64_t fbb, mpz_t lpb, unsigned int t);

/*
  Make the factor bases for the V number fiels. Do not forget to define
   LINESIEVE if you want to set an ideal mod r.

  fb: the V factor bases.
  f: the V polynomials.
  fbb: the V factor base bound.
  t: dimension of the lattice.
  lpb: the V large prime bound.
  V: number of number fields.
*/
void makefb(factor_base_t * fb, mpz_poly_t * f, uint64_t * fbb, unsigned int t,
            mpz_t * lpb, unsigned int V);

#endif /* MAKEFB_H */
