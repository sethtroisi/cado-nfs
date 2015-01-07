#ifndef MAKEFB_H
#define MAKEFB_H

#include <stdint.h>
#include "factor_base.h"

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
