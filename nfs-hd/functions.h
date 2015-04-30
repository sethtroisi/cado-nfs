#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <stdint.h>
#include "cado.h"
#include "utils.h"

/*
  f0: the first polynomial.
  f1: the second polynomial.
  g: a polynomial to build f0 and f1, must have the degree n (extension of the
   field).
  a: the a to build f1. (MNFS)
  b: the b to build f1.
  p: the prime number that define the field.
  h: a poylnomial to build f0 and f1.
  coeff0: lower bound of the coefficients of g.
  coeff1: upper bound of the coefficients of g.
  q: a special-q.
  t: dimension of the lattice.
*/
void function_special_q(mpz_poly_ptr f0, mpz_poly_ptr f1, mpz_poly_ptr g,
    mpz_ptr a, mpz_ptr b, mpz_ptr c, mpz_srcptr p, mpz_poly_srcptr h,
    int coeff0, int coeff1, unsigned int q, unsigned int t,
    unsigned int nb_times, double weight_0, double weight_1);

void function_classical(mpz_poly_ptr f0, mpz_poly_ptr f1, mpz_srcptr p,
    unsigned int n, int coeff0, int coeff1, unsigned int nb_times,
    double weight_0, double weight_1);

#endif // FUNCTIONS_H
