#ifndef INT64_POLY_H
#define INT64_POLY_H

/*
  This file try to copy the mpz_poly way.
  For instance, it can produce many strange thing and I need to debug it. We do
  not want to deal with mpz_poly, but sometimes, we can not.
*/

#include <stdint.h>
#include "sieving_interval.h"
#include "mpz_poly.h"

typedef struct {
  int alloc;
  int deg;
  int64_t * coeff;
} s_int64_poly_t;

typedef s_int64_poly_t int64_poly_t[1];
typedef s_int64_poly_t * int64_poly_ptr;
typedef const s_int64_poly_t * int64_poly_srcptr;

/*
  Initialize a polynomial.
  IN:
    f: int64_poly_ptr, the polynomial.
    d: int, degree of the polyniomial.
*/
void int64_poly_init(int64_poly_ptr f, int d);

/*
  Realloc f to (at least) nc coefficients.
  IN:
    f: int64_poly_ptr, the polynomial.
    nc: int, number of coefficients.
 */
void int64_poly_realloc (int64_poly_ptr f, int nc);

/*
  Copy f to g.
  IN:
    g: mpz_poly_ptr, the receiver polynomial.
    f: mpz_poly_srcptr, the source polynomial
*/
void int64_poly_copy(int64_poly_ptr g, int64_poly_srcptr f);

/*
  Swap f and g. g must just be initialized.
  IN:
     g: mpz_poly_ptr, a polynomial.
     f: mpz_poly_ptr, a polynomial
*/
void int64_poly_swap (int64_poly_ptr f, int64_poly_ptr g);

/*
  Delete the polynomial.
  IN:
    f: int64_poly_t *, the polynomial.
*/
void int64_poly_clear(int64_poly_ptr f);

/*
  Find polynomial degree, deg is an upper bound for the degree.
  IN:
    f: int64_poly_ptr, the polynomial.
    deg: int, upper bound for the degree of f.

*/
void int64_poly_cleandeg(int64_poly_ptr f, int deg);

/*
  Set f to a zero polynomial.
  IN:
    f: int64_poly_ptr, the polynomial.
*/
void int64_poly_set_zero(int64_poly_ptr f);

/*
  Set the ith coefficient of a polynomial.
  IN:
    f: int64_poly_ptr, the polynomial.
    i: int, index of the coefficient.
    coeff: int64_t, new value of the ith coefficient.

*/
void int64_poly_setcoeff(int64_poly_ptr f, int i, int64_t coeff);

/*
  Get the ith coefficient of a polynomial.
  IN:
    i: int, index of the coefficient.
    f: int64_poly_t, the polynomial.
*/
void int64_poly_get_coeff(int64_t * res, int i, int64_poly_srcptr f);

/*
  Set f to x^i.
  IN:
    f: int64_poly_ptr, the polynomial.
    i: int, degree of the polynomial.
 */
void int64_poly_set_xi(int64_poly_ptr f, int i);

/*
  Set f to b*x^i.
  IN:
    f: int64_poly_ptr, the polynomial.
    i: int, degree of the polynomial.
    b: int64_t, multiplicative coefficient.
*/
void int64_poly_set_bxi(int64_poly_ptr f, int i, int64_t b);

/*
  To print a polynomial.
  IN:
    f: int64_poly_srcptr, the polynomial.
*/
void int64_poly_printf(int64_poly_srcptr f);

/*
  To write a polynomial in a file.

  file: the file.
  f: the polynomial.
*/
void int64_poly_fprintf(FILE * file, int64_poly_srcptr f);

/*
  TODO: explain what the function do.
*/
void int64_poly_add_one_sieving_region(int64_poly_ptr f,
                                       sieving_interval_srcptr H);

/*
  Find the maximum of the coefficients of the polynomial.
  IN:
    max: int64_t *, the maximum.
    f: int64_poly_srcptr, the polynomial.
*/
void int64_poly_max(int64_t * max, int64_poly_srcptr f);

/*
  Find the infinite norm of f and set in in.
  IN:
    in: int64_t * , the infinite norm.
    f: int64_poly_srcptr, the polynomial.
*/
void int64_poly_infinite_norm(uint64_t * in, int64_poly_srcptr f);

/* /\* */
/*   TODO: explain what do this function. */
/* *\/ */
/* void int64_poly_addoneptr(int64_poly_t * f, sieving_interval_t H); */

/* /\* */
/*   Return 1 if the polynomial is a zero polynomial, 0 othewise. */
/*   IN: */
/*     f: int64_poly_t, the polynomial. */
/* *\/ */
/* int int64_poly_is_zero(int64_poly_t f); */

/* /\* */
/*   Delete the unsignifcant 0 coefficient and change the degree. */
/*   EX: */
/*     1 + 1*x^1 + 0*x^2 + 2*x^3 + 0*x^4 -> 1 + 1*x^1 + 0*x^2 + 2*x^3. */
/*   IN: */
/*     f: int64_poly_t, the polynomial. */
/* *\/ */
/* void int64_poly_delete_unsignificant(int64_poly_t * f); */

/* /\* */
/*   Compute the resultant between f and g. */
/*   WARNING: it is a int64_t in output, use the mpz version for safety reason. */
/*   IN: */
/*     f: int64_poly_t, a polynomial. */
/*     g: int64_poly_t, a polynomial. */
/*   OUT: */
/*     the value of the resultant. */
/* *\/ */
/* int64_t int64_poly_resultant(int64_poly_t f, int64_poly_t g); */

/* /\* */
/*   Return the leading coefficient of a polynomial. */
/*   IN: */
/*     f: int64_poly_t, the polynomial. */
/*   OUT: */
/*     leading coefficient of f. */
/* *\/ */
/* int64_t int64_poly_leading_coefficient(int64_poly_t f); */

/* /\* */
  
/*  *\/ */
/* void int64_poly_multiply_constant(int64_t a, int64_poly_t * f); */

/* void int64_poly_divide_constant(int64_t a, int64_poly_t * f); */

/* void int64_poly_swap(int64_poly_t * f, int64_poly_t * g); */

/* void int64_poly_add_poly(int64_poly_t a, int64_poly_t b, int64_poly_t * c); */

/* void int64_poly_mult_xi(int64_poly_t *  a, int i); */

/* void int64_poly_mult_bxi(int64_poly_t * a, int64_t b, int i); */

/* void int64_poly_mult_poly(int64_poly_t a, int64_poly_t b, int64_poly_t * */
/* 				c); */

/* void int64_poly_division_with_remainder(int64_poly_t a, int64_poly_t b, */
/* 					int64_poly_t * quo, int64_poly_t * rem); */

/* int64_t int64_poly_content(int64_poly_t f); */

/* void int64_poly_primitive_part(int64_poly_t * f); */

/* void int64_poly_pseudo_division(int64_poly_t a, int64_poly_t b, int64_poly_t * */
/* 				quo, int64_poly_t * rem); */

/* void int64_poly_pseudo_remainder(int64_poly_t a, int64_poly_t b, int64_poly_t * */
/* 				rem); */

/* //TODO: link with mpz_poly.c */
/* /\* void int64_poly_to_mpz_poly(int64_poly_t a, mpz_poly_t * b); *\/ */

void int64_poly_to_mpz_poly(mpz_poly_ptr a, int64_poly_srcptr b);


#endif /* INT64_POLY_H */
