#ifndef DOUBLE_POLY_H_
#define DOUBLE_POLY_H_

#include <stdio.h>
#include <limits.h>
#include "mpz_poly.h"

/* floating point polynomials */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  unsigned int deg;
  double *coeff;         /* array of deg+1 entries */
} double_poly_struct_t;

typedef double_poly_struct_t double_poly_t[1];
typedef double_poly_struct_t * double_poly_ptr;
typedef const double_poly_struct_t * double_poly_srcptr;

/* double_poly.c */
void double_poly_init (double_poly_ptr, unsigned int);
void double_poly_clear (double_poly_ptr);
void double_poly_set (double_poly_ptr, double_poly_srcptr);
double double_poly_eval (double_poly_srcptr, const double);
double double_poly_dichotomy (double_poly_srcptr, double, double, double,
                              unsigned int);
void double_poly_derivative(double_poly_ptr, double_poly_srcptr);
void double_poly_product(double_poly_ptr, double_poly_srcptr, double_poly_srcptr);
void double_poly_sum(double_poly_ptr, double_poly_srcptr, double_poly_srcptr);
void double_poly_subtract(double_poly_ptr, double_poly_srcptr, double_poly_srcptr);
void double_poly_scale(double_poly_ptr, double_poly_srcptr, double);
void double_poly_revert(double_poly_ptr);
unsigned int double_poly_compute_roots(double *, double_poly_ptr, double);
void double_poly_print (FILE *, double_poly_srcptr, char *name);
void double_poly_set_mpz_poly (double_poly_ptr p, mpz_poly_ptr q);

#ifdef __cplusplus
}
#endif

#endif	/* DOUBLE_POLY_H_ */
