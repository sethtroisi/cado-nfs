#ifndef FPOLY_H_
#define FPOLY_H_

#include <stdio.h>

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
double double_poly_eval (double_poly_srcptr, const double);
double double_poly_dichotomy (double_poly_srcptr, double, double, double,
                              unsigned int);
void double_poly_print (FILE *, double_poly_srcptr, char *name);

static inline void
double_poly_scale (double *u, const double *t, unsigned int d, double h)
{
  double hpow;
  u[d] = t[d];
  for (hpow = h; --d != UINT_MAX; hpow *= h) u[d] = t[d] * hpow;
}

#ifdef __cplusplus
}
#endif

#endif	/* FPOLY_H_ */
