#ifndef CADO_POLY_H_
#define CADO_POLY_H_

#include <stdio.h>
#include <gmp.h>

#include "params.h"

/* The maximum degree of polynomials supported. Used for statically 
   allocating storage (i.e. "mpz_t poly[MAXDEGREE]") */
#define MAXDEGREE 10

#define RATIONAL_SIDE   0
#define ALGEBRAIC_SIDE   1

struct cado_poly_side_s {
  mpz_t *f;          /* rational coefficients */
  int degree;        /* degree of polynomial g */
  // unsigned long lim; /* rational factor base bound */
  // int lpb;           /* rational large prime bound is 2^lpbr */
  // int mfb;           /* bound for rational residuals is 2^mfbr */
  // double lambda;     /* lambda sieve parameter on the rational side */
};
typedef struct cado_poly_side_s cado_poly_side[1];
typedef struct cado_poly_side_s * cado_poly_side_ptr;

struct cado_poly_s {
  mpz_t n;        /* number to factor */
  mpz_t m;        /* common root of f and g mod n */
  double skew;    /* skewness */

  cado_poly_side_ptr rat, alg;
  cado_poly_side pols[2];
} cado_poly_struct;
typedef struct cado_poly_s cado_poly[1];
typedef struct cado_poly_s * cado_poly_ptr;

extern const char * sidenames[2];

#ifdef __cplusplus
extern "C" {
#endif

extern void fprint_polynomial (FILE *, mpz_t *, const int);

// This reads a file created by polyselect and fill in the structure
// accordingly. Return 1 if success, 0 if failure (and diagnostic on
// stderr)
extern int cado_poly_read (cado_poly_ptr, const char *filename);
extern int cado_poly_read_stream (cado_poly_ptr, FILE *);
extern void cado_poly_set (cado_poly_ptr p, cado_poly_ptr q);

extern void cado_poly_init (cado_poly_ptr);
extern void cado_poly_clear (cado_poly_ptr);

/* sanity check */
extern void cado_poly_check (cado_poly_ptr);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_POLY_H_ */
