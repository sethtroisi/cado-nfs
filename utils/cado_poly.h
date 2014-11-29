#ifndef CADO_POLY_H_
#define CADO_POLY_H_

#include <stdio.h>
#include <gmp.h>

#include "params.h"
#include "mpz_poly.h"

/* The maximum degree of polynomials supported. Used for statically 
   allocating storage (i.e. "mpz_t poly[MAXDEGREE]") */
#define MAXDEGREE 10

#define RATIONAL_SIDE   0
#define ALGEBRAIC_SIDE   1

#define NB_POLYS_MAX 8 /* maximal number of polynomials in multiple fields */

struct cado_poly_s {
  mpz_t n;        /* number to factor */
  double skew;    /* skewness from poly file, if given, otherwise 0. */

  int nb_polys;   /* number of polynomials used, 2 in most cases */
  mpz_poly_ptr rat, alg;
  mpz_poly_t pols[NB_POLYS_MAX];
};
typedef struct cado_poly_s cado_poly[1];
typedef struct cado_poly_s * cado_poly_ptr;

/* we used to have a true global variable declared here (not extern).
 * This is horrible. */
extern struct cado_poly_s cado_poly_struct;

extern const char * sidenames[2]; // FIXME: 2 or NB_POLYS_MAX?

#ifdef __cplusplus
extern "C" {
#endif

// This reads a file created by polyselect and fill in the structure
// accordingly. Return 1 if success, 0 if failure (and diagnostic on
// stderr)
extern int cado_poly_read (cado_poly_ptr, const char *filename);
extern int cado_poly_read_stream (cado_poly_ptr, FILE *);
extern void cado_poly_set (cado_poly_ptr p, cado_poly_ptr q);

extern void cado_poly_init (cado_poly_ptr);
extern void cado_poly_clear (cado_poly_ptr);

// Compute m as the common root of f and g mod N.
// N is taken as a third argument; it can be a strict factor of the N
// stored in the polynomial.
// If this fails, then in most of the cases we have found a factor of N
// that is given instead of m.
// The return value tells whether it worked (and then m is the common
// root) or it failed (and then m is a factor of N).
extern int cado_poly_getm(mpz_ptr, cado_poly_ptr, mpz_ptr);

/* Return the rational side or -1 if two algebraic side */
extern int cado_poly_get_ratside (cado_poly_ptr);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_POLY_H_ */
