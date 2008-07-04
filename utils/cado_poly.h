#ifndef CADO_POLY_H_
#define CADO_POLY_H_

#include <stdio.h>
#include <gmp.h>

/* cado_poly is defined in cado.h */
#include <cado.h>

#ifdef __cplusplus
extern "C" {
#endif

extern void fprint_polynomial (FILE *, mpz_t *, const int);

// This reads a file created by polyselect and fill in the structure
// accordingly. Return 1 if success, 0 if failure (and diagnostic on
// stderr)
extern int cado_poly_read (cado_poly, char *filename);

extern void cado_poly_init (cado_poly);
extern void cado_poly_clear (cado_poly);

/* sanity check */
extern void check_polynomials (cado_poly);

/* legacy -- to be removed */
static inline int read_polynomial(cado_poly c, char *f)
    __attribute__((deprecated));
static inline int read_polynomial(cado_poly c, char *f)
{ return cado_poly_read(c,f); }


#ifdef __cplusplus
}
#endif

#endif	/* CADO_POLY_H_ */
