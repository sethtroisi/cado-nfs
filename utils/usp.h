#ifndef USP_H_
#define USP_H_

#include <gmp.h> /* for mpz_t */
#include "mpz_poly.h"

/* this structure represents the interval [a/2^ka, b/2^kb] */
typedef struct {
  mpz_t a;
  int ka;
  mpz_t b;
  int kb;
} root_struct;

#ifdef __cplusplus
extern "C" {
#endif

int numberOfRealRoots (mpz_t *p, int n, double T, int verbose, root_struct *R);
double rootRefine (root_struct *r, mpz_t *p, int n, double precision);
void root_struct_init (root_struct *R);
void root_struct_clear (root_struct *R);
int mpz_poly_mpz_roots (mpz_t *r, mpz_poly_t p);

#ifdef __cplusplus
}
#endif

#endif  /* USP_H_ */
