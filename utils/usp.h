#ifndef USP_H_
#define USP_H_

#include <gmp.h> /* for mpz_t */

/* this structure represents the interval [a/2^ka, b/2^kb] */
typedef struct {
  mpz_t a;
  int ka;
  mpz_t b;
  int kb;
} root_struct;

int numberOfRealRoots (mpz_t *p, int n, double T, int verbose, root_struct *R);
double rootRefine (root_struct *r, mpz_t *p, int n);

#endif  /* USP_H_ */
