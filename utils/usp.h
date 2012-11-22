#ifndef USP_H_
#define USP_H_

#include <gmp.h> /* for mpz_t */

typedef struct {
  mpz_t a;
  int ka;
  mpz_t b;
  int kb;
} root_struct;

int numberOfRealRoots (mpz_t *p, int n, double T, int verbose, root_struct *R);

#endif  /* USP_H_ */
