#ifndef MPZ_POLY_H_
#define MPZ_POLY_H_

#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

int lift_root(mpz_t * f, int d, unsigned long pk, unsigned long * r);
int lift_rootz(mpz_t * f, int d, mpz_t pk, mpz_t r);
int mp_poly_cmp (mpz_t *f, mpz_t *g, int d);

#ifdef __cplusplus
}
#endif

#endif	/* MPZ_POLY_H_ */
