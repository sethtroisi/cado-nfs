#ifndef IMPLICIT_MPZ_POLY_H_
#define IMPLICIT_MPZ_POLY_H_

#include <stdint.h>
#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

int lift_root(mpz_t * f, int d, unsigned long pk, unsigned long * r);
int lift_rootz(mpz_t * f, int d, mpz_t pk, mpz_t r);
void mp_poly_homogeneous_eval_siui (mpz_t v, mpz_t *f, const unsigned int d, const int64_t i, const uint64_t j);
void mp_poly_content (mpz_t c, mpz_t *f, const int d);

#ifdef __cplusplus
}
#endif

#endif	/* IMPLICIT_MPZ_POLY_H_ */
