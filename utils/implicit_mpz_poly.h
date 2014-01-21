#ifndef IMPLICIT_MPZ_POLY_H_
#define IMPLICIT_MPZ_POLY_H_

#include <stdint.h>
#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

void mp_poly_content (mpz_t c, mpz_t *f, const int d);

#ifdef __cplusplus
}
#endif

#endif	/* IMPLICIT_MPZ_POLY_H_ */
