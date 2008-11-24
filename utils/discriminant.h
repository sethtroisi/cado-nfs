#ifndef CADO_UTILS_DISCRIMINANT_H
#define CADO_UTILS_DISCRIMINANT_H

#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Computes discriminant of a polynomial */
void discriminant (mpz_t, mpz_t *, const int);

#ifdef __cplusplus
}
#endif

#endif /* CADO_UTILS_DISCRIMINANT_H */
