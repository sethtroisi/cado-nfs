#ifndef CADO_UTILS_LLL_H_
#define CADO_UTILS_LLL_H_

#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  mpz_t **coeff;
  int NumRows, NumCols;
} mat_Z;

long LLL (mpz_t det, mat_Z B, mat_Z* U, mpz_t a, mpz_t b, long verbose);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_UTILS_LLL_H_ */

