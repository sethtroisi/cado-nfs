#ifndef PR_H_
#define PR_H_

#include <stdint.h>
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

uint32_t poly_roots_double(double *poly, uint32_t degree, double complex *roots);
uint32_t poly_roots_longdouble(double *poly, uint32_t degree, long double complex *roots);

#ifdef __cplusplus
}
#endif

#endif	/* PR_H_ */
