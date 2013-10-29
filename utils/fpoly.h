#ifndef FPOLY_H_
#define FPOLY_H_

#include <stdio.h>

/* floating point polynomials */

#ifdef __cplusplus
extern "C" {
#endif

/* fpoly.c */
double fpoly_eval (const double *, const unsigned int, const double);
double fpoly_dichotomy (const double *, const unsigned int, double, double, double, unsigned int);
void   fpoly_print (FILE *, const double *f, const unsigned int deg, char *name);


#ifdef __cplusplus
}
#endif

#endif	/* FPOLY_H_ */
