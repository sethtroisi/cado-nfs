#ifndef __QLAT_H__
#define __QLAT_H__

#include "types.h"

int skewGauss(qlat_t qlat, unsigned int skewness);
void print_qlat_info(qlat_t qlat);
int is_valid_sq(qlat_t qlat, ffspol_srcptr F);

// Compute lambda for an element of the factor base.
// If the result is projective, the set lambda to p.
void compute_lambda(fbprime_ptr lambda,
        fbprime_srcptr p, fbprime_srcptr r, qlat_srcptr qlat);


void ab2ij(ij_t i, ij_t j, fppol_t a, fppol_t b, qlat_t qlat);
void ij2ab(fppol_t a, fppol_t b, ij_t i, ij_t j, qlat_t qlat);

#endif   /* __QLAT_H__ */
