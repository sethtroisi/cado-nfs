#ifndef __QLAT_H__
#define __QLAT_H__

#include "types.h"

int skewGauss(qlat_t qlat, unsigned int skewness);
void print_qlat_info(qlat_t qlat);



void ab2ij(ij_t i, ij_t j, fppol_t a, fppol_t b, qlat_t qlat);
void ij2ab(fppol_t a, fppol_t b, ij_t i, ij_t j, qlat_t qlat);

#endif   /* __QLAT_H__ */
