#ifndef UTILS_MPZ_POLY_H
#define UTILS_MPZ_POLY_H 

#include "cado.h"
#include "utils.h"

/*
 * Factorise a polynomial in characteristic 2. If l is empty, f % 2 is a
 *  polynomial of degree 0 or 0.
 * Not optimal for large polynomials.
 *
 * list: list of factors of f.
 * f: the polynomial we want to factorise.
 */
void mpz_poly_factor2(mpz_poly_factor_list_ptr list, mpz_poly_srcptr f);

#endif /* UTILS_MPZ_POLY_H */
