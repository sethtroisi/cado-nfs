#ifndef UTILS_COFACTORISATION_H
#define UTILS_COFACTORISATION_H

#include "cado.h"
#include "utils.h"
#include "uint64_array.h"
#include "mat_Z.h"
#include "sieving_bound.h"
#include "ideal.h"

//Define an array with all the factors.
typedef struct
{
  mpz_t * factorization;
  unsigned int number;
  unsigned int alloc;
} s_factor_t;

typedef s_factor_t factor_t[1];
typedef s_factor_t * factor_ptr;
typedef const s_factor_t * factor_srcptr;

/*
 * To find the relations.
 *
 * indices: for each number field, positions where the norm is less than the
 *  threshold.
 * number_element: number of elements in the sieving region.
 * lpb: large prime bounds.
 * matrix: MqLLL.
 * f: polynomials that define the number fields.
 * H: sieving bounds.
 * V: number of number fields.
 */
unsigned int find_relations(uint64_array_t * indices, uint64_t number_element,
    unsigned int * lpb, mat_Z_srcptr matrix, mpz_poly_t * f,
    sieving_bound_srcptr H, unsigned int V, ideal_spq_srcptr special_q,
    unsigned int q_side, int main, FILE * outstd);

#endif /* UTILS_COFACTORISATION_H */
