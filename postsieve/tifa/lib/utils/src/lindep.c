//
// Copyright (C) 2006, 2007 INRIA (French National Institute for Research in
// Computer Science and Control)
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
//

/**
 * \file    lindep.c
 * \author  Jerome Milan
 * \date    Mon Mar 13 2006
 * \version 1.0
 */

 /*
  *  Copyright (C) 2006, 2007 INRIA
  *  License: GNU Lesser General Public License (LGPL)
  *  History:
  */

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include "lindep.h"
#include "gauss_elim.h"
#include "first_primes.h"
#include "print_error.h"

//-----------------------------------------------------------------------------
void fill_matrix_trial_div(binary_matrix_t* const matrix,
                           mpz_array_t* const partially_factored,
                           const mpz_array_t* const to_factor,
                           const uint32_array_t* const factor_base) {

    bool flip_bit = false;

    uint32_t nrows        = matrix->nrows;
    uint32_t nrows_to_add = to_factor->length;

    if ((nrows_to_add + nrows) > matrix->nrows_alloced) {
        nrows_to_add = matrix->nrows_alloced - nrows;
    }

    for (uint32_t i = 0; i < nrows_to_add; i++) {

        mpz_init_set(partially_factored->data[i], to_factor->data[i]);

        if (mpz_sgn(partially_factored->data[i]) == -1) {
            flip_matrix_bit(nrows, 0, matrix);
            mpz_neg(partially_factored->data[i], partially_factored->data[i]);
        }
        for (uint32_t j = 0; j < factor_base->length; j++) {

            while (0 != mpz_divisible_ui_p(partially_factored->data[i],
                                           factor_base->data[j]) ) {

                mpz_divexact_ui(
                    partially_factored->data[i],
                    partially_factored->data[i],
                    factor_base->data[j]
                );
                flip_bit = !flip_bit;
            }
            if (flip_bit) {
                flip_matrix_bit(nrows, j + 1, matrix);
                flip_bit = false;
            }
        }
        nrows++;
    }
    partially_factored->length = to_factor->length;
    matrix->nrows = nrows;
}
//-----------------------------------------------------------------------------
void fill_matrix_from_list(binary_matrix_t* const matrix,
                           const mpz_array_t* const smooth_array,
                           const uint32_array_list_t* const list,
                           const uint32_array_t* const factor_base) {
    mpz_t tmp;
    mpz_init(tmp);
    uint32_t ind_p = 0;

    bool flip_bit = false;

    for (uint32_t i = 0; i < smooth_array->length; i++) {

        mpz_set(tmp, smooth_array->data[i]);
        if (list->data[i] != NULL) {
            //
            // Simple trial divisions by the primes given by the
            // uint32_array_t list->data[i].
            //
            for (uint32_t j = 0; j < list->data[i]->length; j++) {

                if (0 != mpz_divisible_ui_p(tmp, list->data[i]->data[j])) {
                    //
                    // We need the index of this prime in the factor base to
                    // properly flip the matrix's bits...
                    //
                    ind_p = index_in_sorted_uint32_array(
                                    list->data[i]->data[j],
                                    factor_base,
                                    0,
                                    factor_base->length);

                    mpz_divexact_ui(tmp, tmp, factor_base->data[ind_p]);

                    flip_bit = !flip_bit;
                }
                while (0 != mpz_divisible_ui_p(tmp, list->data[i]->data[j])) {
                    mpz_divexact_ui(tmp, tmp, factor_base->data[ind_p]);
                    flip_bit = !flip_bit;
                }
                if (flip_bit) {
                    flip_matrix_bit(i, ind_p+1, matrix);
                    flip_bit = false;
                }
            }
        }
    }
    mpz_clear(tmp);
}
//-----------------------------------------------------------------------------
uint32_array_list_t*
find_dependencies(binary_matrix_t* const matrix, ls_res_method_t method) {

    uint32_array_list_t* relations;
    relations = alloc_uint32_array_list(matrix->nrows - matrix->ncols);

    switch (method) {
        case SMART_GAUSS_ELIM:
            gaussian_elim(relations, matrix);
            break;

        default:
            gaussian_elim(relations, matrix);
    }
    return relations;
}
//-----------------------------------------------------------------------------
ecode_t find_factors(mpz_array_t*  const factors,
                     const mpz_t n,
                     const mpz_array_t* const xi_array,
                     const mpz_array_t* const yi_array,
                     const uint32_array_list_t* const dependencies) {

    ecode_t ecode = NO_FACTOR_FOUND;

    mpz_t x;
    mpz_init(x);

    mpz_t y;
    mpz_init(y);

    mpz_t gcd;
    mpz_init(gcd);

    uint32_t irow = 0;

    for (uint32_t idep = 0; idep < dependencies->length; idep++) {
        if (factors->length == factors->alloced) {
            break;
        }
        mpz_set_ui(x, 1);
        mpz_set_ui(y, 1);
        //
        // Compute y   = Sqrt(Prod(yi))
        // Compute x   = Prod(xi) + y = Prod(xi) + Sqrt(Prod(yi))
        // Compute gcd = gcd(x, n)
        //
        for (uint32_t irel = 0;
             irel < dependencies->data[idep]->length;
             irel++) {

            irow = dependencies->data[idep]->data[irel];
            mpz_mul(x, x, xi_array->data[irow]);
            mpz_mod(x, x, n);
            mpz_mul(y, y, yi_array->data[irow]);
        }
        if (0 == mpz_perfect_square_p(y)) {
           //
           // At this stage, if we made no mistake, y should be a product
           // of squares. If that's not the case, something went terribly
           // wrong and it's no use to proceed any further...
           //
           PRINTF_STDERR("\n");
           PRINT_ERROR("%Zd is not a perfect square!\n", y);
           PRINT_ERROR("This should not happen - Factorization aborted!\n");
           return FATAL_INTERNAL_ERROR;
        }
        mpz_sqrt(y, y);
        mpz_add(x, x, y);
        mpz_gcd(gcd, x, n);

        if ( (0 != mpz_cmp(gcd, n)) && (0 != mpz_cmp_ui(gcd, 1)) ) {
            //
            // gcd is a non trivial factor of n. Add it to our list if it's
            // not already there...
            //
            if (NOT_IN_ARRAY == index_in_mpz_array(gcd, factors)) {
                append_mpz_to_array(factors, gcd);
                ecode = SOME_FACTORS_FOUND;
            }
        }
    }
    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(gcd);
    
    return ecode;
}
//-----------------------------------------------------------------------------

