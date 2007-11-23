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
 * \file    res_tdiv.c
 * \author  Jerome Milan
 * \date    Thu Nov 15 2007
 * \version 1.0
 */

 /*
  *  History:
  *    1.0: Thu Nov 15 2007 by JM:
  *         - Initial version.
  */

#include <stdlib.h>
#include <stdio.h>

#include "tifa_config.h"
#include <gmp.h>
#if TIFA_USE_GMP_INTERNAL_FUNCS
    #include "gmp-impl.h"
#endif
#include "res_tdiv.h"
#include "macros.h"
#include "first_primes.h"
#include "hashtable.h"
#include "gmp_utils.h"

//-----------------------------------------------------------------------------
//                 PROTOTYPES OF NON PUBLIC FUNCTION(S)
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
uint32_t res_tdiv_no_ea(smooth_filter_t* const filter);
//------------------------------------------------------------------------------
uint32_t res_tdiv_step(smooth_filter_t* const filter, unsigned long int step);
//------------------------------------------------------------------------------
uint32_t res_tdiv_first(smooth_filter_t* const filter);
//------------------------------------------------------------------------------
uint32_t res_tdiv_last(smooth_filter_t* const filter);
//------------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                 "PUBLIC" FUNCTION(S) IMPLEMENTATION
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
inline uint32_t res_tdiv(smooth_filter_t* const filter,unsigned long int step) {

    if (filter->method == TDIV) {
        return res_tdiv_no_ea(filter);
    }
    if (step == 0) {
        return res_tdiv_first(filter);
    }
    if (step == filter->nsteps) {
        return res_tdiv_last(filter);
    }
    return res_tdiv_step(filter, step - 1);
}
//------------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                  "PRIVATE" FUNCTION(S) IMPLEMENTATION
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
uint32_t res_tdiv_no_ea(smooth_filter_t* const filter) {

    mpz_array_t* const in_xi  = filter->candidate_xi;
    mpz_array_t* const in_yi  = filter->candidate_yi;
    mpz_array_t* const acc_xi = filter->accepted_xi;
    mpz_array_t* const acc_yi = filter->accepted_yi;

    if (acc_yi->length == acc_yi->alloced) {
        return 0;
    }

    hashtable_t* const htable = filter->htable;
    mpz_ptr n = filter->n;

    uint32_t  nprimes = filter->factor_base[0]->length;
    uint32_t* base    = filter->factor_base[0]->data;

    uint32_t used = 0U;

    mpz_t y;
    mpz_init(y);

    mpz_t nsp_inv;
    mpz_init(nsp_inv);

    uint32_t k = in_yi->length;

#if TIFA_USE_GMP_INTERNAL_FUNCS
    mp_limb_t* num  = 0;
    uint32_t   size = 0;
    uint32_t   pi   = 0;
#endif

    while (k != 0U) {

        k--;
        used++;
        mpz_set(y, in_yi->data[k]);

#if TIFA_USE_GMP_INTERNAL_FUNCS
        //
        // Use some internal GMP functions.
        //
        // _NOTE_: We actually use the mpn_divexact_1 function which is not
        //         documented. See:
        //
        //  http://gmplib.org/list-archives/gmp-bugs/2006-December/000660.html
        //
        //         Apparently this "secret GMP function" (as Torbjorn Granlund
        //         put it) is not supposed to be used by client code... Could it
        //         crash your computer? Burn down your office? Give bad results?
        //         Or, more prosaically, will it be removed in future versions?
        //         In any case, maintainers should be aware of this...
        //
        num  = PTR(y);
        size = ABSIZ(y);
        pi   = 0;
        //
        // Divisions by 2 (always in the factor base)
        //
        if ( (num[0] & 1) == 0) {
            mpn_rshift(num, num, size, 1);
            MPN_NORMALIZE(num, size);
            SIZ(y) = size;

            while ( (num[0] & 1) == 0) {
                mpn_rshift(num, num, size, 1);
                MPN_NORMALIZE(num, size);
                SIZ(y) = size;
            }
        }
        //
        // Divisions by odd primes:
        //
        // _WARNING_: The following loop uses the "internal" mpn_modexact_1_odd
        //            and mpn_divexact_1 functions which could potentially be
        //            removed or changed in future GMP versions. See previous
        //            warning note.
        //
        for (pi = 1; pi < nprimes; pi++) {

            if (mpn_modexact_1_odd(num, size, base[pi]) == 0) {

                mpn_divexact_1(num, num, size, base[pi]);
                MPN_NORMALIZE(num, size);
                SIZ(y) = size;

                while (mpn_modexact_1_odd(num, size, base[pi]) == 0) {
                    mpn_divexact_1(num, num, size, base[pi]);
                    MPN_NORMALIZE(num, size);
                    SIZ(y) = size;
                }
            }
        }
#else
        //
        // GMP's internal functions are not available: revert to a simple
        // MPZ based implementation. (Using "public" MPN function doesn't
        // bring any significant benefits according to my tests, at least on
        // Opteron)
        //
        for (uint32_t i = 0; i < nprimes; i++) {
            if (mpz_divisible_ui_p(y, base[i])) {
                mpz_divexact_ui(y, y, base[i]);
                while (mpz_divisible_ui_p(y, base[i])) {
                    mpz_divexact_ui(y, y, base[i]);
                }
            }
        }
#endif

        if (0 == mpz_cmp_ui(y, 1)) {

            mpz_init_set(acc_yi->data[acc_yi->length], in_yi->data[k]);
            acc_yi->length++;

            mpz_init_set(acc_xi->data[acc_xi->length], in_xi->data[k]);
            acc_xi->length++;

            if (acc_xi->length == acc_xi->alloced) {
                break;
            }

        } else {
            //
    		// y is now the non-smooth part of in_yi->data[i].
    		// Is y a large prime?
    		//
            if (mpz_cmp_ui(y, LARGEST_PRIME) <= 0 ) {

                uint32_t nsp_ui = mpz_get_ui(y);
                //
                // Get the index of y in first_primes_array...
                // If ind != NOT_IN_ARRAY, then y is indeed a prime...
                //
                uint32_t ind = index_in_sorted_uint32_array(
                                    nsp_ui,
                                    &first_primes_array,
                                    0,
                                    first_primes_array.length);

                if (ind != NOT_IN_ARRAY) {
                    //
                    // cand_yi->data[k] is the product of a smooth number with
                    // the prime number first_primes_array->data[ind]
                    //
                    uint32_t *key = malloc(sizeof(uint32_t));
                    *key  = ind;

                    //
                    // _NOTE_: Actually it may be possible to leave the pair
                    //         in the hashtable...
                    //
                    // _TO_DO: Try to keep the entry in the hashtable...
                    //
                //mpz_pair_t *found = remove_entry_in_hashtable(htable, key);

                    mpz_pair_t *found = get_entry_in_hashtable(htable,key);

                    if (found != NULL) {
                        //
                        // The hashtable did contain another entry for this
                        // particular prime number: we can compute a new
                        // pair to add in our arrays xi and smooth_yi
                        //
                        mpz_init(acc_yi->data[acc_yi->length]);
                        mpz_init(acc_xi->data[acc_xi->length]);

                        mpz_mul(
                            acc_yi->data[acc_yi->length],
                            in_yi->data[k],
                            found->y
                        );

                        mpz_mul(
                            acc_xi->data[acc_xi->length],
                            in_xi->data[k],
                            found->x
                        );

                        mpz_divexact_ui(
                            acc_yi->data[acc_yi->length],
                            acc_yi->data[acc_yi->length],
                            nsp_ui
                        );

                        mpz_divexact_ui(
                            acc_yi->data[acc_yi->length],
                            acc_yi->data[acc_yi->length],
                            nsp_ui
                        );

                        if (0 != mpz_invert(nsp_inv, y, n)) {
                            mpz_mul(
                                acc_xi->data[acc_xi->length],
                                acc_xi->data[acc_xi->length],
                                nsp_inv
                            );
                        } else {
                            //
                            // This should not happen: since y is a prime,
                            // it has an inverse in Z/nZ... unless n is a
                            // multiple of y! In such a (rare) case, we have
                            // found a factor of n, but handling it is quite
                            // cumbersome, so we just ignore it...
                            //
                            free(key);
                            //clear_mpz_pair(found);
                            //free(found);
                            continue;
                        }
                        acc_yi->length++;
                        acc_xi->length++;

                        free(key);
                        //clear_mpz_pair(found);
                        //free(found);

                        if (acc_yi->length == acc_yi->alloced) {
                            break;
                        }
                    } else {
                        //
                        // Add this pair in the hashtable with the position
                        // of the prime in first_primes_array as key...
                        //
                        mpz_pair_t* pair = malloc(sizeof(mpz_pair_t));
                        mpz_init_set(pair->x, in_xi->data[k]);
                        mpz_init_set(pair->y, in_yi->data[k] );
                        add_entry_in_hashtable(htable, (void*)key, (void*)pair);
                        //
                        // _WARNING_: Do not free key and pair as they are now
                        //            referenced by the hashtable!
                        //
                    }
                }
            }
        }
    }

    mpz_clear(y);
    mpz_clear(nsp_inv);

    return used;
}
//------------------------------------------------------------------------------
uint32_t res_tdiv_step(smooth_filter_t* const filter, unsigned long int step) {

    mpz_array_t* const in_xi   = filter->filtered_xi[step];
    mpz_array_t* const in_yi   = filter->filtered_yi[step];
    mpz_array_t* const in_cof  = filter->cofactors  [step];
    mpz_array_t* const out_xi  = filter->filtered_xi[step + 1];
    mpz_array_t* const out_yi  = filter->filtered_yi[step + 1];
    mpz_array_t* const out_cof = filter->cofactors  [step + 1];
    mpz_array_t* const acc_xi  = filter->accepted_xi;
    mpz_array_t* const acc_yi  = filter->accepted_yi;

    if (acc_yi->length == acc_yi->alloced) {
        return 0;
    }

    uint32_t next_acc;
    uint32_t next_out;

    hashtable_t* htable = filter->htable;
    mpz_ptr n           = filter->n;
    mpz_ptr bound       = filter->bounds[step + 1];

    uint32_t  nprimes = filter->factor_base[step + 1]->length;
    uint32_t* base    = filter->factor_base[step + 1]->data;

    uint32_t used = 0U;

    mpz_t y;
    mpz_init(y);

    mpz_t nsp_inv;
    mpz_init(nsp_inv);

    uint32_t k = in_yi->length;
    
#if TIFA_USE_GMP_INTERNAL_FUNCS
    mp_limb_t* num  = 0;
    uint32_t   size = 0;
    uint32_t   pi   = 0;
#endif

    while (k != 0U) {

        k--;
        used++;
        mpz_set(y, in_yi->data[k]);

#if TIFA_USE_GMP_INTERNAL_FUNCS
        //
        // Use some internal GMP functions.
        //
        // _NOTE_: We actually use the mpn_divexact_1 function which is not
        //         documented. See:
        //
        //  http://gmplib.org/list-archives/gmp-bugs/2006-December/000660.html
        //
        //         Apparently this "secret GMP function" (as Torbjorn Granlund
        //         put it) is not supposed to be used by client code... Could it
        //         crash your computer? Burn down your office? Give bad results?
        //         Or, more prosaically, will it be removed in future versions?
        //         In any case, maintainers should be aware of this...
        //
        num  = PTR(y);
        size = ABSIZ(y);
        pi   = 0;
        //
        // Divisions by 2 (always in the factor base)
        //
        if ( (num[0] & 1) == 0) {
            mpn_rshift(num, num, size, 1);
            MPN_NORMALIZE(num, size);
            SIZ(y) = size;

            while ( (num[0] & 1) == 0) {
                mpn_rshift(num, num, size, 1);
                MPN_NORMALIZE(num, size);
                SIZ(y) = size;
            }
        }
        //
        // Divisions by odd primes:
        //
        // _WARNING_: The following loop uses the "internal" mpn_modexact_1_odd
        //            and mpn_divexact_1 functions which could potentially be
        //            removed or changed in future GMP versions. See previous
        //            warning note.
        //
        for (pi = 1; pi < nprimes; pi++) {

            if (mpn_modexact_1_odd(num, size, base[pi]) == 0) {

                mpn_divexact_1(num, num, size, base[pi]);
                MPN_NORMALIZE(num, size);
                SIZ(y) = size;

                while (mpn_modexact_1_odd(num, size, base[pi]) == 0) {
                    mpn_divexact_1(num, num, size, base[pi]);
                    MPN_NORMALIZE(num, size);
                    SIZ(y) = size;
                }
            }
        }
#else
        //
        // GMP's internal functions are not available: revert to a simple
        // MPZ based implementation. (Using "public" MPN function doesn't
        // bring any significant benefits according to my tests, at least on
        // Opteron)
        //
        for (uint32_t i = 0; i < nprimes; i++) {
            if (mpz_divisible_ui_p(y, base[i])) {
                mpz_divexact_ui(y, y, base[i]);
                while (mpz_divisible_ui_p(y, base[i])) {
                    mpz_divexact_ui(y, y, base[i]);
                }
            }
        }
#endif

        next_acc = acc_yi->length;
        next_out = out_yi->length;

        if (0 == mpz_cmp_ui(y, 1)) {

            mpz_init_set(acc_yi->data[next_acc], in_yi->data[k]);
            mpz_init_set(acc_xi->data[next_acc], in_xi->data[k]);

            mpz_mul(
                acc_yi->data[next_acc],
                acc_yi->data[next_acc],
                in_cof->data[k]
            );

            acc_yi->length++;
            acc_xi->length++;

            if (acc_xi->length == acc_xi->alloced) {
                break;
            }

        } else {
            //
    		// y is now the non-smooth part of in_yi->data[i].
    		// Is y a large prime?
    		//
            if (mpz_cmp_ui(y, LARGEST_PRIME) <= 0 ) {

                uint32_t nsp_ui = mpz_get_ui(y);
                //
                // Get the index of y in first_primes_array...
                // If ind != NOT_IN_ARRAY, then y is indeed a prime...
                //
                uint32_t ind = index_in_sorted_uint32_array(
                                    nsp_ui,
                                    &first_primes_array,
                                    0,
                                    first_primes_array.length);

                if (ind != NOT_IN_ARRAY) {
                    //
                    // cand_yi->data[k] is the product of a smooth number with
                    // the prime number first_primes_array->data[ind]
                    //
                    uint32_t *key = malloc(sizeof(uint32_t));
                    *key  = ind;

                    //
                    // _NOTE_: Actually it may be possible to leave the pair
                    //         in the hashtable...
                    //
                    // _TO_DO: Try to keep the entry in the hashtable...
                    //
                //mpz_pair_t *found = remove_entry_in_hashtable(htable, key);

                    mpz_pair_t *found = get_entry_in_hashtable(htable,key);

                    if (found != NULL) {
                        //
                        // The hashtable did contain another entry for this
                        // particular prime number: we can compute a new
                        // pair to add in our accepted_* arrays.
                        //
                        mpz_init(acc_yi->data[next_acc]);
                        mpz_init(acc_xi->data[next_acc]);

                        mpz_mul(
                            acc_yi->data[next_acc],
                            in_yi->data[k],
                            found->y
                        );

                        mpz_mul(
                            acc_yi->data[next_acc],
                            acc_yi->data[next_acc],
                            in_cof->data[k]
                        );

                        mpz_mul(
                            acc_xi->data[next_acc],
                            in_xi->data[k],
                            found->x
                        );

                        mpz_divexact_ui(
                            acc_yi->data[next_acc],
                            acc_yi->data[next_acc],
                            nsp_ui
                        );

                        mpz_divexact_ui(
                            acc_yi->data[next_acc],
                            acc_yi->data[next_acc],
                            nsp_ui
                        );

                        if (0 != mpz_invert(nsp_inv, y, n)) {
                            mpz_mul(
                                acc_xi->data[next_acc],
                                acc_xi->data[next_acc],
                                nsp_inv
                            );
                        } else {
                            //
                            // This should not happen: since y is a prime,
                            // it has an inverse in Z/nZ... unless n is a
                            // multiple of y! In such a (rare) case, we have
                            // found a factor of n, but handling it is quite
                            // cumbersome, so we just ignore it...
                            //
                            free(key);
                            //clear_mpz_pair(found);
                            //free(found);

                            mpz_clear(acc_yi->data[next_acc]);
                            mpz_clear(acc_xi->data[next_acc]);
                            continue;
                        }
                        acc_yi->length++;
                        acc_xi->length++;

                        free(key);
                        //clear_mpz_pair(found);
                        //free(found);

                        if (acc_yi->length == acc_yi->alloced) {
                            break;
                        }
                    } else {
                        //
                        // Add this pair in the hashtable with the position
                        // of the prime in first_primes_array as key...
                        //
                        mpz_pair_t* pair = malloc(sizeof(mpz_pair_t));
                        mpz_init_set(pair->x, in_xi->data[k]);
                        mpz_init_set(pair->y, in_yi->data[k]);
                        mpz_mul(pair->y, pair->y, in_cof->data[k]);
                        add_entry_in_hashtable(htable, (void*)key, (void*)pair);
                        //
                        // _WARNING_: Do not free key and pair as they are now
                        //            referenced by the hashtable!
                        //
                    }
                }
            } else {
                //
                // Here comes the early abort strategy!
                //
                // in_y->data[i] is not the product of a smooth number by a
                // large prime. If its non smooth part is greater than a certain
                // bound, just discard it. Otherwise it qualifies for the next
                // round! Congratulations!
                //
                if (mpz_cmpabs(y, bound) > 0) {
                    continue;
                }
                mpz_set(out_yi->data[next_out], y);
                mpz_set(out_xi->data[next_out], in_xi->data[k]);

                mpz_mul(
                    out_cof->data[next_out],
                    in_yi->data[k],
                    in_cof->data[k]
                );
                mpz_divexact(
                    out_cof->data[next_out],
                    out_cof->data[next_out],
                    y
                );

                out_yi->length++;
                out_xi->length++;
                out_cof->length++;

                if (out_yi->length == out_yi->alloced) {
                    break;
                }
            }
        }
    }

    mpz_clear(y);
    mpz_clear(nsp_inv);

    return used;
}
//------------------------------------------------------------------------------
uint32_t res_tdiv_first(smooth_filter_t* const filter) {

    mpz_array_t* const in_xi  = filter->candidate_xi;
    mpz_array_t* const in_yi  = filter->candidate_yi;
    mpz_array_t* const out_xi = filter->filtered_xi[0];
    mpz_array_t* const out_yi = filter->filtered_yi[0];
    mpz_array_t* const cofact = filter->cofactors[0];
    mpz_array_t* const acc_xi = filter->accepted_xi;
    mpz_array_t* const acc_yi = filter->accepted_yi;

    if (acc_yi->length == acc_yi->alloced) {
        return 0;
    }

    uint32_t next_acc;    
    uint32_t next_out;

    hashtable_t* const htable = filter->htable;
    mpz_ptr n     = filter->n;
    mpz_ptr bound = filter->bounds[0];

    uint32_t  nprimes = filter->factor_base[0]->length;
    uint32_t* base    = filter->factor_base[0]->data;

    uint32_t used = 0U;

    mpz_t y;
    mpz_init(y);

    mpz_t nsp_inv;
    mpz_init(nsp_inv);

    uint32_t k = in_yi->length;

#if TIFA_USE_GMP_INTERNAL_FUNCS
    mp_limb_t* num  = 0;
    uint32_t   size = 0;
    uint32_t   pi   = 0;
#endif

    while (k != 0U) {

        k--;
        used++;
        mpz_set(y, in_yi->data[k]);

#if TIFA_USE_GMP_INTERNAL_FUNCS
        //
        // Use some internal GMP functions.
        //
        // _NOTE_: We actually use the mpn_divexact_1 function which is not
        //         documented. See:
        //
        //  http://gmplib.org/list-archives/gmp-bugs/2006-December/000660.html
        //
        //         Apparently this "secret GMP function" (as Torbjorn Granlund
        //         put it) is not supposed to be used by client code... Could it
        //         crash your computer? Burn down your office? Give bad results?
        //         Or, more prosaically, will it be removed in future versions?
        //         In any case, maintainers should be aware of this...
        //
        num  = PTR(y);
        size = ABSIZ(y);
        pi   = 0;
        //
        // Divisions by 2 (always in the factor base)
        //
        if ( (num[0] & 1) == 0) {
            mpn_rshift(num, num, size, 1);
            MPN_NORMALIZE(num, size);
            SIZ(y) = size;

            while ( (num[0] & 1) == 0) {
                mpn_rshift(num, num, size, 1);
                MPN_NORMALIZE(num, size);
                SIZ(y) = size;
            }
        }
        //
        // Divisions by odd primes:
        //
        // _WARNING_: The following loop uses the "internal" mpn_modexact_1_odd
        //            and mpn_divexact_1 functions which could potentially be
        //            removed or changed in future GMP versions. See previous
        //            warning note.
        //
        for (pi = 1; pi < nprimes; pi++) {

            if (mpn_modexact_1_odd(num, size, base[pi]) == 0) {

                mpn_divexact_1(num, num, size, base[pi]);
                MPN_NORMALIZE(num, size);
                SIZ(y) = size;

                while (mpn_modexact_1_odd(num, size, base[pi]) == 0) {
                    mpn_divexact_1(num, num, size, base[pi]);
                    MPN_NORMALIZE(num, size);
                    SIZ(y) = size;
                }
            }
        }
#else
        //
        // GMP's internal functions are not available: revert to a simple
        // MPZ based implementation. (Using "public" MPN function doesn't
        // bring any significant benefits according to my tests, at least on
        // Opteron)
        //
        for (uint32_t i = 0; i < nprimes; i++) {
            if (mpz_divisible_ui_p(y, base[i])) {
                mpz_divexact_ui(y, y, base[i]);
                while (mpz_divisible_ui_p(y, base[i])) {
                    mpz_divexact_ui(y, y, base[i]);
                }
            }
        }
#endif

        next_acc = acc_yi->length;
        next_out = out_yi->length;

        if (0 == mpz_cmp_ui(y, 1)) {

            mpz_init_set(acc_yi->data[next_acc], in_yi->data[k]);
            mpz_init_set(acc_xi->data[next_acc], in_xi->data[k]);

            acc_yi->length++;
            acc_xi->length++;

            if (acc_xi->length == acc_xi->alloced) {
                break;
            }

        } else {
            //
    		// y is now the non-smooth part of in_yi->data[i].
    		// Is y a large prime?
    		//
            if (mpz_cmp_ui(y, LARGEST_PRIME) <= 0 ) {

                uint32_t nsp_ui = mpz_get_ui(y);
                //
                // Get the index of y in first_primes_array...
                // If ind != NOT_IN_ARRAY, then y is indeed a prime...
                //
                uint32_t ind = index_in_sorted_uint32_array(
                                    nsp_ui,
                                    &first_primes_array,
                                    0,
                                    first_primes_array.length);

                if (ind != NOT_IN_ARRAY) {
                    //
                    // cand_yi->data[k] is the product of a smooth number with
                    // the prime number first_primes_array->data[ind]
                    //
                    uint32_t *key = malloc(sizeof(uint32_t));
                    *key  = ind;

                    //
                    // _NOTE_: Actually it may be possible to leave the pair
                    //         in the hashtable...
                    //
                    // _TO_DO: Try to keep the entry in the hashtable...
                    //
                //mpz_pair_t *found = remove_entry_in_hashtable(htable, key);

                    mpz_pair_t *found = get_entry_in_hashtable(htable,key);

                    if (found != NULL) {
                        //
                        // The hashtable did contain another entry for this
                        // particular prime number: we can compute a new
                        // pair to add in our accepted_* arrays.
                        //
                        mpz_init(acc_yi->data[next_acc]);
                        mpz_init(acc_xi->data[next_acc]);
                        mpz_mul(
                            acc_yi->data[next_acc],
                            in_yi->data[k],
                            found->y
                        );
                        mpz_mul(
                            acc_xi->data[next_acc],
                            in_xi->data[k],
                            found->x
                        );
                        mpz_divexact_ui(
                            acc_yi->data[next_acc],
                            acc_yi->data[next_acc],
                            nsp_ui
                        );
                        mpz_divexact_ui(
                            acc_yi->data[next_acc],
                            acc_yi->data[next_acc],
                            nsp_ui
                        );

                        if (0 != mpz_invert(nsp_inv, y, n)) {
                            mpz_mul(
                                acc_xi->data[next_acc],
                                acc_xi->data[next_acc],
                                nsp_inv
                            );
                        } else {
                            //
                            // This should not happen: since y is a prime,
                            // it has an inverse in Z/nZ... unless n is a
                            // multiple of y! In such a (rare) case, we have
                            // found a factor of n, but handling it is quite
                            // cumbersome, so we just ignore it...
                            //
                            free(key);
                            //clear_mpz_pair(found);
                            //free(found);

                            mpz_clear(acc_yi->data[next_acc]);
                            mpz_clear(acc_xi->data[next_acc]);
                            continue;
                        }
                        acc_yi->length++;
                        acc_xi->length++;

                        free(key);
                        //clear_mpz_pair(found);
                        //free(found);

                        if (acc_yi->length == acc_yi->alloced) {
                            break;
                        }
                    } else {
                        //
                        // Add this pair in the hashtable with the position
                        // of the prime in first_primes_array as key...
                        //
                        mpz_pair_t* pair = malloc(sizeof(mpz_pair_t));
                        mpz_init_set(pair->x, in_xi->data[k]);
                        mpz_init_set(pair->y, in_yi->data[k]);
                        add_entry_in_hashtable(htable, (void*)key, (void*)pair);
                        //
                        // _WARNING_: Do not free key and pair as they are now
                        //            referenced by the hashtable!
                        //
                    }
                }
            } else {
                //
                // Here comes the early abort strategy!
                //
                // in_y->data[i] is not the product of a smooth number by a
                // large prime. If its non smooth part is greater than a certain
                // bound, just discard it. Otherwise it qualifies for the next
                // round! Congratulations!
                //
                if (mpz_cmpabs(y, bound) > 0) {
                    continue;
                }
                mpz_set(out_yi->data[next_out], y);
                mpz_set(out_xi->data[next_out], in_xi->data[k]);
                mpz_set(cofact->data[next_out], in_yi->data[k]);
                mpz_divexact(cofact->data[next_out], cofact->data[next_out], y);

                out_yi->length++;
                out_xi->length++;
                cofact->length++;

                if (out_yi->length == out_yi->alloced) {
                    break;
                }
            }
        }
    }

    mpz_clear(y);
    mpz_clear(nsp_inv);

    return used;
}
//------------------------------------------------------------------------------
uint32_t res_tdiv_last(smooth_filter_t* const filter) {

    unsigned long int step = filter->nsteps - 1;

    mpz_array_t* const in_xi  = filter->filtered_xi[step];
    mpz_array_t* const in_yi  = filter->filtered_yi[step];
    mpz_array_t* const cofact = filter->cofactors[step];
    mpz_array_t* const acc_xi = filter->accepted_xi;
    mpz_array_t* const acc_yi = filter->accepted_yi;

    if (acc_yi->length == acc_yi->alloced) {
        return 0;
    }

    uint32_t next;

    hashtable_t* htable = filter->htable;
    mpz_ptr n           = filter->n;

    uint32_t  nprimes = filter->factor_base[step + 1]->length;
    uint32_t* base    = filter->factor_base[step + 1]->data;

    uint32_t used = 0U;

    mpz_t y;
    mpz_init(y);

    mpz_t nsp_inv;
    mpz_init(nsp_inv);

    uint32_t k = in_yi->length;

#if TIFA_USE_GMP_INTERNAL_FUNCS
    mp_limb_t* num  = 0;
    uint32_t   size = 0;
    uint32_t   pi   = 0;
#endif

    while (k != 0U) {

        k--;
        used++;
        mpz_set(y, in_yi->data[k]);

#if TIFA_USE_GMP_INTERNAL_FUNCS
        //
        // Use some internal GMP functions.
        //
        // _NOTE_: We actually use the mpn_divexact_1 function which is not
        //         documented. See:
        //
        //  http://gmplib.org/list-archives/gmp-bugs/2006-December/000660.html
        //
        //         Apparently this "secret GMP function" (as Torbjorn Granlund
        //         put it) is not supposed to be used by client code... Could it
        //         crash your computer? Burn down your office? Give bad results?
        //         Or, more prosaically, will it be removed in future versions?
        //         In any case, maintainers should be aware of this...
        //
        num  = PTR(y);
        size = ABSIZ(y);
        pi   = 0;
        //
        // Divisions by 2 (always in the factor base)
        //
        if ( (num[0] & 1) == 0) {
            mpn_rshift(num, num, size, 1);
            MPN_NORMALIZE(num, size);
            SIZ(y) = size;

            while ( (num[0] & 1) == 0) {
                mpn_rshift(num, num, size, 1);
                MPN_NORMALIZE(num, size);
                SIZ(y) = size;
            }
        }
        //
        // Divisions by odd primes:
        //
        // _WARNING_: The following loop uses the "internal" mpn_modexact_1_odd
        //            and mpn_divexact_1 functions which could potentially be
        //            removed or changed in future GMP versions. See previous
        //            warning note.
        //
        for (pi = 1; pi < nprimes; pi++) {

            if (mpn_modexact_1_odd(num, size, base[pi]) == 0) {

                mpn_divexact_1(num, num, size, base[pi]);
                MPN_NORMALIZE(num, size);
                SIZ(y) = size;

                while (mpn_modexact_1_odd(num, size, base[pi]) == 0) {
                    mpn_divexact_1(num, num, size, base[pi]);
                    MPN_NORMALIZE(num, size);
                    SIZ(y) = size;
                }
            }
        }
#else
        //
        // GMP's internal functions are not available: revert to a simple
        // MPZ based implementation. (Using "public" MPN function doesn't
        // bring any significant benefits according to my tests, at least on
        // Opteron)
        //
        for (uint32_t i = 0; i < nprimes; i++) {
            if (mpz_divisible_ui_p(y, base[i])) {
                mpz_divexact_ui(y, y, base[i]);
                while (mpz_divisible_ui_p(y, base[i])) {
                    mpz_divexact_ui(y, y, base[i]);
                }
            }
        }
#endif

        next = acc_yi->length;

        if (0 == mpz_cmp_ui(y, 1)) {

            mpz_init_set(acc_yi->data[next], in_yi->data[k]);
            mpz_init_set(acc_xi->data[next], in_xi->data[k]);

            mpz_mul(acc_yi->data[next], acc_yi->data[next], cofact->data[k]);

            acc_xi->length++;
            acc_yi->length++;

            if (acc_xi->length == acc_xi->alloced) {
                break;
            }

        } else {
            //
    		// y is now the non-smooth part of in_yi->data[i].
    		// Is y a large prime?
    		//
            if (mpz_cmp_ui(y, LARGEST_PRIME) <= 0 ) {

                uint32_t nsp_ui = mpz_get_ui(y);
                //
                // Get the index of y in first_primes_array...
                // If ind != NOT_IN_ARRAY, then y is indeed a prime...
                //
                uint32_t ind = index_in_sorted_uint32_array(
                                    nsp_ui,
                                    &first_primes_array,
                                    0,
                                    first_primes_array.length);

                if (ind != NOT_IN_ARRAY) {
                    //
                    // cand_yi->data[k] is the product of a smooth number with
                    // the prime number first_primes_array->data[ind]
                    //
                    uint32_t *key = malloc(sizeof(uint32_t));
                    *key  = ind;

                    //
                    // _NOTE_: Actually it may be possible to leave the pair
                    //         in the hashtable...
                    //
                    // _TO_DO: Try to keep the entry in the hashtable...
                    //
                //mpz_pair_t *found = remove_entry_in_hashtable(htable, key);

                    mpz_pair_t *found = get_entry_in_hashtable(htable,key);

                    if (found != NULL) {
                        //
                        // The hashtable did contain another entry for this
                        // particular prime number: we can compute a new
                        // pair to add in our accepted_* arrays.
                        //
                        mpz_init(acc_yi->data[next]);
                        mpz_init(acc_xi->data[next]);

                        mpz_mul(acc_yi->data[next], in_yi->data[k], found->y);
                        mpz_mul(acc_xi->data[next], in_xi->data[k], found->x);

                        mpz_mul(
                            acc_yi->data[next],
                            acc_yi->data[next],
                            cofact->data[k]
                        );

                        mpz_divexact_ui(
                            acc_yi->data[next],
                            acc_yi->data[next],
                            nsp_ui
                        );

                        mpz_divexact_ui(
                            acc_yi->data[next],
                            acc_yi->data[next],
                            nsp_ui
                        );

                        if (0 != mpz_invert(nsp_inv, y, n)) {
                            mpz_mul(
                                acc_xi->data[next],
                                acc_xi->data[next],
                                nsp_inv
                            );
                        } else {
                            //
                            // This should not happen: since y is a prime,
                            // it has an inverse in Z/nZ... unless n is a
                            // multiple of y! In such a (rare) case, we have
                            // found a factor of n, but handling it is quite
                            // cumbersome, so we just ignore it...
                            //
                            free(key);
                            //clear_mpz_pair(found);
                            //free(found);

                            mpz_clear(acc_yi->data[next]);
                            mpz_clear(acc_xi->data[next]);
                            continue;
                        }
                        acc_yi->length++;
                        acc_xi->length++;

                        free(key);
                        //clear_mpz_pair(found);
                        //free(found);

                        if (acc_yi->length == acc_yi->alloced) {
                            break;
                        }
                    } else {
                        //
                        // Add this pair in the hashtable with the position
                        // of the prime in first_primes_array as key...
                        //
                        mpz_pair_t* pair = malloc(sizeof(mpz_pair_t));
                        mpz_init_set(pair->x, in_xi->data[k]);
                        mpz_init_set(pair->y, in_yi->data[k]);
                        mpz_mul(pair->y, pair->y, cofact->data[k]);
                        add_entry_in_hashtable(htable, (void*)key, (void*)pair);
                        //
                        // _WARNING_: Do not free key and pair as they are now
                        //            referenced by the hashtable!
                        //
                    }
                }
            }
        }
    }
    mpz_clear(y);
    mpz_clear(nsp_inv);

    return used;
}
//------------------------------------------------------------------------------

