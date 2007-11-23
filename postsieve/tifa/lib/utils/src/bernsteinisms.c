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
 * \file    bernsteinisms.c
 * \author  Jerome Milan
 * \date    Fri Oct 12 2007
 * \version 1.1
 */

 /*
  *  History:
  *
  *  1.1:   Fri Oct 12 2007 by JM:
  *         - Added multi-step early abort strategy (see smooth_filter.h)
  *
  *  1.0.2: Wed Apr 18 2007 by JM:
  *         - Shorten function names (nothing else!) even though we were
  *           already safe since the limit for unique identification is 63
  *           characters for C99...
  *
  *  1.0.1: Fri Dec 15 2006 by JM:
  *         - Bug fixed in bern_51 (c was sometimes improperly computed, duh!)
  *
  *  1.0.0: Fri Mar 3 2006 by JM:
  *         - Initial version.
  */

#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "tifa_config.h"
#include "gmp_utils.h"
#include "bernsteinisms.h"
#include "x_tree.h"
#include "funcs.h"

#include "first_primes.h"
#include "hashtable.h"
#include "macros.h"

#define DEBUG 0
#if DEBUG
    #include <stdio.h>
    #define PRINTF(...)     printf(__VA_ARGS__); fflush(stdout); 
    #define GMP_PRINTF(...) gmp_printf(__VA_ARGS__); fflush(stdout); 
#else
    #define PRINTF(...)     /* */
    #define GMP_PRINTF(...) /* */
#endif

//-----------------------------------------------------------------------------
//                 PROTOTYPES OF NON PUBLIC FUNCTION(S)
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
void bern_63_rec(
    uint32_array_t*,
    const mpz_t,
    const mpz_t[],
    uint32_t nb_nodes,
    uint32_t curnode
);
//------------------------------------------------------------------------------
void bern_71_rec(
    uint32_array_list_t* const,
    uint32_t* const,
    uint32_array_t*,
    const mpz_array_t* const,
    uint32_t
);
//------------------------------------------------------------------------------
uint32_t djb_batch_rt_step(smooth_filter_t* const, unsigned long int);
//------------------------------------------------------------------------------
uint32_t djb_batch_rt_first(smooth_filter_t* const filter);
//------------------------------------------------------------------------------
uint32_t djb_batch_rt_last(smooth_filter_t* const filter);
//------------------------------------------------------------------------------
uint32_t djb_batch_rt_no_ea(smooth_filter_t* const filter);
//------------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                 "PUBLIC" FUNCTION(S) IMPLEMENTATION
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Ref. "How to find small factors of integers", Daniel J. Bernstein
//       http://cr.yp.to/papers/sf.pdf
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
mpz_t* bern_51(uint32_t b, const mpz_t u) {
    //
    // _NOTE_: This is the algorithm 5.1 as explained in D. J. Bernstein's
    //         paper "How to find small factors of integers".
    //
    //         Variable names and notations are the ones used in the paper.
    //
    mpz_t* v = malloc(sizeof(mpz_t));
    mpz_init(*v);
    //
    // Step 1
    //
    if (1 == b) {
        mpz_set_ui(*v, 1);
        return v;
    }
    //
    // Step 2
    //
    uint32_t c;
    c = (b & 1) + (b >> 1);

    //
    // Step 3
    //
    mpz_t* v0 = bern_51(c, u);
    //
    // Step 4
    //
    mpz_t u0;
    mpz_t u1;
    mpz_init(u0);
    mpz_init(u1);

    mpz_tdiv_r_2exp(u0, u, c);
    mpz_tdiv_q_2exp(u1, u, c);
    mpz_tdiv_r_2exp(u1, u1, c);
    //
    // Step 5
    //
    mpz_t z;
    mpz_init_set_ui(z, 1);

    mpz_addmul(z, u0, *v0);
    mpz_tdiv_q_2exp(z, z, c);
    mpz_addmul(z, u1, *v0);
    mpz_tdiv_r_2exp(z, z, c);
    //
    // Step 6
    //
    mpz_mul(u0, *v0, z);
    mpz_mul_2exp(u0, u0, c);
    mpz_tdiv_r_2exp(u1, u0, b);
    mpz_add(*v, *v0, u1);
    //
    // Clean up
    //
    mpz_clear(*v0);
    free(v0);
    mpz_clear(u0);
    mpz_clear(u1);
    mpz_clear(z);

    return v;
}
//------------------------------------------------------------------------------
mpz_t* bern_53(uint32_t b, const mpz_t u, const mpz_t x) {
    //
    // _NOTE_: This is the algorithm 5.3 as explained in D. J. Bernstein's
    //         paper "How to find small factors of integers".
    //
    //         Variable names and notations are the ones used in the paper.
    //
    mpz_t *r = malloc(sizeof(mpz_t));
    mpz_init(*r);
    //
    // Step 1
    //
    mpz_t *v = bern_51(b, u);
    //
    // Step 2
    //
    mpz_t x0;
    mpz_t x1;
    mpz_init(x0);
    mpz_init(x1);

    mpz_tdiv_q_2exp(x1, x, b);
    mpz_tdiv_r_2exp(x0, x, b);
    //
    // Step 3
    //
    mpz_t q;
    mpz_init(q);

    mpz_mul(q, *v, x0);
    mpz_tdiv_r_2exp(q, q, b);
    //
    // Step 4
    //
    mpz_set(*r, x0);
    mpz_addmul(*r, u, q);
    mpz_tdiv_q_2exp(*r, *r, b);
    mpz_add(*r, *r, x1);
    //
    // Clean up
    //
    mpz_clear(*v);
    free(v);
    mpz_clear(x0);
    mpz_clear(x1);
    mpz_clear(q);

    return r;
}
//------------------------------------------------------------------------------
uint32_array_t* bern_63(const mpz_t x, const mpz_tree_t* const tree) {
    //
    // _NOTE_: This is the algorithm 6.3 as explained in D. J. Bernstein's
    //         paper "How to find small factors of integers".
    //
    //         Variable names and notations are (mostly) the ones used
    //         in the paper.
    //

    //
    // Bootstrap...
    //
    uint32_array_t* result = malloc(sizeof(mpz_array_t));

    if (0 == tree->length) {
        result->alloced = 0;
        result->length = 0;
        result->data = NULL;
        return result;
    }
    result->alloced = (tree->length + 1) / 2;
    result->length  = 0;
    result->data    = malloc(result->alloced*sizeof(uint32_t));
    //
    // Recursion: steps 1 to 7
    //
    bern_63_rec(result, x, (const mpz_t*)tree->data, tree->length, 0);

    return result;
}
//------------------------------------------------------------------------------
void bern_63_rec(uint32_array_t* result, const mpz_t x,
                 const mpz_t tree[],
                 uint32_t nb_nodes, uint32_t curnode) {
    //
    // _NOTE_: This is the algorithm 6.3 as explained in D. J. Bernstein's
    //         paper "How to find small factors of integers".
    //
    //         Variable names and notations are (mostly) the ones used
    //         in the paper.
    //

    //
    // Step 1
    //
    const mpz_t* const u = &(tree[curnode]);
    //
    // Step 2
    //
    uint32_t c = mpz_sizeinbase(*u, 2);
    uint32_t d = mpz_sizeinbase(x, 2);

    mpz_t *r;
    if (d > (c + 1)) {
        //
        // Step 3
        //
        r = bern_53(d - c, *u, x);

    } else {
        //
        // Step 4
        //
        r = malloc(sizeof(mpz_t));
        mpz_init_set(*r, x);
    }
    //
    // Step 5
    //
    uint8_t flag = 0;

    if ((1 == nb_nodes)) {
        //
        // The root of the tree has no children
        //
        mpz_t uc;
        mpz_init_set(uc, *u);
        //
        // If (u != 1) : check if r is in {0, u, 2*u, 3*u}
        // If (u == 1) : useless trivial case...
        //
        if (mpz_cmp_ui(*u, 1) != 0) {

            if ((mpz_cmp_ui(*r, 0) == 0) || (mpz_cmp(*r, uc) == 0)) {
                flag++;
            } else {
                mpz_add(uc, uc, *u);
                if ((mpz_cmp(*r, uc) == 0)) {
                    flag++;
                } else {
                    mpz_add(uc, uc, *u);
                    if ((mpz_cmp(*r, uc) == 0)) {
                        flag++;
                    }
                }
            }
        }
        if (flag != 0) {
            result->data[result->length] = mpz_get_ui(*u);
            result->length++;
        }
        mpz_clear(uc);

    } else {
        //
        // Step 6 and 7
        //
        bern_63_rec(result, x, tree, (nb_nodes - 1)/2, 2*curnode + 1);
        bern_63_rec(result, x, tree, (nb_nodes - 1)/2, 2*curnode + 2);
    }
    mpz_clear(*r);
    free(r);
}
//------------------------------------------------------------------------------
void bern_71(uint32_array_list_t* const decomp_list,
             const mpz_array_t* const to_be_factored,
             const uint32_array_t* const odd_primes) {
    //
    // _NOTE_: This is the algorithm 7.1 as explained in D. J. Bernstein's
    //         paper "How to find small factors of integers".
    //
    //         Variable names and notations are (mostly) the ones used
    //         in the paper.
    //

    //
    // Bootstrap. Step 1
    //
    if (0U == to_be_factored->length) {
        for (uint32_t i = 0; i < to_be_factored->length; i++) {
            add_entry_in_uint32_array_list(NULL, decomp_list);
        }
        return;
    }
    //
    // Step 2 - processed only once.
    //
    mpz_tree_t* x = prod_tree(to_be_factored);
    uint32_t nb_int_to_factor = to_be_factored->length;
    //
    // Recursion: steps 3 to 8
    //
    bern_71_rec(decomp_list, &nb_int_to_factor, (uint32_array_t*)odd_primes,
                x, 0);

    clear_mpz_tree(x);
    free(x);
}
//------------------------------------------------------------------------------
void bern_71_rec(uint32_array_list_t* const decomp_list,
                 uint32_t * const nb_int_to_factor,
                 uint32_array_t*  odd_primes,
                 const mpz_array_t* const x, uint32_t curnode) {

    //
    // _TO_DO_: Lots of redundant computations here! Fix this!
    //

    //
    // Step 3
    //
    mpz_tree_t *T = prod_tree_ui((const uint32_array_t*)odd_primes);
    //
    // Step 4
    //
    uint32_array_t *P = bern_63(x->data[curnode], (const mpz_array_t*)T);

    if (0 == P->length) {
        //
        // Step 1 (Yes, indeed! Why call bern_71_rec another time since we
        //         can handle it now...)
        //
        // The list of primes such that x->data[curnode] mod p = 0 is empty.
        // Thus, for every node i in the subtree of root curnode, the list
        // such that subtree[i] mod p = 0 will be empty. Consequently, we have
        // to add in our list as many null lists as there is nodes in the
        // subtree of root curnode.
        //
        // Instead of invoking recursively bern_71_rec, compute this number of
        // nodes by finding the depth of curnode is the global x tree...
        //
        uint8_t depth_curnode = 0;
        while ((1U<<depth_curnode) <= (curnode+1)) {
            depth_curnode++;
        }
        depth_curnode--;

        uint32_t nb_leafs = (x->length + 1)/2;
        uint32_t nb_leafs_subtree = nb_leafs/(1<<depth_curnode);
        //
        // ... and add as many null list as required...
        //
        for (uint32_t i = 0; i < nb_leafs_subtree; i++) {
            if (*nb_int_to_factor != 0U) {
                add_entry_in_uint32_array_list(NULL, decomp_list);
                (*nb_int_to_factor)--;
            }
        }
        clear_uint32_array(P);
        free(P);
        clear_mpz_tree(T);
        free(T);

        return;
    }
    if (curnode >= (x->length/2)) {
        //
        // i.e. curnode is a leaf...
        //
        // Now, check *nb_int_to_factor as there may be useless nodes (i.e.
        // nodes with value 1) in the product tree of all the integers to
        // factor...
        //
        if (*nb_int_to_factor != 0U) {
            //
            // Step 5.
            //
            add_entry_in_uint32_array_list(P, decomp_list);
            (*nb_int_to_factor)--;
        } else {
            clear_uint32_array(P);
            free(P);
        }
    } else {
        //
        // Step 6, 7 and 8
        //
        bern_71_rec(decomp_list, nb_int_to_factor, P, x, 2*curnode+1);
        bern_71_rec(decomp_list, nb_int_to_factor, P, x, 2*curnode+2);

        clear_uint32_array(P);
        free(P);
    }
    clear_mpz_tree(T);
    free(T);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Ref. "How to find smooth part of integers", Daniel J. Bernstein
//       http://cr.yp.to/factorization/smoothparts-20040510.pdf
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
uint32_t bern_21_rt(mpz_array_t* const smooth,
                    const mpz_array_t* const xi,
                    const mpz_t z) {
    //
    // _NOTE_: This is the algorithm 2.1 as explained in D. J. Bernstein's
    //         paper "How to find smooth part of integers".
    //
    //         Variable names and notations are (mostly) the ones used
    //         in the paper.
    //
    if (smooth->length == smooth->alloced) {
        return 0;
    }
    uint32_t retval = xi->length;
    //
    // Step 1 is already precomputed (z = root(prod_tree_ui(pj)))
    //
    // Step 2
    //
    mpz_tree_t* pxitree = prod_tree(xi);
    mpz_tree_t* rtree   = rem_tree(z, pxitree);

    clear_mpz_tree(pxitree);
    free(pxitree);
    //
    // Step 3
    //
    uint32_t e   = 0U;
    uint32_t msb = 0U;

    mpz_t yk;
    mpz_init2(yk, mpz_sizeinbase(xi->data[0], 2) + 8);

    for (uint32_t k = 0U; k < xi->length; k++) {
        msb = mpz_sizeinbase(xi->data[k], 2) - 1;
        e   = most_significant_bit(msb) + 1;

        mpz_set(yk, rtree->data[rtree->length/2 + k]);
        mpz_powm_ui(yk, yk, 1<<e, xi->data[k]);

        //
        // _NOTE_: xi->data[k] is smooth iff (yk == 0) (no need to compute
        //         gcd if we just need to know whether xi->data[k] is smooth
        //         or not).
        //
        if (0 == mpz_cmp_ui(yk, 0)) {
            mpz_init_set(smooth->data[smooth->length], xi->data[k]);
            smooth->length++;
            if (smooth->length == smooth->alloced) {
                retval = k + 1;
                break;
            }
        }
        //
        // _NOTE_: Comment previous if-branch and uncomment the following
        //         if-branch to compute the smooth part of each xi->data[k]
        //
        //mpz_t gcd;
        //mpz_init(gcd);
        //mpz_gcd(gcd, yk, xi->data[k]);
        //if (0 == mpz_cmp(xi->data[k], gcd)) {
        //    mpz_init_set(smooth->data[smooth->length], xi->data[k]);
        //    smooth->length++;
        //    if (smooth->length == smooth->alloced) {
        //        mpz_clear(gcd);
        //        break;
        //    }
        //}
        //mpz_clear(gcd);
    }
    mpz_clear(yk);

    clear_mpz_tree(rtree);
    free(rtree);

    return retval;
}
//------------------------------------------------------------------------------
uint32_t bern_21(mpz_array_t* const smooth,
                 const mpz_array_t* const xi,
                 const mpz_t z) {
    //
    // _NOTE_: Same as bern_21_rt except that no remainder tree is computed
    //         (modulos are computed individually).
    //
    if (smooth->length == smooth->alloced) {
        return 0;
    }
    uint32_t retval = xi->length;
    //
    // Step 1 is already precomputed (z = prod_tree_ui(pj))
    //
    // Step 3 (no step 2 as we don't use a remainder tree)
    //
    uint32_t e = 0U;
    uint32_t msb = 0U;
    mpz_t yk;
    mpz_init(yk);
    for (uint32_t k = 0U; k < xi->length; k++) {
        msb = mpz_sizeinbase(xi->data[k], 2) - 1;
        e   = most_significant_bit(msb) + 1;

        mpz_set(yk, z);
        mpz_mod(yk, yk, xi->data[k]);
        mpz_powm_ui(yk, yk, 1<<e, xi->data[k]);

        //
        // _NOTE_: xi->data[k] is smooth iff (yk == 0) (no need to compute
        //         gcd if we just need to know whether xi->data[k] is smooth
        //         or not).
        //
        if (0 == mpz_cmp_ui(yk, 0)) {
            mpz_init_set(smooth->data[smooth->length], xi->data[k]);
            smooth->length++;
            if (smooth->length == smooth->alloced) {
                retval = k + 1;
                break;
            }
        }
        //
        // _NOTE_: Comment previous if-branch and uncomment the following
        //         if-branch to compute the smooth part of each xi->data[k]
        //
        //mpz_t gcd;
        //mpz_init(gcd);
        //mpz_gcd(gcd, yk, xi->data[k]);
        //if (0 == mpz_cmp(xi->data[k], gcd)) {
        //    mpz_init_set(smooth->data[smooth->length], xi->data[k]);
        //    smooth->length++;
        //    if (smooth->length == smooth->alloced) {
        //        mpz_clear(gcd);
        //        break;
        //    }
        //}
        //mpz_clear(gcd);
    }
    mpz_clear(yk);

    return retval;
}
//------------------------------------------------------------------------------
uint32_t bern_21_rt_pairs(mpz_array_t* const xi,
                          mpz_array_t* const smooth_yi,
                          const mpz_array_t* const cand_xi,
                          const mpz_array_t* const cand_yi,
                          const mpz_t z) {
    //
    // _NOTE_: Same as bern_21_rt but tailored to better suit our
    //         integer factorization problem. Stores the smooth numbers from
    //         cand_yi together with their corresponding numbers in cand_xi
    //         respectively in smooth_yi and xi.
    //
    if (smooth_yi->length == smooth_yi->alloced) {
        return 0;
    }
    uint32_t retval = cand_xi->length;
    //
    // Step 1 is already precomputed (z = root(prod_tree_ui(pj)))
    //
    // Step 2
    //
    mpz_tree_t* pyitree = prod_tree(cand_yi);
    mpz_tree_t* rtree   = rem_tree(z, pyitree);

    clear_mpz_tree(pyitree);
    free(pyitree);
    //
    // Step 3
    //
    uint32_t e = 0U;
    uint32_t msb = 0U;

    mpz_t yk;
    mpz_init(yk);

    for (uint32_t k = 0U; k < cand_yi->length; k++) {
        msb = mpz_sizeinbase(cand_yi->data[k], 2) - 1;
        e   = most_significant_bit(msb) + 1;

        mpz_set(yk, rtree->data[rtree->length/2 + k]);
        mpz_powm_ui(yk, yk, 1<<e, cand_yi->data[k]);
        //
        // _NOTE_: cand_yi->data[k] is smooth iff (yk == 0) (no need to compute
        //         gcd if we just need to know whether cand_yi->data[k] is
        //         smooth or not).
        //
        if (0 == mpz_sgn(yk)) {
            mpz_init_set(smooth_yi->data[smooth_yi->length], cand_yi->data[k]);
            mpz_init_set(xi->data[xi->length], cand_xi->data[k]);

            smooth_yi->length++;
            xi->length++;

            if (smooth_yi->length == smooth_yi->alloced) {
                retval = k + 1;
                break;
            }
        }
    }
    mpz_clear(yk);

    clear_mpz_tree(rtree);
    free(rtree);

    return retval;
}
//------------------------------------------------------------------------------
uint32_t bern_21_pairs(mpz_array_t* const xi,
                       mpz_array_t* const smooth_yi,
                       const mpz_array_t* const cand_xi,
                       const mpz_array_t* const cand_yi,
                       const mpz_t z) {
    //
    // _NOTE_: Same as bern_21_rt_pairs except that no remainder tree is
    //         computed (modulos are computed individually).
    //
    // _KLUDGE_: The major part of the code in the bern_21_rt_pairs and
    //           bern_21_pairs is identical. Such code reuse on a
    //           wide scale is certainly a very bad thing and one can
    //           legitimately argue that this is a huge design mistake.
    //
    // _TO_DO_: Try to refactor these two function to minimize code
    //          copy/pasting while avoiding numerous function calls.
    //          The easiest may be just to add yet another parameter and
    //          perform a test at runtime (performance penalty will be just
    //          lost in the noise anyway).
    //
    if (smooth_yi->length == smooth_yi->alloced) {
        return 0;
    }
    uint32_t retval = cand_xi->length;
    //
    // Step 1 is already precomputed (z = root(prod_tree_ui(pj)))
    //
    // Step 3 (no step 2 as we don't use a remainder tree)
    //
    uint32_t e = 0U;
    uint32_t msb = 0U;

    mpz_t yk;
    mpz_init(yk);

    for (uint32_t k = 0U; k < cand_yi->length; k++) {
        msb = mpz_sizeinbase(cand_yi->data[k], 2) - 1;
        e   = most_significant_bit(msb) + 1;

        mpz_set(yk, z);
        mpz_mod(yk, yk, cand_yi->data[k]);
        mpz_powm_ui(yk, yk, 1<<e, cand_yi->data[k]);
        //
        // _NOTE_: cand_yi->data[k] is smooth iff (yk == 0) (no need to compute
        //         the gcd if we just need to know whether cand_yi->data[k] is
        //         smooth or not).
        //
        if (0 == mpz_sgn(yk)) {
            mpz_init_set(smooth_yi->data[smooth_yi->length], cand_yi->data[k]);
            mpz_init_set(xi->data[xi->length], cand_xi->data[k]);

            smooth_yi->length++;
            xi->length++;

            if (smooth_yi->length == smooth_yi->alloced) {
                retval = k + 1;
                break;
            }
        }
    }
    mpz_clear(yk);

    return retval;
}
//------------------------------------------------------------------------------
uint32_t bern_21_rt_pairs_lp(const mpz_t n,
                             hashtable_t* const htable,
                             mpz_array_t* const xi,
                             mpz_array_t* const smooth_yi,
                             const mpz_array_t* const cand_xi,
                             const mpz_array_t* const cand_yi,
                             const mpz_t z) {
    //
    // _NOTE_: Same as bern_21_rt_pairs but uses the large prime variation.
    //
    if (smooth_yi->length == smooth_yi->alloced) {
        return 0;
    }
    uint32_t retval = cand_xi->length;
    //
    // Step 1 is already precomputed (z = root(prod_tree_ui(pj)))
    //
    // Step 2
    //
    mpz_tree_t* pyitree = prod_tree(cand_yi);
    mpz_tree_t* rtree   = rem_tree(z, pyitree);

    clear_mpz_tree(pyitree);
    free(pyitree);
    //
    // Step 3
    //
    uint32_t e = 0U;
    uint32_t msb = 0U;

    mpz_t yk;
    mpz_init(yk);

    mpz_t gcd;
    mpz_init(gcd);

    mpz_t nsp;
    mpz_init(nsp);

    mpz_t nsp_inv;
    mpz_init(nsp_inv);

    for (uint32_t k = 0U; k < cand_yi->length; k++) {
        msb = mpz_sizeinbase(cand_yi->data[k], 2) - 1;
        e   = most_significant_bit(msb) + 1;

        mpz_set(yk, rtree->data[rtree->length/2 + k]);
        mpz_powm_ui(yk, yk, 1<<e, cand_yi->data[k]);
        //
        // _NOTE_: cand_yi->data[k] is smooth iff (yk == 0)
        //
        if (0 == mpz_sgn(yk)) {
			//
            // We found a good pair (cand_yi->data[k], cand_xi->data[k])!
            // Add cand_yi->data[k] is the smooth_yi array and add
            // cand_xi->data[k] is the xi array.
            //
            mpz_init_set(smooth_yi->data[smooth_yi->length], cand_yi->data[k]);
            mpz_init_set(xi->data[xi->length], cand_xi->data[k]);

            smooth_yi->length++;
            xi->length++;

            if (smooth_yi->length == smooth_yi->alloced) {
                retval = k + 1;
                break;
            }
        } else {
            //
            // Here comes the large prime variation...
            //
            mpz_gcd(gcd, yk, cand_yi->data[k]);
			//
			// gcd is now the smooth part of cand_yi->data[k].
			// Is cand_yi->data[k] the product of gcd by a large prime?
			//
            mpz_divexact(nsp, cand_yi->data[k], gcd);

            if (mpz_cmpabs_ui(nsp, LARGEST_PRIME) <= 0 ) {

                uint32_t nsp_ui = mpz_get_ui(nsp);
                //
                // Get the index of (cand_yi->data[k] / gcd) in
                // first_primes_array... If ind != NOT_IN_ARRAY, then
                // (cand_yi->data[k] / gcd) is indeed a prime...
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
                    mpz_pair_t *found = remove_entry_in_hashtable(htable, key);

                    if (found != NULL) {
                        //
                        // The hashtable did contain another entry for this
                        // particular prime number: we can compute a new
                        // pair to add in our arrays xi and smooth_yi
                        //
                        mpz_init(smooth_yi->data[smooth_yi->length]);
                        mpz_init(xi->data[xi->length]);

                        mpz_mul(smooth_yi->data[smooth_yi->length],
                                cand_yi->data[k], found->y);

                        mpz_mul(xi->data[xi->length],
                                cand_xi->data[k], found->x);

                        mpz_divexact_ui(smooth_yi->data[smooth_yi->length],
                                        smooth_yi->data[smooth_yi->length],
                                        nsp_ui);
                        mpz_divexact_ui(smooth_yi->data[smooth_yi->length],
                                        smooth_yi->data[smooth_yi->length],
                                        nsp_ui);

                        if (0 != mpz_invert(nsp_inv, nsp, n)) {
                            mpz_mul(xi->data[xi->length], xi->data[xi->length],
                                    nsp_inv);
                        } else {
                            //
                            // This should not happen: since nsp is a prime,
                            // it has an inverse in Z/nZ... unless n is a
                            // multiple of nsp! In such a (rare) case, we have
                            // found a factor of n, but handling it is quite
                            // cumbersome, so we just ignore it...
                            //
                            free(key);
                            clear_mpz_pair(found);
                            free(found);

                            continue;
                        }
                        smooth_yi->length++;
                        xi->length++;

                        free(key);
                        clear_mpz_pair(found);
                        free(found);

                        if (smooth_yi->length == smooth_yi->alloced) {
                            retval = k + 1;
                            break;
                        }
                    } else {
                        //
                        // Add this pair in the hashtable with the position
                        // of the prime in first_primes_array as key...
                        //
                        mpz_pair_t* pair = malloc(sizeof(mpz_pair_t));
                        mpz_init_set(pair->x, cand_xi->data[k]);
                        mpz_init_set(pair->y, cand_yi->data[k]);
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
    mpz_clear(yk);
    mpz_clear(gcd);
    mpz_clear(nsp);
    mpz_clear(nsp_inv);

    clear_mpz_tree(rtree);
    free(rtree);

    return retval;
}
//------------------------------------------------------------------------------
uint32_t bern_21_pairs_lp(const mpz_t n,
                          hashtable_t* const htable,
                          mpz_array_t* const xi,
                          mpz_array_t* const smooth_yi,
                          const mpz_array_t* const cand_xi,
                          const mpz_array_t* const cand_yi,
                          const mpz_t z) {
    //
    // _NOTE_: Same as bern_21_pairs but uses the large prime variation.
    //
    // _KLUDGE_: The major part of the code in bern_21_rt_pairs_lp and in
    //           bern_21_pairs_lp is identical. Such code reuse on a wide scale
    //           is certainly a very bad thing and one can legitimately argue
    //           that this is a huge design mistake.
    //
    // _TO_DO_: Try to refactor these two function to minimize code
    //          copy/pasting while avoiding numerous function calls.
    //          The easiest may be just to add yet another parameter and
    //          perform a test at runtime (performance penalty will be just
    //          lost in the noise anyway).
    //
    if (smooth_yi->length == smooth_yi->alloced) {
        return 0;
    }
    uint32_t retval = cand_xi->length;
    //
    // Step 1 is already precomputed (z = root(prod_tree_ui(pj)))
    //
    // Step 3 (no step 2 as we don't use a remainder tree)
    //
    uint32_t e = 0U;
    uint32_t msb = 0U;

    mpz_t yk;
    mpz_init(yk);

    mpz_t gcd;
    mpz_init(gcd);

    mpz_t nsp;
    mpz_init(nsp);

    mpz_t nsp_inv;
    mpz_init(nsp_inv);

    for (uint32_t k = 0U; k < cand_yi->length; k++) {
        msb = mpz_sizeinbase(cand_yi->data[k], 2) - 1;
        e   = most_significant_bit(msb) + 1;

        mpz_set(yk, z);
        mpz_mod(yk, yk, cand_yi->data[k]);
        mpz_powm_ui(yk, yk, 1<<e, cand_yi->data[k]);

        //
        // _NOTE_: cand_yi->data[k] is smooth iff (yk == 0)
        //
        if (0 == mpz_cmp_ui(yk, 0)) {
			//
            // We found a good pair (cand_yi->data[k], cand_xi->data[k])!
            // Add cand_yi->data[k] is the smooth_yi array and add
            // cand_xi->data[k] is the xi array.
            //
            mpz_init_set(smooth_yi->data[smooth_yi->length], cand_yi->data[k]);
            mpz_init_set(xi->data[xi->length], cand_xi->data[k]);

            smooth_yi->length++;
            xi->length++;

            if (smooth_yi->length == smooth_yi->alloced) {
                retval = k + 1;
                break;
            }
        } else {
            //
            // Here comes the large prime variation...
            //
            mpz_gcd(gcd, yk, cand_yi->data[k]);
			//
			// gcd is now the smooth part of cand_yi->data[k].
			// Is cand_yi->data[k] the product of gcd by a large prime?
			//
            mpz_divexact(nsp, cand_yi->data[k], gcd);

            if (mpz_cmpabs_ui(nsp, LARGEST_PRIME) <= 0 ) {

                uint32_t nsp_ui = mpz_get_ui(nsp);
                //
                // Get the index of (cand_yi->data[k] / gcd) in
                // first_primes_array... If ind != NOT_IN_ARRAY, then
                // (cand_yi->data[k] / gcd) is indeed a prime...
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
                    mpz_pair_t *found = remove_entry_in_hashtable(htable, key);

                    if (found != NULL) {
                        //
                        // The hashtable did contain another entry for this
                        // particular prime number: we can compute a new
                        // pair to add in our arrays xi and smooth_yi
                        //
                        mpz_init(smooth_yi->data[smooth_yi->length]);
                        mpz_init(xi->data[xi->length]);

                        mpz_mul(smooth_yi->data[smooth_yi->length],
                                cand_yi->data[k], found->y);

                        mpz_mul(xi->data[xi->length],
                                cand_xi->data[k], found->x);

                        mpz_divexact_ui(smooth_yi->data[smooth_yi->length],
                                        smooth_yi->data[smooth_yi->length],
                                        nsp_ui);
                        mpz_divexact_ui(smooth_yi->data[smooth_yi->length],
                                        smooth_yi->data[smooth_yi->length],
                                        nsp_ui);

                        if (0 != mpz_invert(nsp_inv, nsp, n)) {
                            mpz_mul(xi->data[xi->length], xi->data[xi->length],
                                    nsp_inv);
                        } else {
                            //
                            // This should not happen: since nsp is a prime,
                            // it has an inverse in Z/nZ... unless n is a
                            // multiple of nsp! In such a (rare) case, we have
                            // found a factor of n, but handling it is quite
                            // cumbersome, so we just ignore it...
                            //
                            free(key);
                            clear_mpz_pair(found);
                            free(found);
                            continue;
                        }
                        smooth_yi->length++;
                        xi->length++;

                        free(key);
                        clear_mpz_pair(found);
                        free(found);

                        if (smooth_yi->length == smooth_yi->alloced) {
                            retval = k + 1;
                            break;
                        }
                    } else {
                        //
                        // Add this pair in the hashtable with the position
                        // of the prime in first_primes_array as key...
                        //
                        mpz_pair_t* pair = malloc(sizeof(mpz_pair_t));
                        mpz_init_set(pair->x, cand_xi->data[k]);
                        mpz_init_set(pair->y, cand_yi->data[k]);
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
    mpz_clear(yk);
    mpz_clear(gcd);
    mpz_clear(nsp);
    mpz_clear(nsp_inv);

    return retval;
}
//------------------------------------------------------------------------------
uint32_t bern_21_rt_pairs_siqs(mpz_array_t* const xi,
                               mpz_array_t* const smooth_yi,
                               mpz_array_t* const a_for_smooth_gx,
                               const mpz_array_t* const cand_xi,
                               const mpz_array_t* const cand_yi,
                               const mpz_array_t* const cand_a,
                               const mpz_t z) {
    //
    // _NOTE_: Same as bern_21_rt_pairs but specifically tailored to be used
    //         by the SIQS algorithm.
    //
    if (smooth_yi->length == smooth_yi->alloced) {
        return 0;
    }
    uint32_t retval = cand_xi->length;
    //
    // Step 1 is already precomputed (z = root(prod_tree_ui(pj)))
    //
    // Step 2
    //
    mpz_tree_t* pyitree = prod_tree(cand_yi);
    mpz_tree_t* rtree   = rem_tree(z, pyitree);

    clear_mpz_tree(pyitree);
    free(pyitree);
    //
    // Step 3
    //
    uint32_t e = 0U;
    uint32_t msb = 0U;

    mpz_t yk;
    mpz_init(yk);

    for (uint32_t k = 0U; k < cand_yi->length; k++) {
        msb = mpz_sizeinbase(cand_yi->data[k], 2) - 1;
        e   = most_significant_bit(msb) + 1;

        mpz_set(yk, rtree->data[rtree->length/2 + k]);
        mpz_powm_ui(yk, yk, 1<<e, cand_yi->data[k]);
        //
        // _NOTE_: cand_yi->data[k] is smooth iff (yk == 0)
        //
        if (0 == mpz_sgn(yk)) {
			//
            // We found a good pair (cand_yi->data[k], cand_xi->data[k])!
            // Add cand_yi->data[k] is the smooth_yi array and add
            // cand_xi->data[k] is the xi array.
            //
            mpz_init_set(smooth_yi->data[smooth_yi->length], cand_yi->data[k]);
            mpz_init_set(xi->data[xi->length], cand_xi->data[k]);
            mpz_init_set(
                a_for_smooth_gx->data[a_for_smooth_gx->length],
                cand_a->data[k]
            );
            smooth_yi->length++;
            xi->length++;
            a_for_smooth_gx->length++;

            if (smooth_yi->length == smooth_yi->alloced) {
                retval = k + 1;
                break;
            }
        }
    }
    mpz_clear(yk);
    
    clear_mpz_tree(rtree);
    free(rtree);
    
    return retval;
}
//------------------------------------------------------------------------------
uint32_t bern_21_rt_pairs_lp_siqs(const mpz_t n,
                                  hashtable_t* const htable,
                                  mpz_array_t* const xi,
                                  mpz_array_t* const smooth_yi,
                                  mpz_array_t* const a_for_smooth_gx,
                                  const mpz_array_t* const cand_xi,
                                  const mpz_array_t* const cand_yi,
                                  const mpz_array_t* const cand_a,
                                  const mpz_t z) {

    //
    // _NOTE_: Same as bern_21_rt_pairs_lp but specifically tailored to be used
    //         by the SIQS algorithm as we stores the coefficient "a" of the
    //         SIQS polynomial used.
    //
    if (smooth_yi->length == smooth_yi->alloced) {
        return 0;
    }
    uint32_t retval = cand_xi->length;
    //
    // Step 1 is already precomputed (z = root(prod_tree_ui(pj)))
    //
    // Step 2
    //
    mpz_tree_t* pyitree = prod_tree(cand_yi);
    mpz_tree_t* rtree   = rem_tree(z, pyitree);

    clear_mpz_tree(pyitree);
    free(pyitree);
    //
    // Step 3
    //
    uint32_t e = 0U;
    uint32_t msb = 0U;

    mpz_t yk;
    mpz_init(yk);

    mpz_t gcd;
    mpz_init(gcd);

    mpz_t nsp;
    mpz_init(nsp);

    mpz_t nsp_inv;
    mpz_init(nsp_inv);

    for (uint32_t k = 0U; k < cand_yi->length; k++) {
        msb = mpz_sizeinbase(cand_yi->data[k], 2) - 1;
        e   = most_significant_bit(msb) + 1;

        mpz_set(yk, rtree->data[rtree->length/2 + k]);
        mpz_powm_ui(yk, yk, 1<<e, cand_yi->data[k]);
        //
        // _NOTE_: cand_yi->data[k] is smooth iff (yk == 0)
        //
        if (0 == mpz_sgn(yk)) {
			//
            // We found a good pair (cand_yi->data[k], cand_xi->data[k])!
            // Add cand_yi->data[k] is the smooth_yi array and add
            // cand_xi->data[k] is the xi array.
            //
            mpz_init_set(smooth_yi->data[smooth_yi->length], cand_yi->data[k]);
            mpz_init_set(xi->data[xi->length], cand_xi->data[k]);
            mpz_init_set(
                a_for_smooth_gx->data[a_for_smooth_gx->length],
                cand_a->data[k]
            );
            smooth_yi->length++;
            xi->length++;
            a_for_smooth_gx->length++;

            if (smooth_yi->length == smooth_yi->alloced) {
                retval = k + 1;
                break;
            }
        } else {
            //
            // Here comes the large prime variation...
            //
            mpz_gcd(gcd, yk, cand_yi->data[k]);
			//
			// gcd is now the smooth part of cand_yi->data[k].
			// Is cand_yi->data[k] the product of gcd by a large prime?
			//
            mpz_divexact(nsp, cand_yi->data[k], gcd);
            mpz_abs(nsp, nsp);

            if (mpz_cmp_ui(nsp, LARGEST_PRIME) <= 0 ) {

                uint32_t nsp_ui = mpz_get_ui(nsp);
                //
                // Get the index of (cand_yi->data[k] / gcd) in
                // first_primes_array... If ind != NOT_IN_ARRAY, then
                // (cand_yi->data[k] / gcd) is indeed a prime...
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
                        mpz_init(smooth_yi->data[smooth_yi->length]);
                        mpz_init(xi->data[xi->length]);

                        mpz_init_set(
                            a_for_smooth_gx->data[a_for_smooth_gx->length],
                            cand_a->data[k]
                        );

                        mpz_mul(smooth_yi->data[smooth_yi->length],
                                cand_yi->data[k], found->y);

                        mpz_mul(xi->data[xi->length],
                                cand_xi->data[k], found->x);

                        mpz_divexact_ui(smooth_yi->data[smooth_yi->length],
                                        smooth_yi->data[smooth_yi->length],
                                        nsp_ui);
                        mpz_divexact_ui(smooth_yi->data[smooth_yi->length],
                                        smooth_yi->data[smooth_yi->length],
                                        nsp_ui);

                        if (0 != mpz_invert(nsp_inv, nsp, n)) {
                            mpz_mul(xi->data[xi->length], xi->data[xi->length],
                                    nsp_inv);
                        } else {
                            //
                            // This should not happen: since nsp is a prime,
                            // it has an inverse in Z/nZ... unless n is a
                            // multiple of nsp! In such a (rare) case, we have
                            // found a factor of n, but handling it is quite
                            // cumbersome, so we just ignore it...
                            //
                            free(key);
                            //clear_mpz_pair(found);
                            //free(found);
                            continue;
                        }
                        smooth_yi->length++;
                        xi->length++;
                        a_for_smooth_gx->length++;

                        free(key);
                        //clear_mpz_pair(found);
                        //free(found);

                        if (smooth_yi->length == smooth_yi->alloced) {
                            retval = k + 1;
                            break;
                        }
                    } else {
                        //
                        // Add this pair in the hashtable with the position
                        // of the prime in first_primes_array as key...
                        //
                        mpz_pair_t* pair = malloc(sizeof(mpz_pair_t));
                        mpz_init_set(pair->x, cand_xi->data[k]);
                        mpz_init_set(pair->y, cand_yi->data[k]);
                        mpz_mul(pair->y, pair->y, cand_a->data[k]);
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
    mpz_clear(yk);
    mpz_clear(gcd);
    mpz_clear(nsp);
    mpz_clear(nsp_inv);

    clear_mpz_tree(rtree);
    free(rtree);

    return retval;
}
//------------------------------------------------------------------------------
inline uint32_t djb_batch_rt(smooth_filter_t* const filter,
                             unsigned long int step) {
    //
    // Combining smoothness batches with early abort-like strategies
    // certainly makes no sense. Right now it is not used (filter->nsteps
    // will always be zero). The code's here though should one wants to play 
    // around with it...
    //
    if (filter->nsteps == 0) {
        return djb_batch_rt_no_ea(filter);
    }
    if (step == 0) {
        return djb_batch_rt_first(filter);
    }
    if (step == filter->nsteps) {
        return djb_batch_rt_last(filter);
    }
    return djb_batch_rt_step(filter, step - 1);    
}
//------------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                  "PRIVATE" FUNCTION(S) IMPLEMENTATION
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
uint32_t djb_batch_rt_no_ea(smooth_filter_t* const filter) {
        
    mpz_array_t* const in_xi  = filter->candidate_xi;
    mpz_array_t* const in_yi  = filter->candidate_yi;
    mpz_array_t* const acc_xi = filter->accepted_xi;
    mpz_array_t* const acc_yi = filter->accepted_yi;

    uint32_t next;
        
    hashtable_t* const htable = filter->htable;
    mpz_ptr n = filter->n;
    
    if (acc_yi->length == acc_yi->alloced) {
        return 0;
    }
    uint32_t retval = in_xi->length;
    //
    // Step 2
    //
    mpz_tree_t* const pyitree = prod_tree(in_yi);
    mpz_tree_t* const rtree   = rem_tree(filter->prod_pj[0], pyitree);
    
    clear_mpz_tree(pyitree);
    free(pyitree);
    //
    // Step 3
    //
    uint32_t e = 0U;
    uint32_t msb = 0U;

    mpz_t yk;
    mpz_init(yk);

    mpz_t gcd;
    mpz_init(gcd);

    mpz_t nsp;
    mpz_init(nsp);

    mpz_t nsp_inv;
    mpz_init(nsp_inv);

    for (uint32_t k = 0U; k < in_yi->length; k++) {
        msb = mpz_sizeinbase(in_yi->data[k], 2) - 1;
        e   = most_significant_bit(msb) + 1;
            
        mpz_set(yk, rtree->data[rtree->length/2 + k]);
        mpz_powm_ui(yk, yk, 1<<e, in_yi->data[k]);
        
        next = acc_yi->length;    
        //
        // _NOTE_: in_yi->data[k] is smooth iff (yk == 0) (no need to compute
        //         gcd if we just need to know whether in_yi->data[k] is
        //         smooth or not).
        //
        if (0 == mpz_sgn(yk)) {
            mpz_init_set(acc_yi->data[next], in_yi->data[k]);
            mpz_init_set(acc_xi->data[next], in_xi->data[k]);
        
            acc_xi->length++;
            acc_yi->length++;
                    
            if (acc_yi->length == acc_yi->alloced) {
                retval = k + 1;
                break;
            }
        } else {            
            //
            // Here comes the large prime variation...
            //
            mpz_gcd(gcd, yk, in_yi->data[k]);
    	    //
    	    // gcd is now the smooth part of in_yi->data[k].
    	    // Is in_yi->data[k] the product of gcd by a large prime?
    	    //
            mpz_divexact(nsp, in_yi->data[k], gcd);
            
            if (mpz_cmpabs_ui(nsp, LARGEST_PRIME) <= 0 ) {
            
                uint32_t nsp_ui = mpz_get_ui(nsp);
                //
                // Get the index of (in_yi->data[k] / gcd) in
                // first_primes_array... If ind != NOT_IN_ARRAY, then
                // (in_yi->data[k] / gcd) is indeed a prime...
                //
                uint32_t ind = index_in_sorted_uint32_array(
                                    nsp_ui,
                                    &first_primes_array,
                                    0,
                                    first_primes_array.length
                                );
                
                if (ind != NOT_IN_ARRAY) {
                    //
                    // in_yi->data[k] is the product of a smooth number with
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
            
                        if (0 != mpz_invert(nsp_inv, nsp, n)) {
                            mpz_mul(
                                acc_xi->data[next], 
                                acc_xi->data[next],
                                nsp_inv
                            );
                        
                        } else {
                            //
                            // This should not happen: since nsp is a prime,
                            // it has an inverse in Z/nZ... unless n is a
                            // multiple of nsp! In such a (rare) case, we have
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
                            retval = k + 1;
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
            }
        }
    }
    mpz_clear(yk);
    mpz_clear(gcd);
    mpz_clear(nsp);
    mpz_clear(nsp_inv);

    clear_mpz_tree(rtree);
    free(rtree);

    return retval;   
}
//------------------------------------------------------------------------------
uint32_t djb_batch_rt_step(smooth_filter_t* const filter,
                            unsigned long int step) {

    mpz_array_t* const in_xi   = filter->filtered_xi[step];
    mpz_array_t* const in_yi   = filter->filtered_yi[step];
    mpz_array_t* const in_cof  = filter->cofactors  [step];
    mpz_array_t* const out_xi  = filter->filtered_xi[step + 1];
    mpz_array_t* const out_yi  = filter->filtered_yi[step + 1];
    mpz_array_t* const out_cof = filter->cofactors  [step + 1];
    mpz_array_t* const acc_xi  = filter->accepted_xi;
    mpz_array_t* const acc_yi  = filter->accepted_yi;
    
    uint32_t next_acc;
    uint32_t next_out;
    
    hashtable_t* htable = filter->htable;
    mpz_ptr n           = filter->n;
    mpz_ptr bound       = filter->bounds[step + 1];
    
    if (out_yi->length == out_yi->alloced) {
        return 0;
    }
    uint32_t retval = in_xi->length;

    //
    // Step 2
    //
    mpz_tree_t* pyitree = prod_tree(in_yi);
    mpz_tree_t* rtree   = rem_tree(filter->prod_pj[step + 1], pyitree);

    clear_mpz_tree(pyitree);
    free(pyitree);
    //
    // Step 3
    //
    uint32_t e = 0U;
    uint32_t msb = 0U;

    mpz_t yk;
    mpz_init(yk);

    mpz_t gcd;
    mpz_init(gcd);

    mpz_t nsp;
    mpz_init(nsp);

    mpz_t nsp_inv;
    mpz_init(nsp_inv);

    for (uint32_t k = 0U; k < in_yi->length; k++) {
        msb = mpz_sizeinbase(in_yi->data[k], 2) - 1;
        e   = most_significant_bit(msb) + 1;

        mpz_set(yk, rtree->data[rtree->length/2 + k]);
        mpz_powm_ui(yk, yk, 1<<e, in_yi->data[k]);
        
        next_acc = acc_yi->length;
        next_out = out_yi->length;   
        //
        // _NOTE_: in_y->data[k] is smooth iff (yk == 0) (no need to compute
        //         gcd if we just need to know whether in_y->data[k] is
        //         smooth or not).
        //
        if (0 == mpz_sgn(yk)) {

            mpz_init_set(acc_yi->data[next_acc], in_yi->data[k]);
            mpz_init_set(acc_xi->data[next_acc], in_xi->data[k]);
        
            mpz_mul(
                acc_yi->data[next_acc],
                acc_yi->data[next_acc],
                in_cof->data[k]
            );
        
            acc_xi->length++;
            acc_yi->length++;
        
            if (acc_yi->length == acc_yi->alloced) {
                retval = k + 1;
                break;
            }
        } else {
            //
            // Here comes the large prime variation...
            //
            mpz_gcd(gcd, yk, in_yi->data[k]);
    	    //
    	    // gcd is now the smooth part of in_y->data[k].
    	    // Is in_y->data[k] the product of gcd by a large prime?
    	    //
            mpz_divexact(nsp, in_yi->data[k], gcd);
            
            if (mpz_cmpabs_ui(nsp, LARGEST_PRIME) <= 0) {
            
                uint32_t nsp_ui = mpz_get_ui(nsp);
                //
                // Get the index of (in_yi->data[k] / gcd) in
                // first_primes_array... If ind != NOT_IN_ARRAY, then
                // (in_yi->data[k] / gcd) is indeed a prime...
                //
                uint32_t ind = index_in_sorted_uint32_array(
                                    nsp_ui,
                                    &first_primes_array,
                                    0,
                                    first_primes_array.length);
                
                if (ind != NOT_IN_ARRAY) {
                    //
                    // in_y->data[k] is the product of a smooth number with
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
                        // pair to add in candidate arrays.
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
            
                        if (0 != mpz_invert(nsp_inv, nsp, n)) {
                            mpz_mul(
                                acc_xi->data[next_acc], 
                                acc_xi->data[next_acc],
                                nsp_inv
                            );
                        } else {
                            //
                            // This should not happen: since nsp is a prime,
                            // it has an inverse in Z/nZ... unless n is a
                            // multiple of nsp! In such a (rare) case, we have
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
                            retval = k + 1;
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
                // Here comes the early abort strategy that makes no sense here!
                //
                // The non smooth part of in_y->data[k] is not the product of
                // a smooth number by a large prime. If the non smooth part
                // is greater than a certain bound, just discard it. Otherwise
                // it qualifies for the next round of smoothness batch!
                //
                if (mpz_cmpabs(nsp, bound) > 0) {
                    continue;
                }
                mpz_set(out_yi->data[next_out], nsp);
                mpz_set(out_xi->data[next_out], in_xi->data[k]);
                mpz_mul(out_cof->data[next_out], gcd, in_cof->data[k]);

                out_yi->length++;
                out_xi->length++;
                out_cof->length++;

                if (out_yi->length == out_yi->alloced) {
                    retval = k + 1;
                    break;
                }
            }
        }
    }
    mpz_clear(yk);
    mpz_clear(gcd);
    mpz_clear(nsp);
    mpz_clear(nsp_inv);

    clear_mpz_tree(rtree);
    free(rtree);

    return retval;
}
//------------------------------------------------------------------------------
uint32_t djb_batch_rt_first(smooth_filter_t* const filter) {
    
    mpz_array_t* const in_xi  = filter->candidate_xi;
    mpz_array_t* const in_yi  = filter->candidate_yi;
    mpz_array_t* const out_xi = filter->filtered_xi[0];
    mpz_array_t* const out_yi = filter->filtered_yi[0];
    mpz_array_t* const cofact = filter->cofactors[0];
    mpz_array_t* const acc_xi = filter->accepted_xi;
    mpz_array_t* const acc_yi = filter->accepted_yi;

    uint32_t next_acc;
    uint32_t next_out;
            
    hashtable_t* htable = filter->htable;
    mpz_ptr n           = filter->n;
    mpz_ptr bound       = filter->bounds[0];
    
    if (out_yi->length == out_yi->alloced) {
        return 0;
    }
    uint32_t retval = in_xi->length;

    //
    // Step 2
    //
    mpz_tree_t* pyitree = prod_tree(in_yi);
    mpz_tree_t* rtree   = rem_tree(filter->prod_pj[0], pyitree);

    clear_mpz_tree(pyitree);
    free(pyitree);
    //
    // Step 3
    //
    uint32_t e = 0U;
    uint32_t msb = 0U;

    mpz_t yk;
    mpz_init(yk);

    mpz_t gcd;
    mpz_init(gcd);

    mpz_t nsp;
    mpz_init(nsp);

    mpz_t nsp_inv;
    mpz_init(nsp_inv);

    for (uint32_t k = 0U; k < in_yi->length; k++) {
        msb = mpz_sizeinbase(in_yi->data[k], 2) - 1;
        e   = most_significant_bit(msb) + 1;

        mpz_set(yk, rtree->data[rtree->length/2 + k]);
        mpz_powm_ui(yk, yk, 1<<e, in_yi->data[k]);
        
        next_acc = acc_yi->length;
        next_out = out_yi->length;
        //
        // _NOTE_: in_yi->data[k] is smooth iff (yk == 0) (no need to compute
        //         gcd if we just need to know whether in_yi->data[k] is
        //         smooth or not).
        //
        if (0 == mpz_sgn(yk)) {

            mpz_init_set(acc_yi->data[next_acc], in_yi->data[k]);
            mpz_init_set(acc_xi->data[next_acc], in_xi->data[k]);
        
            acc_xi->length++;
            acc_yi->length++;
                    
            if (acc_yi->length == acc_yi->alloced) {
                retval = k + 1;
                break;
            }
        } else {
            //
            // Here comes the large prime variation...
            //
            mpz_gcd(gcd, yk, in_yi->data[k]);
    	    //
    	    // gcd is now the smooth part of in_yi->data[k].
    	    // Is in_yi->data[k] the product of gcd by a large prime?
    	    //
            mpz_divexact(nsp, in_yi->data[k], gcd);
            
            if (mpz_cmpabs_ui(nsp, LARGEST_PRIME) <= 0 ) {
            
                uint32_t nsp_ui = mpz_get_ui(nsp);
                //
                // Get the index of (in_yi->data[k] / gcd) in
                // first_primes_array... If ind != NOT_IN_ARRAY, then
                // (in_yi->data[k] / gcd) is indeed a prime...
                //
                uint32_t ind = index_in_sorted_uint32_array(
                                    nsp_ui,
                                    &first_primes_array,
                                    0,
                                    first_primes_array.length
                                );
                
                if (ind != NOT_IN_ARRAY) {
                    //
                    // in_yi->data[k] is the product of a smooth number with
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
            
                        if (0 != mpz_invert(nsp_inv, nsp, n)) {
                            mpz_mul(
                                acc_xi->data[next_acc], 
                                acc_xi->data[next_acc],
                                nsp_inv
                            );
                        } else {
                            //
                            // This should not happen: since nsp is a prime,
                            // it has an inverse in Z/nZ... unless n is a
                            // multiple of nsp! In such a (rare) case, we have
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
                            retval = k + 1;
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
                // in_y->data[k] is not the product of a smooth number by a
                // large prime. If its non smooth part is greater than a certain
                // bound, just discard it. Otherwise it qualifies for the next
                // round of smoothness batch! Congratulations!
                //
                if (mpz_cmpabs(nsp, bound) > 0) {
                    continue;
                }
                mpz_set(out_yi->data[next_out], nsp);
                mpz_set(out_xi->data[next_out], in_xi->data[k]);
                mpz_set(cofact->data[next_out], gcd);
                
                out_yi->length++;
                out_xi->length++;
                cofact->length++;
                
                if (out_yi->length == out_yi->alloced) {
                    retval = k + 1;
                    break;
                }
            }
        }
    }
    mpz_clear(yk);
    mpz_clear(gcd);
    mpz_clear(nsp);
    mpz_clear(nsp_inv);

    clear_mpz_tree(rtree);
    free(rtree);

    return retval;   
}
//------------------------------------------------------------------------------
uint32_t djb_batch_rt_last(smooth_filter_t* const filter) {

    unsigned long int step = filter->nsteps - 1;

    mpz_array_t* const in_xi  = filter->filtered_xi[step];
    mpz_array_t* const in_yi  = filter->filtered_yi[step];
    mpz_array_t* const cofact = filter->cofactors[step];
    mpz_array_t* const acc_xi = filter->accepted_xi;
    mpz_array_t* const acc_yi = filter->accepted_yi;

    uint32_t next;
    
    hashtable_t* htable = filter->htable;
    mpz_ptr n           = filter->n;
    
    if (acc_yi->length == acc_yi->alloced) {
        return 0;
    }
    uint32_t retval = in_xi->length;

    //
    // Step 2
    //
    mpz_tree_t* pyitree = prod_tree(in_yi);
    mpz_tree_t* rtree   = rem_tree(filter->prod_pj[step + 1], pyitree);

    clear_mpz_tree(pyitree);
    free(pyitree);
    //
    // Step 3
    //
    uint32_t e = 0U;
    uint32_t msb = 0U;

    mpz_t yk;
    mpz_init(yk);

    mpz_t gcd;
    mpz_init(gcd);

    mpz_t nsp;
    mpz_init(nsp);

    mpz_t nsp_inv;
    mpz_init(nsp_inv);

    for (uint32_t k = 0U; k < in_yi->length; k++) {
        msb = mpz_sizeinbase(in_yi->data[k], 2) - 1;
        e   = most_significant_bit(msb) + 1;

        mpz_set(yk, rtree->data[rtree->length/2 + k]);
        mpz_powm_ui(yk, yk, 1<<e, in_yi->data[k]);
        
        next = acc_yi->length;
        //
        // _NOTE_: in_y->data[k] is smooth iff (yk == 0) (no need to compute
        //         gcd if we just need to know whether in_y->data[k] is
        //         smooth or not).
        //
        if (0 == mpz_sgn(yk)) {
            
            mpz_init_set(acc_yi->data[next], in_yi->data[k]);
            mpz_init_set(acc_xi->data[next], in_xi->data[k]);
        
            mpz_mul(acc_yi->data[next], acc_yi->data[next], cofact->data[k]);
        
            acc_xi->length++;
            acc_yi->length++;
        
            if (acc_yi->length == acc_yi->alloced) {
                retval = k + 1;
                break;
            }
        } else {
            //
            // Here comes the large prime variation...
            //
            mpz_gcd(gcd, yk, in_yi->data[k]);
    	    //
    	    // gcd is now the smooth part of in_y->data[k].
    	    // Is in_y->data[k] the product of gcd by a large prime?
    	    //
            mpz_divexact(nsp, in_yi->data[k], gcd);
            
            if (mpz_cmpabs_ui(nsp, LARGEST_PRIME) <= 0 ) {
            
                uint32_t nsp_ui = mpz_get_ui(nsp);
                //
                // Get the index of (in_y->data[k] / gcd) in
                // first_primes_array... If ind != NOT_IN_ARRAY, then
                // (in_y->data[k] / gcd) is indeed a prime...
                //
                uint32_t ind = index_in_sorted_uint32_array(
                                   nsp_ui,
                                   &first_primes_array,
                                   0,
                                   first_primes_array.length
                               );
                
                if (ind != NOT_IN_ARRAY) {
                    //
                    // in_y->data[k] is the product of a smooth number with
                    // the prime number first_primes_array->data[ind]
                    //
                    uint32_t *key = malloc(sizeof(uint32_t));
                    *key = ind;
            
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
                        // pair to add in our accepted arrays.
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
            
                        if (0 != mpz_invert(nsp_inv, nsp, n)) {
                            mpz_mul(
                                acc_xi->data[next], 
                                acc_xi->data[next],
                                nsp_inv
                            );
                        } else {
                            //
                            // This should not happen: since nsp is a prime,
                            // it has an inverse in Z/nZ... unless n is a
                            // multiple of nsp! In such a (rare) case, we have
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
                            retval = k + 1;
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
    mpz_clear(yk);
    mpz_clear(gcd);
    mpz_clear(nsp);
    mpz_clear(nsp_inv);

    clear_mpz_tree(rtree);
    free(rtree);

    return retval;
}
//------------------------------------------------------------------------------
