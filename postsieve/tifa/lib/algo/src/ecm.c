//
// Copyright (C) 2006, 2007, 2008 INRIA (French National Institute for Research
// in Computer Science and Control)
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
 * \file    ecm.c
 * \author  Jerome Milan
 * \date    Thu Feb 7 2008
 * \version 1.0
 */

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <time.h>

#include "tifa_config.h"
#include "first_primes.h"
#include "array.h"
#include "funcs.h"
#include "x_tree.h"
#include "macros.h"
#include "factoring_machine.h"
#include "ecm.h"
#include "tifa_factor.h"

#define __PREFIX__   "ecm: "
#define __VERBOSE__  TIFA_VERBOSE_ECM
#define __TIMING__   TIFA_TIMING_ECM

#include "messages.h"

//------------------------------------------------------------------------------
//                   NON PUBLIC STRUCTURE(S) AND TYPEDEF(S)
//------------------------------------------------------------------------------

//
// An ad-hoc structure holding the data variables used by the ECM
// implementation.
//
struct struct_ecm_context_t {
    //
    // The number to factor
    //
    mpz_t n;
    //
    // k <- Prod_{pi<=B1} pi^{log(B1)/log(pi)}
    //
    mpz_t k;
    //
    // (x0 : z0) is the starting point on the curve using the Montgomery
    // representation.
    // 
    mpz_t x0;
    mpz_t z0;
    //
    // (xk : zk) is the point obtained at the en d of the first phase.
    // (xk : zk) = [k] (x0 : z0) (uses the Montgomery representation).
    //
    mpz_t xk;
    mpz_t zk;
    //
    // Parameters of a Suyama curve: b.z.y^2 = x^3 + a.z.x^2 + x.z^2
    //
    mpz_t a;
    mpz_t b;
    //
    // Pointer to the found factors. Could be accessed differently, but
    // it's convenient to have a copy of the pointer here...
    //
    mpz_array_t* factors;
    //
    // Suyama's curve parameter
    //
    unsigned long int sigma;
    //
    // Index of first "large prime" i.e. smallest prime greater than B1
    //
    uint32_t  first_lp_index;
    //
    // Temporary variable used in computations...
    //
    mpz_t tmp_add_1;
    mpz_t tmp_add_2;
    mpz_t tmp_add_3;
    mpz_t tmp_add_4;
    mpz_t tmp_add_5;
    mpz_t tmp_add_6;
    mpz_t tmp_mul_1;
    mpz_t tmp_mul_2;
    mpz_t tmp_mul_3;
    mpz_t tmp_mul_4;
    //
    // A pointer to the ECM parameters. Could be accessed differently, but
    // it's convenient to have a copy of the pointer here...
    //
    ecm_params_t*  params;
    //
    // true if the second phase precomputations have been performed.
    //
    bool precomputations_done;
    //
    // Bounds used in phase 1 and 2.
    //
    unsigned long int B1;
    unsigned long int B2;
    //
    // Used in precomputations. See reference article for details. 
    //
    unsigned long int D;
    unsigned long int DO2;
    unsigned long int m_min;
    unsigned long int m_max;
    unsigned int*  Js;
    unsigned int   njs;
    unsigned int*  gcd_table;
    unsigned int** prime_table;
    mpz_t* SX;
    mpz_t* SZ;
    //
    // How many times did we update the context to try more curves (with
    // the same bounds)?
    //
    unsigned int nupdated_more_curves;
    //
    // How many times did we update the context to try larger bounds?
    //
    unsigned int nupdated_larger_bounds;
};
//------------------------------------------------------------------------------
typedef struct struct_ecm_context_t ecm_context_t;
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                      CONSTANT DATA, SYMBOLS AND MACROS
//------------------------------------------------------------------------------
//
// Number of entries in default parameters tables.
//
#define NDFLT_ENTRIES 31
//
// If default_bitsize[i] <= size(n_to_factor) < default_bitsize[i+1], uses
// parameters' values given by default_<param>[i]. 
//
// _WARNING_: These values were obtained from benchmarks for composites
//            numbers p1 * p2 in the range 50 bits - 110 bits, where parameters
//            were choosen to maximize the ratio (time / probabiblity of
//            success). Values from the benchmark were then fitted and 
//            extrapolated to cover the whole 50-200 bits range.
//
//            We caution that this best parameters determination is _very_
//            crude! Needless to say, such extrapolations are known for being
//            "creatively accurate"...
//
// _NOTE_: In order to restrict the range of values combination to bench, the
//         "best" values given in the paper "A Practical Analysis of the 
//         Elliptic Curve Factoring Algorithm" were used as a starting point.
//
// _SEE_: "A Practical Analysis of the Elliptic Curve Factoring Algorithm",
//        R. D. Silverman and S. S. Wagstaff Jr., Mathematics of Computation,
//        Vol. 61, No. 203, (July 1993), pp. 445-462.
//
static const unsigned int default_bitsizes[NDFLT_ENTRIES] = {
    50,
    55,
    60,
    65,
    70,
    75,
    80,
    85,
    90,
    95,
    100,
    105,
    110,
    115,
    120,
    125,
    130,
    135,
    140,
    145,
    150,
    155,
    160,
    165,
    170,
    175,
    180,
    185,
    190,
    200
};

//
// Default values for bound in first phase.
//
static const unsigned long int default_b1s[NDFLT_ENTRIES] = {
    178,
    201,
    233,
    280,
    346,
    441,
    576,
    769,
    1046,
    1441,
    2005,
    2811,
    3963,
    5609,  
    7960,
    11321,
    16121,
    22981,
    32783,
    46787,
    66799,
    95389,
    136243,
    194613,
    278017,
    397187,
    567461,
    810753,
    1158377,
    1655073,
    2364770
};

//
// Default values for ratio (bound_second_phase / bound_first_phase).
//
static const unsigned int default_ratios[NDFLT_ENTRIES] = {
    45,
    44,
    43,
    42,
    41,
    41,
    40,
    39,
    39,
    38,
    37,
    37,
    36,
    35,
    35,
    34,
    34,
    33,
    33,
    32,
    31,
    31,
    30,
    30,
    29,
    29,
    29,
    28,
    28,
    27,
    27
};

//
// Default values for number of curves to use before giving up (if SINGLE_RUN
// factoring mode is used).
//
static const unsigned int default_ncurves[NDFLT_ENTRIES] = {
    6,
    9,
    12,
    18,
    22,
    30,
    36,
    45,
    54,
    65,
    78,
    92,
    110,
    128,
    150,
    175,
    204,
    237,
    275,
    318,
    368,
    425,
    491,
    566,
    653,
    752,
    866,
    997,
    1146,
    1318,
    1516
};
//
// Maximum number of times to try more curves to find a factor before definitely
// giving up and starting again with greater bounds (only if the factorization
// mode is different than SINGLE_RUN).
//
#define MAX_NTRY_MORE_CURVES 4

//
// Multiply the b1 and b2 bounds by BOUND_MULTIPLIER if the previous
// MAX_NTRY_MORE_CURVES tentatives failed to find a factor (only if the 
// factorization mode is different than SINGLE_RUN).
//
#define BOUND_MULTIPLIER 2

//
// Maximum number of times to try more curves with larger bounds (only if the 
// factorization mode is different than SINGLE_RUN).
//
#define MAX_NTRY_LARGER_BOUNDS 4

//------------------------------------------------------------------------------
//                    PROTOTYPES OF NON PUBLIC FUNCTION(S)
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
static ecode_t init_ecm_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t clear_ecm_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t update_ecm_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t switch_to_larger_bounds(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t perform_ecm(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t phase_1(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t phase_2(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t choose_curve(ecm_context_t* const context);
//------------------------------------------------------------------------------
static ecode_t perform_precomputations(ecm_context_t* const context);
//------------------------------------------------------------------------------
static ecode_t recurse(
    mpz_array_t* const,
    uint32_array_t* const,
    const mpz_t,
    factoring_mode_t
);
//------------------------------------------------------------------------------
static void multiply_point(
    mpz_t xk, mpz_t zk,  const mpz_t x0, const mpz_t z0,
    const mpz_t k, ecm_context_t* const context
);
//------------------------------------------------------------------------------
static void double_point(
    mpz_t x2a, mpz_t z2a, const mpz_t xa, const mpz_t za,
    ecm_context_t* const context
);
//------------------------------------------------------------------------------
static void add_points(
    mpz_t xc, mpz_t zc,
    const mpz_t xa, const mpz_t za,
    const mpz_t xb, const mpz_t zb,
    const mpz_t x0, const mpz_t z0,
    ecm_context_t* const context
);
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                 "PUBLIC" FUNCTION(S) IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
ecode_t ecm(mpz_array_t* const factors, uint32_array_t* const multis,
            const mpz_t n, const ecm_params_t* const params,
            const factoring_mode_t mode) {

    PRINT_INTRO_MSG("elliptic curve method");
    INIT_TIMER;
    START_TIMER;

    ecode_t ecode = NO_FACTOR_FOUND;

    factoring_machine_t machine;

    mpz_init_set(machine.n, n);
    machine.mode                = mode;
    machine.params              = (void*) params;
    machine.init_context_func   = init_ecm_context;
    machine.perform_algo_func   = perform_ecm;
    machine.update_context_func = update_ecm_context;
    machine.clear_context_func  = clear_ecm_context;
    machine.recurse_func        = recurse;
    machine.factors             = factors;
    machine.multis              = multis;

    ecode = run_machine(&machine);

    mpz_clear(machine.n);

    STOP_TIMER;
    PRINT_STATUS(machine.success, ecode);

    return ecode;
}
//------------------------------------------------------------------------------
void set_ecm_params_to_default(const mpz_t n __attribute__ ((unused)),
                               ecm_params_t* const params) {
    //
    // _WARNING_: See above note about default parameters. In particular they
    //            were obtained from benchmarking factorization of RSA moduli.
    //
    // _TO_DO_: Modify this function to pass as argument the size of the
    //          factor we expect to find instead of the number to split.
    //
    unsigned int sizen = mpz_sizeinbase(n, 2);
    unsigned int i = 0;
    
    if (sizen < default_bitsizes[0]) {
        i = 0;
    } else {
        i = NDFLT_ENTRIES - 1;
        while (sizen < default_bitsizes[i]) {
            i--;
        }
    }
    params->b1      = default_b1s[i];
    params->b2      = default_ratios[i] * params->b1;
    params->ncurves = default_ncurves[i];
    return;
}
//------------------------------------------------------------------------------
static ecode_t init_ecm_context(factoring_machine_t* const machine) {

    PRINT_INIT_MSG;
    INIT_TIMER;
    START_TIMER;

    srand(time(NULL));

    machine->context             = (void*) malloc(sizeof(ecm_context_t));
    ecm_context_t* const context = (ecm_context_t* const) machine->context;
    ecm_params_t*  const params  = (ecm_params_t* const)  machine->params;

    context->factors = machine->factors;
    context->params  = params;

    mpz_init(context->a);
    mpz_init(context->b);
    mpz_init(context->x0);
    mpz_init(context->z0);
    mpz_init(context->xk);
    mpz_init(context->zk);
    mpz_init(context->k);
    mpz_init_set(context->n, machine->n);

    mpz_init(context->tmp_add_1);
    mpz_init(context->tmp_add_2);
    mpz_init(context->tmp_add_3);
    mpz_init(context->tmp_add_4);
    mpz_init(context->tmp_add_5);
    mpz_init(context->tmp_add_6);
    
    mpz_init(context->tmp_mul_1);
    mpz_init(context->tmp_mul_2);
    mpz_init(context->tmp_mul_3);
    mpz_init(context->tmp_mul_4);
    
    context->sigma  = 5;

    context->precomputations_done = false;
    context->nupdated_more_curves = 0;
    
    //
    // Compute context->k = Prod_{pi<=B1} pi^{log(B1)/log(pi)} via a product
    // tree...
    //
    unsigned int ip    = 0;
    unsigned int power = 0;
    double       logb1 = log(params->b1);
    
    while (first_primes[ip] <= params->b1) {
        ip++;
    }
    context->first_lp_index = ip;

    uint32_array_t* ppowers = alloc_uint32_array(ip);

    for (uint32_t i = 0; i < ip; i++) {
        power = (unsigned int) (logb1 / log(first_primes[i]));
        ppowers->data[i] = pow(first_primes[i], power);
    }
    ppowers->length = ip;

    mpz_tree_t* ptree = prod_tree_ui(ppowers);

    mpz_set(context->k, ptree->data[0]);

    clear_mpz_tree(ptree);
    clear_uint32_array(ppowers);
    free(ptree);
    free(ppowers);

    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//------------------------------------------------------------------------------
static ecode_t clear_ecm_context(factoring_machine_t* const machine) {

    PRINT_CLEAN_MSG;
    INIT_TIMER;
    START_TIMER;

    ecm_context_t* context = (ecm_context_t*) machine->context;

    mpz_clear(context->n);
    mpz_clear(context->a);
    mpz_clear(context->b);
    mpz_clear(context->x0);
    mpz_clear(context->z0);
    mpz_clear(context->xk);
    mpz_clear(context->zk);
    mpz_clear(context->k);

    mpz_clear(context->tmp_add_1);
    mpz_clear(context->tmp_add_2);
    mpz_clear(context->tmp_add_3);
    mpz_clear(context->tmp_add_4);
    mpz_clear(context->tmp_add_5);
    mpz_clear(context->tmp_add_6);

    mpz_clear(context->tmp_mul_1);
    mpz_clear(context->tmp_mul_2);
    mpz_clear(context->tmp_mul_3);
    mpz_clear(context->tmp_mul_4);

    if (context->precomputations_done) {
        for (unsigned int ij = 0; ij < context->njs; ij++) {
            unsigned int j = context->Js[ij];
            mpz_clear(context->SX[j]);
            mpz_clear(context->SZ[j]);
        }
        free(context->SX);
        free(context->SZ);

        unsigned int nrows = context->m_max - context->m_min + 1;
        for (unsigned int i = 0U; i < nrows; i++) {
            free(context->prime_table[i]);
        }
        free(context->prime_table);
        free(context->Js);
        free(context->gcd_table);
    }
    free(context);

    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//------------------------------------------------------------------------------
static ecode_t update_ecm_context(factoring_machine_t* const machine
                                  __attribute__ ((unused))) {

    //
    // Strategy:
    //
    //     a) We try again with the same number of curves using the same
    //        bounds.
    //     b) If we have already tried MAX_NTRY_MORE_CURVES times with the
    //        same bounds, multiply the bounds by BOUND_MULTIPLIER and try
    //        again. If we have changed the bounds for larger ones
    //        MAX_NTRY_LARGER_BOUNDS times, then just admit defeat.
    //
    INIT_TIMER;
    START_TIMER;
    
    ecm_context_t* context = (ecm_context_t*) machine->context;
    ecode_t ecode = SUCCESS;
    
    if (context->nupdated_larger_bounds == MAX_NTRY_LARGER_BOUNDS) {
        PRINT_UPDATE_GIVEUP_MSG;
        ecode = GIVING_UP;

    } else {
        if (context->nupdated_more_curves == MAX_NTRY_MORE_CURVES) {
            //
            // Try with larger b1 and b2 bounds.
            //
            context->nupdated_larger_bounds++;
            context->nupdated_more_curves = 0;
            PRINT_UPDATE_LARGER_BOUNDS_MSG;
            
            ecode = switch_to_larger_bounds(machine); 
            
        } else {
            //
            // Try again with the same number of curves using the same bounds.
            //
            context->nupdated_more_curves++;
            PRINT_UPDATE_MORE_CURVES_MSG;
        }
    }
    STOP_TIMER;
    PRINT_TIMING;

    return ecode;
}
//------------------------------------------------------------------------------
static ecode_t switch_to_larger_bounds(factoring_machine_t* const machine) {
    //
    // Switch to larger bounds and recompute what needs to be recomputed...
    // Returns SUCCESS upon... success. Returns GIVING_UP if the bounds are 
    // already too large, thus aborting the factorization.
    //
    ecm_context_t* context = (ecm_context_t*) machine->context;
    ecm_params_t*  params  = (ecm_params_t*)  machine->params;

    unsigned long int max_bound = TIFA_ULONG_MAX / BOUND_MULTIPLIER;

    if (params->b1 >= max_bound) {
        return GIVING_UP;
    }
    if (params->b2 >= max_bound) {
        return GIVING_UP;
    }

    params->b1      *= BOUND_MULTIPLIER;
    params->b2      *= BOUND_MULTIPLIER;
    params->ncurves *= BOUND_MULTIPLIER;

    if (context->precomputations_done) {
        //
        // _TO_DO_: Actually these clear_* are not needed here (except for
        //          prime_table). We can bypass this memory cleaning by
        //          adding a little bit of logic in perform_precomputations.
        //
        for (unsigned int ij = 0; ij < context->njs; ij++) {
            unsigned int j = context->Js[ij];
            mpz_clear(context->SX[j]);
            mpz_clear(context->SZ[j]);
        }
        free(context->SX);
        free(context->SZ);

        unsigned int nrows = context->m_max - context->m_min + 1;
        for (unsigned int i = 0U; i < nrows; i++) {
            free(context->prime_table[i]);
        }
        free(context->prime_table);
        free(context->Js);
        free(context->gcd_table);
    }
    context->precomputations_done = false;

    //
    // Compute context->k = Prod_{pi<=B1} pi^{log(B1)/log(pi)} via a product
    // tree...
    //
    unsigned int ip    = 0;
    unsigned int power = 0;
    double       logb1 = log(params->b1);
    
    while (first_primes[ip] <= params->b1) {
        ip++;
    }
    context->first_lp_index = ip;

    uint32_array_t* ppowers = alloc_uint32_array(ip);

    for (uint32_t i = 0; i < ip; i++) {
        power = (unsigned int) (logb1 / log(first_primes[i]));
        ppowers->data[i] = pow(first_primes[i], power);
    }
    ppowers->length = ip;

    mpz_tree_t* ptree = prod_tree_ui(ppowers);

    mpz_set(context->k, ptree->data[0]);

    clear_mpz_tree(ptree);
    clear_uint32_array(ppowers);
    free(ptree);
    free(ppowers);

    return SUCCESS;
}
//------------------------------------------------------------------------------
static ecode_t perform_ecm(factoring_machine_t* const machine) {

    INIT_TIMER;
    START_TIMER;

    INIT_NAMED_TIMER(phase_1);
    INIT_NAMED_TIMER(phase_2);
    INIT_NAMED_TIMER(choose);
    
    ecm_context_t* context = (ecm_context_t*) machine->context;
    ecm_params_t*  params  = (ecm_params_t*)  machine->params;

    ecode_t exit_code = NO_FACTOR_FOUND;

    uint32_t ncurves = 0;

    while (ncurves < params->ncurves) {

        ncurves++;
        //
        // Choose a random curve using Suyama's parametrization...
        //
        START_NAMED_TIMER(choose);
        exit_code = choose_curve(context);
        STOP_NAMED_TIMER(choose);

        if (exit_code == SOME_FACTORS_FOUND) {
            //
            // Extremely unlikely!
            //
            break;
        }
        //
        // Perform first phase of ECM...
        //
        START_NAMED_TIMER(phase_1);
        exit_code = phase_1(machine);
        STOP_NAMED_TIMER(phase_1);

        if (exit_code == SOME_FACTORS_FOUND) {
            break;
        }
        if (params->b2 == 0) {
            //
            // No second phase is performed if b2 is set to 0.
            //
            continue;
        }
        if (!context->precomputations_done) {
            perform_precomputations(context);
        }
        //
        // Perform second phase of ECM...
        //
        START_NAMED_TIMER(phase_2);
        exit_code = phase_2(machine);
        STOP_NAMED_TIMER(phase_2);

        if (exit_code == SOME_FACTORS_FOUND) {
            break;
        }
    }
    STOP_TIMER;
    
    PRINT_NCURVES_USED(ncurves);
    
    PRINT_CHOOSED_CURVES_IN;
    PRINT_NAMED_TIMING(choose);
    PRINT_PHASE_1_IN;
    PRINT_NAMED_TIMING(phase_1);
    PRINT_PHASE_2_IN;
    PRINT_NAMED_TIMING(phase_2);
        
    return exit_code;
}
//------------------------------------------------------------------------------
static ecode_t recurse(mpz_array_t* const factors, uint32_array_t* const multis,
                       const mpz_t n, factoring_mode_t mode) {

    return tifa_factor(factors, multis, n, mode);
}
//------------------------------------------------------------------------------
static ecode_t phase_1(factoring_machine_t* const machine) {

    ecm_context_t* const context = (ecm_context_t* const) machine->context;

    multiply_point(
        context->xk, context->zk, context->x0, context->z0, context->k, context
    );

    ecode_t exit_code = NO_FACTOR_FOUND;

    mpz_t g;
    mpz_init(g);

    mpz_gcd(g, context->zk, machine->n);

    if ((mpz_cmp_ui(g, 1) != 0) && (mpz_cmp(g, machine->n) != 0)) {
        append_mpz_to_array(machine->factors, g);
        mpz_divexact(g, machine->n, g);
        append_mpz_to_array(machine->factors, g);
        exit_code = SOME_FACTORS_FOUND;
    }
    mpz_clear(g);

    return exit_code;
}
//------------------------------------------------------------------------------
static ecode_t phase_2(factoring_machine_t* const machine) {
    //
    // This second phase implementation follows the standard continuation
    // and is implemented in a way reminiscent of the description
    // given in the article "Implementing the Elliptic Curve Method of Factoring
    // in Reconfigurable Hardware" by Kris Gaj et al.
    //
    ecm_context_t* const context = (ecm_context_t* const) machine->context;
    ecode_t exit_code = NO_FACTOR_FOUND;

    mpz_t g;
    mpz_t d;

    mpz_init_set_ui(d, 1);
    mpz_init(g);

    unsigned long int DO2   = context->DO2;
    unsigned long int D     = context->D;
    unsigned long int m_min = context->m_min;
    unsigned long int m_max = context->m_max;
    mpz_t*        SX        = context->SX;
    mpz_t*        SZ        = context->SZ;
    unsigned int* Js        = context->Js;
    unsigned int  njs       = context->njs;
    unsigned int** prime_table = context->prime_table;

    mpz_t xq0;
    mpz_t zq0;

    mpz_t xq;
    mpz_t zq;

    mpz_t xr;
    mpz_t zr;

    mpz_t xdiff;
    mpz_t zdiff;

    mpz_t xq_old;
    mpz_t zq_old;

    mpz_init(xq0);
    mpz_init(zq0);

    mpz_init(xq);
    mpz_init(zq);

    mpz_init(xr);
    mpz_init(zr);

    mpz_init(xq_old);
    mpz_init(zq_old);

    mpz_init(xdiff);
    mpz_init(zdiff);

    mpz_t tmp;
    mpz_init(tmp);

    mpz_set(xq0, context->xk);
    mpz_set(zq0, context->zk);

    mpz_set(xq, xq0);
    mpz_set(zq, zq0);

    if (gcd_ulint(1, D) == 1) {
        mpz_set(context->SX[1], xq);
        mpz_set(context->SZ[1], zq);
    }
    double_point(xq, zq, xq, zq, context); // Q = [2] Q0

    mpz_set(xdiff, xq0);
    mpz_set(zdiff, zq0);

    for (unsigned long int j = 2; j <= DO2; j++) {
        if (gcd_ulint(j, D) == 1) {
            mpz_set(context->SX[j], xq);
            mpz_set(context->SZ[j], zq);
        }
        mpz_set(xq_old, xq);
        mpz_set(zq_old, zq);

        add_points(xq, zq, xq, zq, xq0, zq0, xdiff, zdiff, context);

        mpz_set(xdiff, xq_old);
        mpz_set(zdiff, zq_old);
    }
    mpz_set(xq, xq0);
    mpz_set(zq, zq0);

    mpz_set_ui(tmp, D);
    multiply_point(xq, zq, xq0, zq0, tmp, context);     // Q = [D] Q0

    mpz_set_ui(tmp, m_min-1);
    multiply_point(xdiff, zdiff, xq, zq, tmp, context); // DIFF = [m_min-1] Q

    mpz_set_ui(tmp, m_min);
    multiply_point(xr, zr, xq, zq, tmp, context);       // R = [m_min] Q

    for (unsigned int im = 0; im <= (m_max - m_min); im++) {

        unsigned int m = m_min + im;

        for (unsigned int ij = 0; ij < njs; ij++) {
            unsigned int j = Js[ij];
            if (prime_table[im][j] == 1) {
                mpz_mul(tmp, xr, SZ[j]);
                mpz_submul(tmp, zr, SX[j]);
                mpz_mul(d, d, tmp);
                mpz_mod(d, d, context->n);
            }
        }
        mpz_set(xq_old, xr);
        mpz_set(zq_old, zr);

        if (m == 1) {
            double_point(xr, zr, xr, zr, context);
        } else {
            add_points(xr, zr, xr, zr, xq, zq, xdiff, zdiff, context);
        }
        mpz_set(xdiff, xq_old);
        mpz_set(zdiff, zq_old);
    }
    mpz_gcd(g, d, context->n);

    if ((mpz_cmp_ui(g, 1) != 0) && (mpz_cmp(g, machine->n) != 0)) {
        append_mpz_to_array(machine->factors, g);
        mpz_divexact(g, machine->n, g);
        append_mpz_to_array(machine->factors, g);
        exit_code = SOME_FACTORS_FOUND;
    }
    mpz_clear(d);
    mpz_clear(g);
    mpz_clear(xq0);
    mpz_clear(zq0);
    mpz_clear(xq);
    mpz_clear(zq);
    mpz_clear(xr);
    mpz_clear(zr);
    mpz_clear(tmp);

    mpz_clear(xdiff);
    mpz_clear(zdiff);
    mpz_clear(xq_old);
    mpz_clear(zq_old);

    return exit_code;
}
//------------------------------------------------------------------------------
static ecode_t choose_curve(ecm_context_t* const context) {
    //
    // Choose a curve and a starting point using Suyama's parametrization.
    //
    ecode_t exit_code = SUCCESS;

    context->sigma = (rand() & 16777215) + 6;
    
    mpz_ptr n = context->n;

    mpz_ptr v       = context->tmp_add_1;
    mpz_ptr u       = context->tmp_add_2;
    mpz_ptr u3v4    = context->tmp_add_3;
    mpz_ptr u3v4inv = context->tmp_add_4;
    mpz_ptr tmp     = context->tmp_add_5;
    mpz_ptr tmp2    = context->tmp_add_6;

    mpz_set_ui(u, context->sigma);
    mpz_mul(u, u, u);
    mpz_sub_ui(u, u, 5); // u = sigma^2 - 5

    mpz_set_ui(v, context->sigma);
    mpz_mul_ui(v, v, 4); // v = 4 * sigma

    mpz_mul(u3v4, u, u);
    mpz_mul(u3v4, u3v4, u);
    mpz_mul(u3v4, u3v4, v);
    mpz_mul_ui(u3v4, u3v4, 4); // u3v4 = 4 * u^3 * v

    mpz_gcd(tmp, u3v4, n);

    if (mpz_cmp_ui(tmp, 1) != 0) {
        //
        // We have a factor (very unlikely)!
        // Abort everything and save the found factor!
        //
        exit_code = SOME_FACTORS_FOUND;
        append_mpz_to_array(context->factors, tmp);
        mpz_divexact(tmp, n, tmp);
        append_mpz_to_array(context->factors, tmp);

        return exit_code;
    }

    mpz_invert(u3v4inv, u3v4, n); // inv(4 * u^3 * v) mod n

    mpz_sub(tmp, v, u);
    mpz_pow_ui(tmp, tmp, 3); // tmp = (v - u)^3

    mpz_mul_ui(tmp2, u, 3);
    mpz_add(tmp2, tmp2, v); // tmp2 = 3 * u + v

    mpz_mul(context->a, tmp, tmp2);
    mpz_mul(context->a, context->a, u3v4inv);
    mpz_sub_ui(context->a, context->a, 2);

    mpz_set_ui(context->z0, 1);

    mpz_pow_ui(context->x0, u, 4);
    mpz_mul_ui(context->x0, context->x0, 4);
    mpz_mul(context->x0, context->x0, u3v4inv);

    mpz_powm_ui(context->x0, context->x0, 3, n);

    mpz_pow_ui(v, v, 3); // v = v^3

    mpz_gcd(context->b, v, n);

    if (mpz_cmp_ui(context->b, 1) != 0) {
        //
        // We have a factor (very unlikely)!
        // Abort everything and save the found factor!
        //
        exit_code = SOME_FACTORS_FOUND;
        append_mpz_to_array(context->factors, context->b);
        mpz_divexact(context->b, n, tmp);
        append_mpz_to_array(context->factors, context->b);
        return exit_code;
    }
    mpz_invert(context->b, v, n);       // inv(v^3) mod n
    mpz_mul(context->b, context->b, u); // u * inv(v^3) mod n

    mpz_mod(context->a, context->a, n);
    mpz_mod(context->b, context->b, n);

    return exit_code;
}
//------------------------------------------------------------------------------
static void add_points(mpz_t xc, mpz_t zc,
                       const mpz_t xa, const mpz_t za,
                       const mpz_t xb, const mpz_t zb,
                       const mpz_t x0, const mpz_t z0,
                       ecm_context_t* const context) {
    //
    // Set (xc : zc) to (xa : za) + (xb : zb) knowing
    // (x0 : z0) = (xa : za) - (xb : zb).
    // Uses the Montgomery representation.
    //
    // _WARNING_: xc and x0 should point to different memory spaces!
    //
    mpz_ptr n  = context->n;

    mpz_ptr r = context->tmp_add_1;
    mpz_ptr s = context->tmp_add_2;
    mpz_ptr t = context->tmp_add_3;
    mpz_ptr u = context->tmp_add_4;
    mpz_ptr v = context->tmp_add_5;
    mpz_ptr w = context->tmp_add_6;

    mpz_sub(r, xa, za);
    mpz_sub(s, xb, zb);

    mpz_add(t, xa, za);
    mpz_add(u, xb, zb);

    mpz_mul(v, r, u);
    mpz_mul(w, s, t);
    
    mpz_add(xc, v, w);
    mpz_mul(xc, xc, xc);
    mpz_mul(xc, xc, z0);

    mpz_sub(zc, v, w);
    mpz_mul(zc, zc, zc);
    mpz_mul(zc, zc, x0);

    mpz_mod(xc, xc, n);
    mpz_mod(zc, zc, n);

    return;
}
//------------------------------------------------------------------------------
static void double_point(mpz_t x2a, mpz_t z2a,
                         const mpz_t xa, const mpz_t za,
                         ecm_context_t* const context) {
    //
    // Set (x2a : z2a) to [2] (xa : za).
    // Uses the Montgomery representation.
    //
    mpz_ptr n = context->n;
    mpz_ptr a = context->a;

    mpz_ptr plus  = context->tmp_add_1;
    mpz_ptr minus = context->tmp_add_2;
    mpz_ptr xz4   = context->tmp_add_3;
    mpz_ptr tmp   = context->tmp_add_4;

    mpz_add(plus, xa, za);
    mpz_sub(minus, xa, za);

    mpz_mul(plus, plus, plus);
    mpz_mul(minus, minus, minus);

    mpz_sub(xz4, plus, minus);
    mpz_mul(x2a, plus, minus);

    mpz_tdiv_q_2exp(plus, xz4, 2);
    mpz_add_ui(tmp, a, 2);
    mpz_mul(plus, plus, tmp);

    mpz_add(z2a, minus, plus);
    mpz_mul(z2a, z2a, xz4);

    mpz_mod(x2a, x2a, n);
    mpz_mod(z2a, z2a, n);

    return;
}
//------------------------------------------------------------------------------
static void multiply_point(mpz_t xk, mpz_t zk,
                           const mpz_t x0, const mpz_t z0,
                           const mpz_t k,
                           ecm_context_t* const context) {
    //
    // Set (xk : zk) to [k] (x0 : z0) using Montgomery's ladder algorithm.
    // Uses the Montgomery representation.
    //
    mpz_ptr qx = context->tmp_mul_1;
    mpz_ptr qz = context->tmp_mul_2;
    mpz_ptr px = context->tmp_mul_3;
    mpz_ptr pz = context->tmp_mul_4;

    mpz_set(qx, x0);
    mpz_set(qz, z0);

    double_point(px, pz, x0, z0, context);

    int s = mpz_sizeinbase(k, 2) - 2;

    while (s >= 0) {
        if (mpz_tstbit(k, s) == 1) {
            add_points(qx, qz, px, pz, qx, qz, x0, z0, context);
            double_point(px, pz, px, pz, context);
        } else {
            add_points(px, pz, px, pz, qx, qz, x0, z0, context);
            double_point(qx, qz, qx, qz, context);
        }
        s--;
    }
    mpz_set(xk, qx);
    mpz_set(zk, qz);

    return;
}
//------------------------------------------------------------------------------
static ecode_t perform_precomputations(ecm_context_t* const context) {
    //
    // Fills the tables needed for the second phase.
    //
    // These tables depends essentially on the bounds b1 and b2, so as soon
    // as the "optimal" parameters' values are determined, these tables can
    // be hardcoded.
    //
    // Ref: "Implementing the Elliptic Curve Method of Factoring in
    // Reconfigurable Hardware", K. Gaj et al., Cryptographic Hardware and
    // Embedded Systems - CHES 2006.
    //
    INIT_TIMER;
    PRINT_PERFORMING_P2_PRECOMP;
    START_TIMER;

    context->B2 = context->params->b2;
    context->B1 = context->params->b1;
    //
    // _WARNING_: This D value (representing a speed/memory trade-off) is from
    //            the aforementionned paper. Some tests are needed to assess
    //            the better value in our case.
    //
    context->D   = 210;
    context->DO2 = 105;

    context->m_min = (context->B1 + context->DO2) / context->D;
    context->m_max = (context->B2 - context->DO2) / context->D + 1;

    unsigned long int D   = context->D;
    unsigned long int DO2 = context->DO2;
    unsigned long int m_min = context->m_min;
    unsigned long int m_max = context->m_max;
    //
    // _NOTE_: More space than strictly required is allocated in order to keep
    //         a simple implementation without translating indices.
    //
    context->gcd_table = malloc((DO2 + 1) * sizeof(unsigned int));
    context->SX        = malloc((DO2 + 1) * sizeof(mpz_t));
    context->SZ        = malloc((DO2 + 1) * sizeof(mpz_t));
    context->Js        = malloc((DO2 + 1) * sizeof(unsigned int));
    context->njs       = 0;

    for (unsigned long int j = 1; j <= DO2; j++) {
        if (gcd_ulint(j, D) == 1) {
            context->gcd_table[j] = 1;
            mpz_init(context->SX[j]);
            mpz_init(context->SZ[j]);
            context->Js[context->njs] = j;
            context->njs++;
        } else {
            context->gcd_table[j] = 0;
        }
    }
    unsigned int nrows = m_max - m_min + 1;
    unsigned int ncols = DO2 + 1;

    context->prime_table = malloc(nrows * sizeof(unsigned int*));

    for (unsigned int i = 0U; i < nrows; i++) {
        context->prime_table[i] = malloc(ncols * sizeof(unsigned int));
    }
    unsigned int deltam = m_max - m_min;

    for (unsigned int im = 0; im <= deltam; im++) {
        unsigned int m = m_min + im;
        for (unsigned int j = 1; j <= DO2; j++) {
            if (is_prime(m * D + j) || is_prime(abs(m * D - j))) {
                context->prime_table[im][j] = 1;
                continue;
            }
            context->prime_table[im][j] = 0;
        }
    }
    context->precomputations_done = true;

    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//------------------------------------------------------------------------------
