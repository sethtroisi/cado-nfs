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
 * \file    fermat.c
 * \author  Jerome Milan
 * \date    Tue Aug 28 2007
 * \version 1.0
 */
 
 /*
  * History:
  *   1.1: Mon Sep 24 2007 by JM
  *        - Added mostly single precision implementation + multipliers +
  *          roll back to initial version if needed.
  *   1.0: Tue Aug 28 2007 by JM
  *        - Initial version.
  */
  
  //------------------------------------------------------------------------
  // References:
  //------------------------------------------------------------------------
  // - "Speeding Fermat's Factoring Method", James McKee,
  //   Mathematics of Computation,
  //   Volume 68, Number 228, pages 1729-1737.
  //------------------------------------------------------------------------

#include "tifa_config.h"

#include <stdlib.h>
#include <stdbool.h>
#include <inttypes.h>
#include <math.h>
#include <gmp.h>

#include "first_primes.h"
#include "fermat.h"
#include "funcs.h"
#include "macros.h"
#include "linked_list.h"

#define FERMAT_DEBUG 1
#if FERMAT_DEBUG
  #include <stdio.h>
  #define GMP_PRINT(...) gmp_printf(__VA_ARGS__)
  #define PRINT(...)     printf(__VA_ARGS__)
#else
  #define GMP_PRINT(...) /* intentionally left empty */
  #define PRINT(...)     /* intentionally left empty */
#endif

#define __PREFIX__  "fermat: "
#define __VERBOSE__ TIFA_VERBOSE_FERMAT
#define __TIMING__  TIFA_TIMING_FERMAT

#include "messages.h"

//-----------------------------------------------------------------------------
//                         NON PUBLIC DEFINE(S)
//-----------------------------------------------------------------------------
//
// Maximum number of primes to use. Abort the factorization if using MAX_NPRIMES
// primes is not enough to find a factor.
//
#define MAX_NPRIMES (1<<24)
//
// Multiplier used in McKee's "greedy" variant. If n is the number to factor,
// the greedy phase will stop as soon as yi >= MAX_Y_BOUND_MULTIPLIER * n^1/4.
//
#define MAX_Y_BOUND_MULTIPLIER 1
//
// Number of multipliers used.
//
#define NMULTIPLIERS 2

//
// Use a special way to "tag" local variables used by the SQRTM_P2_UI macro to
// avoid name clashes.
//
#define VSP(X) __ ## X ## _sqrtm_p2_ui_var__

//
// Declares and initializes local variables used by the SQRTM_P2_UI macro.
//
#define DECL_SQRTM_P2_UI_VARS() \
    mpz_t VSP(pz);              \
    mpz_t VSP(pz2);             \
    mpz_t VSP(inv);             \
    mpz_t VSP(np);              \
                                \
    mpz_init(VSP(pz));          \
    mpz_init(VSP(pz2));         \
    mpz_init(VSP(inv));         \
    mpz_init(VSP(np));

//
// Clears local variables used by the SQRTM_P2_UI macro.
//
#define CLEAR_SQRTM_P2_UI_VARS()    \
    mpz_clear(VSP(pz));             \
    mpz_clear(VSP(pz2));            \
    mpz_clear(VSP(inv));            \
    mpz_clear(VSP(np));

//
// Computes the modular squareroot of N mod P^2, where P is a prime. Stores
// the modular squareroot (if it exists) in S and stores a return code in ECODE.
//
// Parameters:
//  ECODE: unsigned long int
//      S: mpz_t
//      N: const mpz_t
//      P: unsigned long int
//
#define SQRTM_P2_UI(ECODE, S, N, P)                         \
    mpz_set_ui(VSP(pz), P);                                 \
    mpz_mul(VSP(pz2), VSP(pz), VSP(pz));                    \
                                                            \
    mpz_mod(VSP(np), N, VSP(pz2));                          \
                                                            \
    unsigned long int nmodp = mpz_tdiv_ui(VSP(np), P);      \
    unsigned long int snp   = sqrtm(nmodp, P);              \
                                                            \
    if (snp == UINT32_MAX) {                                \
        /*                                                  \
         * No solution                                      \
         */                                                 \
        ECODE = 1;                                          \
    } else {                                                \
        mpz_set_ui(S, snp);                                 \
        mpz_mul_ui(S, S, snp);                              \
        mpz_neg(S, S);                                      \
        mpz_add(S, S, VSP(np));                             \
                                                            \
        mpz_set_ui(VSP(inv), snp);                          \
        mpz_mul_2exp(VSP(inv), VSP(inv), 1);                \
                                                            \
        mpz_invert(VSP(inv), VSP(inv), VSP(pz));            \
                                                            \
        if (mpz_sgn(VSP(inv)) == 0) {                       \
            /*                                              \
             * Not inversible                               \
             */                                             \
            ECODE = 2;                                      \
        } else {                                            \
            mpz_divexact_ui(S, S, P);                       \
            mpz_mul(S, VSP(inv), S);                        \
            mpz_mod_ui(S, S, P);                            \
            mpz_mul_ui(S, S, P);                            \
            mpz_add_ui(S, S, snp);                          \
            ECODE = 0;                                      \
        }                                                   \
    }

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                   NON PUBLIC STRUCTURE(S) AND TYPEDEF(S)
//-----------------------------------------------------------------------------

//
// An ad-hoc structure holding the data variables used by Fermat's' algorithm.
//
struct struct_fermat_context_t {
    //
    // The number to factor.
    //
    mpz_ptr n;
    //
    // Whether the context has been updated or not
    //
    bool updated;
    //
    // Latest prime used in "simpler-precision" (sic!) version before failure.
    // Full multi-precision implementation will start using this prime.
    //
    unsigned long int latest_p; 
    //
    // Whether the number to factor is too large
    //
    bool integer_too_large;
};
//-----------------------------------------------------------------------------
typedef struct struct_fermat_context_t fermat_context_t;
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                           NON PUBLIC DATA
//-----------------------------------------------------------------------------

//
// List of multipliers used.
//
static unsigned int multipliers[NMULTIPLIERS] = {1, 3};

//-----------------------------------------------------------------------------
//                 PROTOTYPES OF NON PUBLIC FUNCTION(S)
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
static ecode_t init_fermat_context(factoring_machine_t* const machine);
//-----------------------------------------------------------------------------
static ecode_t clear_fermat_context(factoring_machine_t* const machine);
//-----------------------------------------------------------------------------
static ecode_t update_fermat_context(factoring_machine_t* const machine);
//-----------------------------------------------------------------------------
static ecode_t perform_fermat(factoring_machine_t* const machine);
//-----------------------------------------------------------------------------
static ecode_t perform_fermat_simple(factoring_machine_t* const machine);
//-----------------------------------------------------------------------------
static ecode_t perform_fermat_mpz(factoring_machine_t* const machine);
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                 "PUBLIC" FUNCTION(S) IMPLEMENTATION
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void set_fermat_params_to_default(fermat_params_t* const params
                                  __attribute__ ((unused))) {
    //
    // While there is no user parameters for this Fermat implementation, this
    // function is kept as a placeholder: it might be needed in future code
    // revisions.
    //
    return;
}
//-----------------------------------------------------------------------------
static ecode_t recurse(mpz_array_t* const factors, uint32_array_t* const multis,
                       const mpz_t n, factoring_mode_t mode) {

    fermat_params_t params;
    set_fermat_params_to_default(&params);

    return fermat(factors, multis, n, &params, mode);
}
//-----------------------------------------------------------------------------
ecode_t fermat(mpz_array_t* const factors, uint32_array_t* const multis,
               const mpz_t n, const fermat_params_t* const params,
               const factoring_mode_t mode) {

    PRINT_INTRO_MSG("Fermat");
    INIT_TIMER;
    START_TIMER;

    ecode_t ecode = NO_FACTOR_FOUND;

    factoring_machine_t machine;

    mpz_init_set(machine.n, n);

    machine.mode                = mode;
    machine.params              = (void*) params;
    machine.init_context_func   = init_fermat_context;
    machine.perform_algo_func   = perform_fermat;
    machine.update_context_func = update_fermat_context;
    machine.clear_context_func  = clear_fermat_context;
    machine.recurse_func        = recurse;
    machine.factors             = factors;
    machine.multis              = multis;

    ecode = run_machine(&machine);

    mpz_clear(machine.n);

    STOP_TIMER;
    PRINT_STATUS(machine.success, ecode);

    return ecode;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                  "PRIVATE" FUNCTION(S) IMPLEMENTATION
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                  Functions used by factoring machine
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
static ecode_t init_fermat_context(factoring_machine_t* const machine) {
    //
    // Init the Fermat's algorithm implementation specific variables
    //
    PRINT_INIT_MSG;
    INIT_TIMER;
    START_TIMER;

    machine->context          = (void*) malloc(sizeof(fermat_context_t));
    fermat_context_t* context = (fermat_context_t*) machine->context;

    context->n       = machine->n;
    context->updated = false;
    context->integer_too_large = false;

    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//-----------------------------------------------------------------------------
static ecode_t clear_fermat_context(factoring_machine_t* const machine) {

    PRINT_CLEAN_MSG;
    INIT_TIMER;
    START_TIMER;

    fermat_context_t* context = (fermat_context_t*) machine->context;
    free(context);

    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//-----------------------------------------------------------------------------
static ecode_t update_fermat_context(factoring_machine_t* const machine
                                     __attribute__ ((unused))) {

    INIT_TIMER;
    START_TIMER;
    
    fermat_context_t* context = (fermat_context_t*) machine->context;
    
    if (context->integer_too_large) {
        PRINT_UPDATE_GIVEUP_MSG;
        STOP_TIMER;
        PRINT_TIMING;
        return INTEGER_TOO_LARGE;
    }
    ecode_t ecode = SUCCESS;

    if (!context->updated) {
        //
        // One could also try the "basic" variant which, according to Pollard, 
        // is faster... on a 7 Mhz Psion handheld computer with 256K of 
        // memory...
        //
        // For the time being, we switch to a mostly multi-precision'ed
        // (sic) implementation...
        //
        PRINT_UPDATE_MULTIPREC_MSG;
        context->updated = true;
    } else {
        //
        // Actually this should not happen since we're sure we'll eventually
        // fund a factor...
        //
        PRINT_UPDATE_GIVEUP_MSG;
        ecode = GIVING_UP;
    }

    STOP_TIMER;
    PRINT_TIMING;

    return ecode;
}
//-----------------------------------------------------------------------------
static ecode_t perform_fermat(factoring_machine_t* const machine) {
    fermat_context_t* const context = (fermat_context_t*) machine->context;
    //
    // Begin by performing McKee's Fermat using simple precision when possible.
    // If this fails, switch to a mostly multi-precision version. Note
    // than the failure rate with the simpler version rises dramatically with
    // the bitlength of the number to factor.
    //
    if (!context->updated) {
        return perform_fermat_simple(machine);
    } else {
        return perform_fermat_mpz(machine);
    }
}
//-----------------------------------------------------------------------------
static ecode_t perform_fermat_simple(factoring_machine_t* const machine) {

    INIT_TIMER;
    INIT_NAMED_TIMER(sqrtm);
    INIT_NAMED_TIMER(greedy);

    PRINT_FERMAT_FACT_MSG;
    START_TIMER;

    fermat_context_t* const context = (fermat_context_t*) machine->context;
    //
    // Perform a race between the factorization of n and mult*n as suggested by
    // J. McKee in the reference paper "Speeding Fermat's Factoring Method".
    // There are not many comments in the code, but the algorithm is
    // straitforwardly implemented from the description given in the paper.
    //
    mpz_ptr n = context->n;
    
    mpz_t ni[NMULTIPLIERS];
    for (unsigned int i = 0; i < NMULTIPLIERS; i++) {
        mpz_init_set(ni[i], n);
        mpz_mul_ui(ni[i], ni[i], multipliers[i]);
    }
    
    mpz_t sqrtni[NMULTIPLIERS];
    for (unsigned int i = 0; i < NMULTIPLIERS; i++) {
        mpz_init(sqrtni[i]);
        mpz_sqrt(sqrtni[i], ni[i]);
    }
    
    if (0 == mpz_fits_ulong_p(sqrtni[NMULTIPLIERS-1])) {
        //
        // Actually we could try with to factorize without multiplier and so
        // check for (0 == mpz_fits_ulong_p(sqrtn)) to gain a few bits... but
        // this algorithm is so rapidly impractical that it doesn't matter in
        // practice...
        //
        for (unsigned int i = 0; i < NMULTIPLIERS; i++) {
            mpz_clear(ni[i]);
            mpz_clear(sqrtni[i]);
        }
        context->latest_p = 3; 
	    context->integer_too_large = true;
        return INTEGER_TOO_LARGE;
    }
    unsigned long int b[NMULTIPLIERS];
    unsigned long int max_y[NMULTIPLIERS];
    
    for (unsigned int i = 0; i < NMULTIPLIERS; i++) {
        b[i] = mpz_get_ui(sqrtni[i]) + 1;
        max_y[i] = sqrt(b[i]) * MAX_Y_BOUND_MULTIPLIER;
    }
    
    mpz_t Q;
    mpz_t gcd;
    mpz_t ny2;

    mpz_init(Q);
    mpz_init(gcd);
    mpz_init(ny2);
    
    unsigned long int ip = 1;
    unsigned long int p  = first_primes[ip];
    unsigned long int p2 = 0;
    
    long int x0 = 0;
                
    unsigned long int np = 0;
    unsigned long int s  = 0;
    unsigned long int ri = 0;
    unsigned long int xi = 0;
    unsigned long int yi = 0;
    
    ecode_t ecode = NO_FACTOR_FOUND;
            
    while (ip <  NFIRST_PRIMES) {
        
        p = first_primes[ip];
        ip++;
    
        if (p > TIFA_SQRT_LONG_MAX) {
            context->latest_p = p;
            goto clean_and_return;
        }
        p2 = p* p;
        
        for (int i = 0; i < NMULTIPLIERS; i++) {
            //
            // Step 1: Compute x1 and x2, the solutions to the equation
            //         Q(x, 1) = 0 (mod p^2) where Q(X, Y) = (X + b.Y)^2-n.Y^2
            //
            START_NAMED_TIMER(sqrtm);
        
            np = mpz_mod_ui(gcd, ni[i], p2);     
            s  = sqrtm_p2(np, p);
            
            if (s == NO_SQRT_MOD_P2) {
                //
                // n doesn't have modular squareroot mod p^2. Try with the next
                // multiplier or prime.
                //
                continue;
            }
            STOP_NAMED_TIMER(sqrtm);
            START_NAMED_TIMER(greedy);
                        
            //
            // Step 2a: Try to find a factor using x1
            //
            // Compute x0 = x1 mod p^2
            //            
            x0 = s - b[i];
            if (x0 < 0) {
                x0 = p2 - ((-x0) % p2);
            } else {
                x0 = x0 % p2;
            }
            
            //
            // Compute Q(x0, 1) mod n
            //
            mpz_set_ui(Q, x0);
            mpz_add_ui(Q, Q, b[i]);
            mpz_mul(Q, Q, Q);
            mpz_sub(Q, Q, ni[i]);
            
            if (MPZ_IS_SQUARE(Q)) {
                mpz_sqrt(Q, Q);
                mpz_set_ui(gcd, x0);
                mpz_add_ui(gcd, gcd, b[i]);
                mpz_sub(gcd, gcd, Q);
                mpz_gcd(gcd, gcd, ni[i]);
                
                if ((mpz_cmp_ui(gcd, 1) != 0) && (mpz_cmp(gcd, n) != 0) ) {   
                    ecode = SOME_FACTORS_FOUND;
                    goto compute_factors_and_return;
                }
            }
            if (0 == x0) {
                continue;
            }
            
            //
            // McKee's "Greedy" variant
            //
            ri = 1;
            xi = x0;
            yi = 1;
        
            while ((yi < max_y[i]) && (xi != 0)) {            
                //
                // Update ri, xi and yi
                //
                ri  = p2 / xi;
                
                if ((ri * xi) != p2) {
                    ri++;
                }
                xi  = xi * ri - p2;
                yi *= ri;
            
                //
                // Compute Q(xi, yi) mod n
                //
                mpz_set_ui(Q, yi);
                mpz_mul_ui(Q, Q, b[i]);
                mpz_add_ui(Q, Q, xi);
                mpz_mul(Q, Q, Q);
                
                mpz_set_ui(ny2, yi);
                mpz_mul_ui(ny2, ny2, yi);
                mpz_mul(ny2, ny2, ni[i]);
                mpz_sub(Q, Q, ny2);
            
                if (MPZ_IS_SQUARE(Q)) {
                    mpz_sqrt(Q, Q);
                    mpz_set_ui(gcd, yi);
                    mpz_mul_ui(gcd, gcd, b[i]);
                    mpz_add_ui(gcd, gcd, xi);
                    mpz_sub(gcd, gcd, Q);
                    mpz_gcd(gcd, gcd, n);
                    
                    if ((mpz_cmp_ui(gcd, 1) != 0) && (mpz_cmp(gcd, n) != 0) ) {
                        ecode = SOME_FACTORS_FOUND;
                        goto compute_factors_and_return;
                    }
                }
            }
            //
            // Step 2b: Try to find a factor using x2
            //
            // Compute x2 mod p^2
            //
            x0 = - (b[i] + s);
            if (x0 < 0) {
                x0 = p2 - ((-x0) % p2);
            } else {
                x0 = x0 % p2;
            }
            
            //
            // Compute Q(x2, 1) mod n
            //
            mpz_set_ui(Q, x0);
            mpz_add_ui(Q, Q, b[i]);
            mpz_mul(Q, Q, Q);
            mpz_sub(Q, Q, n);
            
            if (MPZ_IS_SQUARE(Q)) {
                mpz_sqrt(Q, Q);
                mpz_set_ui(gcd, x0);
                mpz_add_ui(gcd, gcd, b[i]);
                mpz_sub(gcd, gcd, Q);
                mpz_gcd(gcd, gcd, n);
                
                if ((mpz_cmp_ui(gcd, 1) != 0) && (mpz_cmp(gcd, n) != 0) ) {
                    ecode = SOME_FACTORS_FOUND;
                    goto compute_factors_and_return;
                }
            }
            if (0 == x0) {
                continue;
            }
            //
            // McKee's "Greedy" variant
            //
            ri = 1;
            xi = x0;
            yi = 1;
            
            while ((yi < max_y[i]) && (xi != 0)) {
                //
                // Update ri, xi and yi
                //
                ri  = p2 / xi;
                
                if ((ri * xi) != p2) {
                    ri++;
                }
                xi  = xi * ri - p2;
                yi *= ri;
                            
                //
                // Compute Q(xi, yi) mod n
                //
                mpz_set_ui(Q, yi);
                mpz_mul_ui(Q, Q, b[i]);
                mpz_add_ui(Q, Q, xi);
                mpz_mul(Q, Q, Q);
                
                mpz_set_ui(ny2, yi);
                mpz_mul_ui(ny2, ny2, yi);
                mpz_mul(ny2, ny2, ni[i]);
                mpz_sub(Q, Q, ny2);
                
                if (MPZ_IS_SQUARE(Q)) {
                    mpz_sqrt(Q, Q);
                    mpz_set_ui(gcd, yi);
                    mpz_mul_ui(gcd, gcd, b[i]);
                    mpz_add_ui(gcd, gcd, xi);
                    mpz_sub(gcd, gcd, Q);
                    mpz_gcd(gcd, gcd, n);
                    
                    if ((mpz_cmp_ui(gcd, 1) != 0) && (mpz_cmp(gcd, n) != 0) ) {
                        ecode = SOME_FACTORS_FOUND;
                        goto compute_factors_and_return;
                    }
                }
            }
            STOP_NAMED_TIMER(greedy);  
        }
    }

  compute_factors_and_return:

    PRINT_SQRTM_MSG(GET_NAMED_TIMING(sqrtm));
    PRINT_GREEDY_MSG(GET_NAMED_TIMING(greedy));

    if (ecode == SOME_FACTORS_FOUND) {
        append_mpz_to_array(machine->factors, gcd);
        mpz_divexact(gcd, context->n, gcd);
        append_mpz_to_array(machine->factors, gcd);
    }

  clean_and_return:
    
    for (unsigned int i = 0; i < NMULTIPLIERS; i++) {
        mpz_clear(ni[i]);
        mpz_clear(sqrtni[i]);
    }
        
    mpz_clear(Q);
    mpz_clear(gcd);
    mpz_clear(ny2);

    STOP_TIMER;
    PRINT_FERMAT_FACT_DONE_MSG;
    PRINT_TIMING;

    return ecode;
}
//-----------------------------------------------------------------------------
static ecode_t perform_fermat_mpz(factoring_machine_t* const machine) {

    INIT_TIMER;
    INIT_NAMED_TIMER(nextprime);
    INIT_NAMED_TIMER(sqrtm);
    INIT_NAMED_TIMER(greedy);

    fermat_context_t* const context = (fermat_context_t*) machine->context;

    mpz_ptr n = context->n;

    mpz_t ni[NMULTIPLIERS];
    for (unsigned int i = 0; i < NMULTIPLIERS; i++) {
        mpz_init_set(ni[i], n);
        mpz_mul_ui(ni[i], ni[i], multipliers[i]);
    }
    
    mpz_t sqrtni[NMULTIPLIERS];
    for (unsigned int i = 0; i < NMULTIPLIERS; i++) {
        mpz_init(sqrtni[i]);
        mpz_sqrt(sqrtni[i], ni[i]);
    }

    if (0 == mpz_fits_ulong_p(sqrtni[NMULTIPLIERS-1])) {
        //
        // Actually we could try with to factorize without multiplier... See
        // similar comment in perform_fermat_simple()
        //
        for (unsigned int i = 0; i < NMULTIPLIERS; i++) {
            mpz_clear(ni[i]);
            mpz_clear(sqrtni[i]);
        }
	    context->integer_too_large = true;
        return INTEGER_TOO_LARGE;
    }

    unsigned long int p;
    unsigned long int retval = 0;

    unsigned long int b[NMULTIPLIERS];
    unsigned long int max_y[NMULTIPLIERS];
    
    for (unsigned int i = 0; i < NMULTIPLIERS; i++) {
        b[i] = mpz_get_ui(sqrtni[i]) + 1;
        max_y[i] = sqrt(b[i]) * MAX_Y_BOUND_MULTIPLIER;
    }
    
    mpz_t pz;
    mpz_t pz2;
    mpz_t s;
    mpz_t x0;
    mpz_t xi;
    mpz_t ri;
    mpz_t yi;
    mpz_t Q;
    mpz_t gcd;
    mpz_t ny2;
    
    mpz_init_set_ui(pz, context->latest_p);
    mpz_init(pz2);
    mpz_init(s);
    mpz_init(x0);
    mpz_init(xi);
    mpz_init(ri);
    mpz_init(yi);
    mpz_init(Q);
    mpz_init_set(gcd, n);
    mpz_init(ny2);

    DECL_SQRTM_P2_UI_VARS();

    ecode_t ecode = NO_FACTOR_FOUND;

    PRINT_FERMAT_FACT_MSG;
    START_TIMER;
    //
    // Starts with the least prime greater than 2.n^1/4 and cycle though the
    // prime until we find a factor or the abort limit given by MAX_NPRIMES
    // is reached.
    //
    for (unsigned long int ip = 1; ip <= MAX_NPRIMES; ip++) {

        START_NAMED_TIMER(nextprime);

        mpz_nextprime(pz, pz);

        STOP_NAMED_TIMER(nextprime);

        if (0 == mpz_fits_ulong_p(pz)) {
            //
            // The prime is getting too large. Just stop and admit defeat...
            //
            ecode = FAILURE;
            goto clean_and_return;
        }
        p = mpz_get_ui(pz);
        mpz_mul(pz2, pz, pz);
        
        for (int i = 0; i < NMULTIPLIERS; i++) {
            //
            // Step 1: Compute x1 and x2, the solutions to the equation
            //         Q(x, 1) = 0 (mod p^2) where Q(X, Y) = (X + b.Y)^2 - n.Y^2
            //
            START_NAMED_TIMER(sqrtm);
            
            SQRTM_P2_UI(retval, s, ni[i], p);
            
            STOP_NAMED_TIMER(sqrtm);
            
            if (retval != 0) {
                //
                // n doesn't have modular squareroot mod p^2. Try with the next
                // prime.
                //
                continue;
            }        
            START_NAMED_TIMER(greedy);
            
            //
            // Step 2a: Try to find a factor using x1
            //
            // Compute x0 = x1 mod p^2
            //
            mpz_sub_ui(x0, s, b[i]);
            mpz_mod(x0, x0, pz2);
            
            //
            // Compute Q(x0, 1) mod n
            //
            mpz_add_ui(Q, x0, b[i]);
            mpz_mul(Q, Q, Q);
            mpz_sub(Q, Q, ni[i]);
            
            if (MPZ_IS_SQUARE(Q)) {
                mpz_sqrt(Q, Q);
                mpz_add_ui(x0, x0, b[i]);
                mpz_sub(x0, x0, Q);
                mpz_gcd(gcd, x0, n);
                if ((mpz_cmp_ui(gcd, 1) != 0) && (mpz_cmp(gcd, n) != 0) ) {
                    ecode = SOME_FACTORS_FOUND;
                    goto compute_factors_and_return;
                }
            }
            if (0 == mpz_sgn(x0)) {
                continue;
            }
            
            //
            // McKee's "Greedy" variant
            //
            mpz_set(xi, x0);
            mpz_set_ui(yi, 1);
            mpz_set_ui(ri, 1);
            
            while ((mpz_cmp_ui(yi, max_y[i]) < 0) && (mpz_sgn(xi) != 0)) {
                //
                // Update ri, xi and yi
                //
                mpz_cdiv_qr(ri, xi, pz2, xi);
                mpz_neg(xi, xi);
                mpz_mul(yi, yi, ri);            
                //
                // Compute Q(xi, yi) mod n
                //
                mpz_mul_ui(Q, yi, b[i]);
                mpz_add(Q, Q, xi);
                mpz_mul(Q, Q, Q);
                mpz_mul(ny2, yi, yi);
                mpz_mul(ny2, ny2, ni[i]);
                mpz_sub(Q, Q, ny2);
            
                if (MPZ_IS_SQUARE(Q)) {
                    mpz_sqrt(Q, Q);
                    mpz_mul_ui(gcd, yi, b[i]);
                    mpz_add(gcd, gcd, xi);
                    mpz_sub(gcd, gcd, Q);
                    mpz_gcd(gcd, gcd, n);
                    if ((mpz_cmp_ui(gcd, 1) != 0) && (mpz_cmp(gcd, n) != 0) ) {
                        ecode = SOME_FACTORS_FOUND;
                        goto compute_factors_and_return;
                    }
                }
            }
            //
            // Step 2b: Try to find a factor using x2
            //
            // Compute x2 mod p^2
            //
            mpz_add_ui(x0, s, b[i]);
            mpz_neg(x0, x0);
            mpz_mod(x0, x0, pz2);
            
            //
            // Compute Q(x2, 1) mod n
            //
            mpz_add_ui(Q, x0, b[i]);
            mpz_mul(Q, Q, Q);
            mpz_sub(Q, Q, n);
            
            if (MPZ_IS_SQUARE(Q)) {
                mpz_sqrt(Q, Q);
                mpz_add_ui(x0, x0, b[i]);
                mpz_sub(x0, x0, Q);
                mpz_gcd(gcd, x0, n);
                if ((mpz_cmp_ui(gcd, 1) != 0) && (mpz_cmp(gcd, n) != 0) ) {
                    ecode = SOME_FACTORS_FOUND;
                    goto compute_factors_and_return;
                }
            }
            if (0 == mpz_sgn(x0)) {
                continue;
            }
            //
            // McKee's "Greedy" variant
            //
            mpz_set(xi, x0);
            mpz_set_ui(yi, 1);
            mpz_set_ui(ri, 1);
            
            while ((mpz_cmp_ui(yi, max_y[i]) < 0) && (mpz_sgn(xi) != 0)) {
                //
                // Update ri, xi and yi
                //
                mpz_cdiv_qr(ri, xi, pz2, xi);
                mpz_neg(xi, xi);
                mpz_mul(yi, yi, ri);            
                //
                // Compute Q(xi, yi) mod n
                //
                mpz_mul_ui(Q, yi, b[i]);
                mpz_add(Q, Q, xi);
                mpz_mul(Q, Q, Q);
                mpz_mul(ny2, yi, yi);
                mpz_mul(ny2, ny2, ni[i]);
                mpz_sub(Q, Q, ny2);
                if (MPZ_IS_SQUARE(Q)) {
                    mpz_sqrt(Q, Q);
                    mpz_mul_ui(gcd, yi, b[i]);
                    mpz_add(gcd, gcd, xi);
                    mpz_sub(gcd, gcd, Q);
                    mpz_gcd(gcd, gcd, n);
                    if ((mpz_cmp_ui(gcd, 1) != 0) && (mpz_cmp(gcd, n) != 0) ) {
                        ecode = SOME_FACTORS_FOUND;
                        goto compute_factors_and_return;
                    }
                }
            }
            STOP_NAMED_TIMER(greedy);
        }
    }

  compute_factors_and_return:
  
    PRINT_NEXTPRIME_MSG(GET_NAMED_TIMING(nextprime));
    PRINT_SQRTM_MSG(GET_NAMED_TIMING(sqrtm));
    PRINT_GREEDY_MSG(GET_NAMED_TIMING(greedy));
    
    if (ecode == SOME_FACTORS_FOUND) {
        append_mpz_to_array(machine->factors, gcd);
        mpz_divexact(gcd, context->n, gcd);
        append_mpz_to_array(machine->factors, gcd);
    }

  clean_and_return:

    CLEAR_SQRTM_P2_UI_VARS();

    for (unsigned int i = 0; i < NMULTIPLIERS; i++) {
        mpz_clear(ni[i]);
        mpz_clear(sqrtni[i]);
    }
    
    mpz_clear(pz);
    mpz_clear(pz2);
    mpz_clear(s);
    mpz_clear(x0);
    mpz_clear(xi);
    mpz_clear(ri);
    mpz_clear(yi);
    mpz_clear(Q);
    mpz_clear(gcd);
    mpz_clear(ny2);

    STOP_TIMER;
    PRINT_FERMAT_FACT_DONE_MSG;
    PRINT_TIMING;

    return ecode;
}
//-----------------------------------------------------------------------------
