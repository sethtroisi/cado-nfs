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
 * \file    factoring_program.c
 * \author  Jerome Milan
 * \date    Wed Mar 5 2008
 * \version 1.1
 */

/*
 * History:
 * --------
 *   1.1: Wed Mar 5 2008 by JM:
 *        - Rewritten and simplified since tdiv(...) changed.
 *   1.0: Circa February (March?) 2007 by JM:
 *        - Initial version.
 */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <inttypes.h>
#include <stdbool.h>
#include <string.h>

#include "tifa_config.h"
#include "funcs.h"
#include "macros.h"
#include "tdiv.h"
#include "factoring_program.h"

#include "tool_utils.h"
#include "common_funcs.h"

#define __TIMING__ 1
#include "../../lib/utils/include/timer.h"

//-----------------------------------------------------------------------------
#define PROG_DFLT_ARRAY_LENGTH 12
//-----------------------------------------------------------------------------
void print_factors(mpz_array_t* const, uint32_array_t* const);
//-----------------------------------------------------------------------------
ecode_t run_program(factoring_program_t* const program) {
    
    ecode_t rcode = NO_FACTOR_FOUND;
    //
    // Reads arguments given on the command line (if any) and initializes
    // all the needed parameters.
    //
    program->process_args_func(program);

    //
    // Of course, check if program->n is prime before proceeding.
    //
    if (MPZ_IS_PRIME(program->n)) {
        gmp_printf("\nOOPS: %Zd is (probably) prime!\n", program->n);
        return 0;
    }
    
    //
    // Prints program's parameters.
    //
    printf("Parameters used:\n");
    printf("----------------\n");
    program->print_params_func(program);

    printf("\nInteger to factor:\n");
    printf("------------------\n");
    gmp_printf("\t%Zd\n\n", program->n);

    //
    // Factors & multiplicities for program->n.
    //
    mpz_array_t*    factors = alloc_mpz_array(PROG_DFLT_ARRAY_LENGTH);
    uint32_array_t* multis  = alloc_uint32_array(PROG_DFLT_ARRAY_LENGTH);

    //
    // Proceed with the trial divisions...
    //
    ecode_t ecode = tdiv(factors, multis, program->n, program->nprimes_tdiv);

    if (ecode == COMPLETE_FACTORIZATION_FOUND) {
        printf("\n");
        print_factors(factors, multis);
        rcode = COMPLETE_FACTORIZATION_FOUND;
        goto clear_tdiv_and_return;
    }
    
    //
    // Some factors were found but they cannot account for the complete
    // factorization of program->n.
    //
    mpz_t unfactored;
    mpz_init(unfactored);
    
    if (ecode == NO_FACTOR_FOUND) {
        mpz_init_set(unfactored, program->n);
    } else {
        mpz_init_set(unfactored, factors->data[factors->length - 1]);
        factors->length--;
        multis->length--;
        rcode = PARTIAL_FACTORIZATION_FOUND;
    }
    
    printf("\n");
    printf("Integer to factor after trial division:\n");
    printf("---------------------------------------\n");
    gmp_printf("\t%Zd\n\n", unfactored);
    
    unsigned int shift = 0;
    //
    // First, check if the unfactored part is not a perfect square.
    //
    while (MPZ_IS_SQUARE(unfactored)) {
        mpz_sqrt(unfactored, unfactored);
        shift++;
        if (mpz_cmp_ui(unfactored, 1) == 0) {
            break;
        }
    }
    
    if (MPZ_IS_PRIME(unfactored)) {
        //
        // We have found the complete factorization!
        //
        append_mpz_to_array(factors, unfactored);
        append_uint32_to_array(multis, 1 << shift);
        
        print_factors(factors, multis);
        rcode = COMPLETE_FACTORIZATION_FOUND;
        goto clear_sqt_and_return;
    }
    //
    // Factors & multiplicities for the unfactored part of program->n.
    //
    mpz_array_t*    progfa = alloc_mpz_array(PROG_DFLT_ARRAY_LENGTH);
    uint32_array_t* progmu = alloc_uint32_array(PROG_DFLT_ARRAY_LENGTH);
    
    if (program->verbose || program->timing) {
        printf("%s trace:\n", program->algo_name);
        int len = strlen(program->algo_name);
        for (int i=0; i < len; i++) {
            printf("-");
        }
        printf("-------\n\n");
    }
    
    ecode = program->factoring_algo_func(
                progfa,
                progmu,
                unfactored,
                program->params,
                program->mode
            );
    
    if (program->verbose || program->timing) {
        printf("\n");
    }
        
    switch (ecode) {
        case COMPLETE_FACTORIZATION_FOUND:
            append_mpz_array(factors, progfa);
            append_uint32_array(multis, progmu);
            print_factors(factors, multis);
            rcode = COMPLETE_FACTORIZATION_FOUND;
            break;
        case SOME_FACTORS_FOUND:
        case SOME_COPRIME_FACTORS_FOUND:
        case PARTIAL_FACTORIZATION_FOUND: {
            //
            // Simplify the list of found factors by computing and
            // keeping only coprime factors.
            //
            uint32_t     blen = progfa->length;
            mpz_array_t* base = alloc_mpz_array(blen * (blen + 1));

            find_coprime_base(base, unfactored, progfa);
            ins_sort_mpz_array(base);

            append_mpz_array(factors, base);
            
            print_factors(factors, multis);
            
            clear_mpz_array(base);
            free(base);
            rcode = PARTIAL_FACTORIZATION_FOUND;
            break;
        }
        default:
            break;
    }
    
    clear_mpz_array(progfa);
    clear_uint32_array(progmu);
    
  clear_sqt_and_return:
    
    mpz_clear(unfactored);
    
  clear_tdiv_and_return:
    
    clear_mpz_array(factors);
    clear_uint32_array(multis);
    
    switch (rcode) {
        case COMPLETE_FACTORIZATION_FOUND:
            printf("Factorization is complete.\n\n");
            break;
        case PARTIAL_FACTORIZATION_FOUND:
            printf("Factorization is (potentially) not complete.\n\n");
            break;
        case NO_FACTOR_FOUND:
        default:
            printf("Factorization failed.\n\n");
            break;
    }
    
    return rcode;   
}
//-----------------------------------------------------------------------------
void print_factors(mpz_array_t* const factors, uint32_array_t* const multis) {
    printf("Found %2u factor(s):\n", (unsigned int)factors->length); 
    printf("-------------------\n");
    
    uint32_t nmult = 0;
    if (multis == NULL) {
        nmult = 0;
    } else {
        nmult = multis->length;
    }
    
    for (uint32_t i = 0; i != nmult; i++) {
        //
        // _WARNING_: The multis array can be shorter than the factors array 
        //            since the factors obtained by trial divisions always come 
        //            with their corresponding multiplicities. This assumes
        //            that these factors are given at the beginning of the
        //            factors array.
        //
        gmp_printf("\t%Zd ^ %"PRIu32"\n", factors->data[i], multis->data[i]);
        fflush(stdout);
    }
    for (uint32_t i = nmult; i != factors->length; i++) {
        gmp_printf("\t%Zd\n", factors->data[i]);
        fflush(stdout);
    }
    printf("\n");
}
//-----------------------------------------------------------------------------
