//
// Copyright (C) 2008 INRIA (French National Institute for Research
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
 * \file    ecm_factors.c
 * \author  Jerome Milan
 * \date    Jan 14 2008
 * \version 1.0
 */

#include <tifa_config.h>

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <stdbool.h>

#include "tool_utils.h"
#include "tdiv.h"
#include "ecm.h"
#include "factoring_program.h"
#include "common_funcs.h"

//------------------------------------------------------------------------------
static void print_usage(factoring_program_t* const program) {
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "%15s\n", program->argv[0]);
    fprintf(stderr, "or:\n");
    fprintf(stderr, "%15s\t<number_to_factor>\n", program->argv[0]);
    fprintf(stderr, "or:\n");
    fprintf(stderr, "%15s\t<nprimes_tdiv> <number_to_factor>\n", 
            program->argv[0]);
    fprintf(stderr, "or:\n");
    fprintf(stderr, "%15s <nprimes_tdiv> <b1> <b2> <ncurves> "       
            "<number_to_factor>\n\n", program->argv[0]);

    PRINT_USAGE_WARNING_MSG();
}
//------------------------------------------------------------------------------
static void process_args(factoring_program_t* const program) {
    
    ecm_params_t* params = (ecm_params_t*) program->params;
    
    int    argc = program->argc;
    char** argv = program->argv;
    
    uint32_t* nprimes_tdiv = &(program->nprimes_tdiv);
    uint32_t* nfactors     = &(program->nfactors);
    
    print_hello_msg(program->algo_name);
    
    switch (argc) {
    case 1: {
        //
        // No argument provided: proceed in interactive mode, i.e. get the
        // number to factor and then use ECM optimal default values.
        //
        char str_buffer[MAX_NDIGITS];
        PRINT_ENTER_NUMBER_MSG();
        char* str_factor_me = fgets(str_buffer, MAX_NDIGITS, stdin);
        printf("\n");
        chomp(str_factor_me, MAX_NDIGITS);
    
        if (!is_a_number(str_factor_me, MAX_NDIGITS)) {
            PRINT_NAN_ERROR(str_factor_me);
            exit(-1);
        }
        mpz_init_set_str(program->n, str_factor_me, 10);
        set_ecm_params_to_default(program->n, params);
        *nprimes_tdiv = NPRIMES_TRIAL_DIV;
        *nfactors     = 2;
    
        break;
    }
    case 2: {
        //
        // Only one argument provided: the number to factor. Let the
        // program choose the default CFRAC parameters.
        //
        if (!is_a_number(argv[1], MAX_NDIGITS)) {
            PRINT_NAN_ERROR(argv[1]);
            exit(-1);
        }
        mpz_init_set_str(program->n, argv[1], 10);
        set_ecm_params_to_default(program->n, params);
        *nprimes_tdiv = NPRIMES_TRIAL_DIV;
        *nfactors     = 2;
    
        break;
    }
    case 3: {
        //
        // Two arguments provided: the <nprimes_tdiv> parameter
        // and the number to factor. Let the program choose the default
        // ECM parameters.
        //
        if (!is_a_number(argv[1], MAX_NDIGITS)) {
            PRINT_NAN_ERROR(argv[1]);
            exit(-1);
        }
        *nprimes_tdiv = strtoul(argv[1], NULL, 10);
    
        if (!is_a_number(argv[2], MAX_NDIGITS)) {
            PRINT_NAN_ERROR(argv[2]);
            exit(-1);
        }
        mpz_init_set_str(program->n, argv[2], 10);
        set_ecm_params_to_default(program->n, params);
        *nfactors = 2;
        break;
    }
    case 6: {
        //
        // Fives arguments provided: the <nprimes_tdiv> parameter,
        // the two bounds <b1> and <b2>, the number <ncurves> of curves to try 
        // before giving up and finally the number to factor.
        //
        if (!is_a_number(argv[1], MAX_NDIGITS)) {
            PRINT_NAN_ERROR(argv[1]);
            exit(-1);
        }
        *nprimes_tdiv = strtoul(argv[1], NULL, 10);
        if (!is_a_number(argv[4], MAX_NDIGITS)) {
            PRINT_NAN_ERROR(argv[4]);
            exit(-1);
        }
        mpz_init_set_str(program->n, argv[5], 10);
                
        if (!is_a_number(argv[2], MAX_NDIGITS)) {
            PRINT_NAN_ERROR(argv[2]);
            exit(-1);
        }
        params->b1 = strtoul(argv[2], NULL, 10);
        if (!is_a_number(argv[3], MAX_NDIGITS)) {
            PRINT_NAN_ERROR(argv[3]);
            exit(-1);
        }
        params->b2 = strtoul(argv[3], NULL, 10);
        if (!is_a_number(argv[4], MAX_NDIGITS)) {
            PRINT_NAN_ERROR(argv[4]);
            exit(-1);
        }
        params->ncurves = strtoul(argv[4], NULL, 10);
        
        *nfactors = 2;
        break;
    }
    default:
        PRINT_BAD_ARGC_ERROR();
        print_usage(program);
        exit(-1);
    }    
    return;
}
//------------------------------------------------------------------------------
static void set_params_to_default(factoring_program_t* const program) {
    set_ecm_params_to_default(program->n, (ecm_params_t*) program->params);
}
//------------------------------------------------------------------------------
static void print_params(factoring_program_t* const program) {
    ecm_params_t* params = (ecm_params_t*) program->params;
    printf("\tb1      : %"PRIu32"\n", params->b1);
    printf("\tb2      : %"PRIu32"\n", params->b2);
    printf("\tncurves : %"PRIu32"\n", params->ncurves);
}
//------------------------------------------------------------------------------
static ecode_t ecm_func(mpz_array_t* const factors,
                          uint32_array_t* const multis, const mpz_t n,
                          const void* const params, factoring_mode_t mode) {
    return ecm(factors, multis, n, (const ecm_params_t*) params, mode);
}
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

    ecm_params_t params;
    factoring_program_t program;
        
    program.argc = argc;
    program.argv = argv;

    program.verbose = TIFA_VERBOSE_ECM;
    program.timing  = TIFA_TIMING_ECM;

    program.algo_name = "ECM";
    program.params    = (void*) &params;
    program.mode      = FIND_SOME_FACTORS;
    //program.mode      = FIND_COMPLETE_FACTORIZATION;
    //program.mode      = SINGLE_RUN;

    program.print_usage_func           = print_usage;
    program.print_params_func          = print_params;
    program.process_args_func          = process_args;
    program.factoring_algo_func        = ecm_func;
    program.set_params_to_default_func = set_params_to_default;

    ecode_t ecode = run_program(&program);
        
    return ecode;
}
//------------------------------------------------------------------------------
