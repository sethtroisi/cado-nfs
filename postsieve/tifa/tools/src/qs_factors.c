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
 * \file    qs_factors.c
 * \author  Jerome Milan
 * \date    Tue Mar 13 2007
 * \version 1.1
 */

/*
 * History:
 *   1.1: Tue Mar 13 2007 by JM
 *        - Rewritten to use a factoring_program_t.
 *   1.0: Thu Jul 6 2006 by JM
 *        - Initial version.
 */

#include <tifa_config.h>

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <stdbool.h>

#include "tool_utils.h"
#include "tdiv.h"
#include "qs.h"
#include "factoring_program.h"
#include "common_funcs.h"

//
// The number of arguments accepted is either 0, 1, 2, or (QS_F_MAX_ARGC-1)
//
#define QS_F_MAX_ARGC 8

//------------------------------------------------------------------------------
static void print_usage(factoring_program_t* const program) {
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "%15s\t<nprimes_in_base>\n", program->argv[0]);
    fprintf(stderr, "%15s\t<nprimes_tdiv_smooth_residues>\n", "");
    fprintf(stderr, "%15s\t<nrelations>\n", "");
    fprintf(stderr, "%15s\t<linalg_method>\n", "");
    fprintf(stderr, "%15s\t<use_large_primes>\n", "");
    fprintf(stderr, "%15s\t<nprimes_tdiv>\n", "");
    fprintf(stderr, "%15s\t<number_to_factor>\n", "");
    fprintf(stderr, "or:\n");
    fprintf(stderr, "%15s\n", program->argv[0]);
    fprintf(stderr, "or:\n");
    fprintf(stderr, "%15s\t<number_to_factor>\n", program->argv[0]);
    fprintf(stderr, "or:\n");
    fprintf(stderr, "%15s\t<nprimes_tdiv> <number_to_factor>\n\n",
            program->argv[0]);
    PRINT_USAGE_WARNING_MSG();
}
//------------------------------------------------------------------------------
static void process_args(factoring_program_t* const program) {

    qs_params_t* params = (qs_params_t*) program->params;

    int    argc = program->argc;
    char** argv = program->argv;

    uint32_t* nprimes_tdiv = &(program->nprimes_tdiv);
    uint32_t* nfactors     = &(program->nfactors);

    print_hello_msg(program->algo_name);

    switch (argc) {
    case 1: {
        //
        // No argument provided: proceed in interactive mode, i.e. get the
        // number to factor and then use CFRAC optimal default values.
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
        set_qs_params_to_default(program->n, params);
        *nprimes_tdiv = NPRIMES_TRIAL_DIV;
        *nfactors     = params->nrelations;

        break;
    }
    case 2: {
        //
        // Only one argument provided: the number to factor. Let the
        // program choose the default QS parameters.
        //
        if (!is_a_number(argv[1], MAX_NDIGITS)) {
            PRINT_NAN_ERROR(argv[1]);
            exit(-1);
        }
        mpz_init_set_str(program->n, argv[1], 10);
        set_qs_params_to_default(program->n, params);
        *nprimes_tdiv = NPRIMES_TRIAL_DIV;
        *nfactors     = params->nrelations;
        break;
    }
    case 3: {
        //
        // Two arguments provided: the <nprimes_tdiv> parameter
        // and the number to factor. Let the program choose the default
        // QS parameters.
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
        set_qs_params_to_default(program->n, params);
        *nfactors = params->nrelations;
        break;
    }
    case QS_F_MAX_ARGC: {
        //
        // Read parameters on the command line as they are provided by
        // the factorize.pl script... The script already checks for the validity
        // of the parameters, but let's check one more time while we're at it.
        //
        for (int i = 1; i < argc; i++) {
            if (!is_a_number(argv[i], MAX_NDIGITS)) {
                PRINT_NAN_ERROR(argv[i]);
                print_usage(program);
                exit(-1);
            }
        }
        char** endptr = NULL;
        uint32_t use_lp_variation = 0;

        params->nprimes_in_base = strtoul(argv[1], endptr, 10);
        params->nprimes_tdiv    = strtoul(argv[2], endptr, 10);
        params->nrelations      = strtoul(argv[3], endptr, 10);
        params->linalg_method   = strtoul(argv[4], endptr, 10);
        use_lp_variation        = strtoul(argv[5], endptr, 10);
        *nprimes_tdiv           = strtoul(argv[6], endptr, 10);
        *nfactors               = params->nrelations;

        mpz_init_set_str(program->n, argv[7], 10);

        if (0 == use_lp_variation) {
            params->use_large_primes = false;
        } else {
            params->use_large_primes = true;
        }

        if (params->nprimes_tdiv == 0) {
            params->nprimes_tdiv = 1;
        }
        if (params->nprimes_tdiv > params->nprimes_in_base) {
            params->nprimes_tdiv = params->nprimes_in_base;
        }

        break;
    }
    default:
        PRINT_BAD_ARGC_ERROR();
        print_usage(program);
        exit(-1);
    }
}
//------------------------------------------------------------------------------
static void set_params_to_default(factoring_program_t* const program) {
    set_qs_params_to_default(program->n, (qs_params_t*) program->params);
}
//------------------------------------------------------------------------------
static void print_params(factoring_program_t* const program) {
    qs_params_t* params = (qs_params_t*) program->params;

    printf("\tnprimes_in_base      : %u\n", params->nprimes_in_base);
    printf("\tnprimes_tdiv_residues: %u\n", params->nprimes_tdiv);
    printf("\tnrelations           : %u\n", params->nrelations);
    printf("\tlinalg_method        : %u\n", params->linalg_method);
    printf("\tnprimes_tdiv         : %u\n", program->nprimes_tdiv);
    if (params->use_large_primes) {
        printf("\tuse_large_primes     : yes\n");
    } else {
        printf("\tuse_large_primes     : no\n");
    }
}
//------------------------------------------------------------------------------
static ecode_t qs_func(mpz_array_t* const factors,
                       uint32_array_t* const multis, const mpz_t n,
                       const void* const params, factoring_mode_t mode) {
    return qs(factors, multis, n, (const qs_params_t*) params, mode);
}
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

    qs_params_t params;
    factoring_program_t program;

    program.argc = argc;
    program.argv = argv;

    program.verbose = TIFA_VERBOSE_QS;
    program.timing  = TIFA_TIMING_QS;

    program.algo_name = "QS";
    program.params    = (void*) &params;
    program.mode      = SINGLE_RUN;

    program.print_usage_func           = print_usage;
    program.print_params_func          = print_params;
    program.process_args_func          = process_args;
    program.factoring_algo_func        = qs_func;
    program.set_params_to_default_func = set_params_to_default;

    ecode_t ecode = run_program(&program);

    return ecode;
}
//------------------------------------------------------------------------------
