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
 * \file    tdiv_factors.c
 * \author  Jerome Milan
 * \date    Wed Mar 14 2007
 * \version 1.1
 */

/*
 * History:
 *   1.1: Wed Mar 14 2007 by JM
 *        - Rewritten to use a factoring_program_t.
 *   1.0: Thu Mar 16 2006 by JM
 *        - Initial version.
 */

#include "tifa_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <stdbool.h>

#include "tool_utils.h"
#include "tdiv.h"
#include "factoring_program.h"
#include "common_funcs.h"

//
// The number of arguments accepted is either 0, or TDIV_F_MAX_ARGC
//
#define TDIV_F_MAX_ARGC 3

//------------------------------------------------------------------------------
static void print_usage(factoring_program_t* const program) {
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "%15s <nprimes_to_trial_divide_by> <number_to_factor>\n",
            program->argv[0]);
    fprintf(stderr, "or:\n");
    fprintf(stderr, "%15s\n\n", program->argv[0]);

    PRINT_USAGE_WARNING_MSG();
}
//------------------------------------------------------------------------------
static void process_args(factoring_program_t* const program) {
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
        *nprimes_tdiv = NPRIMES_TRIAL_DIV;
        *nfactors     = 0;
        
        break;
    }
    case TDIV_F_MAX_ARGC: {
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
        
        *nprimes_tdiv = strtoul(argv[1], endptr, 10);
        *nfactors     = 0;
        
        mpz_init_set_str(program->n, argv[2], 10);
        
        break;
    }
    default:
        PRINT_BAD_ARGC_ERROR();
        print_usage(program);
        exit(-1);
    }
}
//------------------------------------------------------------------------------
static void set_params_to_default(factoring_program_t* const program
                                  __attribute__ ((unused))) {
    return;
}
//------------------------------------------------------------------------------
static void print_params(factoring_program_t* const program) {
    printf("\tnprimes_tdiv: %u\n", program->nprimes_tdiv);
}
//------------------------------------------------------------------------------
static ecode_t tdiv_func(mpz_array_t* const factors   __attribute__ ((unused)),
                         uint32_array_t* const multis __attribute__ ((unused)),
                         const mpz_t n              __attribute__ ((unused)),
                         const void* const params   __attribute__ ((unused)),
                         factoring_mode_t mode      __attribute__ ((unused))) {

#if TIFA_VERBOSE_TDIV || TIFA_TIMING_TDIV
    printf("See above trace\n");
#endif
    return SUCCESS;
}
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

    factoring_program_t program;

    program.argc = argc;
    program.argv = argv;

    program.verbose = TIFA_VERBOSE_TDIV;
    program.timing  = TIFA_TIMING_TDIV;

    program.algo_name = "trial division";
    program.params    = NULL;
    program.mode      = FIND_SOME_FACTORS;

    program.print_usage_func           = print_usage;
    program.print_params_func          = print_params;
    program.process_args_func          = process_args;
    program.factoring_algo_func        = tdiv_func;
    program.set_params_to_default_func = set_params_to_default;

    ecode_t ecode = run_program(&program);

    return ecode;
}
//------------------------------------------------------------------------------

