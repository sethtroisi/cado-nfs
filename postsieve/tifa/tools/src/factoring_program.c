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
 * \date    Circa February (March?) 2007
 * \version 1.0
 */

#include <tifa_config.h>

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <inttypes.h>
#include <stdbool.h>
#include <string.h>

#include "tool_utils.h"
#include "tdiv.h"
#include "factoring_program.h"
#include "funcs.h"
#include "common_funcs.h"

//------------------------------------------------------------------------------
ecode_t run_program(factoring_program_t* const program) {
    //
    // Reads arguments given on the command line (if any) and initializes
    // all the needed parameters.
    //
    program->process_args_func(program);

    //
    // Of course, check if program->n is prime before proceeding.
    //
    uint32_t is_prime = mpz_probab_prime_p(program->n, NTRIES_MILLER_RABIN);
    if (2 == is_prime) {
        gmp_printf("\nOOPS: %Zd is prime!\n", program->n);
        return 0;
    }
    if (1 == is_prime) {
        gmp_printf("\nOOPS: %Zd is probably prime!\n", program->n);
        return 0;
    }

    printf("Parameters used:\n");
    printf("----------------\n");
    program->print_params_func(program);

    printf("\n");
    printf("Integer to factor:\n");
    printf("------------------\n");
    gmp_printf("\t%Zd\n", program->n);
    printf("\n");

    //
    // The remaining part of program->n (a)fter the (t)rial (d)ivisions.
    //
    mpz_t n_atd;
    mpz_init(n_atd);

    //
    // mpz_array_t holding the factors to find via trial divisions.
    //
    mpz_array_t* tdivfactors = alloc_mpz_array(program->nprimes_tdiv);

    //
    // uint32_array_t holding the multiplicities to find via trial divisions.
    //
    uint32_array_t* tdivmultis = alloc_uint32_array(program->nprimes_tdiv);

    //
    // Proceed with the trial divisions...
    //
    tdiv(n_atd, tdivfactors, tdivmultis, program->n, program->nprimes_tdiv);

    ecode_t ecode = NO_FACTOR_FOUND;

    if (   (0 != mpz_cmp_ui(n_atd, 1))
        && (0 == mpz_probab_prime_p(n_atd, NTRIES_MILLER_RABIN)) ) {
        //
        // After trial division, n_atd is still not prime.
        // First, check if n_atd is not a perfect square.
        //
        if (0 != mpz_perfect_square_p(n_atd)) {
            mpz_t sqroot;
            mpz_init(sqroot);
            mpz_sqrt(sqroot, n_atd);

            printf("\n");
            printf("Found %2u Factors:\n", tdivfactors->length + 1);
            printf("-----------------\n");
            gmp_printf("%Zd = \n", program->n);
            if (tdivfactors->length != 0) {
                gmp_printf("\t  %Zd",        tdivfactors->data[0]);
                gmp_printf(" ^ %"PRIu32"\n", tdivmultis->data[0]);
                for (uint32_t i = 1; i < tdivfactors->length; i++) {
                    gmp_printf("\t* %Zd",        tdivfactors->data[i]);
                    gmp_printf(" ^ %"PRIu32"\n", tdivmultis->data[i]);
                }
                gmp_printf("\t* %Zd ^ 2\n\n", sqroot);
            } else {
                gmp_printf("\t %Zd ^ 2\n\n", sqroot);
            }
            mpz_clear(sqroot);

        } else {
            //
            // Not a square: continue factorization using the given algorithm.
            //
            printf("\n");
            printf("Integer to factor after trial division:\n");
            printf("---------------------------------------\n");
            gmp_printf("\t%Zd\n\n", n_atd);

            if (program->verbose || program->timing) {
                printf("%s trace:\n", program->algo_name);
                int len = strlen(program->algo_name);
                for (int i=0; i < len; i++) {
                    printf("-");
                }
                printf("-------\n\n");
            }
            mpz_array_t*    progfactors = alloc_mpz_array(program->nfactors);
            uint32_array_t* progmultis  = alloc_uint32_array(program->nfactors);
            //
            // Proceed with the factoring algorithm.
            //
            ecode = program->factoring_algo_func(
                        progfactors,
                        progmultis,
                        n_atd,
                        program->params,
                        program->mode
                    );
            
            if (program->mode != FIND_COMPLETE_FACTORIZATION) {
                //
                // Simplify the list of found factors by computing and
                // keeping only coprime factors.
                //
                uint32_t aflen    = (tdivfactors->length + progfactors->length);
                mpz_array_t* base = alloc_mpz_array(aflen * (aflen + 1));

                find_coprime_base(base, program->n, tdivfactors);
                find_coprime_base(base, program->n, progfactors);
                ins_sort_mpz_array(base);

                if (program->verbose || program->timing) {
                    printf("\n");
                }

                printf("Found %2u Factors:\n", base->length);
                printf("-----------------\n");
                print_mpz_array(base);
                printf("\n");

                clear_mpz_array(base);
                free(base);

            } else {
                if (program->verbose || program->timing) {
                    printf("\n");
                }

                printf("Found %2u Factors:\n",
                       tdivfactors->length + progfactors->length);
                printf("-----------------\n");

                if (tdivfactors->length > 0) {

                    gmp_printf("%Zd = \n", program->n);
                    gmp_printf("\t  %Zd",        tdivfactors->data[0]);
                    gmp_printf(" ^ %"PRIu32"\n", tdivmultis->data[0]);

                    for (uint32_t i = 1; i < tdivfactors->length; i++) {
                        gmp_printf("\t* %Zd",        tdivfactors->data[i]);
                        gmp_printf(" ^ %"PRIu32"\n", tdivmultis->data[i]);
                    }
                    for (uint32_t i = 0; i < progfactors->length; i++) {
                        gmp_printf("\t* %Zd",        progfactors->data[i]);
                        gmp_printf(" ^ %"PRIu32"\n", progmultis->data[i]);
                    }
                } else {
                    gmp_printf("%Zd = \n", program->n);
                    gmp_printf("\t  %Zd",        progfactors->data[0]);
                    gmp_printf(" ^ %"PRIu32"\n", progmultis->data[0]);
                    for (uint32_t i = 1; i < progfactors->length; i++) {
                        gmp_printf("\t* %Zd",        progfactors->data[i]);
                        gmp_printf(" ^ %"PRIu32"\n", progmultis->data[i]);
                    }
                }
                printf("\n");
            }
            clear_mpz_array(progfactors);
            free(progfactors);

            clear_uint32_array(progmultis);
            free(progmultis);
        }
    } else {
        //
        // Trial divisions were (probably!) enough to completely factor the
        // number.
        //
        printf("\n");
        printf("Found %2u Factors:\n", tdivfactors->length + 1);
        printf("-----------------\n");

        gmp_printf("%Zd = \n", program->n);
        gmp_printf("\t  %Zd",        tdivfactors->data[0]);
        gmp_printf(" ^ %"PRIu32"\n", tdivmultis->data[0]);

        for (uint32_t i = 1; i < tdivfactors->length; i++) {
            gmp_printf("\t* %Zd",        tdivfactors->data[i]);
            gmp_printf(" ^ %"PRIu32"\n", tdivmultis->data[i]);
        }
        gmp_printf("\t* %Zd (probably prime)\n", n_atd);
        printf("\n");
    }

    clear_mpz_array(tdivfactors);
    free(tdivfactors);

    clear_uint32_array(tdivmultis);
    free(tdivmultis);

    mpz_clear(program->n);
    mpz_clear(n_atd);

    return ecode;
}
//------------------------------------------------------------------------------

