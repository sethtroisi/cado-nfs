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
 * \file    factoring_machine.c
 * \author  Jerome Milan
 * \date    March 2007
 * \version 0.1
 */

#include <stdlib.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>

#include "factoring_machine.h"
#include "funcs.h"

//------------------------------------------------------------------------------
//                         NON PUBLIC DEFINE(S)
//------------------------------------------------------------------------------

//
// Number of Miller-Rabin iterations to perform for each compositeness test.
//
#define NMILLER_RABIN       32
//
// Mere shortcut to GMP's mpz_probab_prime_p
//
#define IS_PRIME(X) (0 != mpz_probab_prime_p((X), NMILLER_RABIN))
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                 PROTOTYPES OF NON PUBLIC FUNCTION(S)
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
ecode_t single_run(factoring_machine_t*);
//------------------------------------------------------------------------------
ecode_t find_some_factors(factoring_machine_t*);
//------------------------------------------------------------------------------
ecode_t find_some_coprime_factors(factoring_machine_t*);
//------------------------------------------------------------------------------
ecode_t find_some_prime_factors(factoring_machine_t*);
//------------------------------------------------------------------------------
ecode_t find_complete_factorization(factoring_machine_t*);
//-----------------------------------------------------------------------------
bool are_all_prime(const mpz_array_t* const);
//------------------------------------------------------------------------------
void compute_cofactor(
    mpz_t, uint32_array_t* const,
    const mpz_t, const mpz_array_t* const
);
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                 "PUBLIC" FUNCTION(S) IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
inline ecode_t run_machine(factoring_machine_t* machine) {
    switch (machine->mode) {
        case SINGLE_RUN:
            return single_run(machine);
            break;
        case FIND_SOME_FACTORS:
            return find_some_factors(machine);
            break;
        case FIND_SOME_PRIME_FACTORS:
            return find_some_prime_factors(machine);
            break;
        case FIND_SOME_COPRIME_FACTORS:
            return find_some_coprime_factors(machine);
            break;
        case FIND_COMPLETE_FACTORIZATION:
            return find_complete_factorization(machine);
            break;
        default:
            return UNKNOWN_FACTORING_MODE;
    }
}
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                  "PRIVATE" FUNCTION(S) IMPLEMENTATION
//------------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                  Functions directly used by run_machine
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
ecode_t single_run(factoring_machine_t* machine) {
    //
    // Perform a single run of the factoring algorithm.
    // Don't try to do anything "smartish" if the algorithm fails to find a
    // factor...
    //
    ecode_t exit_code;

    exit_code = machine->init_context_func(machine);
    exit_code = machine->perform_algo_func(machine);

    if (exit_code == SOME_FACTORS_FOUND) {
        machine->success = true;
    } else {
        machine->success = false;
    }
    (void) machine->clear_context_func(machine);

    return exit_code;
}
//-----------------------------------------------------------------------------
ecode_t find_some_factors(factoring_machine_t* machine) {
    //
    // Try to find a least "some factor(s)" even if they are not prime.
    //
    ecode_t exit_code = NO_FACTOR_FOUND;

    exit_code = machine->init_context_func(machine);
    exit_code = machine->perform_algo_func(machine);

    while (true) {
        //
        // Try again and again until we either:
        //    - find some factor(s)
        //    - reach the abort limit
        //    - or end up with a really bad error...
        //
        if (exit_code == SOME_FACTORS_FOUND) {
            break;
        }
        if (exit_code == FATAL_INTERNAL_ERROR) {
            break;
        }
        exit_code = machine->update_context_func(machine);

        if (exit_code != SUCCESS) {
            break;
        }
        exit_code = machine->perform_algo_func(machine);
    }
    machine->clear_context_func(machine);

    if (exit_code == SOME_FACTORS_FOUND) {
        machine->success = true;
    } else {
        machine->success = false;
    }
    return exit_code;
}
//-----------------------------------------------------------------------------
ecode_t find_some_coprime_factors(factoring_machine_t* machine) {
    //
    // Try to find "some coprime factors" but don't try to find _all_ prime
    // factors.
    //
    // Compared to find_some_factors(), this involves some processing on the
    // found factors to obtain a coprime base. Actually there is no guarantee
    // that the found factors are prime.
    //
    ecode_t exit_code = NO_FACTOR_FOUND;

    exit_code = machine->init_context_func(machine);
    exit_code = machine->perform_algo_func(machine);

    while (true) {
        //
        // Try again and again until we either:
        //    - find some factor(s)
        //    - reach the abort limit
        //    - or end up with a really bad error...
        //
        if (exit_code == SOME_FACTORS_FOUND) {
            break;
        }
        if (exit_code == FATAL_INTERNAL_ERROR) {
            break;
        }
        exit_code = machine->update_context_func(machine);

        if (exit_code != SUCCESS) {
            break;
        }
        exit_code = machine->perform_algo_func(machine);
    }
    machine->clear_context_func(machine);

    if (exit_code == SOME_FACTORS_FOUND) {
        //
        // Compute a coprime base for all the found factors.
        //
        uint32_t     flen = machine->factors->length;
        mpz_array_t* base = alloc_mpz_array(flen * (flen + 1));

        find_coprime_base(base, machine->n, machine->factors);
        swap_mpz_array(machine->factors, base);

        clear_mpz_array(base);
        free(base);

        exit_code = SOME_COPRIME_FACTORS_FOUND;
        machine->success = true;
    } else {
        machine->success = false;
    }
    return exit_code;
}
//-----------------------------------------------------------------------------
ecode_t find_some_prime_factors(factoring_machine_t* machine) {
    //
    // Try to find "some prime factors" but don't try to find _all_ prime
    // factors.
    //
    // Compared to find_some_factors(), this involves some processing on the
    // found factors (to obtain a coprime base) and some compositeness checks.
    // Note that we return FAILURE as soon as one of the found factor is NOT
    // prime. SUCCESS is therefore restricted to the case where all found
    // factors (after processing) are prime. It is debatable whether or not
    // this is a suitable behaviour. Indeed we could return SUCCESS as soon as
    // _one_ of the found factors is prime, which would actually better match
    // the name of the function. But then, what would be the point of such a
    // factoring mode? Heck, I'm "not even sure" (to put it mildly) that it
    // is useful at all the way it is...
    //
    ecode_t exit_code = NO_FACTOR_FOUND;

    exit_code = machine->init_context_func(machine);
    exit_code = machine->perform_algo_func(machine);

    while (true) {
        //
        // Try again and again until we either:
        //    - find some factor(s)
        //    - reach the abort limit
        //    - or end up with a really bad error...
        //
        if (exit_code == SOME_FACTORS_FOUND) {
            break;
        }
        if (exit_code == FATAL_INTERNAL_ERROR) {
            break;
        }
        exit_code = machine->update_context_func(machine);

        if (exit_code != SUCCESS) {
            break;
        }
        exit_code = machine->perform_algo_func(machine);
    }
    machine->clear_context_func(machine);

    if (exit_code == SOME_FACTORS_FOUND) {
        //
        // Compute a coprime base for all the found factors.
        //
        uint32_t     flen = machine->factors->length;
        mpz_array_t* base = alloc_mpz_array(flen * (flen + 1));

        find_coprime_base(base, machine->n, machine->factors);
        swap_mpz_array(machine->factors, base);

        clear_mpz_array(base);
        free(base);

        if (are_all_prime(machine->factors)) {
            exit_code = SOME_PRIME_FACTORS_FOUND;
            machine->success = true;
        } else {
            //
            // _NOTE_: We fail as soon as one element in the coprime base is
            //         not prime. Is that a sensible behavior?
            //
            exit_code = SOME_COPRIME_FACTORS_FOUND;
            machine->success = false;
        }
    } else {
        machine->success = false;
    }
    return exit_code;
}
//-----------------------------------------------------------------------------
ecode_t find_complete_factorization(factoring_machine_t* machine) {
    //
    // Try to find a "complete factorization", that is try to find _all_ prime
    // factors together with their multiplicities.
    //
    // _NOTE_: Things are getting quite ugly here. The code looks opaque
    //         and the overall strategy certainly could be improved.
    //
    // _NOTE_: This function is actually written in a finate state machine
    //         style because earlier (read "even buggier"!) versions followed
    //         a much more complicated code flow. After some clean-ups (and
    //         concessions!) things turned out to be simpler than previously
    //         thought so using the FSM nomenclature is now completely overkill
    //         and certainly degrades the readability of the function...
    //
    // _TO_DO_: Just read the previous note to infer what's to do! :-)
    //
    machine_state_t state = INITIAL_STATE;
    ecode_t exit_code     = NO_FACTOR_FOUND;
    ecode_t retcode       = FAILURE;

    while (state != FINAL_STATE) {

        switch (state) {

        case INITIAL_STATE:
            retcode = machine->init_context_func(machine);
            state = RUN_ALGO_STATE;
            break;

        case RUN_ALGO_STATE:
            retcode = machine->perform_algo_func(machine);
            
            switch (retcode) {
            
            case SOME_FACTORS_FOUND:
                state = TEST_PRIMALITY_STATE;
                break;
            
            case FATAL_INTERNAL_ERROR:
                state     = FAILURE_STATE;
                exit_code = FATAL_INTERNAL_ERROR;
                break;
            
            default:
                state = UPDATE_CONTEXT_STATE;
            }
            break;

        case UPDATE_CONTEXT_STATE:
            retcode = machine->update_context_func(machine);
            if (retcode != SUCCESS) {
                state = FAILURE_STATE;
            } else {
                state = RUN_ALGO_STATE;
            }
            break;

        case TEST_PRIMALITY_STATE: {
            //
            // Compute a coprime base for all the found factors, then check
            // the primality of each element in said base.
            //
            uint32_t     flen = machine->factors->length;
            mpz_array_t* base = alloc_mpz_array(flen * (flen + 1));

            find_coprime_base(base, machine->n, machine->factors);

            swap_mpz_array(machine->factors, base);

            clear_mpz_array(base);
            free(base);

            if (are_all_prime(machine->factors)) {
                state = TEST_FACTORIZATION_STATE;
            } else {
                state = RECURSIVE_FACTORIZATION_STATE;
            }
            break;
        }

        case TEST_FACTORIZATION_STATE: {
            //
            // The found factors stored in machine->factors are all prime!
            // Check if these factors are enough to obtain the complete
            // factorization.
            //
            mpz_t cofactor;
            mpz_init(cofactor);

            compute_cofactor(
                cofactor,
                machine->multis,
                machine->n,
                machine->factors
            );

            if (mpz_cmp_ui(cofactor, 1) == 0) {
                state     = SUCCESS_STATE;
                exit_code = COMPLETE_FACTORIZATION_FOUND;

                mpz_clear(cofactor);
                break;
            }
            if (IS_PRIME(cofactor)) {

                append_mpz_to_array(machine->factors, cofactor);
                append_uint32_to_array(machine->multis, 1);

                state     = SUCCESS_STATE;
                exit_code = COMPLETE_FACTORIZATION_FOUND;

                mpz_clear(cofactor);
                break;
            }
            //
            // If we reach this point then cofactor is different from 1 and
            // it is not a prime. So let's factorize it.
            //
            mpz_array_t*    morefactors = alloc_mpz_array(ELONGATION);
            uint32_array_t* moremultis  = alloc_uint32_array(ELONGATION);

            retcode = machine->recurse_func(
                          morefactors,
                          moremultis,
                          cofactor,
                          machine->mode
                      );

            switch (retcode) {

            case COMPLETE_FACTORIZATION_FOUND:
                //
                // Since machine->factors holds only prime factors and since
                // we know that these factors cannot divide cofactor, we
                // can just concatenate our result arrays to obtain the
                // complete factorization.
                //
                append_mpz_array(machine->factors, morefactors);
                append_uint32_array(machine->multis, moremultis);
                state     = SUCCESS_STATE;
                exit_code = COMPLETE_FACTORIZATION_FOUND;
                break;

            case PARTIAL_FACTORIZATION_FOUND:
                //
                // Since machine->factors holds only prime factors and since
                // we know that these factors cannot divide cofactor, we
                // can just concatenate our result arrays. However, we do
                // not have the complete factorization...
                //
                append_mpz_array(machine->factors, morefactors);
                append_uint32_array(machine->multis, moremultis);
                state     = FAILURE_STATE;
                exit_code = PARTIAL_FACTORIZATION_FOUND;
                break;
                
            case FATAL_INTERNAL_ERROR:
            case NO_FACTOR_FOUND:
            case GIVING_UP:
                //
                // What should we do here? For the time being we just return
                // with the found factors and give up. Ideally, we could
                // try different factoring algorithms here before resigning.
                //
                append_mpz_to_array(machine->factors, cofactor);
                append_uint32_to_array(machine->multis, 1);
                state     = FAILURE_STATE;
                exit_code = PARTIAL_FACTORIZATION_FOUND;
                break;

            default:
                //
                // This switch case cannot be accessed... presumably...
                //
                append_mpz_to_array(machine->factors, cofactor);
                append_uint32_to_array(machine->multis, 1);
                state     = FAILURE_STATE;
                exit_code = PARTIAL_FACTORIZATION_FOUND;
                break;
            }

            clear_mpz_array(morefactors);
            clear_uint32_array(moremultis);

            free(morefactors);
            free(moremultis);

            mpz_clear(cofactor);
            break;
        }

        case RECURSIVE_FACTORIZATION_STATE: {
            //
            // Some coprime factors were found. They are, however, not
            // necessarily prime, so we should factor the non prime ones.
            //
            uint32_t alloced = 2 * machine->factors->length;

            mpz_array_t*    all_factors = alloc_mpz_array(alloced);
            uint32_array_t* all_multis  = alloc_uint32_array(alloced);

            bool factorization_is_partial = false;

            for (uint32_t i = 0; i < machine->factors->length; i++) {

                if (IS_PRIME(machine->factors->data[i])) {

                    append_mpz_to_array(all_factors, machine->factors->data[i]);
                    append_uint32_to_array(
                        all_multis,
                        machine->multis->data[i]
                    );

                } else {
                    //
                    // Proceed to factorize the non prime factor...
                    //
                    mpz_array_t*   morefactors = alloc_mpz_array(ELONGATION);
                    uint32_array_t* moremultis = alloc_uint32_array(ELONGATION);

                    retcode = machine->recurse_func(
                                  morefactors,
                                  moremultis,
                                  machine->factors->data[i],
                                  machine->mode
                              );

                    switch (retcode) {

                    case COMPLETE_FACTORIZATION_FOUND:
                        append_mpz_array(machine->factors, morefactors);
                        append_uint32_array(machine->multis, moremultis);
                        break;

                    case PARTIAL_FACTORIZATION_FOUND:
                        factorization_is_partial = true;
                        append_mpz_array(machine->factors, morefactors);
                        append_uint32_array(machine->multis, moremultis);
                        break;

                    case FATAL_INTERNAL_ERROR:
                    case NO_FACTOR_FOUND:
                    case GIVING_UP:
                        //
                        // Same remark as in the test factorization state:
                        //
                        // For the time being we just return with the found
                        // factors and give up. Ideally, we could try different
                        // factoring algorithms here before resigning.
                        //
                    default:
                        factorization_is_partial = true;
                        append_mpz_to_array(
                            machine->factors,
                            machine->factors->data[i]
                        );
                        append_uint32_to_array(machine->multis, 1);

                        break;
                    }
                    clear_mpz_array(morefactors);
                    clear_uint32_array(moremultis);

                    free(morefactors);
                    free(moremultis);
                }
            }
            swap_mpz_array(machine->factors, all_factors);
            swap_uint32_array(machine->multis, all_multis);

            clear_mpz_array(all_factors);
            clear_uint32_array(all_multis);

            free(all_factors);
            free(all_multis);

            if (factorization_is_partial) {
                state     = FAILURE_STATE;
                exit_code = PARTIAL_FACTORIZATION_FOUND;
            } else {
                state     = SUCCESS_STATE;
                exit_code = COMPLETE_FACTORIZATION_FOUND;
            }
            break;
        }

        case SUCCESS_STATE:
            machine->success = true;
            state = CLEAN_STATE;
            break;

        case FAILURE_STATE:
            machine->success = false;
            state = CLEAN_STATE;
            break;

        case CLEAN_STATE:
            machine->clear_context_func(machine);
            state = FINAL_STATE;
            break;

        default:
            return FATAL_INTERNAL_ERROR;
        }
    }
    return exit_code;
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                       Various utility functions
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
bool are_all_prime(const mpz_array_t* const array) {
    //
    // Returns true if all elements in array are (probably) prime.
    // Returns false otherwise.
    //
    bool all_prime = true;
    for (uint32_t i = 0; i < array->length; i++) {
        if (!IS_PRIME(array->data[i])) {
            all_prime = false;
            break;
        }
    }
    return all_prime;
}
//-----------------------------------------------------------------------------
void compute_cofactor(mpz_t cofactor, uint32_array_t* const multis,
                      const mpz_t n, const mpz_array_t* const factors) {
    //
    // Computes the cofactor of 'n' obtained after trial division by its known
    // factors stored in the coprime base 'factors' and stores it in, you
    // guessed it, 'cofactor'.
    //
    mpz_set(cofactor, n);

    multis->length = 0;

    for (uint32_t i = 0; i < factors->length; i++) {
        //
        // _NOTE_: We already know that 'n' is divisible by each element in
        //         the array 'factors' _and_ we known that these factors are
        //         coprime with each other (i.e. 'factors' is a coprime base).
        //
        mpz_divexact(cofactor, cofactor, factors->data[i]);

        append_uint32_to_array(multis, 1);
    }
    if (mpz_cmp_ui(cofactor, 1) != 0) {
        for (uint32_t i = 0; i < factors->length; i++) {
            while (mpz_divisible_p(cofactor, factors->data[i])) {
                mpz_divexact(cofactor, cofactor, factors->data[i]);
                multis->data[i]++;
            }
            if (mpz_cmp_ui(cofactor, 1) == 0) {
                break;
            }
        }
    }
}
//-----------------------------------------------------------------------------
