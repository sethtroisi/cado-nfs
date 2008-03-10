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
 * \file    factoring_machine.c
 * \author  Jerome Milan
 * \date    Fri Jan 11 2008
 * \version 0.2
 */

 /*
  *  History:
  *    0.2: Fri Jan 11 2008 by JM:
  *         - Rewrote function find_complete_factorization(...) (was 
  *           complicated, obscure and completely buggy). 
  *    0.1: March 2007 by JM:
  *         - Initial version.
  */

#include <stdlib.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>

#include "macros.h"
#include "funcs.h"
#include "factoring_machine.h"

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
ecode_t compute_multiplicities(
    uint32_array_t* const multis,
    const mpz_t n,
    const mpz_array_t* const factors
);
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                 "PUBLIC" FUNCTION(S) IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
ecode_t run_machine(factoring_machine_t* machine) {
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
    ecode_t exit_code = NO_FACTOR_FOUND;
    machine->success  = false;
    
    exit_code = machine->init_context_func(machine);
    exit_code = machine->perform_algo_func(machine);
    
    while (true) {
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

    if (exit_code != SOME_FACTORS_FOUND) {
        return exit_code;
    }
    //
    // We have found at least one factor. Compute a coprime base for all the
    // found factors, then compute the multiplicity of each element in said
    // base (even if they are not all prime).
    //
    uint32_t     flen = machine->factors->length;
    mpz_array_t* base = alloc_mpz_array(flen * (flen + 1));
        
    find_coprime_base(base, machine->n, machine->factors);
    swap_mpz_array(machine->factors, base);
    
    clear_mpz_array(base);
    
    exit_code = compute_multiplicities(
                    machine->multis,
                    machine->n,
                    machine->factors
                );
    
    if (exit_code != SUCCESS) {
        return exit_code;
    }
    
    if (are_all_prime(machine->factors)) {
        //
        // The found factors stored in machine->factors are all prime!
        // This means that we have found the complete factorization of
        // machine->n.
        //
        machine->success = true;
        return COMPLETE_FACTORIZATION_FOUND;
    }
    
    //
    // We have a coprime base together with the multiplicity of each of its
    // element. However, some of these factors are not prime so we should
    // factorize them.
    //
    uint32_t alloced = 2 * machine->factors->length;
    
    //
    // The all_* arrays will hold all the (hopefully prime) factors of
    // machine->n and their multiplicities.
    //
    mpz_array_t*    all_factors = alloc_mpz_array(alloced);
    uint32_array_t* all_multis  = alloc_uint32_array(alloced);

    bool factorization_is_partial = false;

    for (uint32_t i = 0; i < machine->factors->length; i++) {
        
        const mpz_ptr curfactor = machine->factors->data[i];
        const uint32_t curmulti = machine->multis->data[i];
        
        if (MPZ_IS_PRIME(curfactor)) {
            
            append_mpz_to_array(all_factors, curfactor);
            append_uint32_to_array(all_multis, curmulti);

        } else {
            //
            // Proceed to factorize this non prime factor...
            //
            mpz_array_t*   morefactors = alloc_mpz_array(ELONGATION);
            uint32_array_t* moremultis = alloc_uint32_array(ELONGATION);

            ecode_t retcode = machine->recurse_func(
                                  morefactors,
                                  moremultis,
                                  curfactor,
                                  FIND_COMPLETE_FACTORIZATION
                              );

            switch (retcode) {

            case COMPLETE_FACTORIZATION_FOUND:
                for (uint32_t j = 0; j < moremultis->length; j++) {
                    moremultis->data[j] *= curmulti;
                }
                append_mpz_array(all_factors, morefactors);
                append_uint32_array(all_multis, moremultis);
                break;

            case PARTIAL_FACTORIZATION_FOUND:
                factorization_is_partial = true;
                
                for (uint32_t j = 0; j < moremultis->length; j++) {
                    moremultis->data[j] *= machine->multis->data[i];
                }
                append_mpz_array(all_factors, morefactors);
                append_uint32_array(all_multis, moremultis);
                break;

            case FATAL_INTERNAL_ERROR:
            case NO_FACTOR_FOUND:
            case GIVING_UP:
            default:
                //
                // For the time being we just return with the found
                // factors and give up. Ideally, we could try different
                // factoring algorithms here before resigning.
                //
                factorization_is_partial = true;
                append_mpz_to_array(all_factors, curfactor);
                append_uint32_to_array(all_multis, curmulti);
                break;
            }
            clear_mpz_array(morefactors);
            clear_uint32_array(moremultis);
        }
    }
    swap_mpz_array(machine->factors, all_factors);
    swap_uint32_array(machine->multis, all_multis);

    clear_mpz_array(all_factors);
    clear_uint32_array(all_multis);

    if (factorization_is_partial) {
        return PARTIAL_FACTORIZATION_FOUND;
    } else {
        machine->success = true;
        return COMPLETE_FACTORIZATION_FOUND;
    }
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
        if (!MPZ_IS_PRIME(array->data[i])) {
            all_prime = false;
            break;
        }
    }
    return all_prime;
}
//-----------------------------------------------------------------------------
ecode_t compute_multiplicities(uint32_array_t* const multis,
                               const mpz_t n,
                               const mpz_array_t* const factors) {
    //
    // Computes the multiplicities of the factors of 'n' given by the coprime
    // base 'factors' and stores them in the 'multis' array.
    //
    mpz_t cofactor;
    mpz_init_set(cofactor, n);

    multis->length = 0;

    for (uint32_t i = 0; i < factors->length; i++) {
        //
        // _NOTE_: We already know that 'n' is divisible by each element in
        //         the array 'factors' _and_ we known that these factors are
        //         coprime with each other (i.e. 'factors' is a coprime base)
        //         so no divisibility test is required for the first division.
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
    if (mpz_cmp_ui(cofactor, 1) != 0) {
        //
        // This should not happen since 'n' is supposed to be
        // smooth of the coprime base given by 'factors'! We could return
        // the found factors and multiplities but that would mask a potential 
        // nasty bug in the code, so let's be convervative and return an error.
        //
        mpz_clear(cofactor);
        return FATAL_INTERNAL_ERROR;
    }
    mpz_clear(cofactor);
    return SUCCESS;
}
//-----------------------------------------------------------------------------
