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
 * \file    trial_div.c
 * \author  Jerome Milan
 * \date    Thu Mar 15 2007
 * \version 1.2
 */

 /*
  *  History:
  *
  *  1.2: Thu Mar 15 2007 by JM
  *       - Completely refactored to use a factoring_machine_t.
  *  1.1: Fri Jan 26 2007 by JM
  *       - Added MPN implementation with TIFA_USE_GMP_INTERNAL_FUNCS switch.
  *  1.0: Sat Mar 18 2006
  *       - Initial version by JM.
  */

#include <stdlib.h>

#include "tifa_config.h"
#include <gmp.h>
#if TIFA_USE_GMP_INTERNAL_FUNCS
    #include "gmp-impl.h"
#endif

#include "tdiv.h"
#include "first_primes.h"

#define __PREFIX__  "tdiv: "
#define __VERBOSE__ TIFA_VERBOSE_TDIV
#define __TIMING__  TIFA_TIMING_TDIV

#include <messages.h>

//------------------------------------------------------------------------------
//                   NON PUBLIC STRUCTURE(S) AND TYPEDEF(S)
//------------------------------------------------------------------------------

//
// An ad-hoc structure holding the data variables used by the trial division
// implementation.
//
struct struct_tdiv_context_t {
    //
    // Pointer to the unfactored part of the number to factor.
    //
    mpz_ptr n_atd;

    //
    // Number of primes used for trial divisions.
    //
    uint32_t nprimes;
};
//------------------------------------------------------------------------------
typedef struct struct_tdiv_context_t tdiv_context_t;
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                  "PRIVATE" FUNCTION(S) IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
static ecode_t init_tdiv_context(factoring_machine_t* const machine
                                 __attribute__ ((unused))) {
    PRINT_INIT_MSG;
    INIT_TIMER;
    START_TIMER;
    //
    // There is indeed nothing to do here. Timing measurements and output
    // messages are only kept to be consistent with the other factorization
    // functions.
    //
    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//------------------------------------------------------------------------------
static ecode_t clear_tdiv_context(factoring_machine_t* const machine
                                  __attribute__ ((unused))) {
    PRINT_CLEAN_MSG;
    INIT_TIMER;
    START_TIMER;
    //
    // Again, there is nothing to do here. Timing measurements and output
    // messages are only kept to be consistent with the other factorization
    // functions.
    //
    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//------------------------------------------------------------------------------
static ecode_t update_tdiv_context(factoring_machine_t* const machine
                                   __attribute__ ((unused))) {

    INIT_TIMER;
    START_TIMER;
    //
    // For the time being, just give up the factorization. Of course, we
    // could go on dividing by more primes, but would that be really helpful?
    //
    PRINT_UPDATE_GIVEUP_MSG;
    STOP_TIMER;
    PRINT_TIMING;

    return GIVING_UP;
}
//------------------------------------------------------------------------------
static ecode_t perform_tdiv(factoring_machine_t* const machine) {

    tdiv_context_t* context = (tdiv_context_t*) machine->context;

    if (context->nprimes == 0) {
        return NO_FACTOR_FOUND;
    }

    if (context->nprimes > NFIRST_PRIMES) {
        //
        // We are limited to the NFIRST_PRIMES precomputed primes and
        // cowardly refuse to generate more primes "on the fly".
        //
        context->nprimes = NFIRST_PRIMES;
    }
    mpz_t factor;
    mpz_init(factor);
    uint32_t nprimes = context->nprimes;

    INIT_TIMER;
    START_TIMER;
    PRINT_DIVIDING_MSG(nprimes);

#if TIFA_USE_GMP_INTERNAL_FUNCS
    //
    // Use some internal GMP functions.
    //
    // _NOTE_: We actually use the mpn_divexact_1 function which is not
    //         documented. See:
    //
    //      http://gmplib.org/list-archives/gmp-bugs/2006-December/000660.html
    //
    //         Apparently this "secret GMP function" (as Torbjorn Granlund put
    //         it) is not supposed to be used by client code... Could it crash
    //         your computer? Burn down your office? Give bad results? Or, more
    //         prosaically, will it be removed in future versions? In any case,
    //         maintainers should be aware of this...
    //
    mp_limb_t* num  = PTR(context->n_atd);
    uint32_t   size = ABSIZ(context->n_atd);
    uint32_t   pi   = 0;

    //
    // Divisions by 2
    //
    if ( (num[0] & 1) == 0) {

        mpz_setbit(factor, 1);

        append_mpz_to_array(machine->factors, factor);
        append_uint32_to_array(machine->multis, 1);

        mpn_rshift(num, num, size, 1);
        MPN_NORMALIZE(num, size);
        SIZ(context->n_atd) = size;

        while ( (num[0] & 1) == 0) {
            mpn_rshift(num, num, size, 1);
            MPN_NORMALIZE(num, size);
            SIZ(context->n_atd) = size;

            machine->multis->data[machine->multis->length - 1]++;
        }
    }
    //
    // Divisions by odd primes:
    //
    // _WARNING_: The following loop uses the "internal" mpn_modexact_1_odd
    //            and mpn_divexact_1 functions which could potentially be
    //            removed or changed in future GMP versions. See previous
    //            warning note.
    //
    for (pi = 1; pi < nprimes; pi++) {

        if (mpn_modexact_1_odd(num, size, first_primes[pi]) == 0) {

            mpn_divexact_1(num, num, size, first_primes[pi]);
            MPN_NORMALIZE(num, size);
            SIZ(context->n_atd) = size;

            mpz_set_ui(factor, first_primes[pi]);

            append_mpz_to_array(machine->factors, factor);
            append_uint32_to_array(machine->multis, 1);

            while (mpn_modexact_1_odd(num, size, first_primes[pi]) == 0) {
                mpn_divexact_1(num, num, size, first_primes[pi]);
                MPN_NORMALIZE(num, size);
                SIZ(context->n_atd) = size;
                machine->multis->data[machine->multis->length - 1]++;
            }
        }
    }

#else

    //
    // GMP's internal functions are not available: revert to a simple
    // MPZ based implementation. (Using "public" MPN function doesn't bring any
    // significant benefits according to my tests, at least on Opteron)
    //
    for (uint32_t i = 0; i < nprimes; i++) {
        if (mpz_divisible_ui_p(context->n_atd, first_primes[i])) {

            mpz_set_ui(factor, first_primes[i]);
            mpz_divexact_ui(context->n_atd, context->n_atd, first_primes[i]);

            append_mpz_to_array(machine->factors, factor);
            append_uint32_to_array(machine->multis, 1);

            while (mpz_divisible_ui_p(context->n_atd, first_primes[i])) {
                mpz_divexact_ui(context->n_atd, context->n_atd,first_primes[i]);
                machine->multis->data[machine->multis->length - 1]++;
            }
        }
    }
#endif

    mpz_clear(factor);

    STOP_TIMER;
    PRINT_TIMING;

    ecode_t ecode;
    if (machine->factors->length != 0) {
        ecode = SOME_FACTORS_FOUND;
    } else {
        ecode = NO_FACTOR_FOUND;
    }
    return ecode;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                 "PUBLIC" FUNCTION(S) IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
ecode_t tdiv(mpz_t n_atd, mpz_array_t* const factors,
             uint32_array_t* const multis, const mpz_t n,
             const uint32_t nprimes) {

    PRINT_INTRO_MSG("trial division");
    INIT_TIMER;
    START_TIMER;

    ecode_t ecode = NO_FACTOR_FOUND;

    factoring_machine_t machine;

    tdiv_context_t context;

    context.n_atd   = n_atd;
    context.nprimes = nprimes;

    mpz_init_set(machine.n, n);
    mpz_set(context.n_atd, machine.n);

    machine.mode                = FIND_SOME_FACTORS;
    machine.params              = NULL;
    machine.context             = (void*)&context;
    machine.init_context_func   = init_tdiv_context;
    machine.perform_algo_func   = perform_tdiv;
    machine.update_context_func = update_tdiv_context;
    machine.clear_context_func  = clear_tdiv_context;
    machine.factors             = factors;
    machine.multis              = multis;

    ecode = run_machine(&machine);

    mpz_clear(machine.n);

    STOP_TIMER;
    PRINT_STATUS(machine.success, ecode);

    return ecode;
}
//------------------------------------------------------------------------------

