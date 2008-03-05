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
 * \file    trial_div.c
 * \author  Jerome Milan
 * \date    Tue Mar 4 2008
 * \version 1.3
 */

 /*
  * History:
  * --------
  *  1.3: Tue Mar  4 2008 by JM:
  *       - Rewritten and simplified to get rid of the factoring_machine_t.
  *  1.2: Thu Mar 15 2007 by JM
  *       - Completely refactored to use a factoring_machine_t.
  *  1.1: Fri Jan 26 2007 by JM
  *       - Added MPN implementation with TIFA_USE_GMP_INTERNAL_FUNCS switch.
  *  1.0: Sat Mar 18 2006
  *       - Initial version by JM.
  */

#include <stdlib.h>
#include <stdbool.h>
#include <gmp.h>

#include "tifa_config.h"
#if TIFA_USE_GMP_INTERNAL_FUNCS
    #include "gmp-impl.h"
#endif

#include "first_primes.h"
#include "factoring_machine.h"
#include "macros.h"
#include "tdiv.h"

#define __PREFIX__  "tdiv: "
#define __VERBOSE__ TIFA_VERBOSE_TDIV
#define __TIMING__  TIFA_TIMING_TDIV

#include "messages.h"

//------------------------------------------------------------------------------
//                 "PUBLIC" FUNCTION(S) IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
ecode_t tdiv(mpz_array_t* const factors, uint32_array_t* const multis,
             const mpz_t n, const uint32_t nprimes_tdiv) {
    
    //
    // _WARNING_: The tdiv function does not have the same signature than the 
    //            other factoring functions.
    //
    PRINT_INTRO_MSG("trial division");
    INIT_TIMER;
    START_TIMER;
    
    bool have_found_factors = false;
    
    uint32_t nprimes = nprimes_tdiv;
    
    if (nprimes == 0) {
        nprimes = TDIV_DFLT_NPRIMES_TDIV;
    }
    if (nprimes > NFIRST_PRIMES) {
        //
        // We are limited to the NFIRST_PRIMES precomputed primes and
        // cowardly refuse to generate more primes "on the fly".
        //
        nprimes = NFIRST_PRIMES;
    }
    PRINT_DIVIDING_MSG(nprimes);
    
    mpz_t factor;
    mpz_init(factor);
    
    mpz_t n_atd;            // stands for "n After Trial Division"
    mpz_init_set(n_atd, n);
    
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
    mp_limb_t* num  = PTR(n_atd);
    uint32_t   size = ABSIZ(n_atd);
    uint32_t   pi   = 0;

    //
    // Divisions by 2
    //
    if (IS_EVEN(num[0])) {

        mpz_setbit(factor, 1);

        append_mpz_to_array(factors, factor);
        append_uint32_to_array(multis, 1);
        
        mpn_rshift(num, num, size, 1);
        MPN_NORMALIZE(num, size);
        SIZ(n_atd) = size;

        while (IS_EVEN(num[0])) {
            mpn_rshift(num, num, size, 1);
            MPN_NORMALIZE(num, size);
            SIZ(n_atd) = size;

            multis->data[multis->length - 1]++;
        }
        have_found_factors = true;
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

        if (mpz_cmp_ui(n_atd, 1) == 0) {
            break;
        }
        
        if (mpn_modexact_1_odd(num, size, first_primes[pi]) == 0) {

            mpn_divexact_1(num, num, size, first_primes[pi]);
            MPN_NORMALIZE(num, size);
            SIZ(n_atd) = size;

            mpz_set_ui(factor, first_primes[pi]);

            append_mpz_to_array(factors, factor);
            append_uint32_to_array(multis, 1);
            
            while (mpn_modexact_1_odd(num, size, first_primes[pi]) == 0) {
                mpn_divexact_1(num, num, size, first_primes[pi]);
                MPN_NORMALIZE(num, size);
                SIZ(n_atd) = size;
                multis->data[multis->length - 1]++;
            }
            have_found_factors = true;
        }
    }
#else
    //
    // GMP's internal functions are not available: revert to a simple
    // MPZ based implementation. (Using "public" MPN function doesn't bring any
    // significant benefits according to my tests, at least on Opteron)
    //
    for (uint32_t i = 0; i < nprimes; i++) {
        if (mpz_cmp_ui(n_atd, 1) == 0) {
            break;
        }
        if (mpz_divisible_ui_p(n_atd, first_primes[i])) {

            mpz_set_ui(factor, first_primes[i]);
            mpz_divexact_ui(n_atd, n_atd, first_primes[i]);

            append_mpz_to_array(factors, factor);
            append_uint32_to_array(multis, 1);
            
            while (mpz_divisible_ui_p(n_atd, first_primes[i])) {
                mpz_divexact_ui(n_atd, n_atd,first_primes[i]);
                multis->data[multis->length - 1]++;
            }
            have_found_factors = true;
        }
    }
#endif
   
    PRINT_TIMING;
    
    ecode_t ecode   = NO_FACTOR_FOUND;
    bool    success = false;
    
    if (have_found_factors) {
        if (mpz_cmp_ui(n_atd, 1) == 0) {
            ecode   = COMPLETE_FACTORIZATION_FOUND;
            success = true;
        } else {
            append_mpz_to_array(factors, n_atd);
            append_uint32_to_array(multis, 1);
            if (MPZ_IS_PRIME(n_atd)) {
                ecode   = COMPLETE_FACTORIZATION_FOUND;
                success = true;
            } else {
                ecode = SOME_FACTORS_FOUND;
            }
        }
    }
    mpz_clear(factor);
    mpz_clear(n_atd);
    
    STOP_TIMER;
    PRINT_STATUS(success, ecode);

    return ecode;
}
//------------------------------------------------------------------------------
