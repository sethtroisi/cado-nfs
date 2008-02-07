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
 * \file    tifa_factor.c
 * \author  Jerome Milan
 * \date    Thu Dec 13 2007
 * \version 1.0
 */

#include <stdlib.h>

#include "tifa_config.h"
#include "macros.h"
#include "siqs.h"
#include "squfof.h"
#include "tifa_factor.h"

#include <stdio.h>

//-----------------------------------------------------------------------------
ecode_t tifa_factor(mpz_array_t* const factors, uint32_array_t* const multis,
                    const mpz_t n, const factoring_mode_t mode) {
    
    //
    // First of all, let's check that we're not dealing with a perfect square
    // which can be problematic for some algorithm (i.e. SQUFOF). Note that
    // we can not (conveniently) rely on the client code to perform such a
    // check since tifa_factor(...) is indirectly and recursively called with
    // a factoring_machine_t when using the FIND_*_FACTORIZATION modes.
    //
    if (MPZ_IS_SQUARE(n)) {
        mpz_t rootn;
        mpz_init(rootn);
        mpz_sqrt(rootn, n);

        if (MPZ_IS_PRIME(rootn)) {
            append_mpz_to_array(factors, rootn);
            if (mode == FIND_COMPLETE_FACTORIZATION) {   
                append_uint32_to_array(multis, 2);
            }
            //
            // We have found the complete factorization but the exit code to
            // return actually depends on the factorization mode used.
            //
            return mode_to_outcome[mode];
        }
        uint32_t old_length = 0;
        
        if (mode == FIND_COMPLETE_FACTORIZATION) {
            old_length = multis->length;
        }
        ecode_t ecode = tifa_factor(factors, multis, rootn, mode);
        
        if (mode == FIND_COMPLETE_FACTORIZATION) {
            for (uint32_t i = old_length; i < multis->length; i++) {
                multis->data[i] <<= 1;
            }
        }
        mpz_clear(rootn);
        
        return ecode;
    }
    unsigned long int sizen = mpz_sizeinbase(n, 2);

#if TIFA_WORDSIZE < 64
    //
    // On 32 bit machines TIFA's single precision SQUFOF will start to fail
    // for 59-ish bit numbers.
    //
    unsigned long int squfof_threshold = 58;
#else
    //
    // SIQS starts to be better than SQUFOF for 64 bit numbers.
    //
    unsigned long int squfof_threshold = 64;
#endif

    //
    // _NOTE_: The following choice of algorithm has been obtained by benching
    //         all TIFA's factoring algorithms with RSA moduli (i.e. numbers
    //         that factors in two primes of equal size). For more general
    //         numbers, things can be somewhat different, particularly for 
    //         SQUFOF whose theoretical complexity directly depends on the
    //         number of prime factors.
    //
    // _TO_DO_: For the time being, no trial division step is performed so
    //          depending on the situation, this should be done before calling
    //          tifa_factor. Should we add a trial division step here?
    //
    if (sizen <= squfof_threshold) {
        squfof_params_t params;
        set_squfof_params_to_default(&params);

        return squfof(factors, multis, n, &params, mode);
    }
    siqs_params_t params;
    set_siqs_params_to_default(n, &params);
    
    return siqs(factors, multis, n, &params, mode);
}
//-----------------------------------------------------------------------------

