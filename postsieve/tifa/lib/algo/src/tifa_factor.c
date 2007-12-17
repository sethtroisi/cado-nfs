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

#include "cfrac.h"
#include "qs.h"
#include "siqs.h"
#include "squfof.h"
#include "tdiv.h"

#include "tifa_factor.h"

//-----------------------------------------------------------------------------
ecode_t tifa_factor(mpz_array_t* const factors, uint32_array_t* const multis,
                    const mpz_t n, const factoring_mode_t mode) {
    
    unsigned long int sizen = mpz_sizeinbase(n, 2);
    
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
    //          tifa_factor. To do: add the trial division step here...
    //
    if (sizen <= 64) {
        squfof_params_t params;
        set_squfof_params_to_default(&params);
        
        return squfof(factors, multis, n, &params, mode);
    }
    siqs_params_t params;
    set_siqs_params_to_default(n, &params);
    
    return siqs(factors, multis, n, &params, mode);
}
//-----------------------------------------------------------------------------

