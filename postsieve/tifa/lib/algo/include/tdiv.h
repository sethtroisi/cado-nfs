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
 * \file    tdiv.h
 * \author  Jerome Milan
 * \date    Sat Mar 18 2006
 * \version 1.0
 *
 * \brief The trial division factorization algorithm.
 *
 * Naive (partial) factorization via trial divisions by a few small primes.
 */

 /*
  *  Copyright (C) 2006, 2007  INRIA
  *  Licence: GNU Lesser General Public License (LGPL)
  *  History:
  */

#if !defined(_TIFA_TDIV_H_)
   /**
    * \def _TIFA_TDIV_H_
    * Standard include guard.
    */
#define _TIFA_TDIV_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <inttypes.h>
#include <gmp.h>

#include "array.h"
#include "factoring_machine.h"
#include "exit_codes.h"

   /**
    * \def TDIV_DFLT_NPRIMES_TDIV
    * Default number of the first primes to use for trial division.
    */
#define TDIV_DFLT_NPRIMES_TDIV     (NFIRST_PRIMES/32)
   
   /**
    * \brief Integer factorization via trial division (TDIV).
    *
    * Attempts to factor the integer \c n via trial division by the first
    * \c nprimes primes. Found factors are then stored in the array \c factors 
    * and multiplicities are stored in <tt>multis</tt>.
    * 
    * Returns:    
    * \li \c COMPLETE_FACTORIZATION_FOUND if the complete factorization of
    *                                     \c n was found.
    *
    * \li \c SOME_FACTORS_FOUND if some factors were found but could not
    *                           account for the complete factorization of \c n.
    *                           In that case, the unfactored part of \c n is
    *                           stored in
    *                           \c factors->data[\c factors->lenth - 1].
    *
    * \li \c NO_FACTOR_FOUND if no factor were found.
    *
    * \warning The prime numbers are not computed but read from a table. 
    * Consequently the number of primes \c nprimes should be less than or equal 
    * to \c NFIRST_PRIMES (defined in \link array.h \endlink). If \c nprimes is
    * zero, then the default value \c DFLT_TDIV_NPRIMES will be used instead.
    *
    * \param[out] factors Pointer to the found factors of <tt>n</tt>.
    * \param[out] multis  Pointer to the multiplicities of the found factors.
    * \param[in]  n       The integer to factor.
    * \param[in]  nprimes The number of primes to trial divide \c n by.
    * \return     An exit code.
    */
ecode_t tdiv(
    mpz_array_t* const factors,
    uint32_array_t* const multis,
    const mpz_t n,
    const uint32_t nprimes
);

#ifdef __cplusplus
}
#endif

#endif
