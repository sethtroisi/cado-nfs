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
    * \brief Integer factorization via trial division (TDIV).
    *
    * Attempts to factor the integer \c n via trial division by
    * the first \c nprimes primes. Found factors are then stored in
    * \c factors and the remaining unfactored part of \c n is stored in
    * \c <tt>n_atd</tt>. Additionally, multiplicities are stored in
    * the array <tt>multis</tt>.
    *
    * \warning The prime numbers are not computed but read from a
    * table. Consequently the number of primes \c nprimes should be
    * less than or equal to <tt>NFIRST_PRIMES</tt> defined in
    * \link array.h \endlink .
    *
    * \param[out] n_atd   The unfactored part of \c n.
    * \param[out] factors Pointer to the found factors of <tt>n</tt>.
    * \param[out] multis  Pointer to the multiplicities of the found factors.
    * \param[in]  n       The integer to factor.
    * \param[in]  nprimes The number of primes to trial divide \c n by.
    */
ecode_t tdiv(
    mpz_t n_atd,
    mpz_array_t* const factors,
    uint32_array_t* const multis,
    const mpz_t n,
    const uint32_t nprimes
);

#ifdef __cplusplus
}
#endif

#endif
