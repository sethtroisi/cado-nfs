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
 * \file    tifa_factor.h
 * \author  Jerome Milan
 * \date    Thu Dec 13 2007
 * \version 1.0
 *
 * \brief TIFA's generic factorization function.
 *
 * This is the TIFA library's generic factorization function: it picks the most 
 * suitable factoring algorithm depending on the size of the number to factor.
 */

#if !defined(_TIFA_TIFA_FACTOR_H_)
   /**
    * \def _TIFA_TIFA_FACTOR_H_
    * Standard include guard.
    */
#define _TIFA_TIFA_FACTOR_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <gmp.h>

#include "array.h"
#include "factoring_machine.h"
#include "exit_codes.h"

   /**
    * \brief Generic Integer factorization.
    *
    * Attempts to factor the non perfect square integer \c n with the
    * most suitable algorithm (chosen according to the size of \c n)
    * and with the factoring mode given by <tt>mode</tt>. Found factors are
    * then stored in \c factors. Additionally, if the factoring mode used is
    * set to \c FIND_COMPLETE_FACTORIZATION, factors' multiplicities are stored
    * in the array <tt>multis</tt>. For the time being, no trial divisions
    * are performed so depending on the situation, it could be worthwhile to
    * carry out such a step <i>before</i> calling \c tifa_factor.
    *
    * \note If the factoring mode used is different from
    * \c FIND_COMPLETE_FACTORIZATION, \c multis is allowed to be a NULL
    * pointer. Otherwise, using a NULL pointer will lead to a fatal error.
    *
    * \warning If the \c factors and \c multis arrays have not enough room
    * to store the found factors (and the multiplicities, if any), they will
    * be automatically resized to accommodate the data. This has to be kept
    * in mind when trying to do ingenious stuff with memory management (hint:
    * don't try to be clever here).
    *
    * \param[out] factors Pointer to the found factors of <tt>n</tt>.
    * \param[out] multis Pointer to the multiplicities of the found factors
    *                    (only computed if \c mode is set to
    *                     \c FIND_COMPLETE_FACTORIZATION).
    * \param[in] n The non perfect square integer to factor.
    * \param[in] mode The factoring mode to use.
    * \return An exit code.
    */
ecode_t tifa_factor(
    mpz_array_t* const factors,
    uint32_array_t* const multis,
    const mpz_t n,
    const factoring_mode_t mode
);

#ifdef __cplusplus
}
#endif

#endif

