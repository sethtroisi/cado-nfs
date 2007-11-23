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
 * \file    siqs.h
 * \author  Jerome Milan
 * \date
 * \version 1.0
 *
 * \brief The Self-Initializing Quadratic Sieve factorization algorithm.
 */

#if !defined(_TIFA_SIQS_H_)
   /**
    * \def _TIFA_SIQS_H_
    * Standard include guard.
    */
#define _TIFA_SIQS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <inttypes.h>
#include <gmp.h>

#include "array.h"
#include "lindep.h"

#include "factoring_machine.h"
#include "exit_codes.h"

   /**
    * \def SIQS_DFLT_SIEVE_HALF_WIDTH
    * Default sieving half-width.
    */
#define SIQS_DFLT_SIEVE_HALF_WIDTH 500000

   /**
    * \def SIQS_DFLT_NPRIMES_IN_BASE
    * Default number of prime numbers composing the factor base on which
    * to factor the residues.
    */
#define SIQS_DFLT_NPRIMES_IN_BASE  NFIRST_PRIMES/16
   /**
    * \def SIQS_DFLT_NPRIMES_TDIV
    * Default number of the first primes to use in the trial division
    * of the residues.
    */
#define SIQS_DFLT_NPRIMES_TDIV     NFIRST_PRIMES/16
   /**
    * \def SIQS_DFLT_NRELATIONS
    * Default number of congruence relations to find before attempting the
    * factorization of the large integer.
    */
#define SIQS_DFLT_NRELATIONS       24
   /**
    * \def SIQS_DFLT_LSR_METHOD
    * Default linear system resolution method to use.
    */
#define SIQS_DFLT_LINALG_METHOD    SMART_GAUSS_ELIM
   /**
    * \def SIQS_DFLT_USE_LARGE_PRIMES
    * Use the large prime variation by default.
    */
#define SIQS_DFLT_USE_LARGE_PRIMES true

   /**
    * \struct struct_siqs_params_t siqs.h lib/utils/include/siqs.h
    * \brief  Defines the variable parameters used in the SIQS algorithm.
    *
    * This structure defines the set of the variable parameters used in the
    * SIQS algorithm.
    */
struct struct_siqs_params_t {
       /**
        * Sieve's half-width, i.e. the SIQS will sieve in the interval
        * <tt>[-sieve_half_width, sieve_half_width]</tt>.
        */
    uint32_t sieve_half_width;
       /**
        * Number of prime numbers composing the factor base on which to factor
        * the residues.
        */
    uint32_t nprimes_in_base;
       /**
        * Number of the first primes to use in the trial division
        * of the residues.
        *
        * \warning \c nprimes_tdiv should be greater than or equal
        * to 1.
        */
    uint32_t nprimes_tdiv;
       /**
        * Number of congruence relations to find before attempting the
        * factorization of the large integer.
        */
    uint32_t nrelations;
       /**
        * Linear system resolution method to use.
        */
    linalg_method_t linalg_method;
       /**
        * True if we use the large prime variation.
        * False otherwise.
        */
    bool use_large_primes;
};

   /**
    * \typedef siqs_params_t
    * \brief Equivalent to <tt>struct struct_siqs_params_t</tt>.
    */
typedef struct struct_siqs_params_t siqs_params_t;

   /**
    * \brief Fills a \c siqs_params_t with default values.
    *
    * Fills a \c siqs_params_t with default values.
    *
    * \param[in]  n      The number to factor.
    * \param[out] params A pointer to the \c siqs_params_t structure to fill.
    */
void set_siqs_params_to_default(const mpz_t n, siqs_params_t* const params);

   /**
    * \brief Integer factorization via the self-initializing quadratic sieve
    * (SIQS) algorithm.
    *
    * Attempts to factor the non perfect square integer \c n with the
    * SIQS algorithm, using the set of parameters given by <tt>params</tt>
    * and the factoring mode given by <tt>mode</tt>. Found factors are then
    * stored in \c factors. Additionally, if the factoring mode used is set
    * to FIND_COMPLETE_FACTORIZATION, factors' multiplicities are stored in
    * the array <tt>multis</tt>.
    *
    * \note If the factoring mode used is different from
    * FIND_COMPLETE_FACTORIZATION, \c multis is allowed to be a NULL pointer.
    * Otherwise, using a NULL pointer will lead to a fatal error.
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
    *                     FIND_COMPLETE_FACTORIZATION).
    * \param[in] n The non perfect square integer to factor.
    * \param[in] params Pointer to the values of the parameters used in the
    *                   SIQS algorithm.
    * \param[in] mode The factoring mode to use.
    * \return An exit code.
    */
ecode_t siqs(
    mpz_array_t* const factors,
    uint32_array_t* const multis,
    const mpz_t n,
    const siqs_params_t* const params,
    const factoring_mode_t mode
);

#ifdef __cplusplus
}
#endif

#endif
