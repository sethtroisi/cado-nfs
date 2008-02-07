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
 * \file    cfrac.h
 * \author  Jerome Milan
 * \date    Tue Mar 14 2006
 * \version 1.0
 *
 * \brief The CFRAC factorization algorithm.
 *
 * This is the TIFA library's implementation of the CFRAC factorization
 * algorithm from M. A. Morrison and J. Brillhart, together with the large
 * prime variation.
 *
 * \see "A Method of Factoring and the Factorization of F_7", M. A. Morrison and 
 * J. Brillhart, <i>Mathematics of Computation</i>, vol 29, #129, Jan 1975, 
 * pages 183-205.
 */

#if !defined(_TIFA_CFRAC_H_)
   /**
    * \def _TIFA_CFRAC_H_
    * Standard include guard.
    */
#define _TIFA_CFRAC_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <gmp.h>

#include "first_primes.h"
#include "array.h"
#include "lindep.h"
#include "smooth_filter.h"
#include "factoring_machine.h"
#include "exit_codes.h"

   /**
    * \def CFRAC_DFLT_NPRIMES_IN_BASE
    * Default number of prime numbers composing the factor base on which
    * to factor the residues.
    */
#define CFRAC_DFLT_NPRIMES_IN_BASE  (NFIRST_PRIMES/16)
   /**
    * \def CFRAC_DFLT_NPRIMES_TDIV
    * Default number of the first primes to use in the trial division
    * of the residues.
    */
#define CFRAC_DFLT_NPRIMES_TDIV     (NFIRST_PRIMES/16)
   /**
    * \def CFRAC_DFLT_NRELATIONS
    * Default number of congruence relations to find before attempting the
    * factorization of the large integer.
    */
#define CFRAC_DFLT_NRELATIONS       32
   /**
    * \def CFRAC_DFLT_LSR_METHOD
    * Default linear system resolution method to use.
    */
#define CFRAC_DFLT_LINALG_METHOD    SMART_GAUSS_ELIM
   /**
    * \def CFRAC_DFLT_USE_LARGE_PRIMES
    * Use the large prime variation by default.
    */
#define CFRAC_DFLT_USE_LARGE_PRIMES true

   /**
    * \struct struct_cfrac_params_t cfrac.h lib/algo/include/cfrac.h
    * \brief  Defines the variable parameters used in the CFRAC algorithm.
    *
    * This structure defines the set of the variable parameters used in the
    * CFRAC algorithm.
    */
struct struct_cfrac_params_t {
       /**
        * Number of prime numbers composing the factor base on which to factor
        * the residues.
        */
    uint32_t nprimes_in_base;
       /**
        * Number of the first primes to use in the trial division
        * of the residues known to be smooth.
        *
        * \warning \c nprimes_tdiv should be greater than or equal to 1.
        */
    uint32_t nprimes_tdiv;
       /**
        * Number of linear dependences to find.
        */
    uint32_t nrelations;
       /**
        * Linear system resolution method to use.
        */
    linalg_method_t linalg_method;
       /**
        * True if we use the single large prime variation.
        * False otherwise.
        */
    bool use_large_primes;
       /**
        * Method to use to detect smooth residues and relations.
        */
    smooth_filter_method_t filter_method;
       /**
        * Number of steps in the early abort strategy. If zero, no early
        * abort is performed. Only used is \c linalg_method is set to
        * \c TDIV or \c TDIV_EARLY_ABORT.
        *
        * \note \c nsteps should be less than or equal to \c MAX_NSTEPS,
        * as defined in smooth_filter.h.
        */
    unsigned short int nsteps_early_abort;
};

   /**
    * \typedef cfrac_params_t
    * \brief Equivalent to <tt>struct struct_cfrac_params_t</tt>.
    */
typedef struct struct_cfrac_params_t cfrac_params_t;

   /**
    * \brief Fills a \c cfrac_params_t with "good" default values.
    *
    * Fills a \c cfrac_params_t with "good" default values choosen according
    * to the size of the number n to factor.
    *
    * \warning There is no guarantee that the choosen parameter values will
    * be the best ones for a given number to factor. However, provided that the
    * number to factor is between 40 and 200 bits long, the choosen values
    * should be nearly optimal.
    *
    * \param[in]  n The \c mpz_t integer to factor.
    * \param[out] params A pointer to the \c cfrac_params_t structure to fill.
    */
void set_cfrac_params_to_default(const mpz_t n, cfrac_params_t* const params);

   /**
    * \brief Integer factorization via the continued fraction (CFRAC) algorithm.
    *
    * Attempts to factor the non perfect square integer \c n with the
    * CFRAC algorithm, using the set of parameters given by <tt>params</tt>
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
    * \warning The "no large primes" variant is currently disabled.
    *
    * \param[out] factors Pointer to the found factors of <tt>n</tt>.
    * \param[out] multis Pointer to the multiplicities of the found factors
    *                    (only computed if \c mode is set to
    *                     FIND_COMPLETE_FACTORIZATION).
    * \param[in] n The non perfect square integer to factor.
    * \param[in] params Pointer to the values of the parameters used in the
    *                   CFRAC algorithm.
    * \param[in] mode The factoring mode to use.
    * \return An exit code.
    */
ecode_t cfrac(
    mpz_array_t* const factors,
    uint32_array_t* const multis,
    const mpz_t n,
    const cfrac_params_t* const params,
    const factoring_mode_t mode
);

#ifdef __cplusplus
}
#endif

#endif
