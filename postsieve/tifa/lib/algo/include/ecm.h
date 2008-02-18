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
 * \file    ecm.h
 * \author  Jerome Milan
 * \date    Thu Feb 7 2008
 * \version 1.0
 *
 * \brief The elliptic curve method of integer factorization (ECM).
 *
 * This is the TIFA library's implementation of the 'ECM' factorization
 * algorithm. The second phase of the algorithm follows the standard 
 * continuation and is implemented in a way reminiscent of the description
 * given in the article "Implementing the Elliptic Curve Method of Factoring
 * in Reconfigurable Hardware" by Kris Gaj et al.
 *
 * \see "Implementing the Elliptic Curve Method of Factoring in Reconfigurable
 * Hardware", K. Gaj et al., <i>Cryptographic Hardware and Embedded Systems - 
 * CHES 2006</i>.
 *
 * \warning This is merely a toy-implementation of ECM without any smart
 * optimizations. More work is certainly needed to make it competitive
 * for small numbers. Large numbers are, of course, out of the scope of this
 * library.
 */

#if !defined(_TIFA_ECM_H_)
   /**
    * \def _TIFA_ECM_H_
    * Standard include guard.
    */
#define _TIFA_ECM_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <gmp.h>

#include "array.h"
#include "factoring_machine.h"
#include "exit_codes.h"

   /**
    * \struct struct_ecm_params_t ecm.h lib/algo/include/ecm.h
    * \brief  Defines the variable parameters used in ECM.
    *
    * This structure defines the set of the variable parameters used in the
    * elliptic curve method (ECM).
    */
struct struct_ecm_params_t {
        /**
         * Bound used in the first phase of ECM.
         */
    uint32_t b1;
        /**
         * Bound used in the second phase of ECM. If set to 0, no second phase
         * is performed.
         *
         * \warning Due to a current limitation in the code, it is \c required
         * than \c b2 == 0 (no second phase) or \c b2 > 105. Failure to assess 
         * such a condition will lead to unpredictable behaviour.
         */
    uint32_t b2;
        /**
         * Number of curves to try before giving up the factorization when
         * using the SINGLE_RUN factoring mode.
         */
    uint32_t ncurves;
};
   /**
    * \typedef ecm_params_t
    * \brief Equivalent to <tt>struct struct_ecm_params_t</tt>.
    */
typedef struct struct_ecm_params_t ecm_params_t;

   /**
    * \brief Fills an \c ecm_params_t with "good" default values.
    *
    * Fills an \c ecm_params_t with "good" default values choosen according
    * to the size of the number n to factor.
    *
    * \warning This is, for the time being, a dummy function. Parameters
    * are \e not set to suitable values \e at \e all! Do \e not use it:
    * for the time being, you should choose the parameters by yourself!
    * Shocking!
    *
    * \param[in]  n The \c mpz_t integer to factor.
    * \param[out] params A pointer to the \c ecm_params_t structure to fill.
    */
void set_ecm_params_to_default(const mpz_t n, ecm_params_t* const params);

   /**
    * \brief Integer factorization with the elliptic curve method (ECM).
    *
    * Attempts to factor the non perfect square integer \c n with the
    * ECM, using the set of parameters given by <tt>params</tt> and the
    * factoring mode given by <tt>mode</tt>. Found factors are then
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
    * \param[in] params Pointer to the values of the parameters used in the ECM.
    * \param[in] mode The factoring mode to use.
    * \return An exit code.
    */
ecode_t ecm(
    mpz_array_t* const factors,
    uint32_array_t* const multis,
    const mpz_t n,
    const ecm_params_t* const params,
    const factoring_mode_t mode
);

#ifdef __cplusplus
}
#endif

#endif
