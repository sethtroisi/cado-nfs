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
 * \file    fermat.h
 * \author  Jerome Milan
 * \date    Tue Aug 28 2007
 * \version 1.0
 *
 * \brief McKee's variant of the Fermat factorization algorithm.
 *
 * This is the TIFA library's implementation of James McKee's proposed speedup
 * of the Fermat factorization algorithm (SQUFOF), based on the description
 * given by McKee in its paper "Speeding Fermat's Factoring Method".
 *
 * \note This implementation can only factor numbers whose size is less than
 * twice the size of an \c unsigned \c long \c int.
 *
 * \see "Speeding Fermat's Factoring Method", James McKee.
 * Mathematics of Computation, Volume 68, Number 228, pages 1729-1737.
 *
 */

#if !defined(_TIFA_FERMAT_H_)
   /**
    * \def _TIFA_FERMAT_H_
    * Standard include guard.
    */
#define _TIFA_FERMAT_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <gmp.h>

#include "array.h"
#include "factoring_machine.h"
#include "exit_codes.h"

   /**
    * \struct struct_fermat_params_t fermat.h lib/utils/include/fermat.h
    * \brief  Defines the variable parameters used in Fermat's algorithm
    *         (dummy structure).
    *
    * This structure is intended to define the set of the variable parameters 
    * used in Fermat's algorithm.
    *
    * \warning For the time being, this is a completely unused dummy structure
    *          which is kept only as a placeholder should the need for user
    *          parameters arise in future code revisions.
    */
struct struct_fermat_params_t {
       /**
        * Unused dummy variable.
        */
    unsigned int _dummy_variable_;
};

   /**
    * \typedef fermat_params_t
    * \brief Equivalent to <tt>struct struct_fermat_params_t</tt>.
    */
typedef struct struct_fermat_params_t fermat_params_t;

   /**
    * \brief Fills a \c fermat_params_t with default values (dummy function).
    *
    * This function is intended to fill a \c fermat_params_t with default
    * values.
    *
    * \warning For the time being, this is a dummy function which does 
    *          absolutely nothing at all, but is kept only as a placeholder 
    *          should the need for user parameters arise in future code 
    *          revisions.
    *
    * \param params A pointer to the \c fermat_params_t structure to fill.
    */
void set_fermat_params_to_default(fermat_params_t* const params);

   /**
    * \brief Integer factorization via McKee's speedup of Fermat's factorization
    * algorithm.
    *
    * Attempts to factor the non perfect square integer \c n using James McKee's
    * proposed enhancement of Fermat's algorithm, using the factoring mode given 
    * by <tt>mode</tt>. Found factors are then stored in \c factors. 
    * Additionally, if the factoring mode used is set to 
    * FIND_COMPLETE_FACTORIZATION, factors' multiplicities are stored in the 
    * array <tt>multis</tt>.
    *
    * \warning This implementation can only factor numbers whose sizes in bits
    * are strictly less than twice the size of an <tt>unsigned long int</tt>.
    * This choice was made to maximize the use of single precision operations. 
    * Such a limitation should not be much of a problem since Fermat's algorithm
    * is mostly used to factor very small integers (up to, say, 20 decimal 
    * digits).
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
    * \param[in] params Fermat's algorithm parameters (currently unused).
    * \param[in] mode The factoring mode to use.
    * \return An exit code.
    */
ecode_t fermat(
    mpz_array_t* const factors,
    uint32_array_t* const multis,
    const mpz_t n,
    const fermat_params_t* const params,
    const factoring_mode_t mode
);

#ifdef __cplusplus
}
#endif

#endif
