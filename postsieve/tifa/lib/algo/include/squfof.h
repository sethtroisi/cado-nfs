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
 * \file    squfof.h
 * \author  Jerome Milan
 * \date    Wed Jan 24 2006
 * \version 1.0
 *
 * \brief The SQUFOF factorization algorithm.
 *
 * This is the TIFA library's implementation of Shanks' square form
 * factorization algorithm (SQUFOF), based on the description given by Jason
 * Gower and Samuel Wagstaff in their paper "Square Form Factorization" to be
 * published in Mathematics of Computation.
 *
 * \note This implementation can only factor numbers whose size is less than
 * twice the size of an \c unsigned \c long \c int.
 *
 * \see "Square Form Factorization", Jason E. Gower & Samuel S. Wagstaff Jr.
 * <i>Mathematics of Computation</i>, S 0025-5718(07)02010-8, Article 
 * electronically published on May 14, 2007.
 *
 * \see "Square Form Factorization", Jason E. Gower, PhD thesis, Purdue
 * University, December 2004.
 *
 * \see For a description of the "large step algorithm" used to quickly
 * jump over several forms, see: "On the Parallel Generation of the Residues for 
 * the Continued Fraction Factoring Algorithm", Hugh C. Williams, Marvin C. 
 * Wunderlich, <i>Mathematics of Computation</i>, Volume 48, Number 177,
 * January 1987, pages 405-423
 */

 /*
  *  Copyright (C) 2006, 2007  INRIA
  *  Licence: GNU Lesser General Public License (LGPL)
  */

#if !defined(_TIFA_SQUFOF_H_)
   /**
    * \def _TIFA_SQUFOF_H_
    * Standard include guard.
    */
#define _TIFA_SQUFOF_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <gmp.h>

#include "array.h"
#include "factoring_machine.h"
#include "exit_codes.h"

   /**
    * \struct struct_squfof_params_t squfof.h lib/utils/include/squfof.h
    * \brief  Defines the variable parameters used in the SQUFOF algorithm
    *         (dummy structure).
    *
    * This structure is intended to define the set of the variable parameters 
    * used in the SQUFOF algorithm.
    *
    * \warning For the time being, this is a completely unused dummy structure
    *          which is kept only as a placeholder should the need for user
    *          parameters arise in future code revisions.
    */
struct struct_squfof_params_t {
       /**
        * Unused dummy variable.
        */
    unsigned int _dummy_variable_;
};

   /**
    * \typedef squfof_params_t
    * \brief Equivalent to <tt>struct struct_squfof_params_t</tt>.
    */
typedef struct struct_squfof_params_t squfof_params_t;


   /**
    * \brief Fills a \c squfof_params_t with default values (dummy function).
    *
    * This function is intended to fill a \c squfof_params_t with default
    * values.
    *
    * \warning For the time being, this is a dummy function which does 
    *          absolutely nothing at all, but is kept only as a placeholder 
    *          should the need for user parameters arise in future code 
    *          revisions.
    *
    * \param params A pointer to the \c squfof_params_t structure to fill.
    */
void set_squfof_params_to_default(squfof_params_t* const params);

   /**
    * \brief Integer factorization via the square form factorization
    * (SQUFOF) algorithm.
    *
    * Attempts to factor the non perfect square integer \c n with the
    * SQUFOF algorithm, using the factoring mode given by <tt>mode</tt>.
    * Found factors are then stored in \c factors. Additionally, if the
    * factoring mode used is set to FIND_COMPLETE_FACTORIZATION,
    * factors' multiplicities are stored in the array <tt>multis</tt>.
    *
    * \warning This implementation can only factor numbers whose sizes in bits
    * are strictly less than twice the size of an <tt>unsigned long int</tt>
    * (the exact limit depending on the number to factor and the multiplier
    * used). This choice was made because most of the computations are
    * then performed using only single precision operations. Such a limitation
    * should not be much of a problem since SQUFOF is mostly used to factor
    * very small integers (up to, say, 20 decimal digits).
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
    * \param[in] params SQUFOF's parameters (currently unused).
    * \param[in] mode The factoring mode to use.
    * \return An exit code.
    */
ecode_t squfof(
    mpz_array_t* const factors,
    uint32_array_t* const multis,
    const mpz_t n,
    const squfof_params_t* const params,
    const factoring_mode_t mode
);

#ifdef __cplusplus
}
#endif

#endif
