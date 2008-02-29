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
 * \file    factoring_machine.h
 * \author  Jerome Milan
 * \date    Fri Jan 11 2008
 * \version 1.1
 *
 * \brief Abstraction of an integer factorization algorithm.
 *
 * Implements a machine-like abstraction of an integer factorization algorithm.
 */

 /*
  *  History:
  *    1.1: Fri Jan 11 2008 by JM:
  *         - Removed machine_state_enum and machine_state_t.
  *    1.0: March 2007 by JM:
  *         - Initial version.
  */

#if !defined(_TIFA_FACTORING_MACHINE_H_)
   /**
    * \def _TIFA_FACTORING_MACHINE_H_
    * Standard include guard.
    */
#define _TIFA_FACTORING_MACHINE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <gmp.h>

#include "array.h"
#include "exit_codes.h"

   /**
    * \enum factoring_mode_enum
    *
    * An enumeration of the factoring mode available to the implemented
    * factorization algorithm.
    */
enum factoring_mode_enum {
       /**
        * Perform only a single run of the factorization algorithm.
        */
    SINGLE_RUN,
       /**
        * Run the factorization algorithm until either some factors are found
        * or the abort limit (defined on a per-algorithm basis) is reached.
        */
    FIND_SOME_FACTORS,
       /**
        * Run the factorization algorithm until either some coprime
        * factors are found or the abort limit (defined on a per-algorithm
        * basis) is reached.
        */
    FIND_SOME_COPRIME_FACTORS,
       /**
        * Run the factorization algorithm until either some prime factors are
        * found or the abort limit (defined on a per-algorithm basis) is
        * reached.
        *
        * \note This is probably not useful at all. Why would we discard found
        * factors even if they are not prime? This should better be left to
        * the client application.
        */
    FIND_SOME_PRIME_FACTORS,
       /**
        * Run the factorization algorithm until either the complete
        * factorization (as a product of prime numbers) is found or the abort
        * limit (defined on a per-algorithm basis) is reached.
        */
    FIND_COMPLETE_FACTORIZATION
};
   /**
    * \typedef factoring_mode_t
    * \brief Equivalent to <tt>struct factoring_mode_enum</tt>.
    */
typedef enum factoring_mode_enum factoring_mode_t;

   /**
    * Global constant array mapping factoring modes to their respective best
    * outcome.
    */
static const int mode_to_outcome[5] = {
    SOME_FACTORS_FOUND,
    SOME_FACTORS_FOUND,
    SOME_COPRIME_FACTORS_FOUND,
    SOME_PRIME_FACTORS_FOUND,
    COMPLETE_FACTORIZATION_FOUND
};

   /**
    * \struct factoring_machine_struct factoring_machine.h
    *         tools/include/factoring_machine.h
    *
    * \brief  Defines a structure to represent the logic behind all
    * factorization algorithms.
    *
    * This structure defines a set of data and functions to represent the logic
    * behind all factorization algorithms. The idea is to be able to write down
    * the factoring process' boilerplate once and for all so that actual
    * factorization algorithm can use such a structure by merely "filling-in"
    * the gaps.
    *
    */
struct factoring_machine_struct {
       /**
        * The integer to factor.
        */
    mpz_t n;
       /**
        * The factoring mode to use.
        */
    factoring_mode_t mode;
       /**
        * The context of the factorization algorithm. This is a pointer to
        * an algorithm implementation dependant structure holding all variables,
        * data, and functions needed by the implementation.
        */
    void* context;
       /**
        * The parameters used by the factorization algorithm. This is a pointer
        * to an algorithm implementation dependant structure holding the
        * algorithm parameters needed by the implementation.
        */
    void* params;
       /**
        * A pointer to a function initializing the algorithm context.
        *
        * \param (unnamed) A pointer to the \c factoring_machine_t used by the
        *        actual algorithm implementation.
        */
    ecode_t (*init_context_func)   (struct factoring_machine_struct* const);
       /**
        * A pointer to a function performing the actual factorization stage
        * of the factorization algorithm (by opposition to the initialization
        * stage for example).
        *
        * \param (unnamed) A pointer to the \c factoring_machine_t used by the
        *         actual algorithm implementation.
        */
    ecode_t (*perform_algo_func)   (struct factoring_machine_struct* const);
       /**
        * A pointer to a function updating the context of the factorization
        * algorithm. This function is responsible of the definition of the
        * factorization strategy should something bad happens (e.g. what to
        * do if no factors are found after the first run?).
        *
        * \param (unnamed) A pointer to the \c factoring_machine_t used by the
        *        actual algorithm implementation.
        */
    ecode_t (*update_context_func) (struct factoring_machine_struct* const);
       /**
        * A pointer to a function clearing all memory space used by the
        * context.
        *
        * \param (unnamed) A pointer to the \c factoring_machine_t used by the
        * actual algorithm implementation.
        */
    ecode_t (*clear_context_func)  (struct factoring_machine_struct* const);
       /**
        * A pointer to the function to use to perform recursive factorization
        * (i.e. factorization of the non-prime factors found).
        *
        * \param (unnamed) A pointer to an \c mpz_array_t to hold the found
        *        factors.
        * \param (unnamed) A pointer to an \c uint32_array_t to hold the
        *        multiplicities.
        * \param (unnamed) The non-prime factor to factorize.
        * \param (unnamed) The \c factoring_mode_t to use.
        */
    ecode_t (*recurse_func) (mpz_array_t* const,
                             uint32_array_t* const,
                             const mpz_t,
                             factoring_mode_t);
       /**
        * A pointer to a \c mpz_array_t to hold the found factors.
        */
    mpz_array_t* factors;
       /**
        * A pointer to a \c uint32_array_t to hold the multiplicities of the
        * found factors.
        */
    uint32_array_t* multis;
       /**
        * true if the algorithm succeeds.
        * false otherwise.
        *
        * \note The notion of success is given by the \c factoring_mode_t mode
        * used and not on whether or not some factors are found.
        */
    bool success;
};
   /**
    * \typedef factoring_machine_t
    * \brief Equivalent to <tt>struct factoring_machine_struct</tt>.
    */
typedef struct factoring_machine_struct factoring_machine_t;

   /**
    * \brief Attempt to factor an integer.
    *
    * Attempt to factor an integer with all parameters given by \c machine.
    *
    * \note This function is meant to be a starting point for implementations
    * of factorization algorithms and is obviously not intended to be directly
    * used as a factoring function all by itself.
    *
    * \param machine A pointer to the \c factoring_machine_t to use.
    */
inline ecode_t run_machine(factoring_machine_t* machine);

#ifdef __cplusplus
}
#endif

#endif

