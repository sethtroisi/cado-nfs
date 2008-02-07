//
// Copyright (C) 2006, 2007, 2008 INRIA (French National Institute for Research
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
 * \file    factoring_program.h
 * \author  Jerome Milan
 * \date    Wed Mar 12 2007
 * \version 1.0
 *
 * \brief The logic common to all TIFA's factorization executable programs.
 */

#if !defined(_TIFA_FACTORING_PROGRAM_H_)
   /**
    * \def _TIFA_FACTORING_PROGRAM_H_
    * Standard include guard.
    */
#define _TIFA_FACTORING_PROGRAM_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <gmp.h>

#include "array.h"
#include "exit_codes.h"
#include "factoring_machine.h"

   /**
    * \struct factoring_program_struct factoring_program.h
    *         tools/include/factoring_program.h
    *
    * \brief  Defines a structure to represent the logic behind all
    * factorization programs.
    *
    * This structure defines a set of data and functions to represent the logic
    * behind all factorization programs. The idea is to be able to write down
    * the actual factoring process once and for all so that actual factorization
    * program can use such a structure by merely "filling-in" the gaps.
    *
    */
struct factoring_program_struct {
       /**
        * Argument count as given by the "main" function.
        */
    int argc;
       /**
        * Argument values as given by the "main" function.
        */
    char** argv;
       /**
        * Integer to factor.
        */
    mpz_t n;
       /**
        * Factoring mode to use.
        */
    factoring_mode_t mode;
       /**
        * Verbosity option. Should be either 1 (be verbose) or 0 (stay laconic).
        */
    int verbose;
       /**
        * Timing option. Should be either 1 (proceed with timings) or 0 (do not
        * perform timings).
        */
    int timing;
       /**
        * Number of primes used in the trial division of the number to factor.
        */
    uint32_t nprimes_tdiv;
       /**
        * The maximum number of factors to collect (excluding the factors
        * found in the trial division stage).
        */
    uint32_t nfactors;
       /**
        * Name of the factoring algorithm to use (preferably a short acronym).
        */
    char* algo_name;
       /**
        * A pointer to the parameters of the factoring algorithm to use.
        */
    void* params;
       /**
        * A pointer to a function printing the usage of the factoring
        * program.
        *
        * \param A pointer to the \c factoring_program_t used by the actual
        * factoring program.
        */
    void    (*print_usage_func)     (struct factoring_program_struct* const);
       /**
        * A pointer to a function printing the values of the parameters used
        * by the actual factoring program.
        *
        * \param A pointer to the \c factoring_program_t used by the actual
        * factoring program.
        */
    void    (*print_params_func)    (struct factoring_program_struct* const);
       /**
        * A pointer to a function reading arguments on the command line and
        * setting the parameters of the actual factoring program.
        *
        * \param A pointer to the \c factoring_program_t used by the actual
        * factoring program.
        */
    void    (*process_args_func)    (struct factoring_program_struct* const);
       /**
        * A pointer to a function implementing the factorization algorithm
        * to use.
        *
        * \param factors An \c mpz_array_t to hold the found factors.
        * \param factors A \c uint32_array_t to hold the found multiplicities.
        * \param n       The integer to factor.
        * \param params  A pointer to the parameters to be used in the
        *                factorization algorithm.
        * \param mode The factorization mode to use.
        */
    ecode_t (*factoring_algo_func)  (mpz_array_t* const factors,
                                     uint32_array_t* const multis,
                                     const mpz_t n,
                                     const void* const params,
                                     factoring_mode_t mode);
       /**
        * A pointer to a function setting the algorithm's parameters to some
        * default values.
        *
        * \param A pointer to the \c factoring_program_t used by the actual
        * factoring program.
        */
    void (*set_params_to_default_func) (struct factoring_program_struct* const);

};
   /**
    * \typedef factoring_program_t
    * \brief Equivalent to <tt>struct factoring_program_struct</tt>.
    */
typedef struct factoring_program_struct factoring_program_t;

   /**
    * \brief Run a factoring program.
    *
    * Run an actual factoring program from the command line.
    *
    * \param program The \c factoring_program_t to run.
    */
ecode_t run_program(factoring_program_t* const program);

#ifdef __cplusplus
}
#endif

#endif

