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
 * \file    exit_codes.h
 * \author  Jerome Milan
 * \date    Circa February (March?) 2007
 * \version 1.0
 *
 * \brief Exit codes used by/in some of the TIFA functions
 *
 * Defines several exit codes used by/in some of the TIFA functions together
 * with their string representations.
 */

#if !defined(_TIFA_EXIT_CODES_H_)
#define _TIFA_EXIT_CODES_H_

#ifdef __cplusplus
extern "C" {
#endif

   /**
    * \enum ecode_enum
    *
    * An enumeration of the possible exit codes used by/in some TIFA functions.
    */
enum ecode_enum {
        /**
         * Used by a \c factoring_machine_t to indicate that the
         * \c factoring_mode_t passed as parameter is not valid.
         */
    UNKNOWN_FACTORING_MODE,
        /**
         * Used by the factorization algorithm to indicate that \e some factors
         * were found. In that case, the factors' multiplicities are
         * not computed.
         */
    SOME_FACTORS_FOUND,
        /**
         * Used by the factorization algorithm to indicate that \e some prime
         * factors were found. In that case, the factors' multiplicities are
         * not computed.
         */
    SOME_PRIME_FACTORS_FOUND,
        /**
         * Used by the factorization algorithm to indicate that \e some coprime
         * factors were found. In that case, the factors' multiplicities are
         * not computed.
         */
    SOME_COPRIME_FACTORS_FOUND,
        /**
         * Used by the factorization algorithm to indicate that a partial
         * factorization (in terms of a set of coprimes and multiplicities)
         * was found. The term "partial" refers to the fact that some found
         * factors may not be prime. However the product of the found factors
         * (taking into account their associated multiplicities) does yield
         * the original number to factor.
         */
    PARTIAL_FACTORIZATION_FOUND,
        /**
         * Used by the factorization algorithm to indicate that a complete
         * factorization (in terms of primes and multiplicities) was found.
         */
    COMPLETE_FACTORIZATION_FOUND,
        /**
         * Used by the factorization algorithm to indicate that no factor was
         * found.
         */
    NO_FACTOR_FOUND,
        /**
         * Generic exit code used to indicate a serious internal error, possibly
         * leading to an unpredictable behavior.
         */
    FATAL_INTERNAL_ERROR,
        /**
         * Used by the SQUFOF algorithm to indicate that the queue overflowed,
         * thus leading to give up the factorization process.
         */
    QUEUE_OVERFLOW,
        /**
         * Used by the SQUFOF algorithm to indicate that no proper form was
         * found, thus leading to give up the factorization process.
         */
    NO_PROPER_FORM_FOUND,
        /**
         * Used to indicate that an abort limit has been reached leading to
         * give up the current operation.
         */
    GIVING_UP,
        /**
         * Used by the SQUFOF and Fermat/McKee implementations to indicate that 
         * the integer to factor is too large and cannot be processed.
         */
    INTEGER_TOO_LARGE,
        /**
         * Generic exit code used to indicate that an operation succeeded.
         */
    SUCCESS,
        /**
         * Generic exit code used to indicate that an operation failed.
         */
    FAILURE
};

   /**
    * \typedef ecode_t
    * \brief Equivalent to <tt>struct ecode_enum</tt>.
    */
typedef enum ecode_enum ecode_t;

   /**
    * Global constant array mapping exit codes to their string representations.
    */
static const char* const ecode_to_str[14] = {
    "unknown factoring mode",
    "some factors found",
    "some prime factors found",
    "some coprime factors found",
    "partial factorization found",
    "complete factorization found",
    "no factor found",
    "fatal internal error",
    "queue overflow",
    "no proper form found",
    "giving up",
    "number to factor is too large",
    "success",
    "failure"
};

   /**
    * \def PRINT_ECODE(ECODE)
    * Macro printing the string representation of the exit code \c ECODE
    * on the standard output, followed by a newline.
    */
#define PRINT_ECODE(ECODE) printf("%s\n", ecode_to_str[ECODE]);

#ifdef __cplusplus
}
#endif

#endif
