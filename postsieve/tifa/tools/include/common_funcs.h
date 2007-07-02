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
 * \file    common_funcs.h
 * \author  Jerome Milan
 * \date    Wed Mar 12 2007
 * \version 1.0
 *
 * \brief Miscellaneous functions and macros used by the "tool" programs.
 */

#if !defined(_TIFA_COMMON_FUNCS_H_)
   /**
    * \def _TIFA_COMMON_FUNCS_H_
    * Standard include guard.
    */
#define _TIFA_COMMON_FUNCS_H_

#ifdef __cplusplus
extern "C" {
#endif

    //
    // Various macros used to keep output messages consistent throughout the
    // different factoring programs.
    //

   /**
    * \def PRINT_ABORT_MSG()
    * Macro printing an "program has aborted" message on the standard error.
    */
#define PRINT_ABORT_MSG() fprintf(stderr, "Program aborted\n");

   /**
    * \def PRINT_NAN_ERROR(X)
    * Macro printing a "X is not an integer" error message on the standard
    * error.
    */
#define PRINT_NAN_ERROR(X) do {                             \
    fprintf(stderr, "\nERROR: %s is not an integer!\n", X); \
    PRINT_ABORT_MSG();                                      \
} while (0)

   /**
    * \def PRINT_BAD_ARGC_ERROR()
    * Macro printing a "bad number of argument" error message on the standard
    * error.
    */
#define PRINT_BAD_ARGC_ERROR() do {                           \
    fprintf(stderr, "\nERROR: Bad number of arguments!\n\n"); \
} while (0)

   /**
    * \def PRINT_ENTER_NUMBER_MSG()
    * Macro displaying a prompt asking the user to enter the integer to factor.
    */
#define PRINT_ENTER_NUMBER_MSG() printf("> Enter the integer to factor: ")

   /**
    * \def PRINT_USAGE_WARNING_MSG()
    * Macro printing a boilerplate usage warning on the standard error.
    */
#define PRINT_USAGE_WARNING_MSG() do {                                       \
    fprintf(stderr, "Please, use the perl wrapper factorize.pl instead of"); \
    fprintf(stderr, " a direct invocation of this program.\n");              \
} while (0)

   /**
    * \def MAX_NDIGITS
    * Maximal number of decimal digits of the number to factor.
    */
#define MAX_NDIGITS 256
   /**
    * \def NTRIES_MILLER_RABIN
    * Number of iterations to use in the Miller-Rabin composition tests.
    */
#define NTRIES_MILLER_RABIN 32
   /**
    * \def NPRIMES_TRIAL_DIV
    * Default number of prime numbers used in trial division of number to
    * factor.
    */
#define NPRIMES_TRIAL_DIV 1024

   /**
    * \brief Function used by the "tool" programs to print a greeting message.
    *
    * \param[in] name Name of the factoring program.
    */
void print_hello_msg(char* name);

   /**
    * \brief Function used by the "tool" programs to print a bye-bye message.
    *
    * This function could be used by the "tool" programs to print a bye-bye
    * message (it is not used right now).
    */
void print_bye_msg();

#ifdef __cplusplus
}
#endif

#endif

