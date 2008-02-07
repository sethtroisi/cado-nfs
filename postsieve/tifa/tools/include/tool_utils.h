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
 * \file    tool_utils.h
 * \author  Jerome Milan
 * \date    Wed Mar 22 2006
 * \version 1.0
 *
 * \brief Miscellaneous helpful functions.
 *
 * Miscellaneous functions used by the "tool programs".
 */

#if !defined(_TIFA_TOOL_UTILS_H_)
   /**
    * \def _TIFA_TOOL_UTILS_H_
    * Standard include guard.
    */
#define _TIFA_TOOL_UTILS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <inttypes.h>
#include <stdbool.h>

#include <gmp.h>

#include "array.h"

   /**
    * \brief Does a string \c str_n read as a number?
    *
    * Returns \c true if the string \c str_n represents a number in the
    * decimal base, \c false otherwise.
    *
    * \param[in] str_n The string to check.
    * \param[in] length The maximum length of the string to check.
    */
bool is_a_number(const char* const str_n, uint32_t length);

   /**
    * \brief <tt>NULL</tt> terminates a string.
    *
    * <tt>NULL</tt> terminates the string \c str_n at the first occurence
    * of a newline encountered. If no newline is found \c str_n is left
    * unchanged.
    *
    * \note This function is actually quite different from the Perl builtin
    * \c chomp function. It's name should probably be changed to avoid
    * possible confusion.
    *
    * \param[in] str_n The string to \c NULL terminate.
    * \param[in] length The maximum length of the string to check.
    */
void chomp(char* const str_n, uint32_t length);

#ifdef __cplusplus
}
#endif

#endif
