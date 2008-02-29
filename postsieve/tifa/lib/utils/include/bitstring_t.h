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
 * \file    bitstring_t.h
 * \author  Jerome Milan
 * \date    Fri Jan 26 2006
 * \version 1.0
 *
 * \brief Preprocessor defines for 'string of bit' type.
 *
 * Defines several preprocessor \c define symbols related to the type used
 * to represent strings of bit. These symbolss are kept separately to avoid
 * too much namespace pollution.
 */

 /*
  *  Copyright (C) 2006, 2007 INRIA
  *  License: GNU Lesser General Public License (LGPL)
  *  History:
  */

#if !defined(_TIFA_BITSTRING_T_H_)
   /**
    * \def _TIFA_BITSTRING_T_H_
    * Standard include guard.
    */
#define _TIFA_BITSTRING_T_H_

#ifdef __cplusplus
extern "C" {
#endif

//
// _NOTE_: This header file defines symbol that should be kept internal to the
//         TIFA code to avoid namespace pollution (or too many TIFA_* symbol
//         that don't need to be external). Do not include this header
//         in another TIFA header!
//

//
// tifa_config.h is included to get the type TIFA_BITSTRING_T. This type
// could alternatively be defined here but it may be of interest to users
// of precompiled packages, so it is kept exposed in the public
// tifa_config.h headers.
//
#include "tifa_config.h"

//
// sizeof and number of bits of the native C type used to represent sets
// of bits. These are merely aliases to TIFA_SIZEOF_BITSTRING_T and
// TIFA_BITSTRING_T_BITSIZE to try to keep macro names not too long...
//
#define SIZEOF_BITSTRING_T  TIFA_SIZEOF_BITSTRING_T
#define BITSTRING_T_BITSIZE TIFA_BITSTRING_T_BITSIZE

//
// If BITSTRING_T_BITSIZE is a power of two (and that will almost always
// be the case unless one chose a weird underlying type to represent sets of
// bits), define BITSTRING_T_SIZE_IS_POW_OF_TWO to 1 and assign that power
// to POW_TWO_BITSTRING_T_SIZE
//
#if (    (BITSTRING_T_BITSIZE == 8)   || (BITSTRING_T_BITSIZE == 16)    \
      || (BITSTRING_T_BITSIZE == 32)  || (BITSTRING_T_BITSIZE == 64)    \
      || (BITSTRING_T_BITSIZE == 128) || (BITSTRING_T_BITSIZE == 256)   \
      || (BITSTRING_T_BITSIZE == 512) || (BITSTRING_T_BITSIZE == 1024)  \
    )
    #define BITSTRING_T_SIZE_IS_POW_OF_TWO 1
#else
    #define BITSTRING_T_SIZE_IS_POW_OF_TWO 0
#endif

#if BITSTRING_T_SIZE_IS_POW_OF_TWO
    //
    // BITSTRING_T_BITSIZE = 2 ^ POW_TWO_BITSTRING_T_SIZE
    //
    #if   (BITSTRING_T_BITSIZE == 8)
        #define POW_TWO_BITSTRING_T_SIZE  3

    #elif (BITSTRING_T_BITSIZE == 16)
        #define POW_TWO_BITSTRING_T_SIZE  4

    #elif (BITSTRING_T_BITSIZE == 32)
        #define POW_TWO_BITSTRING_T_SIZE  5

    #elif (BITSTRING_T_BITSIZE == 64)
        #define POW_TWO_BITSTRING_T_SIZE  6

    #elif (BITSTRING_T_BITSIZE == 128)
        #define POW_TWO_BITSTRING_T_SIZE  7

    #elif (BITSTRING_T_BITSIZE == 256)
        #define POW_TWO_BITSTRING_T_SIZE  8

    #elif (BITSTRING_T_BITSIZE == 512)
        #define POW_TWO_BITSTRING_T_SIZE  9

    #elif (BITSTRING_T_BITSIZE == 1024)
        #define POW_TWO_BITSTRING_T_SIZE  10

    #endif
#endif

#ifdef __cplusplus
}
#endif

#endif

