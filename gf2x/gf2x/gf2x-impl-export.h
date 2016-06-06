/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009, 2010, 2013, 2015
   Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann

   This program is free software; you can redistribute it and/or modify it
   under the terms of either:
    - If the archive contains a file named toom-gpl.c (not a trivial
    placeholder), the GNU General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your
    option) any later version.
    - If the archive contains a file named toom-gpl.c which is a trivial
    placeholder, the GNU Lesser General Public License as published by
    the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.
   
   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the license text for more details.
   
   You should have received a copy of the GNU General Public License as
   well as the GNU Lesser General Public License along with this program;
   see the files COPYING and COPYING.LIB.  If not, write to the Free
   Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
   02110-1301, USA.
*/

#ifndef GF2X_IMPL_EXPORT_H_
#define GF2X_IMPL_EXPORT_H_

#include "gf2x/gf2x-config-export.h"

/* This file contains helper macros which are used in inline functions in
 * gf2x, and therefore potentially exposed to the users which include the
 * gf2x-small.h files, for getting the inlines.
 *
 * We also activate the relevant includes.
 */

#ifndef	GF2X_MAYBE_UNUSED
#if defined(__GNUC__)
#define GF2X_MAYBE_UNUSED __attribute__ ((unused))
#else
#define GF2X_MAYBE_UNUSED
#endif
#endif

#include <stdint.h>

#ifdef  GF2X_HAVE_SSE2_SUPPORT
#include <emmintrin.h>
#if defined(__GNUC__) && __GNUC__ == 4 &&__GNUC_MINOR__ == 1
#define _gf2x_mm_cvtsi64_m64(u) _mm_cvtsi64x_m64((u))
#else
#define _gf2x_mm_cvtsi64_m64(u) _mm_cvtsi64_m64((u))
#endif
/* _m128i from 2 int64_t's */
#define _gf2x_mm_setr_epi64(lo, hi)                      		\
    _mm_setr_epi64(                                      		\
            _gf2x_mm_cvtsi64_m64((int64_t) (lo)),       		\
            _gf2x_mm_cvtsi64_m64((int64_t) (hi))        		\
        )
/* _m128i from 1 int64_t's */
#define _gf2x_mm_set1_epi64(u) _mm_set1_epi64( _gf2x_mm_cvtsi64_m64((int64_t) (u)))
/* _m128i from 2 int64_t CONSTANTS (and try to get suffix right) */
#define _gf2x_mm_setr_epi64_c(lo, hi)                    		\
    _mm_setr_epi64(                                      		\
            _gf2x_mm_cvtsi64_m64(INT64_C(lo)),          		\
            _gf2x_mm_cvtsi64_m64(INT64_C(hi))           		\
        )
/* _m128i from 1 int64_t CONSTANT (and try to get suffix right) */
#define _gf2x_mm_set1_epi64_c(u) _mm_set1_epi64( _gf2x_mm_cvtsi64_m64(INT64_C(u)))
/* and same for 32-bits (which, for some, have SSE-2) */
#define _gf2x_mm_setr_epi32(a0, a1, a2, a3)				\
    _mm_setr_epi32(                                      		\
            (int32_t) (a0),						\
            (int32_t) (a1),						\
            (int32_t) (a2),						\
            (int32_t) (a3)						\
            )
#define _gf2x_mm_set1_epi32(u) _mm_set1_epi32( (int32_t) (u))
#define _gf2x_mm_setr_epi32_c(a0, a1, a2, a3)				\
    _mm_setr_epi32(                                      		\
            (INT32_C(a0)),          					\
            (INT32_C(a1)),           					\
            (INT32_C(a2)),          					\
            (INT32_C(a3))           					\
        )
#define _gf2x_mm_set1_epi32_c(u) _mm_set1_epi32(INT32_C(u))
#endif

#ifdef  GF2X_HAVE_PCLMUL_SUPPORT
#include <wmmintrin.h>
#endif

#endif	/* GF2X_IMPL_EXPORT_H_ */
