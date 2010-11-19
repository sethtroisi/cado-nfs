/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009, 2010
   Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 2.1 of the License, or (at
   your option) any later version.
   
   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
   License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with CADO-NFS; see the file COPYING.  If not, write to
   the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
   Boston, MA 02110-1301, USA.
*/
/* Implements 128x128 -> 256 bit product using pclmulqdq instruction. */

#ifndef GF2X_MUL2_H_
#define GF2X_MUL2_H_

#include "gf2x.h"
#include "gf2x/gf2x-impl.h"
/* All gf2x source files for lowlevel functions must include gf2x-small.h
 * This is mandatory for the tuning mechanism. */
#include "gf2x/gf2x-small.h"

#include <stdint.h>
#include <wmmintrin.h>

#if GF2X_WORDSIZE != 64
#error "This code is for 64-bit only"
#endif

#include "gf2x/gf2x-config.h"

#ifndef HAVE_PCLMUL_SUPPORT
#error "This code needs pclmul support"
#endif

GF2X_STORAGE_CLASS_mul2
#ifdef CARRY
/* {t, 4} <- {s1, 2} * {s2, 2}, and {c, 2} <- {s1+1, 1} * {s2+1, 1} */
void
gf2x_mul2c (unsigned long *t, unsigned long const *s1, unsigned long const *s2,
            unsigned long *c)
#else
#ifdef BORROW
void
/* {t, 4} <- {s1, 2} * {s2, 2}, knowing {c, 2} = {s1+1, 1} * {s2+1, 1} */
gf2x_mul2b (unsigned long *t, unsigned long const *s1, unsigned long const *s2,
            unsigned long const *c)
#else
void gf2x_mul2(unsigned long * t, unsigned long const * s1,
        unsigned long const * s2)
#endif
#endif
{
    typedef union {
        __v2di s;
        unsigned long x[2];
    } __v2di_proxy;

    __v2di ss1, ss2, s1s, s2s;
    __v2di_proxy t00, tk;
#ifndef BORROW
    __v2di_proxy t11;
#endif
    ss1 = _mm_loadu_si128((__v2di *)s1);
    ss2 = _mm_loadu_si128((__v2di *)s2);


    t00.s = _mm_clmulepi64_si128(ss1, ss2, 0);
#ifndef BORROW
    t11.s = _mm_clmulepi64_si128(ss1, ss2, 17);
#endif
    
    s1s = _mm_shuffle_epi32(ss1, 78);
    ss1 ^= s1s;
    s2s = _mm_shuffle_epi32(ss2, 78);
    ss2 ^= s2s;
    
    tk.s = _mm_clmulepi64_si128(ss1, ss2, 0);

#ifndef BORROW
    tk.s ^= t00.s ^ t11.s;
#endif

    /* store result */
    t[0] = t00.x[0];
#ifdef BORROW
    t[1] = t00.x[1] ^ tk.x[0] ^ t00.x[0] ^ c[0];
    t[2] = c[0] ^ tk.x[1] ^ t00.x[1] ^ c[1];
    t[3] = c[1];
#else
    t[1] = t00.x[1] ^ tk.x[0];
    t[2] = t11.x[0] ^ tk.x[1];
    t[3] = t11.x[1];
#endif
#ifdef CARRY
    c[0] = t11.x[0];
    c[1] = t11.x[1];
#endif
}
#endif  /* GF2X_MUL2_H_ */
