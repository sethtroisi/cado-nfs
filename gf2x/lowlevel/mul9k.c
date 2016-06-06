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

#ifndef GF2X_MUL9_H_
#define GF2X_MUL9_H_

#include "gf2x.h"
/* All gf2x source files for lowlevel functions must include gf2x-small.h
 * This is mandatory for the tuning mechanism. */
#include "gf2x/gf2x-small.h"

/* 1 call to mul4 and 2 calls to mul5, i.e., 35 multiplications with
   mul5clk_c */
GF2X_STORAGE_CLASS_mul9
void gf2x_mul9 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
    /* specialized Karatsuba, RPB 20070520 */
    /* slightly faster on bogong than version with loops */
    unsigned long aa[5], bb[5], ab[10], ab5, ab6, ab7, ab8, ab9;
    gf2x_mul4 (c+10, a+5, b+5);
    gf2x_mul5 (c, a, b);
    aa[0] = a[0] ^ a[5];
    aa[1] = a[1] ^ a[6];
    aa[2] = a[2] ^ a[7];
    aa[3] = a[3] ^ a[8];
    aa[4] = a[4];
    bb[0] = b[0] ^ b[5];
    bb[1] = b[1] ^ b[6];
    bb[2] = b[2] ^ b[7];
    bb[3] = b[3] ^ b[8];
    bb[4] = b[4];
    gf2x_mul5 (ab, aa, bb);
    ab5 = ab[5] ^ c[5];
    ab6 = ab[6] ^ c[6];
    ab7 = ab[7] ^ c[7];
    ab8 = ab[8] ^ c[8];
    ab9 = ab[9] ^ c[9];
    c[5] ^= ab[0] ^ c[0] ^ c[10];
    c[6] ^= ab[1] ^ c[1] ^ c[11];
    c[7] ^= ab[2] ^ c[2] ^ c[12];
    c[8] ^= ab[3] ^ c[3] ^ c[13];
    c[9] ^= ab[4] ^ c[4] ^ c[14];
    c[10] ^= ab5 ^ c[15];
    c[11] ^= ab6 ^ c[16];
    c[12] ^= ab7 ^ c[17];
    c[13] ^= ab8;
    c[14] ^= ab9;
}

#endif  /* GF2X_MUL9_H_ */
