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

#ifndef GF2X_MUL8_H_
#define GF2X_MUL8_H_

#include "gf2x.h"
/* All gf2x source files for lowlevel functions must include gf2x-small.h
 * This is mandatory for the tuning mechanism. */
#include "gf2x/gf2x-small.h"

/* 3 calls to mul4, i.e., 27 multiplications */
GF2X_STORAGE_CLASS_mul8
void gf2x_mul8 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
    /* specialized Karatsuba, RPB 20070518 */
    /* slightly faster on bogong than version with loops */
    /* this version uses minimal temporary storage (12 = 3*n/2 words) */
    unsigned long aa[4], bb[4], cc[4];
    gf2x_mul4 (c+8, a+4, b+4);
    gf2x_mul4 (c, a, b);
    cc[0] = c[4] ^ c[8];
    cc[1] = c[5] ^ c[9];
    cc[2] = c[6] ^ c[10];
    cc[3] = c[7] ^ c[11];
    aa[0] = a[0] ^ a[4];
    aa[1] = a[1] ^ a[5];
    aa[2] = a[2] ^ a[6];
    aa[3] = a[3] ^ a[7];
    bb[0] = b[0] ^ b[4];
    bb[1] = b[1] ^ b[5];
    bb[2] = b[2] ^ b[6];
    bb[3] = b[3] ^ b[7];
    gf2x_mul4 (c+4, aa, bb);
    c[4]  ^= c[0]  ^ cc[0];
    c[5]  ^= c[1]  ^ cc[1];
    c[6]  ^= c[2]  ^ cc[2];
    c[7]  ^= c[3]  ^ cc[3];
    c[8]  ^= c[12] ^ cc[0];
    c[9]  ^= c[13] ^ cc[1];
    c[10] ^= c[14] ^ cc[2];
    c[11] ^= c[15] ^ cc[3];
}

#endif  /* GF2X_MUL8_H_ */
