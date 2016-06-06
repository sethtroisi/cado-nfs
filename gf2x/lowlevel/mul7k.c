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

#ifndef GF2X_MUL7_H_
#define GF2X_MUL7_H_

#include "gf2x.h"
/* All gf2x source files for lowlevel functions must include gf2x-small.h
 * This is mandatory for the tuning mechanism. */
#include "gf2x/gf2x-small.h"

/* specialized Karatsuba, RPB 20070518 */
/* slightly faster on bogong than version with loops */
    GF2X_STORAGE_CLASS_mul7
void gf2x_mul7 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
    unsigned long aa[4], bb[4], ab[8], ab4, ab5, ab6, ab7;
    gf2x_mul3 (c+8, a+4, b+4);
    gf2x_mul4 (c, a, b);
    aa[0] = a[0] ^ a[4];
    aa[1] = a[1] ^ a[5];
    aa[2] = a[2] ^ a[6];
    aa[3] = a[3];
    bb[0] = b[0] ^ b[4];
    bb[1] = b[1] ^ b[5];
    bb[2] = b[2] ^ b[6];
    bb[3] = b[3];
    gf2x_mul4 (ab, aa, bb);
    ab4 = ab[4] ^ c[4];
    ab5 = ab[5] ^ c[5];
    ab6 = ab[6] ^ c[6];
    ab7 = ab[7] ^ c[7];
    c[4] ^= ab[0] ^ c[0] ^ c[8];
    c[5] ^= ab[1] ^ c[1] ^ c[9];
    c[6] ^= ab[2] ^ c[2] ^ c[10];
    c[7] ^= ab[3] ^ c[3] ^ c[11];
    c[8] ^= ab4 ^ c[12];
    c[9] ^= ab5 ^ c[13];
    c[10] ^= ab6;
    c[11] ^= ab7;
}

#endif  /* GF2X_MUL7_H_ */
