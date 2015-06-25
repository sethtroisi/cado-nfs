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

#ifndef GF2X_MUL5_H_
#define GF2X_MUL5_H_

#include "gf2x.h"
/* All gf2x source files for lowlevel functions must include gf2x-small.h
 * This is mandatory for the tuning mechanism. */
#include "gf2x/gf2x-small.h"

    GF2X_STORAGE_CLASS_mul5
void gf2x_mul5 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
    /* specialized Karatsuba, RPB 20070518 */
    unsigned long aa[3], bb[3], ab[6], ab3, ab4, ab5;
    gf2x_mul2 (c+6, a+3, b+3);
    gf2x_mul3 (c, a, b);
    aa[0] = a[0] ^ a[3];
    aa[1] = a[1] ^ a[4];
    aa[2] = a[2];
    bb[0] = b[0] ^ b[3];
    bb[1] = b[1] ^ b[4];
    bb[2] = b[2];
    gf2x_mul3 (ab, aa, bb);
    ab3 = ab[3] ^ c[3];
    ab4 = ab[4] ^ c[4];
    ab5 = ab[5] ^ c[5];
    c[3] ^= ab[0] ^ c[0] ^ c[6];
    c[4] ^= ab[1] ^ c[1] ^ c[7];
    c[5] ^= ab[2] ^ c[2] ^ c[8];
    c[6] ^= ab3 ^ c[9];
    c[7] ^= ab4;
    c[8] ^= ab5;
}

#endif  /* GF2X_MUL5_H_ */
