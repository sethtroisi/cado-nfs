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

#ifndef GF2X_MUL4_H_
#define GF2X_MUL4_H_

#include "gf2x.h"
/* All gf2x source files for lowlevel functions must include gf2x-small.h
 * This is mandatory for the tuning mechanism. */
#include "gf2x/gf2x-small.h"

/* specialized Karatsuba */
GF2X_STORAGE_CLASS_mul4
void gf2x_mul4 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
  unsigned long aa[2], bb[2], ab[4];
  unsigned long lo[4], hi[4];
  gf2x_mul2 (lo, a, b);
  gf2x_mul2 (hi, a + 2, b + 2);
  aa[0] = a[0] ^ a[2];
  aa[1] = a[1] ^ a[3];
  bb[0] = b[0] ^ b[2];
  bb[1] = b[1] ^ b[3];
  unsigned long c24 = lo[2] ^ hi[0];
  unsigned long c35 = lo[3] ^ hi[1];
  gf2x_mul2 (ab, aa, bb);
  c[0] = lo[0];
  c[1] = lo[1];
  c[2] = ab[0] ^ lo[0] ^ c24;
  c[3] = ab[1] ^ lo[1] ^ c35;
  c[4] = ab[2] ^ hi[2] ^ c24;
  c[5] = ab[3] ^ hi[3] ^ c35;
  c[6] = hi[2];
  c[7] = hi[3];
}

#endif  /* GF2X_MUL4_H_ */
