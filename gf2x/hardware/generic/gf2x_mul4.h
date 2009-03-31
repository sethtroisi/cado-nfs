/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009
   Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2 of the License, or (at your
   option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
   more details.

   You should have received a copy of the GNU General Public License along
   with this program; see the file COPYING.  If not, write to the Free
   Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
   02111-1307, USA.
*/

#ifndef GF2X_MUL4_H_
#define GF2X_MUL4_H_

#include "gf2x.h"
#include "gf2x/gf2x-small.h"

/* specialized Karatsuba */
GF2X_STORAGE_CLASS_mul4
void gf2x_mul4 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
  unsigned long aa[2], bb[2], ab[4], c24, c35;
  gf2x_mul2 (c, a, b);
  gf2x_mul2 (c + 4, a + 2, b + 2);
  aa[0] = a[0] ^ a[2];
  aa[1] = a[1] ^ a[3];
  bb[0] = b[0] ^ b[2];
  bb[1] = b[1] ^ b[3];
  c24 = c[2] ^ c[4];
  c35 = c[3] ^ c[5];
  gf2x_mul2 (ab, aa, bb);
  c[2] = ab[0] ^ c[0] ^ c24;
  c[3] = ab[1] ^ c[1] ^ c35;
  c[4] = ab[2] ^ c[6] ^ c24;
  c[5] = ab[3] ^ c[7] ^ c35;
}

#endif  /* GF2X_MUL4_H_ */
