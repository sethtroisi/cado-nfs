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

#ifndef tuning_GF2X_MUL2_H_
#define tuning_GF2X_MUL2_H_

#include "gf2x.h"
#include "gf2x/gf2x-small.h"

#ifdef  TUNING
#undef  GF2X_STORAGE_CLASS_mul2
#define GF2X_STORAGE_CLASS_mul2 /**/
#endif

GF2X_STORAGE_CLASS_mul2 void tuning_gf2x_mul2(unsigned long *c, const unsigned long *a, const unsigned long *b)
{
   unsigned long t;
   unsigned long u[2];

   gf2x_mul1 (c, a[0], b[0]);
   gf2x_mul1 (c+2, a[1], b[1]);
   t    = c[1]^c[2];
   gf2x_mul1 (u, a[0]^a[1], b[0]^b[1]);
   c[1] = c[0]^u[0]^t;
   c[2] = c[3]^u[1]^t;
}

#endif  /* tuning_GF2X_MUL2_H_ */
