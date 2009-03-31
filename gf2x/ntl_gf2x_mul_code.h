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

/* Include this from within NTL's GF2X.c in order to replace the NTL
 * multiplication routine for binary polynomials.
 */

#include "gf2x.h"

void mul(GF2X& c, const GF2X& a, const GF2X& b)
{
   long sa = a.xrep.length();
   long sb = b.xrep.length();

   if (sa <= 0 || sb <= 0) {
      clear(c);
      return;
   }
   c.xrep.SetLength(sa+sb);
   gf2x_mul(c.xrep.elts(), a.xrep.elts(), sa, b.xrep.elts(), sb);
   c.normalize();
}
