/* 64-bit multiplication routine over GF(2)[x], using AVX.

  Copyright 2009 Paul Zimmermann.

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

/* assumes ulong is 64 bit */

#include "wmmintrin.h"

/* cf http://www.intel.com/software/products/compilers/docs/clin/main_cls/intref_cls/common/intref_aes_intrinsics.htm */

/* c[0,1] <- a * b */
static ulong
mul1 (ulong *c, ulong a, ulong b)
{
  __m128i v1[1] = {a, 0};
  __m128i v2[1] = {b, 0);
  __m128i v3;
  __m128i *v4 = (__m128i*) c;

  v3 = _mm_clmulepi64_si128(__m128i v1, __m128i v2, 0);
  v4[0] = v4;
}
