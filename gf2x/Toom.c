/* General Toom_Cook multiplication, calls KarMul, Toom3Mul, Toom3WMul
   or Toom4Mul depending on which is expected to be the fastest.

  Copyright 2007 Richard P. Brent.

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

/*   stk should point to a block of sufficient memory for any of these
     routines (toomspace(n) <= 5*n+17 words is enough).
     Output c must not overlap inputs a, b. 
     The output c is a*b (where a, b and c are in GF(2)[x]).
     RPB, 20070510 */
static inline void Toom (_ntl_ulong *c, const _ntl_ulong *a, 
			const _ntl_ulong *b, long n, _ntl_ulong *stk)
{
  while (n && a[n-1] == 0 && b[n-1] == 0)
    {
      c[2*n-1] = 0;
      c[2*n-2] = 0;
      n --;
    }

  switch (BestToom(n))
    {
    case 1:
      KarMul (c, a, b, n, stk);
      return;
    case 2:
      Toom3Mul (c, a, b, n, stk);
      return;
    case 3:
      Toom3WMul (c, a, b, n, stk);
      return;
    default:
      Toom4Mul (c, a, b, n, stk);
      return;      
    }
}
