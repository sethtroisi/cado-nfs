/* Returns the worst-case space (in words) needed by the Toom-Cook routines
   KarMul, Toom3Mul, Toom3WMul, Toom4Mul.
 
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

/* 
 The memory sp(n) necessary for Toom3WMul satisfies
 sp(n) <== (n lt 8) ? 19 : 8*(floor(n/3) + 3) + sp(floor(n/3) + 2),
 sp(7) <= 19.

 It is assumed that KarMul is called for n < 8 <= MUL_TOOMW_THRESHOLD
 and requires space KarMulMem(n) <= 4*ceil(n/2) + KarMulMem(ceil(n/2)),
 KarMulMem(7) <= 19.  The memory for Toom3Mul and Toom4Mul is no larger 
 than that for Toom3WMul.
 
 Note: KarMulMem(7) is now 0, but would increase if MUL_KARA_THRESHOLD
       were reduced. We have not changed ToomSpace as a small overestimate
       in space is not harmful.
*/

#ifndef MUL_KARA_THRESHOLD
#define MUL_KARA_THRESHOLD 10
#endif

#ifndef MUL_TOOM_THRESHOLD
#define MUL_TOOM_THRESHOLD 17
#endif

#ifndef MUL_TOOMW_THRESHOLD
#define MUL_TOOMW_THRESHOLD 10
#endif

#ifndef MUL_TOOMU_THRESHOLD
#define MUL_TOOMU_THRESHOLD 33
#endif

#if (MUL_KARA_THRESHOLD < 5)
#error "MUL_KARA_THRESHOLD assumed to be at least 5"
#endif

#if (MUL_TOOMW_THRESHOLD < 8)
#error "MUL_TOOMW_THRESHOLD assumed to be at least 8"
#endif

static long toomspace (long n)

  {
  long sp;
  long low = (MUL_KARA_THRESHOLD < MUL_TOOMW_THRESHOLD) ?
  	      MUL_KARA_THRESHOLD : MUL_TOOMW_THRESHOLD;
  if (n < low)
    return 0;
  sp = 19;					// KarMulMem (7) <= 19
  while (n >= 8)
    {
    n = n/3 + 2;
    sp += 8*(n+1);
    }
  return sp;
  }

/* Returns upper bound on space required by Toom3uMul (c, a, sa, b, stk): 
   2*sa + 32 + toomspace(sa/4 + 4) */

static long toomuspace (long sa)

  {
  if (sa < MUL_TOOMU_THRESHOLD)
    return 0;
  else  
    return 2*sa + 32 + toomspace(sa/4 + 4);  
  }
