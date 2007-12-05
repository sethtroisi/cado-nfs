/* Generates low-level multiplication routines over GF(2)[x].

  Copyright 2003, 2007 Pierrick Gaudry and Paul Zimmermann.

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

#include <stdio.h>
#include <stdlib.h>

int
main (int argc, char *argv[])
{
  long i;
  unsigned long NTL_BITS_PER_LONG, K, CHOP, REM, MASK, fn, mask2;

  if (argc != 4)
    {
      fprintf (stderr, "Usage: %s NTL_BITS_PER_LONG k fn\n", argv[0]);
      fprintf (stderr, "       where fn is 0 for mul1, 1 for Mul1, and 2 for AddMul1.\n");
      exit (1);
    }
  
  NTL_BITS_PER_LONG = atoi (argv[1]);
  K = atoi (argv[2]);
  fn = atoi (argv[3]);

  printf ("/* This file was generated automatically with\n");
  printf ("   %s %u %u %u. Don't edit it! */\n\n",
          argv[0], NTL_BITS_PER_LONG, K, fn);

  printf ("\n/*\n"
  "  This program is free software; you can redistribute it and/or modify it\n"
  "  under the terms of the GNU General Public License as published by the\n"
  "  Free Software Foundation; either version 2 of the License, or (at your\n"
  "  option) any later version.\n"
  "\n"
  "  This program is distributed in the hope that it will be useful, but WITHOUT\n"
  "  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or\n"
  "  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for\n"
  "  more details.\n"
  "\n"
  "  You should have received a copy of the GNU General Public License along\n"
  "  with this program; see the file COPYING.  If not, write to the Free\n"
  "  Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA\n"
  "  02111-1307, USA.\n"
          "*/\n\n");

  if (fn == 0)
    {
      printf ("static inline void\n");
      printf ("mul1 (_ntl_ulong *c, _ntl_ulong a, _ntl_ulong b)\n");
    }
  else if (fn == 1)
    {
      printf ("static _ntl_ulong\n");
      printf ("Mul1 (_ntl_ulong *cp, const _ntl_ulong *bp, long sb, _ntl_ulong a)\n");
    }
  else if (fn == 2)
    {
      printf ("static _ntl_ulong\n");
      printf ("AddMul1 (_ntl_ulong *dp, const _ntl_ulong *cp, const _ntl_ulong* bp, long sb,\n");
      printf ("        _ntl_ulong a)\n");
    }
  else
    {
      fprintf (stderr, "Error: unexpected fn\n");
      exit (1);
    }
  printf ("{\n");
  if (fn > 0) /* Mul1/AddMul1 */
    {
      printf ("   long i;\n");
      printf ("   _ntl_ulong carry = 0, b;\n");
    }
  printf ("   _ntl_ulong hi, lo, tmp, A[%u];\n", 1<<K);
  printf ("\n");
  printf ("   A[0] = 0;\t\t");
  /* do not truncate: fix non-considered bit products afterwards */
  printf ("A[1] = a;\n");
  for (i = 2; i < (1 << K); i++)
    {
      if (i % 2 == 0)
        printf ("   A[%u] = A[%u] << 1;\t", i, i / 2);
      else
        printf ("A[%u] = A[%u] ^ a;\n", i, i - 1);
    }
  printf ("\n");

  if (fn > 0)
    {
      printf ("   for (i = 0; i < sb; i++)\n");
      printf ("     {\n");
      printf ("       b = bp[i];\n");
    }
  REM = NTL_BITS_PER_LONG - 2 * K;
  MASK = (1UL << K) - 1UL;
  for (CHOP = 0; CHOP + 2 * K < NTL_BITS_PER_LONG; CHOP += 2 * K);
  CHOP = NTL_BITS_PER_LONG - CHOP;
  /* now 1 <= CHOP <= 2 * K, with NTL_BITS_PER_LONG - CHOP multiple of 2K */
  if (CHOP <= K) /* one slice is enough for the upper part */
    printf ("   lo = A[b >> %u];\n", NTL_BITS_PER_LONG - CHOP);
  else /* two slices for the upper part */
    printf ("   lo = (A[b >> %u] << %u) ^ A[(b >> %u) & %u];\n",
            NTL_BITS_PER_LONG - (CHOP - K), K, NTL_BITS_PER_LONG - CHOP, MASK);
  printf ("   hi = lo >> %u;\n", NTL_BITS_PER_LONG - 2 * K);
  for (i = NTL_BITS_PER_LONG - CHOP - 2 * K; i >= 2 * K; i -= 2 * K)
    {
      printf ("   lo = (lo << %u) ^ (A[(b >> %u) & %u] << %u) ^ A[(b >> %u) & %u];\n",
              2 * K, i + K, MASK, K, i, MASK);
      printf ("   hi = (hi << %u) | (lo >> %u);\n", 2 * K, REM);
    }
  /* special case for i=0 since a shift of 0 is undefined */
  printf ("   lo = (lo << %u) ^ (A[(b >> %u) & %u] << %u) ^ A[b & %u];\n",
              2 * K, K, MASK, K, MASK);

  /* now perform the repair step for 2*K */
  MASK = 0;
  for (i = 0; i < NTL_BITS_PER_LONG; i += 2 * K)
    MASK |= 1UL << i;
  MASK = ~MASK;
  mask2 = MASK;
  printf ("\n");
  for (i = 1; i < 2 * K; i++)
    {
      /* bit NTL_BITS_PER_LONG-i from a was not considered in blocks of
         K bits from b for index j >= i */
      /* Gaudry's no-branching optimization */
      printf ("   tmp = -((a >> %u) & 1);\n", NTL_BITS_PER_LONG - i);
      printf ("   tmp &= ((b & 0x%lx) >> %u);\n", mask2, i);
      printf ("   hi = hi ^ tmp;\n");
      mask2 = (mask2 << 1) & MASK;
    }
  printf ("\n");
  if (fn == 0)
    {
      printf ("   c[0] = lo;\n");
      printf ("   c[1] = hi;\n");
    }
  else if (fn == 1)
    {
      printf ("   cp[i] = carry ^ lo;\n");
      printf ("   carry = hi;\n");
    }
  else if (fn == 2)
    {
      printf ("   dp[i] = cp[i] ^ (carry ^ lo);\n");
      printf ("   carry = hi;\n");
    }
  if (fn > 0)
    {
      printf ("    }\n");
      printf ("   return carry;\n");
    }
  printf ("}\n");
}
