/* Program to tune Toom-Cook multiplication over GF(2).
   
  Copyright 2007 Richard Brent and Paul Zimmermann.

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

/* How to use this program:

   1) Back up the files besttoom.h and mparam.h as these will
      be overwritten by tunetoom and tunefft respectively.
   2) compile this program (tunetoom) in the NTL source directory.
      You will need to copy the mulfft-bit.c and HalfGCD.c files to this
      directory.  (And some other files?)
   3) run tunetoom, giving as argument the maximum word size, for example
      ./tunetoom 2000 will tune multiplication of polynomials up to
      degree 128000 on a 64-bit machine. For higher degrees the FFT
      is probably faster - see tunefft.
      tunetoom prints detailed logs on stdout, and writes the tuning
      information in the file besttoom.h.  There are two phases - first
      the "balanced" routines with equal-sized inputs are tuned, then the 
      "unbalanced" routines where one input is about twice as large as the 
      other are tuned.
   4) compile and run tunefft to tune FFT multiplication
      (see instructions in tunefft.c).
*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h> /* for LONG_MAX */
#include <assert.h>
#include <NTL/GF2X.h>

NTL_CLIENT

#define TUNING

#define BESTMAX 2048

short best_tab[BESTMAX];

short best_utab[] = {0};

#include "HalfGCD.c"

#define BESTMIN (MUL_KARA_THRESHOLD-1)

#define BESTMINU (MUL_TOOMU_THRESHOLD-1)

#define MINTIME 0.5 /* time resolution */

FILE *fp;

#define TIME(x, i)				\
  { long j, k = 1;				\
    double s0 = GetTime ();			\
    do {					\
      for (j = 0; j < k; j++) (i);		\
      k = 2 * k;				\
      x = GetTime () - s0;			\
    } while (x < MINTIME);			\
    (x) = (x) / (k - 1);			\
  }

void random (_ntl_ulong *a, long n)
{
  for (long i = 0; i < n; i++)
    a[i] = RandomWord ();
}

void check (const _ntl_ulong *a, const _ntl_ulong *b, long n, long flag)
{
  for (long i = 0; i < n; i++)
    {
    if (a[i] != b[i])
      {
      if (flag == 2)
        printf("Error detected: Toom3Mul and KarMul give different results\n");
      if (flag == 3)
        printf("Error detected: Toom3WMul and KarMul give different results\n");
      if (flag == 4)
        printf("Error detected: Toom4Mul and KarMul give different results\n");
      if ((flag < 2) || (flag > 4))
        printf("Error in call to check, illegal flag %ld\n", flag);
      exit (1);  
      }
    }  
} 

void tunetoom (long tablesz)
{
  extern short best_tab[BESTMAX];
  long high, n, k, i;
  double T3[1], TK[1], TW[1], T4[1];
  double mint;
  _ntl_ulong *a, *b, *c, *d, *t;

  high = tablesz;
  if (high < BESTMIN) high = BESTMIN;
    
  if (high > BESTMAX)
    {
    printf ("Increase constant BESTMAX in tunetoom.c to %ld\n", high);
    exit (1);
    }
    
  for (n = 1; n <= BESTMAX; n++)
    best_tab[n-1] = 1;  
  
  printf ("Generating entries best_tab[%ld..%ld]\n", 1, high);
  fflush (stdout);

  a = (_ntl_ulong*) malloc (high * sizeof (_ntl_ulong));
  b = (_ntl_ulong*) malloc (high * sizeof (_ntl_ulong));
  c = (_ntl_ulong*) malloc (2 * high * sizeof (_ntl_ulong));
  d = (_ntl_ulong*) malloc (2 * high * sizeof (_ntl_ulong));
  t = (_ntl_ulong*) malloc (toomspace(high) * sizeof (_ntl_ulong));

  for (n = BESTMIN+1; n <= high; n++)

      {
      TK[0] = T3[0] = TW[0] = T4[0] = 0.0;
      printf ("%ld ", n);
      fflush (stdout);
      random (a, n);
      random (b, n);
      if (n >= MUL_KARA_THRESHOLD)
        TIME (TK[0], KarMul (c, a, b, n, t));
      if (n >= MUL_TOOM_THRESHOLD)
        {
        TIME (T3[0], Toom3Mul (d, a, b, n, t));
        check (c, d, 2*n, 2);
        }
      if (n >= MUL_TOOMW_THRESHOLD)
        {
        TIME (TW[0], Toom3WMul (d, a, b, n, t));
        check (c, d, 2*n, 3);
        }
      if (n >= MUL_TOOM4_THRESHOLD)   
        {
        TIME (T4[0], Toom4Mul (d, a, b, n, t));
        check (c, d, 2*n, 4);
        }
      printf ("TC2:%1.2e TC3:%1.2e TC3W:%1.2e TC4:%1.2e ", 
      	       TK[0], T3[0], TW[0], T4[0]);
      mint = TK[0];
      k = 1;
      if ((T3[0] < mint) && (n >= MUL_TOOM_THRESHOLD))
        { mint = T3[0]; k = 2; }
      if ((TW[0] < mint) && (n >= MUL_TOOMW_THRESHOLD)) 
        { mint = TW[0]; k = 3; }
      if ((T4[0] < mint) && (n >= MUL_TOOM4_THRESHOLD))
        { mint = T4[0]; k = 4; }
      best_tab[n-1] = (short)k;			// Final value of best_tab[n]  
      printf ("best:%1.2e ", mint);
      if (k == 1) printf("TC2");
      if (k == 2) printf("TC3");
      if (k == 3) printf("TC3W");
      if (k == 4) printf("TC4");
      printf ("\n");
      fflush (stdout);
      }

  free (a);
  free (b);
  free (c);
  free (d);
  free (t);

  while (best_tab[high-1] == 4) 
    high--;				// No need to include final "4" entries

  fprintf (fp, "/* file automatically generated, edit at your peril */\n\n");
  fprintf (fp, "#define BEST_TOOM_TABLE {\\\n");

  for (n = 1; n <= high; n++) 
    {
    if ((n-1)%20 == 0)
      fprintf (fp, "      ");
    fprintf (fp, "%d", best_tab[n-1]);			// Write new values
    if (n == high)
      fprintf (fp, " }\n");
    else
      {
      fprintf (fp, ",");
      if (n%20 == 0)
        fprintf (fp, "\\\n");
      }
    }
  return;
}

/* Forms c := a*b where b has size sb, a has size sa = (sb+1)/2 (words),
   representing polynomials over GF(2).  c needs space for sa+sb words. 
   Needs space 2*sa + toomspace(sa) words in stk[0] ... 
   
   The code is essentially the same as in HalfGCD.c   */

static void mul21 (_ntl_ulong *c, const _ntl_ulong *b, long sb,
                    const _ntl_ulong *a, _ntl_ulong *stk)
              
{
  long i, j;
  long sa = (sb+1)/2;
  long sc = sa + sb;
  _ntl_ulong *v;
  v = stk;
  stk += 2*sa;
  for (i = 0; i < sc; i++)
    c[i] = 0; 
  do 
    {
    if (sa == 0)
      break;

    if (sa == 1)
      {
      c[sb] ^= AddMul1 (c, c, b, sb, a[0]);
      break;
      }
              
    for (i = 0; i+sa <= sb; i += sa)
      {
      Toom (v, a, b + i, sa, stk);	    // Generic Toom-Cook mult.
      for (j = 0; j < 2*sa; j++)
        c[i+j] ^= v[j];
      }  

    { const _ntl_ulong *t; t = a; a = b + i; b = t; }
    { long t; t = sa; sa = sb - i; sb = t; }
    c = c + i;
    }
  while (1);
}

void checku (const _ntl_ulong *a, const _ntl_ulong *b, long n)
{
  for (long i = 0; i < n; i++)
    {
    if (a[i] != b[i])
      {
      printf("Error detected: Toom3uMul and default give different results\n");
      printf("index %ld\n", i);
      exit (1);  
      }
    }  
} 

void tuneutoom (long tabsz)

{
  short best_utab[BESTMAX];
  long high, n, k, i;
  double T3[1], TK[1];
  double mint;
  _ntl_ulong *a, *b, *c, *d, *t;

  high = tabsz;
  if (high < BESTMINU) high = BESTMINU;
    
  if (high > BESTMAX)
    {
    printf ("Increase constant BESTMAX in tuneutoom.c to %ld\n", high);
    exit (1);
    }
    
  for (n = 1; n <= BESTMAX; n++)
    best_utab[n-1] = 0;  
  
  printf ("Generating entries best_utab[%ld..%ld]\n", 1, high);
  fflush (stdout);

  long sa = high;
  long sb = (sa+1)/2;

  long sp1 = toomuspace(sa);		// space for Toom3uMul
  long sp2 = toomspace(sb) + 2*sb;	// space for mul21
  long sp = (sp1 > sp2) ? sp1 : sp2;
      
  a = (_ntl_ulong*) malloc (sa * sizeof (_ntl_ulong));
  b = (_ntl_ulong*) malloc (sb * sizeof (_ntl_ulong));
  c = (_ntl_ulong*) malloc (3 * sb * sizeof (_ntl_ulong));
  d = (_ntl_ulong*) malloc (3 * sb * sizeof (_ntl_ulong));
  t = (_ntl_ulong*) malloc (sp * sizeof (_ntl_ulong));

  
  for (sa = BESTMINU; sa <= high; sa++)
      {
      sb = (sa+1)/2;
      random (a, sa);
      random (b, sb);
      TK[0] = T3[0] = 0.0;
      printf ("%ld ", sa);
      fflush (stdout);
      TIME (TK[0], mul21 (c, a, sa, b, t));
      if (sa >= MUL_TOOMU_THRESHOLD)
        {
        TIME (T3[0], Toom3uMul (d, a, sa, b, t));
        checku (c, d, sa+sb);
        }
      printf ("default:%1.2e TC3U:%1.2e ", TK[0], T3[0]);
      mint = TK[0];
      k = 0;
      if ((T3[0] < mint) && (sa >= MUL_TOOMU_THRESHOLD))
        { mint = T3[0]; k = 1; }
      best_utab[sa-1] = (short)k;		// Final value
      printf ("best:%1.2e ", mint);
      if (k == 0) printf("default");
      if (k == 1) printf("TC3U");
      printf ("\n");
      fflush (stdout);
      }

  free (a);
  free (b);
  free (c);
  free (d);
  free (t);

  while (best_utab[high-1] == 1) 
    high--;				// No need to include final "1" entries

  fprintf (fp, "\n#define BEST_UTOOM_TABLE {\\\n");

  for (n = 1; n <= high; n++) 
    {
    if ((n-1)%20 == 0)
      fprintf (fp, "      ");
    fprintf (fp, "%d", best_utab[n-1]);			// Write new values
    if (n == high)
      fprintf (fp, " }\n");
    else
      {
      fprintf (fp, ",");
      if (n%20 == 0)
        fprintf (fp, "\\\n");
      }
    }

  return;
}

int main (int argc, char *argv[])
{
  long tabsz1, tabsz2;
  if ((argc != 2) && (argc != 3))
    {
      fprintf (stderr, "Usage: tunetoom table-size1 [table-size2]\n");
      fprintf (stderr, " where %ld <= table-size1 <= %ld\n", 
               BESTMIN, BESTMAX);
      fprintf (stderr, " and  %ld <= table-size2 <= %ld\n", 
               BESTMINU, BESTMAX);
      exit (1);
    }

  tabsz1 = tabsz2 = atoi (argv[1]);
  if (argc == 3)
  tabsz2 = atoi (argv[2]);

  fp = fopen ("besttoom.h", "w");

  tunetoom (tabsz1);		// Tune balanced routines

  fflush (fp);
  
  tuneutoom (tabsz2);		// Tune unbalanced routines
  
  fclose (fp);

  printf ("tunetoom finished, output in besttoom.h\n");
  fflush (stdout);

  return 0;
}

