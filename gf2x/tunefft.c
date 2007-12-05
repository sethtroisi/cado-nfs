/* Program to tune the FFT multiplication over GF(2).
   
  Copyright 2007 Paul Zimmermann.

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

   0) First run tunetoom to tune Toom-Cook multiplication (see the
      instructions in tunetoom.c).
   1) compile this program (tunefft) in the NTL source directory.
   2) run tunefft, giving as argument the maximum word size, for example
      ./tunefft 1000000 will tune multiplication of polynomials up to
      degree 64 million on a 64-bit machine.
      This prints detailed logs on stdout, and puts the tuning table in 
      the file mparam.h.
   3) compile factor.c.
*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h> /* for LONG_MAX */
#include <NTL/GF2X.h>

NTL_CLIENT

/* This version of tunefft uses the midpoint of each stair */
   
#define MUL_FFT_THRESHOLD 1000		/* Must be at least 28, but it
					   saves time to set larger */
#define STEPMAX 50

#include "HalfGCD.c"

#define MINTIME 0.5 /* time resolution */

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

/* Return the largest m >= n such that a product of m*m words with an FFT of
   length K leads to the same size of pointwise products as with n words.
   More precisely if f(n) = K/3 * ceil(6*n*NTL_BITS_PER_LONG/K^2) then m is the
   largest integer such that
   ceil(f(n)/NTL_BITS_PER_LONG) = ceil(f(m)/NTL_BITS_PER_LONG).
   In fact we take for m the largest integer such that g(n) = g(m) with
   g(n) = ceil(6*n*NTL_BITS_PER_LONG/K^2).
*/

long end_of_stair (long n, long K)
{
  long g = 6 * n * NTL_BITS_PER_LONG;
  
  g = 1 + (g - 1) / K; /* ceil(g/K) */
  g = 1 + (g - 1) / K; /* ceil(g/K^2) */
  g = g * K * K;
  return g / (6 * NTL_BITS_PER_LONG);
}

long next_step (long n, long K)
{
  long n1 = end_of_stair (n, K);		// Next stair
  long n2 = n + (n + STEPMAX - 1)/STEPMAX;	// ceil(n*(1.0 + 1.0/STEPMAX)
  if (n1 < n2)
    return n1;					// Return the minimum so
  else						// steps are not too large
    return n2;  
}

void random (_ntl_ulong *a, long n)
{
  for (long i = 0; i < n; i++)
    a[i] = RandomWord ();
}

void check (const _ntl_ulong *a, const _ntl_ulong *b, long n, 
	long K, long flag)
{
  for (long i = 0; i < n; i++)
    {
    if (a[i] != b[i])
      {
      printf("\nError detected: FFT%ld (K=%ld) and Toom differ at %ld\n",
        flag, K, n);
      exit (1);  
      }
    }  
} 

int main (int argc, char *argv[])
{
  long maxn, mid, n, n1, n2, ns, i, fpkt;
  long besti; /* 0 for TC, 1, 2, ... for FFT(K0*3^(bestK-1)) */
  long bestK, oldbestK = -1;
  long K, K0 = 3; /* try K0, 3*K0, 9*K0 */
  double T[4]; /* T[0] is for TC, T[1] for K0, T[2] for 3*K0, T[3] for 9*K0 */
  double t1[4], t2[4];
  _ntl_ulong *a, *b, *c, *t, *u;
  FILE *fp;

  if (argc != 2)
    {
      fprintf (stderr, "Usage: tunefft max_word_size\n");
      exit (1);
    }

  maxn = atoi (argv[1]);
  if (maxn <= 0) maxn = 1000000;		// default

  printf ("Tuning FFT multiplication to wordsize %ld\n\n", maxn);

  fp = fopen ("mparam.h", "w");
  fpkt = 0;
  fprintf (fp, "/* file automatically generated, do not edit  */\n");
  fprintf (fp, "/* {n, K} means use FFT(|K|) up from n words, */\n"); 
  fprintf (fp, "/* where K=1 stands for Toom-Cook 3, K < 0 means use FFT2 */\n");
  fprintf (fp, "#define MUL_FFT_TABLE {");
  fflush (fp);

  a = (_ntl_ulong*) malloc (maxn * sizeof (_ntl_ulong));
  b = (_ntl_ulong*) malloc (maxn * sizeof (_ntl_ulong));
  c = (_ntl_ulong*) malloc (2 * maxn * sizeof (_ntl_ulong));
  u = (_ntl_ulong*) malloc (2 * maxn * sizeof (_ntl_ulong));
  t = (_ntl_ulong*) malloc (toomspace(maxn) * sizeof (_ntl_ulong));

  random (a, maxn);
  random (b, maxn);

/* Skip n if (2*n < MUL_FFT_THRESHOLD) as this is too small for the FFT */

  for (n = MUL_FFT_THRESHOLD/2 + 1; n <= maxn;)
    {
      n2 = next_step (n, 3 * K0);			// End of interval
      if (n2 > maxn)					// Only go as far
        n2 = maxn;					// as maxn.
      mid = (n + n2)/2;					// Mid-point
      printf ("%ld..%ld ", n, n2);
      fflush (stdout);
        
      TIME (T[0], Toom (u, a, b, mid, t));		// Time Toom-Cook
      printf ("TC:%1.1e ", T[0]);
      fflush (stdout);
      besti = 0;
      bestK = 1;
      for (K = K0, i = 1; i <= 3; i++)
	{
	  TIME (t1[i], FFTMul (c, a, mid, b, mid, K));
	  check (c, u, 2*mid, K, 1);			// Compare results			
	  TIME (t2[i], FFTMul2 (c, a, mid, b, mid, K));
	  check (c, u, 2*mid, K, 2);			// Compare results
	  if (t1[i] < t2[i])
	    {
	    T[i] = t1[i];
	    printf ("F1(%ld):%1.1e ", K, T[i]);
	    }
	  else
	    {
	    T[i] = t2[i];
	    printf ("F2(%ld):%1.1e ", K, T[i]);
	    }
	  fflush (stdout);
	  if (T[i] < T[besti])
	    {
	    besti = i;
	    bestK = (t2[i] > t1[i]) ? K: -K;	/* -K for FFT2(|K|) */
	    }
	  K *= 3;
	}
//    printf ("best:");
      if (bestK == 1)
	printf ("TC");
      else
        {
        if (bestK > 0)
  	  printf ("F1(%ld)", bestK);
  	else
  	  printf ("F2(%ld)", -bestK);  
	}
      printf ("\n");
      fflush (stdout);

      if (bestK != oldbestK)
	  n1 = (n == MUL_FFT_THRESHOLD/2 + 1) ? 1 : n; 

      if (T[3] < T[1] && T[3] < T[2])
	K0 *= 3;
      else if (T[1] < T[2] && T[1] < T[3] && K0 > 3)
	K0 /= 3;
	
      /* go to next size */
      ns = n;
      n = next_step (n, 3 * K0);    /* middle value of K */
      if (n > n2) n = n2;	    /* end of last stair if K0 increased */
      n++;
      if (n < mid)                  /* revert to former n if K0 decreased */
        n = ns;
      else
        {
        if (bestK != oldbestK)
	  { 
	  if (n1 == 1) bestK = 1;
	  fprintf (fp, "{%ld, %ld}, ", n1, bestK);
	  fpkt++;
	  if ((fpkt%4) == 0)
	    fprintf (fp, "\\\n      ");
	  fflush (fp);
	  oldbestK = bestK;
	  }
	}
    }

  fprintf (fp, "{%ld, 0}}\n", LONG_MAX);
  fclose (fp);

  free (a);
  free (b);
  free (c);
  free (t);
  free (u);

  printf ("tunefft finished, output in mparam.h\n");
  fflush (stdout);
  
  return 0;
}
