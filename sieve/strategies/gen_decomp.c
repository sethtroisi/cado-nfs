/* gen_decomp mfb fbb computes approximations of all decompositions of
   mfb-bit integers into primes larger than fbb. Example:
$ gen_decomp 60 524288
60 20 40 4.700208e+14 # 1
60 21 39 4.596914e+14 # 2
60 22 38 4.502935e+14 # 3
60 23 37 4.420430e+14 # 4
60 24 36 4.351291e+14 # 5
60 25 35 4.297433e+14 # 6
60 26 34 4.254332e+14 # 7
60 27 33 4.214013e+14 # 8
60 28 32 4.197615e+14 # 9
60 29 31 4.183856e+14 # 10
60 30 30 2.088895e+14 # 11
60 20 41 5.930007e+14 # 12
60 21 40 5.788613e+14 # 13
60 22 39 5.664933e+14 # 14
60 23 38 5.560455e+14 # 15
60 24 37 5.460721e+14 # 16
60 25 36 5.393045e+14 # 17
60 26 35 5.322997e+14 # 18
60 27 34 5.281604e+14 # 19
60 28 33 5.249884e+14 # 20
60 29 32 5.234872e+14 # 21
60 30 31 5.197604e+14 # 22
60 20 20 20 2.477904e+12 # 23
60 20 20 21 3.529032e+13 # 24
60 20 20 22 1.034259e+13 # 25
60 20 21 21 1.030657e+13 # 26
*/

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "macros.h"

#define PREC 1000000

/* number of primes <= 2^n */
unsigned long A7053[256] = {0,1,2,4,6,11,18,31,54,97,172,309,564,1028,1900,
                            3512,6542,12251,23000,43390,82025,155611,295947,
                            564163,1077871,2063689,3957809,7603553,14630843,
                            28192750,54400028,105097565,203280221,393615806,
                            762939111,1480206279,2874398515,5586502348,
                            10866266172,21151907950,41203088796,80316571436,
                            156661034233,305761713237,597116381732,
                            1166746786182,2280998753949,4461632979717,
                            8731188863470,17094432576778,33483379603407,
                            65612899915304,128625503610475,0,};

/* return prime_pi(2^i) */
double
prime_pi (unsigned long i)
{
  if (A7053[i] != 0)
    return (double) A7053[i];
  else
    {
      double x = ldexp (1.0, i);
      return x / log (x);
    }
}

#if 0

/* generate a random i-bit integer, uniformly distributed in [2^(i-1), 2^i[ */
double
gen_random (unsigned long i)
{
  double min;

  min = ldexp (1.0, i - 1);
  return min * (1.0 + drand48 ());
}

#else /* better version: takes into account the 1/log(x) bias */

/* generate a random i-bit integer, distributed with density in 1/log(x). */
double
gen_random (unsigned long i)
{
  double n0 = prime_pi (i - 1);
  double n1 = prime_pi (i);
  double a, b, n;

  /* we assume the n-th prime is in a*n*log(n)+b, thus we want:
     a*n0*log(n0) + b = 2^(i-1)
     a*n1*log(n1) + b = 2^i */
  a = ldexp (1.0, i - 1) / (n1 * log (n1) - n0 * log (n0));
  b = ldexp (1.0, i) - a * n1 * log (n1);
  n = n0 + (n1 - n0) * drand48 ();
  return a * n * log (n) + b;
}

#endif

double T[256] = {0,};

/* estimates number of products of primes:
   l[0]-bit * l[1]-bit * ... * l[k-1]-bit
   with product of mfb bits, each prime should be >= lim */
void
doit (unsigned int *l, unsigned int k, unsigned long mfb, unsigned long lim)
{
  unsigned long ok, tot, i, j;
  double S, p, r, rmin, rmax;
  static int count = 0;

  ok = tot = 0;
  for (S = 1.0, i = 0; i < k; S *= T[l[i]], i++);
  if (S == 0.0)
    return;
  rmin = ldexp (1.0, mfb - 1);
  rmax = ldexp (1.0, mfb);
  while (tot < PREC)
    {
      for (r = 1.0, i = 0; i < k; i++)
        {
          p = gen_random (l[i]);
          if (p >= lim)
            r *= p;
          else
            r = 0.0;
        }
      tot ++;
      if (rmin <= r && r <= rmax)
        ok ++;
    }
  printf ("%lu", mfb);
  for (i = 0; i < k; i++)
    {
      for (j = 1; i >= j && l[i-j] == l[i]; j++);
      S /= (double) j;
    }
  for (i = 0; i < k; i++)
    printf (" %d", l[i]);
  printf (" %e # %d\n", (double) ok * S / (double) tot, ++count);
  fflush (stdout);
}

int
main (int argc, char *argv[])
{
  unsigned long mfb, lim, i, j, k, kmax;
  unsigned int l[256], u;
  int imin, t, s;
  double p0, p1;

  ASSERT_ALWAYS (argc == 3);
  mfb = atoi (argv[1]);
  lim = atoi (argv[2]);
  imin = ceil (log2 ((double) (lim + 1)));
  for (i = 1; (i < 256) && i <= mfb; i++)
    {
      p0 = ldexp (1.0, i - 1);
      p0 = (p0 < (double) lim) ? (double) lim / log ((double) lim)
        : prime_pi (i - 1);
      p1 = ldexp (1.0, i);
      p1 = (p1 <= (double) lim) ? (double) lim / log ((double) lim)
        : prime_pi (i);
      T[i] = p1 - p0;
      // if (T[i] != 0) printf ("T[%lu]=%.0f\n", i, T[i]);
    }
  kmax = floor (log (ldexp (1.0, mfb)) / log ((double) lim));
  for (k = 2; k <= kmax; k++)
    {
      /* a product l[0]-bit * l[1]-bit * ... * l[k-1]-bit can have
         l[0] + l[1] + ... + l[k-1] - (k-1) bits */
      for (j = 0; j < k; j++)
        {
          for (i = 0; i < k - 1; i++)
            l[i] = imin;
          l[k-1] = mfb + j - (k - 1) * imin;
          while (1)
            {
              doit (l, k, mfb, lim);
              /* next partition of mfb + j */
              for (t = k - 1; t > 0; t--)
                /* try to increase l[t-1] and decrease l[t] */
                if (l[t-1] + 1 <= l[t])
                  {
                    l[t-1] ++;
                    s = 1;
                    for (u = t; u < k; u++)
                      {                  
                        s -= l[u];
                        l[u] = l[u-1];
                        s += l[u];
                      }
                    /* the total has increased by s */
                    l[k-1] -= s;
                    if (l[k-2] <= l[k-1])
                      break;
                  }
              if (t == 0)
                break;
            }
        }    
    }
  return 0;
}
