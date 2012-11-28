#include <stdio.h>
#include <stdlib.h>
#include <limits.h> /* for ULONG_MAX */
#include <assert.h>

/* some hard-coded values (for now) */
#define B 63
#define DEGREE 6
#define n0 13
#define N0 (1<<n0)
#define n1 11
#define N1 (1<<n1)
#define THRESHOLD -4.0

/* return the number of significant bits of l */
static unsigned long
nbits (unsigned long l)
{
  unsigned long d = 0;

  while (l > 0)
    d++, l >>= 1;
  return d;
}

/* return a * b mod l, where a is reduced mod l */
static unsigned long
mulmod (unsigned long a, unsigned long b, unsigned long l)
{
  unsigned long r, n;

  n = nbits (b);
  r = 0;
  while (n > 0)
    {
      r = r << 1;
      if ((b >> --n) & 1)
        r = r ^ a;
      if ((r ^ l) < r)
        r = r ^ l;
    }
  return r;
}

/* return f(v) mod l where f has degree d */
static unsigned long
evalmod (unsigned long *f, unsigned long d, unsigned long v, unsigned long l)
{
  unsigned long r;

  r = f[d];
  while (d-- > 0)
    {
      r = mulmod (r, v, l);
      r = r ^ f[d];
    }
  return r;
}

/* initialize Gray code G[i] for 0 <= i < n:
   0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, ...
   where G[i] is the index of the highest bit changing between i and i+1 */
void
gray (unsigned long *G, unsigned long n)
{
  unsigned long i, d, g;

  for (i = 0; i < n; i++)
    {
      d = (i + 1) ^ i;
      g = 0;
      while (d > 1)
        d >>= 1, g++;
      G[i] = g;
    }
}

int
main (int argc, char *argv[])
{
  /* list of irreducible polynomials of index <= B */
  unsigned long L[] = {2, /* t */
                       3, /* t + 1 */
                       7, /* t^2 + t + 1 */
                       11, /* t^3 + t + 1 */
                       13, /* t^3 + t^2 + 1 */
                       19, /* t^4 + t + 1 */
                       25, /* t^4 + t^3 + 1 */
                       31, /* t^4 + t^3 + t^2 + t + 1 */
                       37, /* t^5 + t^2 + 1 */
                       41, /* t^5 + t^3 + 1 */
                       47, /* t^5 + t^3 + t^2 + t + 1 */
                       55, /* t^5 + t^4 + t^2 + t + 1 */
                       59, /* t^5 + t^4 + t^3 + t + 1 */
                       61, /* t^5 + t^4 + t^3 + t^2 + 1 */
                       ULONG_MAX};
  unsigned long f[DEGREE + 1] = {0, 0, 0, 0, 0, 0, 1};
  float **s, s0, s1, **sl, sj;
  unsigned long i, j, *l, d, D, f0, fj, r, *G, K, k, hl, m, nsol;

  /* assume B is a power of two - 1 */
  assert ((B & (B + 1)) == 0);
  assert (B <= N0);
  assert (N1 <= N0);
  assert (f[0] == 0);
  assert (f[1] == 0);
  s = (float**) malloc (N1 * sizeof (float*));
  for (j = 0; j < N1; j++)
    {
      s[j] = (float*) malloc (N0 * sizeof (float));
      for (i = 0; i < N0; i++)
        s[j][i] = 0.0;
    }
  sl = (float**) malloc (N1 * sizeof (float*));
  for (j = 0; j < N1; j++)
    sl[j] = (float*) malloc ((B + 1) * sizeof (float));
  G = (unsigned long*) malloc ((N0/2) * sizeof (unsigned long));

  for (l = L;  *l != ULONG_MAX; l++)
    {
      d = nbits (*l) - 1; /* degree of l */
      D = 1UL << d;
      s0 = (float) d / (float) (D - 1);
      s1 = (float) d * (float) D / (float) (D * D - 1);
      gray (G, 1 << (n0 - d));
      for (j = 0; j < N1; j++)
        for (i = 0; i < D; i++)
          sl[j][i] = s0;
      for (r = 0; r < D; r++)
        {
          f0 = evalmod (f, DEGREE, r, *l);
          sl[0][f0] -= s1;
          for (j = 1; j < N1; j++)
            {
              fj = f0 ^ mulmod (r, j, *l);
              sl[j][fj] -= s1;
            }
        }
      K = 1 << (n0 - d);
      for (j = 0; j < N1; j++)
        for (r = 0, hl = 0; r < D; r++)
          {
            sj = sl[j][r];
            for (k = 0, hl = 0; k < K; k++)
              {
                m = hl ^ r;
                s[j][m] += sj;
                hl ^= *l << G[k];
              }
          }
    }

  for (nsol = j = 0; j < N1; j++)
    for (i = 1; i < N0; i++)
      if (s[j][i] < THRESHOLD)
        {
          nsol ++;
          printf ("f1=%lu f0=%lu %f\n", j, i, s[j][i]);
        }
  printf ("Found %lu polynomials below threshold\n", nsol);

  free (G);
  for (j = 0; j < N1; j++)
    free (sl[j]);
  free (sl);
  for (j = 0; j < N1; j++)
    free (s[j]);
  free (s);
  
  return 0;
}
