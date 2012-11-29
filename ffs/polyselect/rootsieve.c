#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h> /* for ULONG_MAX */
#include <assert.h>
#include <float.h>

/* some hard-coded values (for now) */
#define B 63
#define DEGREE 6

float **s, **sl, threshold = -DBL_MAX;
unsigned long nsol = 0;
unsigned long f[DEGREE + 1] = {0, 0, 0, 0, 0, 0, 1};
uint8_t *G;

/* return the number of significant bits of l, 0 for l = 0 */
static unsigned long
nbits (unsigned long l)
{
  unsigned long d = 0;

  while (l > 0)
    d++, l >>= 1;
  return d;
}

static void
print_poly (unsigned long *f, unsigned long d)
{
  int first = 1;
  unsigned long n, fd;

  d = d + 1;
  while (d-- > 0)
    {
      fd = f[d];
      n = nbits (fd);
      if (n > 0)
        {
          if (first == 0)
            printf ("+");
          first = 0;
          if (fd == (1UL << (n - 1)))
            {
              if (n - 1 > 0)
                {
                  if (n - 1 == 1)
                    printf ("t*");
                  else
                    printf ("t^%lu*", n - 1);
                }
            }
          else
            {
              if (d > 0)
                printf ("(");
              if (n - 1 == 1)
                printf ("t");
              else
                printf ("t^%lu", n - 1);
              n--;
              while (n--)
                {
                  if ((fd >> n) & 1)
                    {
                      if (n == 0)
                        printf ("+1");
                      else if (n == 1)
                        printf ("+t");
                      else
                        printf ("+t^%lu", n);
                    }
                }
              if (d > 0)
                printf (")*");
            }
          if (d == 1)
            printf ("y");
          else if (d > 1)
            printf ("y^%lu", d);
        }
    }
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

static unsigned long
reduce (unsigned long r, unsigned long l)
{
  unsigned long n = nbits (l);
  unsigned long m = nbits (r);

  while (m >= n)
    {
      r = r ^ (l << (m - n));
      /* now bit m-1 of r should be 0 */
      m--;
      while (m >= n && ((r >> (m-1)) == 0))
        m--;
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
  r = reduce (r, l);
  return r;
}

/* initialize Gray code G[i] for 0 <= i < n:
   0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, ...
   where G[i] is the index of the highest bit changing between i and i+1 */
void
gray (unsigned long n)
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
is_maybe_irreducible (unsigned long *f, unsigned long d)
{
  /* does x+1 divide f? */
  if (evalmod (f, d, 1, ~0UL) == 0)
    return 0;

  /* does x+t divide f? */
  if (evalmod (f, d, 2, ~0UL) == 0)
    return 0;

  /* does x+t+1 divide f? */
  if (evalmod (f, d, 3, ~0UL) == 0)
    return 0;

  /* does x+t^2 divide f? */
  if (evalmod (f, d, 4, ~0UL) == 0)
    return 0;

  /* does x+t^2+1 divide f? */
  if (evalmod (f, d, 5, ~0UL) == 0)
    return 0;

  /* does x+t^2+t divide f? */
  if (evalmod (f, d, 6, ~0UL) == 0)
    return 0;

  /* does x+t^2+t+1 divide f? */
  if (evalmod (f, d, 7, ~0UL) == 0)
    return 0;

  return 1;
}

static void
all1 (unsigned long *f, unsigned long n1, unsigned long n0)
{
  unsigned long i, j, *l, d, D, r, f0, fj, K, N1 = 1 << n1, N0 = 1 << n0;
  unsigned long hl, k, m;
  float s0, s1, sj;
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

  f[0] = f[1] = 0;
  for (j = 0; j < N1; j++)
    for (i = 0; i < N0; i++)
      s[j][i] = 0.0;
  for (l = L;  *l != ULONG_MAX; l++)
    {
      d = nbits (*l) - 1; /* degree of l */
      D = 1UL << d;
      s0 = (float) d / (float) (D - 1);
      s1 = (float) d * (float) D / (float) (D * D - 1);
      for (j = 0; j < N1; j++)
        for (i = 0; i < D; i++)
          sl[j][i] = s0;
      for (r = 0; r < D; r++)
        {
          f0 = evalmod (f, DEGREE, r, *l);
          if (f0 > B)
            {
              printf ("f0=%lu r=%lu *l=%lu\n", f0, r, *l);
              print_poly (f, DEGREE);
              printf ("\n");
            }
          assert (f0 <= B);
          sl[0][f0] -= s1;
          for (j = 1; j < N1; j++)
            {
              fj = f0 ^ mulmod (r, j, *l);
              assert (j < N1);
              assert (fj <= B);
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
                if (m < N0)
                  s[j][m] += sj;
                hl ^= *l << G[k];
              }
          }
    }

  for (j = 0; j < N1; j++)
    /* we start from i=1 since for i=0, we have f0=0, thus f is divisible
       by y */
    for (i = 1; i < N0; i++)
      if (s[j][i] < threshold)
        {
          f[0] = i;
          f[1] = j;
          if (is_maybe_irreducible (f, DEGREE))
            {
              nsol ++;
              print_poly (f, DEGREE);
              printf (" %f\n", s[j][i]);
              fflush (stdout);
            }
        }
}

int
main (int argc, char *argv[])
{
  unsigned long i, j;
  unsigned long n0 = 0, n1 = 0, n2 = 0, N0, N1, N2, f2;

  while (argc > 1 && argv[1][0] == '-')
    {
      if (argc > 2 && strcmp (argv[1], "-threshold") == 0)
        {
          threshold = atof (argv[2]);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-n0") == 0)
        {
          n0 = atoi (argv[2]);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-n1") == 0)
        {
          n1 = atoi (argv[2]);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-n2") == 0)
        {
          n2 = atoi (argv[2]);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-f5") == 0)
        {
          f[5] = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else
        {
          fprintf (stderr, "Invalid option: %s\n", argv[1]);
          exit (1);
        }
    }

  N0 = 1 << n0;
  N1 = 1 << n1;
  N2 = 1 << n2;

  if (threshold == -DBL_MAX)
    {
      fprintf (stderr, "Missing -threshold option\n");
      exit (1);
    }

  /* assume B is a power of two - 1 */
  assert ((B & (B + 1)) == 0);
  assert (B <= N0);
  assert (N1 <= N0);
  s = (float**) malloc (N1 * sizeof (float*));
  for (j = 0; j < N1; j++)
    s[j] = (float*) malloc (N0 * sizeof (float));
  sl = (float**) malloc (N1 * sizeof (float*));
  for (j = 0; j < N1; j++)
    sl[j] = (float*) malloc ((B + 1) * sizeof (float));
  G = (uint8_t*) malloc ((N0/2) * sizeof (uint8_t));
  gray (N0/2);

  for (f2 = 0; f2 < N2; f2++)
    {
      f[2] = f2;
      all1 (f, n1, n0);
    }
  printf ("Found %lu polynomial(s) below threshold\n", nsol);

  free (G);
  for (j = 0; j < N1; j++)
    free (sl[j]);
  free (sl);
  for (j = 0; j < N1; j++)
    free (s[j]);
  free (s);
  
  return 0;
}
