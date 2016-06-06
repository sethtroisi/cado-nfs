#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h> /* for ULONG_MAX */
#include <assert.h>
#include <float.h>
#include <sys/types.h>
#include <sys/resource.h>

#define USE_FQPOL
#ifdef USE_FQPOL
#include "types.h"
#include "ffspol.h"
#include "fqpol.h"
#include "params.h"
#include "polyfactor.h"
#endif

/* some hard-coded values (for now) */
#define DEGREE 6

float **s, **sl, threshold = -DBL_MAX;
unsigned long nsol = 0, B = 0;
unsigned long f[DEGREE + 1] = {0, 0, 0, 0, 0, 0, 1};
uint8_t *G;
int init_time = 0, sieve_time = 0;

/* return the number of significant bits of l, 0 for l = 0 */
static unsigned long
nbits (unsigned long l)
{
  unsigned long d = 0;

  while (l > 0)
    d++, l >>= 1;
  return d;
}

int
cputime ()
{
  struct rusage rus;

  getrusage (0, &rus);
  return rus.ru_utime.tv_sec * 1000 + rus.ru_utime.tv_usec / 1000;
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
          if (fd == (1UL << (n - 1))) /* only one coefficient t^k */
            {
              if (n - 1 > 0)
                {
                  if (n - 1 == 1)
                    printf ("t");
                  else
                    printf ("t^%lu", n - 1);
                  if (d > 0)
                    printf ("*");
                }
              else if (d == 0)
                printf ("1");
            }
          else /* two or more coefficients */
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

/* return a * b */
static unsigned long
mul (unsigned long a, unsigned long b)
{
  unsigned long r, n;

  n = nbits (b);
  r = 0;
  while (n > 0)
    {
      r = r << 1;
      if ((b >> --n) & 1)
        r = r ^ a;
    }
  return r;
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
is_divisible_by_y2_plus_a1y1_plus_a0 (unsigned long *f, unsigned long d,
                                      unsigned long a1, unsigned long a0)
{
  unsigned long g[DEGREE + 1], i;

  for (i = 0; i <= d; i++)
    g[i] = f[i];
  while (d >= 2)
    {
      /* g -= g[d]*y^d - g[d]*a1*y^(d-1) - g[d]*a0*y^(d-2) */
      g[d-1] ^= mul (g[d], a1);
      g[d-2] ^= mul (g[d], a0);
      d--;
    }
  return g[1] == 0 && g[0] == 0;
}

int
is_maybe_irreducible (unsigned long *f, unsigned long d)
{
  unsigned long l0[] = {1, 2, 3, 4, 5, 6, 7, ULONG_MAX}, i, a0, a1;
  unsigned long l1[][2] = {{1,1},{1,2},{1,3},{1,4},{1,5},{1,6},{1,7},
                      {2,0},{2,1},{2,2},{2,4},{2,5},{2,6},{2,7},
                      {3,0},{3,1},{3,3},{3,4},{3,5},{3,6},{3,7},
                      {4,1},{4,2},{4,3},{4,4},{4,6},{4,7},
                      {5,1},{5,2},{5,3},{5,5},{5,6},{5,7},
                      {6,0},{6,2},{6,3},{6,4},{6,5},{6,6},
                      {7,0},{7,1},{7,2},{7,3},{7,4},{7,5},{7,7},
                      {8,0},{8,1},{8,2},{8,3},{8,4},{8,5},{8,7},
                      {9,0},{9,1},{9,2},{9,3},{9,5},{9,6},{9,7},
                      {10,0},{10,1},{10,2},{10,3},{10,4},{10,6},
                      {11,0},{11,1},{11,2},{11,3},{11,4},{11,5},{11,6},{11,7},
                      {12,0},{12,1},{12,2},{12,3},{12,5},{12,6},
                      {13,0},{13,1},{13,2},{13,3},{13,4},{13,5},{13,6},{13,7},
                      {14,0},{14,1},{14,2},{14,3},{14,4},{14,6},{14,7},
                      {15,0},{15,1},{15,2},{15,3},{15,4},{15,5},{15,7},
                      {16,1},{16,2},{16,3},{16,4},{16,5},{16,6},{16,7},
                      {17,1},{17,2},{17,3},{17,4},{17,5},{17,6},{17,7},
                      {18,0},{18,2},{18,3},{18,4},{18,5},{18,6},{18,7},
                      {19,0},{19,1},{19,2},{19,3},{19,4},{19,5},{19,6},{19,7},
                      {20,2},{20,3},{20,4},{20,5},{20,6},{20,7},
                      {21,1},{21,2},{21,3},{21,4},{21,5},{21,6},{21,7},
                      {22,0},{22,1},{22,2},{22,3},{22,4},{22,5},{22,6},{22,7},
                      {23,0},{23,1},{23,2},{23,3},{23,4},{23,5},{23,6},{23,7},
                      {24,0},{24,1},{24,3},{24,4},{24,5},{24,6},{24,7},
                      {25,0},{25,1},{25,2},{25,3},{25,4},{25,5},{25,6},{25,7},
                      {26,0},{26,1},{26,2},{26,3},{26,4},{26,5},{26,6},{26,7},
                      {27,0},{27,1},{27,3},{27,4},{27,5},{27,6},{27,7},
                      {28,0},{28,1},{28,2},{28,4},{28,5},{28,6},{28,7},
                      {29,0},{29,1},{29,2},{29,3},{29,4},{29,5},{29,6},{29,7},
                      {30,0},{30,1},{30,2},{30,4},{30,5},{30,6},{30,7},
                      {31,0},{31,1},{31,2},{31,3},{31,4},{31,5},{31,6},{31,7},
                      {ULONG_MAX,0}};
  unsigned long l2[][3] = {{4,2,0},{5,3,0},{6,6,2},{16,2,3},{16,4,0},{16,7,3},{17,3,2},{17,5,0},{17,7,2},{19,7,2},{20,2,2},{20,3,3},{20,6,0},{21,7,0},{24,0,2},{30,0,3},{34,10,0},{35,3,2},{40,12,0},{42,14,0},{64,8,0},{65,9,0},{69,11,0},{81,13,0},{84,14,0},{85,15,0},{113,11,2},{114,12,3},{120,8,2},{120,15,3},{ULONG_MAX,0,0}};

  /* divisors of degree 1 in y */
  for (i = 0; (a0 = l0[i]) != ULONG_MAX; i++)
    /* does y+a0 divide f? */
    if (evalmod (f, d, a0, ~0UL) == 0)
      return 0;

  /* divisors of degree 2 in y */
  for (i = 0; (a0 = l1[i][0]) != ULONG_MAX; i++)
    if (is_divisible_by_y2_plus_a1y1_plus_a0 (f, d, l1[i][1], a0))
      return 0;

  return 1;
}

#ifdef USE_FQPOL
/* This function is duplicated from sieve4ffs/sieve.c */
static int
sq_is_irreducible (sq_srcptr p)
{
  fppol_t P;
  int ret;

  fppol_init (P);
  fppol_set_sq (P, p);
  ret = fppol_is_irreducible (P);
  fppol_clear (P);
  return ret;
}

/* cf files testirred.c and fqpol.c */
static int
is_irreducible (unsigned long *f, unsigned long d)
{
  ffspol_t ffspol;
  sq_t ell; /* polynomial used for testing */
  fq_info_t Fq;
  fqpol_t Fbar;
  unsigned int i;
  int ret;
  fppol64_t f64;
  int DEG_ELL = 13;
  int NTEST = 10;

  assert (d == 6);
  ffspol_init2 (ffspol, d + 1);
  for (i = 0; i <= d; i++)
    {
      fppol64_set_ui (f64, f[i], 64, 0);
      fppol_set_64 (ffspol->coeffs[i], f64);
    }
  ffspol->deg = d;

  sq_set_ti (ell, DEG_ELL);
  for (i = ret = 0; i < NTEST && ret == 0; i++)
    {
      do
        sq_monic_set_next (ell, ell, 64);
      while (!sq_is_irreducible (ell));
      fq_info_init (Fq, ell);
      fqpol_init (Fbar);
      fqpol_set_ffspol (Fbar, ffspol, Fq);
      if (Fbar->deg == d && fqpol_is_irreducible (Fbar, Fq))
        ret = 1;
      fqpol_clear (Fbar);
      fq_info_clear (Fq);
    }
  ffspol_clear (ffspol);
  return ret;
}
#endif

static void
all1 (unsigned long *f, unsigned long n1, unsigned long n0)
{
  unsigned long i, j, *l, d, D, r, f0, fj, K, N1 = 1 << n1, N0 = 1 << n0;
  unsigned long hl, k, m, p;
  float s0, s1, sj;
  int is_power;
  /* irreducible polynomials and powers of increasing index */
  unsigned long L[] = {2, 3, 4, 5, 7, 8, 11, 13, 15, 16, 17, 19, 21, 25, 31,
                       32, 37, 41, 47, 51, 55, 59, 61, 64, 67, 69, 73, 81, 85,
                       87, 91, 97, 103, 107, 109, 115, 117, 128, 131, 137, 143,
                       145, 157, 167, 171, 185, 191, 193, 203, 211, 213, 229,
                       239, 241, 247, 253, 255, 256, 257, 261, 273, 283, 285,
                       299, 301, 313, 319, 321, 333, 341, 351, 355, 357, 361,
                       369, 375, 379, 391, 395, 397, 415, 419, 425, 433, 445,
                       451, 463, 471, 477, 487, 499, 501, 505, 512, 515, 529,
                       535, 539, 545, 557, 563, 587, 601, 607, 613, 617, 623,
                       631, 637, 647, 661, 665, 675, 677, 687, 695, 701, 719,
                       721, 731, 743, 757, 761, 769, 771, 787, 789, 799, 803,
                       817, 827, 841, 847, 859, 865, 875, 877, 883, 895, 901,
                       911, 925, 929, 949, 953, 967, 971, 973, 981, 985, 995,
                       1001, 1019, 1024, 1033, 1039, 1041, 1051, 1053, 1063,
                       1069, 1077, 1089, 1095, 1107, 1109, 1123, 1125, 1135,
                       1153, 1163, 1177, 1193, 1199, 1221, 1225, 1239, 1255,
                       1261, 1267, 1279, 1285, 1291, 1293, 1301, 1305, 1311,
                       1315, 1329, 1341, 1347, 1349, 1361, 1367, 1377, 1383,
                       1387, 1413, 1423, 1431, 1435, 1441, 1451, 1465, 1473,
                       1479, 1509, 1527, 1531, 1555, 1557, 1571, 1573, 1585,
                       1591, 1603, 1615, 1617, 1627, 1657, 1663, 1669, 1673,
                       1703, 1709, 1717, 1727, 1729, 1741, 1747, 1759, 1783,
                       1789, 1807, 1809, 1815, 1821, 1825, 1835, 1845, 1849,
                       1863, 1869, 1877, 1881, 1891, 1911, 1915, 1917, 1921,
                       1927, 1933, 1939, 1961, 1969, 1989, 2011, 2027, 2035,
                       2041, 2047, ULONG_MAX};
  /* powers of irreducible polynomials, with corresponding power */
  unsigned long P[][2] = {{4, 2}, {5, 2}, {8, 3}, {15, 3}, {16, 4}, {17, 4},
                          {21, 2}, {32, 5}, {51, 5}, {64, 6}, {69, 2}, {81, 2},
                          {85, 6}, {107, 3}, {128, 7}, {255, 7}, {256, 8},
                          {257, 8}, {261, 2}, {273, 4}, {321, 2}, {341, 2},
                          {512, 9}, {743, 3}, {771, 9}, {925, 3}, {1024, 10},
                          {1041, 2}, {1089, 2}, {1109, 2}, {1285, 10},
                          {1301, 2}, {1349, 2}, {1361, 2}, {1911, 5},
                          {ULONG_MAX, 0}};

  f[0] = f[1] = 0;
  for (j = 0; j < N1; j++)
    for (i = 0; i < N0; i++)
      s[j][i] = 0.0;
  for (l = L, p = 0; *l < B; l++)
    {
      unsigned long d0, D0;
      while (P[p][0] < *l)
        p++;
      if (P[p][0] == *l)
        is_power = P[p][1];
      else
        is_power = 1;
      d = nbits (*l) - 1; /* degree of l */
      D = 1UL << d;
      d0 = d / is_power; /* degree of irreducible */
      D0 = 1UL << d0;
      if (is_power == 1 && (B < D * D)) /* only power is 1 */
        {
          s0 = (float) d / (float) (D - 1);
          s1 = (float) d * (float) D / (float) (D * D - 1);
        }
      else if ((1UL << (d + d0)) <= B) /* not the largest power */
        {
          if (is_power == 1)
            s0 = (float) d / (float) (D - 1);
          else
            s0 = (float) 0.0;
          s1 = (float) d0 / (float) (D0 + 1) / (float) (1UL << (d - d0));
        }
      else
        {
          assert (is_power > 1);
          s0 = (float) 0.0;
          s1 = (float) d0 * (float) D0 / (float) (D0 * D0 - 1)
            / (float) (1UL << (d - d0));
        }
      for (j = 0; j < N1; j++)
        for (i = 0; i < D; i++)
          sl[j][i] = s0;
      init_time -= cputime ();
      for (r = 0; r < D; r++)
        {
          f0 = evalmod (f, DEGREE, r, *l);
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
      init_time += cputime ();
      K = 1 << (n0 - d);
      sieve_time -= cputime ();
//if (D < 4) {
if (1) {
      for (j = 0; j < N1; j++)
        for (r = 0, hl = 0; r < D; r++)
          {
            sj = sl[j][r];
            for (k = 0, hl = 0; k < K; k++)
              {
                m = hl ^ r;
                assert (m < N0);
                s[j][m] += sj;
                hl ^= *l << G[k];
              }
          }
} else { // Experimental version, which seems to be faster, but
         // further tests are needed before putting it in prod.
      for (j = 0; j < N1; j++)
        for (unsigned long r0 = 0, hl = 0; r0 < D; r0+=4)
          {
            unsigned long r1, r2, r3;
            r1=r0+1;
            r2=r0+2;
            r3=r0+3;
            float sj0,sj1,sj2,sj3;
            sj0 = sl[j][r0];
            sj1 = sl[j][r1];
            sj2 = sl[j][r2];
            sj3 = sl[j][r3];
            for (k = 0, hl = 0; k < K; k++)
              {
                m = hl ^ r0;
                assert (m < N0);
                s[j][m] += sj0;
                m = hl ^ r1;
                assert (m < N0);
                s[j][m] += sj1;
                m = hl ^ r2;
                assert (m < N0);
                s[j][m] += sj2;
                m = hl ^ r3;
                assert (m < N0);
                s[j][m] += sj3;
                hl ^= *l << G[k];
              }
          }
}
      sieve_time += cputime ();
    }

  for (j = 0; j < N1; j++)
    /* we start from i=1 since for i=0, we have f0=0, thus f is divisible
       by y */
    for (i = 1; i < N0; i++)
      if (s[j][i] < threshold)
        {
          f[0] = i;
          f[1] = j;
          if (is_maybe_irreducible (f, DEGREE) &&
              is_irreducible (f, DEGREE))
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
  unsigned long j;
  unsigned long n0 = 0, n1 = 0, n2 = 0, n3 = 0, n4 = 0, n5 = 0,
    N0, N1, N2, N3, N4, N5, f5min = 0, f4min = 0, f3min = 0, f2min = 0;

  /* print command line */
  fprintf (stderr, "%s", argv[0]);
  for (j = 1; j < argc; j++)
    fprintf (stderr, " %s", argv[j]);
  fprintf (stderr, "\n");

  while (argc > 1 && argv[1][0] == '-')
    {
      if (argc > 2 && strcmp (argv[1], "-threshold") == 0)
        {
          threshold = atof (argv[2]);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-B") == 0)
        {
          B = atoi (argv[2]);
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
      else if (argc > 2 && strcmp (argv[1], "-n3") == 0)
        {
          n3 = atoi (argv[2]);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-n4") == 0)
        {
          n4 = atoi (argv[2]);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-n5") == 0)
        {
          n5 = atoi (argv[2]);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-f5min") == 0)
        {
          f5min = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-f4min") == 0)
        {
          f4min = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-f3min") == 0)
        {
          f3min = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-f2min") == 0)
        {
          f2min = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-f4") == 0)
        {
          f[4] = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-f3") == 0)
        {
          f[3] = strtoul (argv[2], NULL, 10);
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
  N3 = 1 << n3;
  N4 = 1 << n4;
  N5 = 1 << n5;

  if (threshold == -DBL_MAX)
    {
      fprintf (stderr, "Missing -threshold option\n");
      exit (1);
    }

  if (B == 0)
    {
      fprintf (stderr, "Missing -B option\n");
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

  for (f[5] = f5min; f[5] < f5min + N5; f[5]++)
    for (f[4] = f4min; f[4] < f4min + N4; f[4]++)
      for (f[3] = f3min; f[3] < f3min + N3; f[3]++)
        for (f[2] = f2min; f[2] < f2min + N2; f[2]++)
          all1 (f, n1, n0);
  fprintf (stderr, "Found %lu polynomial(s) below threshold\n", nsol);
  fprintf (stderr, "Init time %dms, sieve time %dms\n", init_time, sieve_time);

  free (G);
  for (j = 0; j < N1; j++)
    free (sl[j]);
  free (sl);
  for (j = 0; j < N1; j++)
    free (s[j]);
  free (s);

  return 0;
}
