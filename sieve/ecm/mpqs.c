/* tiny MPQS implementation, specially tuned for 128-bit input */

// #define TRACE 446179
// #define TRACE_P 937

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/resource.h>
#include "gmp.h"
#include "macros.h"

#define DIM 256

/* For i < 50, isprime_table[i] == 1 iff i is prime */
static unsigned char isprime_table[] = {
0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 
0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 
0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0};

static const size_t isprime_table_size = 
    sizeof(isprime_table) / sizeof(isprime_table[0]);

long
cputime ()
{
  struct rusage rus;

  getrusage (RUSAGE_SELF, &rus);
  /* This overflows a 32 bit signed int after 2147483s = 24.85 days */
  return rus.ru_utime.tv_sec * 1000L + rus.ru_utime.tv_usec / 1000L;
}

int
jacobi (unsigned long a, unsigned long b)
{
  mpz_t aa, bb;
  int res;

  mpz_init_set_ui (aa, a);
  mpz_init_set_ui (bb, b);
  res = mpz_jacobi (aa, bb);
  mpz_clear (aa);
  mpz_clear (bb);
  return res;
}

/* return b^e mod n, assuming n has at most 32 bits */
static uint64_t
mod_pow_uint64 (uint64_t b, uint64_t e, uint64_t n)
{
  uint64_t r = 1, f = 1, y = b;

  while (f <= e)
    {
      if (e & f)
        r = (r * y) % n;
      f <<= 1;
      y = (y * y) % n;
    }
  
  return r;
}

/* Uses Tonelli-Shanks, more precisely Algorithm from Table 1 in [1].
   [1] Adleman-Manders-Miller Root Extraction Method Revisited, Zhengjun Cao,
   Qian Sha, Xiao Fan, 2011, http://arxiv.org/abs/1111.4877.
   Solve x^2 = rr (mod p).
*/
static uint64_t
tonelli_shanks (uint64_t rr, uint64_t p)
{
  uint64_t q, s, i, j, l;
  uint64_t aa, hh, delta, dd, bb, zz;

  if (p == 2)
    return rr;

  assert (p <= 4294967295UL);

  /* write p-1 = q*2^s with q odd */
  for (q = p-1, s = 0; (q&1) == 0; q/=2, s++);
  
  zz = 1;
  i = 1;
  do {
    /* zz is equal to i (mod p) */
    zz += 1;
    i++;
  } while (i < isprime_table_size && (!isprime_table[i] || jacobi (zz, p) != -1));
  assert (i < isprime_table_size);
  aa = mod_pow_uint64 (zz, q, p);
  bb = mod_pow_uint64 (rr, q, p);
  hh = 1;
  for (j = 1; j < s; j++)
    {
      dd = bb;
      for (l = 0; l < s - 1 - j; l++)
        dd = (dd * dd) % p;
      if (dd != 1)
        {
          hh = (hh * aa) % p;
          aa = (aa * aa) % p;
          bb = (bb * aa) % p;
        }
      else
        aa = (aa * aa) % p;
    }
  delta = mod_pow_uint64 (rr, (q + 1) >> 1, p);
  hh = (hh * delta) % p;

  return hh;
}

/* Given k1 such that k1^2 = N (mod p),
   return k1 and k2 which are the roots of (a*x+b)^2 = N (mod p).
   Assume a is odd.
*/
unsigned long
findroot (unsigned long *k2, mpz_t a, mpz_t b, unsigned long p,
          unsigned long k1)
{
  /* special case for p=2: since a is odd, x = k1-b (mod 2) */
  if (p == 2)
    return (k1 - mpz_tstbit (b, 0)) & 1;

  /* the two roots are (k1-b)/a and (-k1-b)/a */
  unsigned long bmodp = mpz_fdiv_ui (b, p), inva;
  mpz_t t;
  mpz_init_set_ui (t, p);
  mpz_invert (t, a, t);
  inva = mpz_get_ui (t);
  mpz_clear (t);
  *k2 = (k1 == 0) ? 0 : p - k1;
  k1 = (k1 >= bmodp) ? k1 - bmodp : k1 + p - bmodp;
  k1 = (k1 * inva) % p;
  *k2 = (*k2 >= bmodp) ? *k2 - bmodp : *k2 + p - bmodp;
  *k2 = (*k2 * inva) % p;
  if (*k2 < k1)
    {
      unsigned long tmp = *k2;
      *k2 = k1;
      return tmp;
    }
  return k1;
}

/* check that ((a*i+b)^2-N)/a is smooth */
int
is_smooth (mpz_t a, unsigned long i, mpz_t b, mpz_t N, unsigned long *P,
           unsigned long ncol, mpz_t row)
{
  mpz_t c;
  int res = 0;

  mpz_set_ui (row, 0);
  mpz_init (c);
  mpz_mul_si (c, a, i);
  mpz_add (c, c, b);
  mpz_mul (c, c, c);
  mpz_sub (c, c, N);
  assert (mpz_divisible_p (c, a));
  mpz_divexact (c, c, a);
  if (mpz_sgn (c) < 0)
    {
      mpz_setbit (row, 0); /* sign is in column 0 */
      mpz_neg (c, c);
    }
  if (mpz_cmp_ui (c, P[ncol-1]) > 0 && mpz_probab_prime_p (c, 1))
    {
      res = 0;
      goto end;
    }
  for (i = 0; i < ncol; i++)
    {
      if (mpz_divisible_ui_p (c, P[i]))
        {
          int e = 0;
          do {
            mpz_divexact_ui (c, c, P[i]);
            e ++;
          } while (mpz_divisible_ui_p (c, P[i]));
          if (e & 1)
            mpz_setbit (row, i + 1);
          if (mpz_cmp_ui (c, 1) == 0)
            {
              res = 1;
              goto end;
            }
          if (mpz_cmp_ui (c, P[ncol-1]) > 0 && mpz_probab_prime_p (c, 1))
            {
              res = 0;
              goto end;
            }
        }
    }
 end:
  mpz_clear (c);
  return res;
}

static inline void
update (unsigned char *S, unsigned long i, unsigned long p MAYBE_UNUSED,
        unsigned char logp, unsigned int M MAYBE_UNUSED)
{
  S[i] += logp;
#ifdef TRACE
  if (i == M + TRACE)
    printf ("%d: %lu %u\n", TRACE, p, S[i]);
#endif
}

/* put factor in z */
void
gauss (mpz_t z, mpz_t *Mat, int nrel, int ncol, mpz_t *X, mpz_t N)
{
  int i, j, k;

  for (j = 0, i = 0; j < ncol; j++)
    {
      /* try to zero all bits in column j, starting from row i */
      for (k = i; k < nrel; k++)
        if (mpz_tstbit (Mat[k], j))
          break;
      if (k == nrel)
        continue; /* go to next column */
      if (i < k) /* swap rows i and k */
        mpz_swap (Mat[i], Mat[k]);
      /* now we have a pivot in (i,j) */
      for (k = i + 1; k < nrel; k++)
        if (mpz_tstbit (Mat[k], j))
          mpz_xor (Mat[k], Mat[k], Mat[i]);
      i++;
    }
  /* all rows from i to nrel-1 are dependencies */
  mpz_t x, y;
  mpz_init (x);
  mpz_init (y);
  while (i < nrel)
    {
      printf ("Trying dependency %d\n", i);
      mpz_set_ui (x, 1);
      mpz_set_ui (y, 1);
      for (j = 0; j < nrel; j++)
        if (mpz_tstbit (Mat[i], ncol + j))
          {
            mpz_mul (x, x, X[j]);
            mpz_mul (z, X[j], X[j]);
            mpz_sub (z, z, N);
            mpz_mul (y, y, z);
          }
      assert (mpz_perfect_square_p (y));
      mpz_sqrt (y, y);
      mpz_sub (z, x, y);
      mpz_gcd (z, z, N);
      if (mpz_cmp_ui (z, 1) > 0 && mpz_cmp (z, N) < 0)
        {
          gmp_printf ("gcd=%Zd\n", z);
          break;
        }
      i++;
    }
  mpz_clear (x);
  mpz_clear (y);
}

/* Put in f a factor of N, using factor base of ncol primes.
   Assume N is odd. */
void
mpqs (mpz_t f, mpz_t N, long ncol)
{
  mpz_t *Mat, *X;
  unsigned char *S, *Logp, logp;
  unsigned long *P, p, q, k;
  long M, i, j, *K, wrel, nrel = 0;
  mpz_t R, a, b, c, sqrta;
  double maxnorm, radix, logradix;
  long st0 = cputime (), st, init_time, sieve_time = 0, check_time = 0;
  long gauss_time, total_time = 0;

  /* assume N is odd */
  assert (mpz_fdiv_ui (N, 2) == 1);

  mpz_init (a);
  mpz_init (sqrta);
  mpz_init (b);
  mpz_init (c);

  /* compute factor base */
  P = malloc (ncol * sizeof (unsigned long));
  K = malloc (ncol * sizeof (unsigned long));
  mpz_set_ui (a, 1);
  /* a prime p can appear in x^2 mod N only if p is a square modulo N */
  for (j = 0; j < ncol;)
    {
      mpz_nextprime (a, a);
      if (mpz_jacobi (a, N) == 1)
        P[j++] = mpz_get_ui (a);
    }
  printf ("largest prime is %lu\n", P[ncol - 1]);
  wrel = ncol + 10; /* wanted number of relations */
  Mat = (mpz_t*) malloc (wrel * sizeof (mpz_t));
  for (i = 0; i < wrel; i++)
    mpz_init (Mat[i]);
  X = (mpz_t*) malloc (wrel * sizeof (mpz_t));
  for (i = 0; i < wrel; i++)
    mpz_init (X[i]);

  M = (1<<16); /* (half) size of sieve array */

#define GUARD 3
#define MAXS (255-GUARD)

  /* initialize sieve area [-M, M-1] */
  S = malloc (2 * M * sizeof (char));
  Logp = malloc (ncol * sizeof (char));

  /* initialize square roots mod p */
  for (j = 0; j < ncol; j++)
    {
      unsigned long Np;
      p = P[j];
      Np = mpz_fdiv_ui (N, p);
      K[j] = tonelli_shanks (Np, p);
    }

  /* we want 'a' near sqrt(2*N)/M */
  mpz_mul_ui (a, N, 2);
  mpz_sqrt (a, a);
  mpz_tdiv_q_ui (a, a, M);

  /* we want 'a' to be a square */
  mpz_sqrt (sqrta, a);

  st = cputime ();
  printf ("init: %ldms\n", st);
  init_time = st - st0;

  while (nrel < wrel) {
  sieve_time -= st;
  do
    mpz_nextprime (sqrta, sqrta);
  while (mpz_jacobi (sqrta, N) != 1);
  unsigned long aui = mpz_get_ui (sqrta);
  mpz_mul (a, sqrta, sqrta);
  /* we want b^2-n divisible by a */
  k = tonelli_shanks (mpz_fdiv_ui (N, aui), aui);
  /* we want b = k + a*t: b^2 = k^2 + 2*k*a*t (mod a^2), thus
     k^2 + 2*k*a*t = N (mod a^2): t = ((N-k^2)/a)/(2*k) mod a */
  mpz_ui_pow_ui (b, k, 2);
  mpz_sub (b, N, b);
  assert (mpz_divisible_p (b, sqrta));
  mpz_divexact (b, b, sqrta);
  mpz_set_ui (c, k);
  mpz_mul_ui (c, c, 2);
  mpz_invert (c, c, sqrta);
  mpz_mul (b, b, c);
  mpz_mod (b, b, sqrta);
  mpz_mul (b, b, sqrta);
  mpz_add_ui (b, b, k);
  // gmp_printf ("a=%Zd\nb=%Zd\n", a, b);

  /* we now have (a*x+b)^2 - N = a*Q(x) where Q(x) = a*x^2 + 2*b*x + c */

  /* initialize radix and Logp[]: we want radix^MAXS = maxnorm ~ N/a */
  maxnorm = mpz_get_d (N) / mpz_get_d (a);
  radix = pow (maxnorm, 1.0 / (double) MAXS);
  logradix = log (radix);
  for (j = 0; j < ncol; j++)
    Logp[j] = (char) (log ((double) P[j]) / logradix + 0.5);

  memset (S, 0, 2 * M * sizeof (char));

  /* sieve */
  mpz_submul_ui (b, a, M);
  for (j = 0; j < ncol; j++)
    {
      unsigned long Np, k2;
      long i2;
      p = P[j];
      k = K[j]; /* k^2 = N (mod p) */
#ifdef TRACE_P
      if (p == TRACE_P)
        printf ("p=%lu k=%lu\n", p, k);
#endif
      k = findroot (&k2, a, b, p, k);
#ifdef TRACE_P
      if (p == TRACE_P)
        printf ("p=%lu k=%lu k2=%lu\n", p, k, k2);
#endif
      logp = Logp[j];
      if (p == 2)
        {
          /* Sieve 2, 4 and 8 simultaneously:
             if N = 1 (mod 8) then 8 divides x^2-N for x = 1, 3, 5, 7
             if N = 3 (mod 8) then 2 divides x^2-N for x = 1, 3, 5, 7
             if N = 5 (mod 8) then 4 divides x^2-N for x = 1, 3, 5, 7
             if N = 7 (mod 8) then 2 divides x^2-N for x = 1, 3, 5, 7.
          */
          Np = mpz_getlimbn (N, 0) & 7; /* N mod 8 */
          if (Np == 1) /* 2^3 always divides */
            logp = (char) (log ((double) 8) / logradix + 0.5);
          else if (Np == 5) /* 2^2 always divides */
            logp = (char) (log ((double) 4) / logradix + 0.5);
          for (i = k; i < 2*M; i += 2)
            update (S, i, p, logp, M);
        }
      else
        {
          for (i = k, i2 = k2; i2 < 2*M; i += p, i2 += p)
            {
              update (S, i, p, logp, M);
              update (S, i2, p, logp, M);
            }
          if (i < 2*M)
            update (S, i, p, logp, M);

          /* sieve squares */
          q = p * p;
          if (q < P[ncol - 1])
            {
              unsigned long Ap, Bp;
              unsigned char logp2;
              logp2 = (char) (log ((double) q) / logradix + 0.5) - logp;
              Np = mpz_fdiv_ui (N, q);
              Ap = mpz_fdiv_ui (a, q);
              Bp = mpz_fdiv_ui (b, q);
              while (k < q)
                {
                  unsigned long y = (Ap * k + Bp) % q;
                  /* check if (a*k+b)^2 = N (mod q) */
                  if ((y * y) % q == Np)
                    for (i = k; i < 2*M; i += q)
                      update (S, i, p, logp2, M);
                  k += p;
                  if (k2 >= q)
                    break;
                  y = (Ap * k2 + Bp) % q;
                  if ((y * y) % q == Np)
                    for (i = k2; i < 2*M; i += q)
                      update (S, i, p, logp2, M);
                  k2 += p;
                }
            }
        }
    }

  st = cputime ();
  sieve_time += st;
  check_time -= st;

#ifdef TRACE
  printf ("%d: S=%d\n", TRACE, S[M + TRACE]);
#endif

  /* find smooth locations: the norm for location i is about 2*sqrt(N)*i */
  int count = 0, count2 = 0;
  double ad, bd, y;
  double Nd = mpz_get_d (N);
  /* we go from about sqrt(N) for i=0 to N/a for i=M */
  unsigned char threshold = (unsigned char) (log (Nd) / logradix);
  ad = mpz_get_d (a);
  bd = mpz_get_d (b);
  for (i = -M; i < M && nrel < wrel; i++)
    {
#ifdef TRACE
      if (i == TRACE)
        printf ("i=%ld: S=%u threshold=%u\n", i, S[M+i], threshold);
#endif
      if (S[M + i] >= threshold)
        {
          count2 ++;
          y = ad * (double) (M+i) + bd;
          y = fabs ((y * y - Nd) / ad);
          y = log (y) / logradix;
#ifdef TRACE
          if (i == TRACE)
            printf ("i=%ld: S=%u y=%f\n", i, S[M+i], y);
#endif
          if (S[M + i] >= (unsigned char) y - 15) {
            count ++;
            if (is_smooth (a, M + i, b, N, P, ncol, Mat[nrel]))
              {
                mpz_setbit (Mat[nrel], ncol + 1 + nrel);
                mpz_set (X[nrel], b);
                mpz_addmul_ui (X[nrel], a, M + i);
                nrel ++;
                // printf ("%ld: i=%ld\n", nrel, i);
              }
          }
        }
    }
  st = cputime ();
  check_time += st;
  gmp_printf ("sqrta=%Zd: %d/%d above threshold, total %ld relations in %ldms (%.2f r/s)\n",
          sqrta, count2, count, nrel, st, (double) nrel / (double) st);
  }

  gauss_time = st;
  gauss (f, Mat, nrel, ncol + 1, X, N);
  st = cputime ();
  gauss_time = st - gauss_time;
  total_time = st - st0;

  free (Logp);
  free (S);
  mpz_clear (R);
  mpz_clear (a);
  mpz_clear (sqrta);
  mpz_clear (b);
  mpz_clear (c);
  for (i = 0; i < wrel; i++)
    mpz_clear (Mat[i]);
  free (Mat);
  for (i = 0; i < wrel; i++)
    mpz_clear (X[i]);
  free (X);
  free (P);
  free (K);

  printf ("Total time: %ldms (init %ld, sieve %ld, check %ld, gauss %ld)\n",
          total_time, init_time, sieve_time, check_time, gauss_time);
}

/* On tarte.loria.fr:
$ ./mpqs 270788552349171139784543548689828248993 900
Total time: 144ms (init 20, sieve 64, check 48, gauss 12)
*/
int
main (int argc, char *argv[])
{
  mpz_t N, f;
  unsigned long ncol;

  assert (argc > 2);
  mpz_init (N);
  mpz_init (f);
  mpz_set_str (N, argv[1], 10);
  ncol = atol (argv[2]);
  mpqs (f, N, ncol);
  mpz_clear (N);
  mpz_clear (f);
}
