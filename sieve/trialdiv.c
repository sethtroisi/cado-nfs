#include <stdio.h>
#include <stdlib.h>
#include "ularith.h"
#include "modredc_ul.h"
#include "trialdiv.h"

static const int verbose = 0;

static void
trialdiv_init_divisor (trialdiv_divisor_t *d, const unsigned long p)
{
  /* TODO: test if p is too large */
  ASSERT (p % 2UL == 1UL);
  d->p = p;
  if (p == 1UL)
    d->w[0] = 0UL;
  else
    ularith_div_2ul_ul_ul_r (&(d->w[0]), 0UL, 1UL, p);
  ularith_div_2ul_ul_ul_r (&(d->w[1]), 0UL, d->w[0], p);
  ularith_div_2ul_ul_ul_r (&(d->w[2]), 0UL, d->w[1], p);
  ularith_div_2ul_ul_ul_r (&(d->w[3]), 0UL, d->w[2], p);
  ularith_div_2ul_ul_ul_r (&(d->w[4]), 0UL, d->w[3], p);
  d->pinv = modredcul_invmodul (p);
  d->plim = (unsigned long) (-1L) / p;
}

/* Trial division for integers with 1 unsigned long */
static inline int
trialdiv_div1 (const unsigned long *n, const trialdiv_divisor_t *d)
{
  return n[0] * d->pinv <= d->plim;
}

/* Trial division for integers with 2 unsigned long */
static inline int
trialdiv_div2 (const unsigned long *n, const trialdiv_divisor_t *d)
{
  unsigned long r0, r1, x0, x1;
  __asm__ ("# trialdiv_div2");
  ularith_mul_ul_ul_2ul (&x0, &x1, n[1], d->w[0]); /* n_1 * (w^1 % p) */
  r0 = n[0];
  r1 = x1;
  ularith_add_ul_2ul (&r0, &r1, x0);
  x0 = r1 * d->w[0]; /* TODO: optimal? */
  x0 += r0;
  if (x0 < r0)
    x0 += d->w[0];
  x0 *= d->pinv;
  return x0 <= d->plim;
}

/* Trial division for integers with 3 unsigned long */
static inline int
trialdiv_div3 (const unsigned long *n, const trialdiv_divisor_t *d)
{
  unsigned long r0, r1, x0, x1;
  if (verbose)
    printf ("p = %lu, n2:n1:n0 = %lu * w^2 + %lu * w + %lu", 
            d->p, n[2], n[1], n[0]);
  ularith_mul_ul_ul_2ul (&x0, &x1, n[1], d->w[0]); /* n_1 * (w^1 % p) */
  r0 = n[0];
  r1 = x1;
  ularith_add_ul_2ul (&r0, &r1, x0);
  ularith_mul_ul_ul_2ul (&x0, &x1, n[2], d->w[1]); /* n_2 * (w^2 % p) */
  ularith_add_2ul_2ul (&r0, &r1, x0, x1);
  if (verbose)
    printf ("  r1:r0 = %lu * w + %lu", r1, r0);
  x0 = r1 * d->w[0];
  x0 += r0;
  if (x0 < r0)
    x0 += d->w[0];
  if (verbose)
    printf ("  x0 = %lu", x0);
  x0 *= d->pinv;
  return x0 <= d->plim;
}

/* Trial division for integers with 4 unsigned long */
static inline int
trialdiv_div4 (const unsigned long *n, const trialdiv_divisor_t *d)
{
  unsigned long r0, r1, x0, x1;
  ularith_mul_ul_ul_2ul (&x0, &x1, n[1], d->w[0]); /* n_1 * (w^1 % p) */
  r0 = n[0];
  r1 = x1;
  ularith_add_ul_2ul (&r0, &r1, x0);
  ularith_mul_ul_ul_2ul (&x0, &x1, n[2], d->w[1]); /* n_2 * (w^2 % p) */
  ularith_add_2ul_2ul (&r0, &r1, x0, x1);
  ularith_mul_ul_ul_2ul (&x0, &x1, n[3], d->w[2]); /* n_3 * (w^3 % p) */
  ularith_add_2ul_2ul (&r0, &r1, x0, x1);
  x0 = r1 * d->w[0];
  x0 += r0;
  if (x0 < r0)
    x0 += d->w[0];
  x0 *= d->pinv;
  return x0 <= d->plim;
}

/* Trial division for integers with 5 unsigned long */
static inline int
trialdiv_div5 (const unsigned long *n, const trialdiv_divisor_t *d)
{
  unsigned long r0, r1, x0, x1;
  ularith_mul_ul_ul_2ul (&x0, &x1, n[1], d->w[0]); /* n_1 * (w^1 % p) */
  r0 = n[0];
  r1 = x1;
  ularith_add_ul_2ul (&r0, &r1, x0);
  ularith_mul_ul_ul_2ul (&x0, &x1, n[2], d->w[1]); /* n_2 * (w^2 % p) */
  ularith_add_2ul_2ul (&r0, &r1, x0, x1);
  ularith_mul_ul_ul_2ul (&x0, &x1, n[3], d->w[2]); /* n_3 * (w^3 % p) */
  ularith_add_2ul_2ul (&r0, &r1, x0, x1);
  ularith_mul_ul_ul_2ul (&x0, &x1, n[4], d->w[3]); /* n_4 * (w^4 % p) */
  ularith_add_2ul_2ul (&r0, &r1, x0, x1);
  x0 = r1 * d->w[0];
  x0 += r0;
  if (x0 < r0)
    x0 += d->w[0];
  x0 *= d->pinv;
  return x0 <= d->plim;
}

/* Trial division for integers with 6 unsigned long */
static inline int
trialdiv_div6 (const unsigned long *n, const trialdiv_divisor_t *d)
{
  unsigned long r0, r1, x0, x1;
  ularith_mul_ul_ul_2ul (&x0, &x1, n[1], d->w[0]); /* n_1 * (w^1 % p) */
  r0 = n[0];
  r1 = x1;
  ularith_add_ul_2ul (&r0, &r1, x0);
  ularith_mul_ul_ul_2ul (&x0, &x1, n[2], d->w[1]); /* n_2 * (w^2 % p) */
  ularith_add_2ul_2ul (&r0, &r1, x0, x1);
  ularith_mul_ul_ul_2ul (&x0, &x1, n[3], d->w[2]); /* n_3 * (w^3 % p) */
  ularith_add_2ul_2ul (&r0, &r1, x0, x1);
  ularith_mul_ul_ul_2ul (&x0, &x1, n[4], d->w[3]); /* n_4 * (w^4 % p) */
  ularith_add_2ul_2ul (&r0, &r1, x0, x1);
  ularith_mul_ul_ul_2ul (&x0, &x1, n[5], d->w[4]); /* n_5 * (w^5 % p) */
  ularith_add_2ul_2ul (&r0, &r1, x0, x1);
  x0 = r1 * d->w[0];
  x0 += r0;
  if (x0 < r0)
    x0 += d->w[0];
  x0 *= d->pinv;
  return x0 <= d->plim;
}


/* Divides primes in d out of N and stores them (with multiplicity) in f */

int 
trialdiv (unsigned long *f, mpz_t N, const trialdiv_divisor_t *d)
{
  int n = 0;
  
  while (mpz_cmp_ui (N, 1UL) > 0)
    {
      size_t s = mpz_size(N);
      if (verbose) gmp_printf ("s = %d, N = %Zd, ", s, N);
      if (s > 6)
        abort ();
      if (s == 1)
        {
          while (!trialdiv_div1 (N[0]._mp_d, d))
            d++;
        }
      else if (s == 2)
        {
          while (!trialdiv_div2 (N[0]._mp_d, d))
            d++;
        }
      else if (s == 3)
        {
          while (!trialdiv_div3 (N[0]._mp_d, d))
            d++;
        }
      else if (s == 4)
        {
          while (!trialdiv_div4 (N[0]._mp_d, d))
            d++;
        }
      else if (s == 5)
        {
          while (!trialdiv_div5 (N[0]._mp_d, d))
            d++;
        }
      else if (s == 6)
        {
          while (!trialdiv_div6 (N[0]._mp_d, d))
            d++;
        }
      if (verbose) printf ("\n");

      if (d->p == 1UL)
        break;

      ASSERT (mpz_divisible_ui_p (N, d->p));

      do {
       f[n++] = d->p;
        mpz_divexact_ui (N, N, d->p);
      } while (mpz_divisible_ui_p (N, d->p));
      d++;
    }
  
  return n;
}


/* Initialise a trialdiv_divisor_t array with the nr primes stored in *f.
   This function allcates memory for the array, inits each entry, and puts
   a sentinel at the end. */
trialdiv_divisor_t *
trialdiv_init (const fbprime_t *f, const unsigned int nr)
{
  trialdiv_divisor_t *d;
  unsigned int i;
  
  d = (trialdiv_divisor_t *) malloc ((nr + 1) * sizeof (trialdiv_divisor_t));
  ASSERT (d != NULL);
  for (i = 0; i < nr; i++)
    trialdiv_init_divisor (&(d[i]), (unsigned long) (f[i]));
  trialdiv_init_divisor (&(d[nr]), 1UL);

  return d;
}

void
trialdiv_clear (trialdiv_divisor_t *d)
{
  free (d);
}

#if 0
int main (int argc, char **argv)
{
  trialdiv_divisor_t d[2];
  unsigned long factors[16];
  int i, r = 0, len = 1;
  mpz_t M, N;
  
  if (argc > 1)
    len = atoi (argv[1]);
  
  mpz_init (M);
  mpz_init (N);
  mpz_set_ui (N, 1UL);
  mpz_mul_2exp (N, N, 8 * sizeof(unsigned long) * len);
  mpz_tdiv_q_ui (N, N, 3UL);
    
  trialdiv_init_divisor (&(d[0]), 101UL);
  trialdiv_init_divisor (&(d[1]), 1UL);

  for (i=0; i < 100000000; i++)
    {
      mpz_set (M, N);
      r += trialdiv (factors, M, d);
      mpz_add_ui (N, N, 1UL);
    }
 
  mpz_clear (M);
  mpz_clear (N);
  printf ("%d\n", r);
  return 0;
}
#endif
