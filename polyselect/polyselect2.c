/* This version implements Kleinjung's original algorithm, with a search in a
   hash table, and M = P^2

   Reference:
   Polynomial selection, Thorsten Kleinjung, talk at the CADO workshop on
   integer factorization, Nancy, October 2008.
   http://cado.gforge.inria.fr/workshop/slides/kleinjung.pdf
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <pthread.h>
#include "utils.h"
#include "auxiliary.h"

// #define DEBUG

#define MAX_THREADS 16

/* hash table structure */
typedef struct
{
  uint32_t *p;              /* contains the primes */
  int64_t *i;               /* contains the values of r such that p^2
                               divides N - (m0 + r)^2 */
  unsigned long alloc;      /* total allocated size */
  unsigned long size;       /* number of entries in hash table */
} __hash_struct;
typedef __hash_struct hash_t[1];

/* thread structure */
typedef struct
{
  mpz_t N;
  unsigned int d;
  unsigned long ad;
  int thread;
} __tab_struct;
typedef __tab_struct tab_t[1];

/* structure to store roots of (m0+x)^d = N (mod p^2) for special-q variant */
typedef struct
{
  unsigned int alloc;   /* allocated size */
  unsigned int size;    /* used size */
  unsigned int *p;      /* primes */
  unsigned int *nr;     /* number of roots of x^d = N (mod p) */
  unsigned long **roots; /* roots of (m0+x)^d = N (mod p^2) */
} _roots_struct;
typedef _roots_struct roots_t[1];

/* read-only global variables */
int verbose = 0, incr = 60;
unsigned int nr = 0; /* minimum number of real roots wanted */
double max_norm = DBL_MAX; /* maximal wanted norm (before rotation) */
uint32_t *Primes;

/* read-write global variables */
pthread_mutex_t lock; /* used as mutual exclusion lock for those variables */
int found = 0;
double potential_collisions = 0.0;

void
roots_init (roots_t R)
{
  R->alloc = 0;
  R->size = 0;
  R->p = NULL;
  R->nr = NULL;
  R->roots = NULL;
}

void
roots_realloc (roots_t R, unsigned long newalloc)
{
  assert (newalloc >= R->size);
  R->alloc = newalloc;
  R->p = realloc (R->p, newalloc * sizeof (unsigned int));
  if (R->p == NULL)
    {
      fprintf (stderr, "Error, cannot reallocate memory in roots_realloc\n");
      exit (1);
    }
  R->nr = realloc (R->nr, newalloc * sizeof (unsigned int));
  if (R->nr == NULL)
    {
      fprintf (stderr, "Error, cannot reallocate memory in roots_realloc\n");
      exit (1);
    }
  R->roots = realloc (R->roots, newalloc * sizeof (unsigned long*));
  if (R->roots == NULL)
    {
      fprintf (stderr, "Error, cannot reallocate memory in roots_realloc\n");
      exit (1);
    }
}

void
roots_add (roots_t R, unsigned int p, unsigned int nr, unsigned long *roots)
{
  unsigned int i;

  if (nr == 0)
    return;
  if (R->size == R->alloc)
    roots_realloc (R, R->alloc + R->alloc / 2 + 1);
  R->p[R->size] = p;
  R->nr[R->size] = nr;
  R->roots[R->size] = malloc (nr * sizeof (unsigned long));
  if (R->roots[R->size] == NULL)
    {
      fprintf (stderr, "Error, cannot allocate memory in roots_add\n");
      exit (1);
    }
  for (i = 0; i < nr; i++)
    R->roots[R->size][i] = roots[i];
  R->size ++;
}

void
roots_clear (roots_t R)
{
  unsigned int i;

  free (R->p);
  free (R->nr);
  for (i = 0; i < R->size; i++)
    free (R->roots[i]);
  free (R->roots);
}

/* rq is a root of N = (m0 + rq)^d mod (q^2) */
static void
match (unsigned long p1, unsigned long p2, int64_t i, mpz_t m0,
       unsigned long ad, unsigned int d, mpz_t N, unsigned long q,
       unsigned long rq)
{
  mpz_t l, mtilde, m, adm1, t, k, *f, g[2];
  unsigned int j, nroots;
  int cmp;
  double skew, logmu, alpha;

#ifdef DEBUG
  printf ("Found match: (%lu,%ld) (%lu,%ld) for ad=%lu, q=%lu, rq=%lu\n",
	  p1, i, p2, i, ad, q, rq);
  gmp_printf ("m0=%Zd\n", m0);
#endif

  mpz_init (l);
  mpz_init (m);
  mpz_init (t);
  mpz_init (k);
  mpz_init (adm1);
  mpz_init (mtilde);
  mpz_init (g[0]);
  mpz_init (g[1]);
  f = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
  if (f == NULL)
    {
      fprintf (stderr, "Error, cannot allocate memory in match\n");
      exit (1);
    }
  for (j = 0; j <= d; j++)
    mpz_init (f[j]);
  /* we have l = p1*p2*q */
  mpz_set_ui (l, p1);
  mpz_mul_ui (l, l, p2);
  mpz_mul_ui (l, l, q);
  /* mtilde = m0 + rq + i*q^2 */
  mpz_set_si (mtilde, i);
  mpz_mul_ui (mtilde, mtilde, q * q);
  mpz_add_ui (mtilde, mtilde, rq);
  mpz_add (mtilde, mtilde, m0);
  /* we want mtilde = d*ad*m + a_{d-1}*l with 0 <= a_{d-1} < d*ad.
     We have a_{d-1} = mtilde/l mod (d*ad). */
  mpz_set_ui (m, ad);
  mpz_mul_ui (m, m, d);
  if (mpz_invert (adm1, l, m) == 0)
    {
      fprintf (stderr, "Error in 1/l mod (d*ad)\n");
      exit (1);
    }
  mpz_mul (adm1, adm1, mtilde);
  mpz_mod (adm1, adm1, m);
  mpz_mul (m, adm1, l);
  mpz_sub (m, mtilde, m);
#ifdef DEBUG
  if (mpz_divisible_ui_p (m, d) == 0)
    {
      fprintf (stderr, "Error: m-a_{d-1}*l not divisible by d\n");
      exit (1);
    }
#endif
  mpz_divexact_ui (m, m, d);
#ifdef DEBUG
  if (mpz_divisible_ui_p (m, ad) == 0)
    {
      fprintf (stderr, "Error: (m-a_{d-1}*l)/d not divisible by ad\n");
      exit (1);
    }
#endif
  mpz_divexact_ui (m, m, ad);
#ifdef DEBUG
  printf ("Raw polynomial:\n");
  gmp_printf ("Y1: %Zd\nY0: -%Zd\n", l, m);
#endif
  mpz_set (g[1], l);
  mpz_neg (g[0], m);
  mpz_set_ui (f[d], ad);
  mpz_pow_ui (t, m, d);
  mpz_mul_ui (t, t, ad);
  mpz_sub (t, N, t);
  mpz_set (f[d-1], adm1);
#ifdef DEBUG
  if (mpz_divisible_p (t, l) == 0)
    {
      fprintf (stderr, "Error: t not divisible by l\n");
      exit (1);
    }
#endif
  mpz_divexact (t, t, l);
  mpz_pow_ui (mtilde, m, d-1);
  mpz_mul (mtilde, mtilde, adm1);
  mpz_sub (t, t, mtilde);
  for (j = d - 2; j > 0; j--)
    {
#ifdef DEBUG
      if (mpz_divisible_p (t, l) == 0)
        {
          fprintf (stderr, "Error: t not divisible by l\n");
          exit (1);
        }
#endif
      mpz_divexact (t, t, l);
      /* t = a_j*m^j + l*R thus a_j = t/m^j mod l */
      mpz_pow_ui (mtilde, m, j);
      mpz_fdiv_q (adm1, t, mtilde); /* t -> adm1 * mtilde + t */
      mpz_invert (k, mtilde, l); /* search adm1 + k such that
                                    t = (adm1 + k) * m^j mod l */
      mpz_mul (k, k, t);
      mpz_sub (k, k, adm1);
      mpz_mod (k, k, l);

      mpz_mul_2exp (k, k, 1);
      cmp = mpz_cmp (k, l);
      mpz_div_2exp (k, k, 1);
      if (cmp >= 0)
        mpz_sub (k, k, l);
      mpz_add (adm1, adm1, k);
      mpz_set (f[j], adm1);
      /* subtract adm1*m^j */
      mpz_submul (t, mtilde, adm1);
    }
#ifdef DEBUG
  if (mpz_divisible_p (t, l) == 0)
    {
      fprintf (stderr, "Error: t not divisible by l\n");
      exit (1);
    }
#endif
  mpz_divexact (t, t, l);
  mpz_set (f[0], t);
#ifdef DEBUG
  for (j = d + 1; j -- != 0; )
    gmp_printf ("c%u: %Zd\n", j, f[j]);
  nroots = numberOfRealRoots (f, d, 0, 0);
  skew = SKEWNESS (f, d, SKEWNESS_DEFAULT_PREC);
  logmu = LOGNORM (f, d, skew);
  alpha = get_alpha (f, d, ALPHA_BOUND);
  printf ("# lognorm %1.2f, alpha %1.2f, E %1.2f, %u rroots\n",
          logmu, alpha, logmu + alpha, nroots);
  gmp_printf ("Optimized polynomial:\n");
#endif

  optimize (f, d, g, 0);
  nroots = numberOfRealRoots (f, d, 0, 0);
  skew = SKEWNESS (f, d, SKEWNESS_DEFAULT_PREC);
  logmu = LOGNORM (f, d, skew);
  if (nroots >= nr && logmu <= max_norm)
    {
      alpha = get_alpha (f, d, ALPHA_BOUND);
      gmp_printf ("n: %Zd\n", N);
      gmp_printf ("Y1: %Zd\nY0: %Zd\n", g[1], g[0]);
      for (j = d + 1; j -- != 0; )
        gmp_printf ("c%u: %Zd\n", j, f[j]);
      printf ("# lognorm %1.2f, alpha %1.2f, E %1.2f, %u rroots\n",
              logmu, alpha, logmu + alpha, nroots);
      printf ("\n");
      fflush (stdout);
      pthread_mutex_lock (&lock);
      found ++;
      pthread_mutex_unlock (&lock);
    }

  mpz_clear (l);
  mpz_clear (m);
  mpz_clear (t);
  mpz_clear (k);
  mpz_clear (adm1);
  mpz_clear (mtilde);
  mpz_clear (g[0]);
  mpz_clear (g[1]);
  for (j = 0; j <= d; j++)
    mpz_clear (f[j]);
  free (f);
}

static void
hash_init (hash_t H)
{
  unsigned long j;

  H->alloc = 1;
  H->p = (uint32_t*) malloc (H->alloc * sizeof (uint32_t));
  if (H->p == NULL)
    {
      fprintf (stderr, "Error, cannot allocate memory in hash_init\n");
      exit (1);
    }
  H->i = (int64_t*) malloc (H->alloc * sizeof (int64_t));
  if (H->i == NULL)
    {
      fprintf (stderr, "Error, cannot allocate memory in hash_init\n");
      exit (1);
    }
  for (j = 0; j < H->alloc; j++)
    H->p[j] = 0;
  H->size = 0;
}

static void hash_grow (hash_t H);

/* rq is a root of N = (m0 + rq)^d mod (q^2) */
static void
hash_add (hash_t H, unsigned long p, int64_t i, mpz_t m0, unsigned long ad,
          unsigned int d, mpz_t N, unsigned long q, unsigned long rq)
{
  unsigned long h;

  if (H->size >= H->alloc)
    hash_grow (H);
  if (i >= 0)
    h = i % H->alloc;
  else
    {
      h = H->alloc - ((-i) % H->alloc);
      if (h == H->alloc)
        h = 0;
    }
  while (H->p[h] != 0)
    {
      if (m0 != NULL && H->i[h] == i && H->p[h] != p)
        match (H->p[h], p, i, m0, ad, d, N, q, rq);
      if (++h == H->alloc)
        h = 0;
    }
  H->p[h] = p;
  H->i[h] = i;
  H->size ++;
}

static void
hash_clear (hash_t H)
{
  free (H->p);
  free (H->i);
}

static void
hash_grow (hash_t H)
{
  unsigned long j, old_alloc;
  uint32_t *old_p;
  int64_t *old_i;

  old_alloc = H->alloc;
  old_p = H->p;
  old_i = H->i;

  H->alloc = 2 * old_alloc;
  H->p = (uint32_t*) malloc (H->alloc * sizeof (uint32_t));
  if (H->p == NULL)
    {
      fprintf (stderr, "Error, cannot allocate memory in hash_grow\n");
      exit (1);
    }
  for (j = 0; j < H->alloc; j++)
    H->p[j] = 0;
  H->i = (int64_t*) malloc (H->alloc * sizeof (int64_t));
  if (H->i == NULL)
    {
      fprintf (stderr, "Error, cannot allocate memory in hash_grow\n");
      exit (1);
    }
  H->size = 0;
  for (j = 0; j < old_alloc; j++)
    if (old_p[j] != 0)
      hash_add (H, old_p[j], old_i[j], NULL, 0, 0, NULL, 0, 0);
  free (old_p);
  free (old_i);
}

#if 0
static double
hash_mean_value (hash_t H)
{
  double s = 0;
  unsigned long j;

  for (j = 0; j < H->alloc; j++)
    if (H->p[j] != 0)
      s += fabs ((double) H->i[j]);
  return s / (double) H->size;
}
#endif

static void
initPrimes (unsigned long P)
{
  unsigned long maxprimes = (unsigned long) 2.0 * (double) P /
    log (2.0 * (double) P) - (double) P / log ((double) P);
  unsigned long nprimes = 0;
  unsigned long p;

  Primes = (uint32_t*) malloc (maxprimes * sizeof (uint32_t));
  if (Primes == NULL)
    {
      fprintf (stderr, "Error, cannot allocate memory in initPrimes\n");
      exit (1);
    }
  for (p = 2; p < P; p = getprime (p));
  while (p <= 2 * P)
    {
      if (nprimes + 1 >= maxprimes)
        {
          maxprimes += maxprimes / 10;
          Primes = (uint32_t*) realloc (Primes, maxprimes * sizeof (uint32_t));
          if (Primes == NULL)
            {
              fprintf (stderr, "Error, cannot reallocate memory in initPrimes\n");
              exit (1);
            }
        }
      Primes[nprimes++] = p;
      p = getprime (p);
    }
  getprime (0); /* free the memory used by getprime */
  Primes[nprimes++] = 0; /* end of list */
  Primes = (uint32_t*) realloc (Primes, nprimes * sizeof (uint32_t));
  if (Primes == NULL)
    {
      fprintf (stderr, "Error, cannot reallocate memory in initPrimes\n");
      exit (1);
    }
}

static void
clearPrimes (void)
{
  free (Primes);
}

/* return 1/a mod p, assume 0 <= a < p */
static unsigned long
invert (unsigned long a, unsigned long p)
{
  modulusul_t q;
  residueul_t b;

  modul_initmod_ul (q, p);
  modul_init (b, q);
  assert (a < p);
  modul_set_ul_reduced (b, a, q);
  modul_inv (b, b, q);
  a = modul_get_ul (b, q);
  modul_clear (b, q);
  modul_clearmod (q);
  return a;
}

/* lift the n roots r[0..n-1] of N = x^d (mod p) to roots of
   N = (m0 + r)^d (mod p^2) */
void
roots_lift (unsigned long *r, mpz_t N, unsigned long d, mpz_t m0,
	    unsigned long p, unsigned int n)
{
  unsigned int j;
  mpz_t tmp, lambda;
  unsigned long inv, pp;

  mpz_init (tmp);
  mpz_init (lambda);
  pp = p * p;
  for (j = 0; j < n; j++)
    {
     /* we have for r=r[j]: r^d = N (mod p), lift mod p^2:
	 (r+lambda*p)^d = N (mod p^2) implies
	 r^d + d*lambda*p*r^(d-1) = N (mod p^2)
	 lambda = (N - r^d)/(p*d*r^(d-1)) mod p */
      mpz_ui_pow_ui (tmp, r[j], d - 1);
      mpz_mul_ui (lambda, tmp, r[j]);    /* lambda = r^d */
      mpz_sub (lambda, N, lambda);
      mpz_divexact_ui (lambda, lambda, p);
      mpz_mul_ui (tmp, tmp, d);         /* tmp = d*r^(d-1) */
      inv = invert (mpz_fdiv_ui (tmp, p), p);
      mpz_mul_ui (lambda, lambda, inv * p); /* inv * p fits in 64 bits if
					       p < 2^32 */
      mpz_add_ui (lambda, lambda, r[j]); /* now lambda^d = N (mod p^2) */

      /* subtract m0 to get roots of (m0+r)^d = N (mod p^2) */
      mpz_sub (lambda, lambda, m0);
      r[j] = mpz_fdiv_ui (lambda, pp);
    }
  mpz_clear (tmp);
  mpz_clear (lambda);
}

static void
newAlgo (mpz_t N, unsigned long d, unsigned long ad)
{
  mpz_t *f, lambda, tmp, m0, Ntilde, ump;
  unsigned long *r, *rq, i, j, p, nprimes, nr, q, nq = 0, qq, qmax;
  hash_t H;
  modulusul_t pp;
  roots_t R;
  int st;
  long M, u, ppl = 0, Mq;
  double dM, pc1, pc2;

  roots_init (R);
  mpz_init (Ntilde);
  mpz_init (m0);
  mpz_init (lambda);
  mpz_init (tmp);
  mpz_init (ump);
  r = (unsigned long*) malloc (d * sizeof (unsigned long));
  if (r == NULL)
    {
      fprintf (stderr, "Error, cannot allocate memory in newAlgo\n");
      exit (1);
    }
  rq = (unsigned long*) malloc (d * sizeof (unsigned long));
  if (rq == NULL)
    {
      fprintf (stderr, "Error, cannot allocate memory in newAlgo\n");
      exit (1);
    }
  f = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
  if (f == NULL)
    {
      fprintf (stderr, "Error, cannot allocate memory in newAlgo\n");
      exit (1);
    }
  for (i = 0; i <= d; i++)
    mpz_init (f[i]);
  mpz_set_ui (f[d], 1);

  /* compute Ntilde, ... from N, ... */
  mpz_set_ui (Ntilde, ad);
  mpz_mul_ui (Ntilde, Ntilde, d);
  mpz_pow_ui (Ntilde, Ntilde, d - 1);
  mpz_mul_ui (Ntilde, Ntilde, d);
  mpz_mul (Ntilde, Ntilde, N); /* d^d * ad^(d-1) * N */
  mpz_root (m0, Ntilde, d);

  qmax = 64;
  /* since the special-q variant gives values of i larger by q^2,
     we take as bound for i P^2*qmax^2 */
  dM = pow ((double) Primes[0] * (double) qmax, 2.0);
  if (dM > (double) LONG_MAX)
    M = LONG_MAX;
  else
    M = (long) dM;
  // printf ("M=%ld\n", M);

  st = cputime ();
  hash_init (H);
  for (nprimes = 0; (p = Primes[nprimes++]) != 0;)
    {
      if ((d * ad) % p == 0)
        continue;
      modul_initmod_ul (pp, p * p);
      /* we want p^2 | N - (m0 + i)^d, thus
         (m0 + i)^d = N (mod p^2) or m0 + i = N^(1/d) mod p^2 */
      mpz_mod_ui (f[0], Ntilde, p);
      mpz_neg (f[0], f[0]); /* f = x^d - N */
      nr = poly_roots_ulong (r, f, d, p);
      roots_lift (r, Ntilde, d, m0, p, nr);
      roots_add (R, p, nr, r);
      for (j = 0; j < nr; j++)
	{
          /* only consider r[j] and r[j] - pp */
	  ppl = (long) modul_getmod_ul (pp);
#if 1
	  hash_add (H, p, r[j], m0, ad, d, N, 1, 0);
	  hash_add (H, p, r[j] - ppl, m0, ad, d, N, 1, 0);
#else
	  u = r[j];
	  do
	    {
	      hash_add (H, p, u, m0, ad, d, N, 1, 0);
	      if (u > M - ppl)
		break;
	      else
		u += ppl; /* u <= M */
	    }
	  while (1);
	  u = r[j] - ppl;
	  do
	    {
	      hash_add (H, p, u, m0, ad, d, N, 1, 0);
	      if (u < -M + ppl)
		break;
	      else
		u -= ppl; /* u <= M */
	    }
	  while (1);
#endif
	}
      modul_clearmod (pp);
    }
  // printf ("q=1: Hsize=%lu, took %dms\n", H->size, cputime () - st);
  pc1 = 0.5 * pow ((double) H->size, 2.0);
  // printf ("%1.0f potential collisions for q=1, mean value=%1.2f\n", pc1, hash_mean_value (H) / (double) M);
  hash_clear (H);

  roots_realloc (R, R->size); /* free unused space */

  /* Special-q variant: we consider q < log(P)^2 */
  for (q = 2, pc2 = 0.0; q <= qmax; q = getprime (q))
    {
      if ((d * ad) % q == 0) /* we need q to be coprime with d and ad */
	continue;
      mpz_mod_ui (f[0], Ntilde, q);
      mpz_neg (f[0], f[0]); /* f = x^d - N */
      nq = poly_roots_ulong (rq, f, d, q);
      /* lift the roots mod q^2 */
      qq = q * q;
      Mq = M / qq; /* since q^2 is implicit in m0 + u*q^2 below */
      roots_lift (rq, Ntilde, d, m0, q, nq);
      for (i = 0; i < nq; i++)
	{
	  st = cputime ();
	  hash_init (H);
	  for (nprimes = 0; nprimes < R->size; nprimes ++)
	    {
	      residueul_t k, inv;
	      unsigned long rj;

	      p = R->p[nprimes]; /* is coprime to d*ad by construction,
				    and to q since q < Primes[0] */
	      modul_initmod_ul (pp, p * p);
	      ppl = (long) modul_getmod_ul (pp);
	      nr = R->nr[nprimes];
	      for (j = 0; j < nr; j++)
		{
		  rj = R->roots[nprimes][j];
		  /* we have p^2 | N - (m0 + rj)^d and
		     q^2 | N - (m0 + rq[i])^d
		     thus (p*q)^2 | N - (m0 + rq[i] + k*q^2)^d
		     for rq[i] + k*q^2 = rj (mod p^2), i.e.,
		     k = (rj-rq[i])/q^2 mod p^2 */
		  modul_init (k, pp);
		  modul_init (inv, pp);
		  /* since q < p and rq[i] < q^2, we have rq[i] < p^2 */
		  u = rj >= rq[i] ? rj - rq[i] : (rj + p * p) - rq[i];
		  modul_set_ul_reduced (k, u, pp);
		  modul_set_ul_reduced (inv, qq, pp);
		  modul_inv (inv, inv, pp);
		  modul_mul (k, k, inv, pp);
	      
		  /* only consider k and k - pp */
		  u = (long) modul_get_ul (k, pp);
#if 1
		  hash_add (H, p, u, m0, ad, d, N, q, rq[i]);
		  hash_add (H, p, u - ppl, m0, ad, d, N, q, rq[i]);
#else
		  do
		    {
		      hash_add (H, p, u, m0, ad, d, N, q, rq[i]);
		      if (u > Mq - ppl)
			break;
		      else
			u += ppl; /* u <= Mq */
		    }
		  while (1);
		  u = (long) modul_get_ul (k, pp) - ppl;
		  do
		    {
		      hash_add (H, p, u, m0, ad, d, N, q, rq[i]);
		      if (u < -Mq + ppl)
			break;
		      else
			u -= ppl; /* u <= Mq */
		    }
		  while (1);
#endif	      
		  modul_clear (k, pp);
		  modul_clear (inv, pp);
		}
	      modul_clearmod (pp);
	    }
	  // printf ("q=%lu, r=%lu: Hsize=%lu took %dms\n", q, rq[i], H->size, cputime () - st);
	  // printf ("q=%lu mean value=%1.2f\n", q, hash_mean_value (H) / (double) Mq);
	  hash_clear (H);
	}
      pc2 += 0.5 * pow ((double) H->size, 2.0);
    }
  // printf ("%1.0f potential collisions for q>1\n", pc2);
  getprime (0); /* free the memory used by getprime */

  pthread_mutex_lock (&lock);
  potential_collisions += pc1 +  pc2;
  pthread_mutex_unlock (&lock);

  for (i = 0; i <= d; i++)
    mpz_clear (f[i]);
  free (f);
  free (r);
  free (rq);
  mpz_clear (lambda);
  mpz_clear (tmp);
  mpz_clear (ump);
  mpz_clear (m0);
  mpz_clear (Ntilde);
  roots_clear (R);
}

void*
one_thread (void* args)
{
  tab_t *tab = (tab_t*) args;

  newAlgo (tab[0]->N, tab[0]->d, tab[0]->ad);
  pthread_exit (NULL);
}

int
main (int argc, char *argv[])
{
  mpz_t N;
  unsigned int d = 0;
  unsigned long P, admin = 0, admax = ULONG_MAX, totP = 0;
  int tries = 0, i, nthreads = 1, nr = 0, st;
  tab_t *T;
#ifdef MAX_THREADS
  pthread_t tid[MAX_THREADS];
#endif

  /* printf command-line */
  printf ("#");
  for (i = 0; i < argc; i++)
    printf (" %s", argv[i]);
  printf ("\n");
  fflush (stdout);

  mpz_init (N);

  while (argc >= 2 && argv[1][0] == '-')
    {
      if (strcmp (argv[1], "-t") == 0)
        {
          nthreads = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
      else if (strcmp (argv[1], "-nr") == 0)
        {
          nr = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
      else if (strcmp (argv[1], "-maxnorm") == 0)
        {
          max_norm = atof (argv[2]);
          argv += 2;
          argc -= 2;
        }
      else if (strcmp (argv[1], "-admin") == 0)
        {
          admin = strtoul (argv[2], NULL, 10);
          argv += 2;
          argc -= 2;
        }
      else if (strcmp (argv[1], "-admax") == 0)
        {
          admax = strtoul (argv[2], NULL, 10);
          argv += 2;
          argc -= 2;
        }
      else if (strcmp (argv[1], "-incr") == 0)
        {
          incr = strtoul (argv[2], NULL, 10);
          argv += 2;
          argc -= 2;
        }
      else if (strcmp (argv[1], "-N") == 0)
        {
          mpz_set_str (N, argv[2], 10);
          argv += 2;
          argc -= 2;
        }
      else if (strcmp (argv[1], "-d") == 0)
        {
	  d = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
      else if (strcmp (argv[1], "-v") == 0)
        {
          verbose ++;
          argv += 1;
          argc -= 1;
        }
      else
        {
          fprintf (stderr, "Invalid option: %s\n", argv[1]);
          exit (1);
        }
    }

  if (argc != 2)
    {
      fprintf (stderr, "Usage: %s [-t nthreads -nr nnn -admin nnn -admax nnn -N nnn -d nnn] P\n", argv[0]);
      exit (1);
    }

#ifdef MAX_THREADS
  if (nthreads > MAX_THREADS)
    {
      fprintf (stderr, "Error, nthreads should be <= %d\n", MAX_THREADS);
      exit (1);
    }
#endif

  if (mpz_cmp_ui (N, 0) <= 0 || d == 0)
    {
      fprintf (stderr, "Error, missing -N number or -d degree\n");
      exit (1);
    }

  P = atoi (argv[1]);
  st = cputime ();
  initPrimes (P);
  printf ("Initializing primes took %dms\n", cputime () - st);
  T = malloc (nthreads * sizeof (tab_t));
  if (T == NULL)
    {
      fprintf (stderr, "Error, cannot allocate memory in main\n");
      exit (1);
    }
  for (i = 0; i < nthreads ; i++)
    {
      mpz_init_set (T[i]->N, N);
      T[i]->d = d;
      T[i]->thread = i;
    }
  if (admin == 0)
    admin = incr;
  while (admin <= admax)
    {
      for (i = 0; i < nthreads ; i++)
        {
          tries ++;
          totP += P;
          if (verbose >= 1)
            {
              gmp_printf ("%d ad=%lu\r", tries, admin);
              fflush (stdout);
            }
          T[i]->ad = admin;
#ifndef MAX_THREADS
          newAlgo (N, d, admin);
#else
          pthread_create (&tid[i], NULL, one_thread, (void *) (T+i));
#endif
          admin += incr;
        }
#ifdef MAX_THREADS
      for (i = 0 ; i < nthreads ; i++)
        pthread_join(tid[i], NULL);
#endif
      if (totP > 100000000)
        {
          printf ("# ad=%lu time=%dms pot. collisions=%1.2e (%1.2e/s)\n",
		  admin, cputime (), potential_collisions,
		  1000.0 * potential_collisions / cputime ());
          fflush (stdout);
          totP = 0;
        }
    }
  printf ("Tried %d ad-value(s), found %d polynomial(s)\n", tries, found);
  printf ("# potential collisions=%1.2e (%1.2e/s)\n",
	  potential_collisions, 1000.0 * potential_collisions
	  / (double) cputime ());
  for (i = 0; i < nthreads ; i++)
    mpz_clear (T[i]->N);
  free (T);
  mpz_clear (N);
  clearPrimes ();

  return 0;
}
