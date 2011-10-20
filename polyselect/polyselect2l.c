/*
  polyselect2l.c is a variant of Paul Zimmermann's polyselect2.c

  [1. Run and parameters]

  The parameters are similar to those in polyselect2.c, except the following,

  "-np xxx" denotes the number of special-q's trials for each ad;

  "-lq xxx" denotes the number of small factors (< 251) in the special-q;

  "-maxnorm xxx" only optimize raw polynomials with size <= xxx.
    If the raw polynomial is not good enough, we will still stream
	it to STDERR for further reference.

  Please report bugs to shi.bai AT anu.edu.au.
*/

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <pthread.h>
#include "utils.h"
#include "auxiliary.h"
#include "murphyE.h"

//#define DEBUG
#define NEW_ROOTSIEVE
#define MAX_THREADS 16
#define MAXQ 256
#define LEN_SPECIAL_Q 54
#define SPECIAL_Q {1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, ULONG_MAX}

extern int MAX_k;
void rootsieve_polyselect ( mpz_t *f, int d, mpz_t m, mpz_t l, mpz_t N, int max_k, int verbose );

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

/* structure to store information on N, d, ad, etc... */
typedef struct
{
	 mpz_t N;
	 unsigned long d;
	 unsigned long ad;
	 mpz_t Ntilde;
	 mpz_t m0;
} _header_struct;
typedef _header_struct header_t[1];

/* read-only global variables */
static int verbose = 0, incr = DEFAULT_INCR, default_MAX_k;
double max_norm = DBL_MAX; /* maximal wanted norm (before rotation) */
uint32_t *Primes;
char *out = NULL; /* output file for msieve input (msieve.dat.m) */
cado_poly best_poly, curr_poly;
double best_E = 0.0; /* Murphy's E (the larger the better) */
int nq = LEN_SPECIAL_Q; /* default has 54 primes as special-q */
int lq = 1; /* default use a single prime as special-q */
double exp_rot[] = {0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 0};

/* read-write global variables */
pthread_mutex_t lock; /* used as mutual exclusion lock for those variables */
int tot_found = 0; /* total number of polynomials */
int found = 0; /* number of polynomials below maxnorm */
double potential_collisions = 0.0, aver_lognorm = 0.0;
unsigned long collisions = 0;
unsigned long collisions_good = 0;
double total_adminus2;
double best_logmu[11];

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
roots_print (roots_t R)
{
	 unsigned int i, j;
	 for (i = 0; i < R->size; i++) {
		  fprintf (stderr, "p: %u, r: ", R->p[i]);
		  for (j = 0; j < R->nr[i]; j ++)
			   fprintf (stderr, "%lu ", R->roots[i][j]);
		  fprintf (stderr, "\n");
	 }
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


/* print poly info */
static void
print_poly_info ( mpz_t *f,
				  unsigned int d,
				  mpz_t g[2] )
{
	 unsigned int i, nroots;
	 double skew, logmu, alpha;

	 gmp_printf ("Y1: %Zd\nY0: %Zd\n", g[1], g[0]);
	 for (i = d + 1; i -- != 0; )
		  gmp_printf ("c%u: %Zd\n", i, f[i]);

	 nroots = numberOfRealRoots (f, d, 0, 0);
	 skew = L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
	 logmu = L2_lognorm (f, d, skew, DEFAULT_L2_METHOD);
	 alpha = get_alpha (f, d, ALPHA_BOUND);
	 printf ("# lognorm %1.2f, skew %1.2f, alpha %1.2f, E %1.2f,  exp_E %1.2f, %u rroots\n",
			 logmu, skew, alpha, logmu + alpha,
			 logmu - sqrt (2.0 * exp_rot[d] * log (skew)),
			 nroots);
}


/* rq is a root of N = (m0 + rq)^d mod (q^2) */
static void
match (unsigned long p1, unsigned long p2, int64_t i, mpz_t m0,
       unsigned long ad, unsigned int d, mpz_t N, unsigned long q,
       mpz_t rq)
{
	 mpz_t l, mtilde, m, adm1, t, k, *f, g[2], *fold, gold[2], qq;
	 unsigned int j;
	 int cmp;
	 double skew, logmu, E;
	 /* the expected rotation space is S^5 for degree 6 */

#ifdef DEBUG
	 gmp_printf ("Found match: (%lu,%ld) (%lu,%ld) for ad=%lu, q=%lu, rq=%Zd\n",
			 p1, i, p2, i, ad, q, rq);
	 gmp_printf ("m0=%Zd\n", m0);
#endif

	 mpz_init (l);
	 mpz_init (m);
	 mpz_init (t);
	 mpz_init (k);
	 mpz_init (qq);
	 mpz_init (adm1);
	 mpz_init (mtilde);
	 mpz_init (g[0]);
	 mpz_init (g[1]);
	 mpz_init (gold[0]);
	 mpz_init (gold[1]);
	 f = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
	 fold = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
	 if (f == NULL || fold == NULL)
	 {
		  fprintf (stderr, "Error, cannot allocate memory in match\n");
		  exit (1);
	 }
	 for (j = 0; j <= d; j++) {
		  mpz_init (f[j]);
		  mpz_init (fold[j]);
	 }
	 /* we have l = p1*p2*q */
	 mpz_set_ui (l, p1);
	 mpz_mul_ui (l, l, p2);
	 mpz_mul_ui (l, l, q);
	 /* mtilde = m0 + rq + i*q^2 */
	 mpz_set_si (mtilde, i);
	 mpz_set_ui (qq, q);
	 mpz_mul_ui (qq, qq, q);
	 mpz_mul (mtilde, mtilde, qq);
	 mpz_add (mtilde, mtilde, rq);
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

	 /* information on all polynomials */
	 pthread_mutex_lock (&lock);
	 total_adminus2 += (double) mpz_sizeinbase (f[d-2], 2);
	 collisions ++;
	 pthread_mutex_unlock (&lock);
	 // gmp_printf ("# a_{d-2}=%Zd\n", f[d-2]);

	 /* save unoptimized polynomial to fold */
	 for (i = d + 1; i -- != 0; )
		  mpz_set (fold[i], f[i]);
	 mpz_set (gold[1], g[1]);
	 mpz_set (gold[0], g[0]);

	 /* old lognorm */
	 skew = L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
	 logmu = L2_lognorm (f, d, skew, DEFAULT_L2_METHOD);

	 /* if the polynomial has norm < "-maxnorm", we optimize it */
	 if (logmu <= max_norm)
	 {

#ifdef DEBUG
		  /* unoptimized poly */
		  printf ("# Raw polynomial:\n");
		  print_poly_info (fold, d, gold);
#endif

		  /* optimize size */
		  optimize (f, d, g, 0, 1);

/* root sieve */
#ifndef NEW_ROOTSIEVE
		  unsigned long alim = 2000;
		  long jmin, kmin;
#endif
		  mpz_neg (m, g[0]);

#ifdef NEW_ROOTSIEVE
		  if (d > 4) {
			   rootsieve_polyselect (f, d, m, g[1], N, MAX_k, 0); // verbose = 2 to see details.
			   mpz_neg (g[0], m);
		  }
		  else {
			   unsigned long alim = 2000;
			   long jmin, kmin;
			   rotate (f, d, alim, m, g[1], &jmin, &kmin, 0, verbose, DEFAULT_L2_METHOD);
			   mpz_neg (g[0], m);
			   /* optimize again, but only translation */
			   optimize_aux (f, d, g, 0, 0, CIRCULAR);
		  }
#else
		  rotate (f, d, alim, m, g[1], &jmin, &kmin, 0, verbose, DEFAULT_L2_METHOD);
		  mpz_neg (g[0], m);
		  /* optimize again, but only translation */
		  optimize_aux (f, d, g, 0, 0, CIRCULAR);
#endif

		  skew = L2_skewness (f, d, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
		  logmu = L2_lognorm (f, d, skew, DEFAULT_L2_METHOD);

		  for (i = 10; i > 0 && logmu < best_logmu[i-1]; i--)
			   best_logmu[i] = best_logmu[i-1];
		  best_logmu[i] = logmu;

		  pthread_mutex_lock (&lock);
		  collisions_good ++;
		  aver_lognorm += logmu;
		  tot_found ++;
		  pthread_mutex_unlock (&lock);

		  /* MurphyE */
		  mpz_set (curr_poly->g[0], g[0]);
		  mpz_set (curr_poly->g[1], g[1]);
		  for (j = d + 1; j -- != 0; )
			   mpz_set (curr_poly->f[j], f[j]);
		  curr_poly->skew = skew;
		  E =  MurphyE (curr_poly, BOUND_F, BOUND_G, AREA, MURPHY_K);

		  mpz_neg (m, g[0]);

		  pthread_mutex_lock (&lock);
		  if (E > best_E)
		  {
			   best_E = E;
			   cado_poly_set (best_poly, curr_poly);
		  }
		  if (out != NULL) /* msieve output */
		  {
			   FILE *fp;
			   fp = fopen (out, (found == 0) ? "w" : "a");
			   if (fp == NULL)
			   {
					fprintf (stderr, "Error, cannot open file %s\n", out);
					exit (1);
			   }
			   fprintf (fp, "0");
			   for (j = d + 1; j -- != 0; )
					gmp_fprintf (fp, " %Zd", f[j]);
			   mpz_neg (m, g[0]);
			   gmp_fprintf (fp, " %Zd %Zd\n", g[1], m);
			   fclose (fp);
		  }
		  found ++;
		  pthread_mutex_unlock (&lock);

		  /* print optimized polynomial */
#ifdef DEBUG
		  gmp_printf ("# Optimized polynomial:\n");
#endif

		  gmp_printf ("n: %Zd\n", N);
		  print_poly_info ( f, d, g );
		  printf ("# Murphy's E(Bf=%.0f,Bg=%.0f,area=%.2e)=%1.2e (best so far %1.2e)\n",
				  BOUND_F, BOUND_G, AREA, E, best_E);
		  printf ("\n");
		  fflush (stdout);
	 }
	 else {
		  gmp_fprintf (stderr, "# Skip polynomial: %.2f, ad: %lu, l: %Zd, m: %Zd\n", logmu, ad, l, m);
	 }

	 mpz_clear (l);
	 mpz_clear (m);
	 mpz_clear (t);
	 mpz_clear (k);
	 mpz_clear (qq);
	 mpz_clear (adm1);
	 mpz_clear (mtilde);
	 mpz_clear (g[0]);
	 mpz_clear (g[1]);
	 mpz_clear (gold[0]);
	 mpz_clear (gold[1]);
	 for (j = 0; j <= d; j++) {
		  mpz_clear (f[j]);
		  mpz_clear (fold[j]);
	 }
	 free (f);
	 free (fold);
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
          unsigned int d, mpz_t N, unsigned long q, mpz_t rq)
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
	 mpz_t tmp;
	 mpz_init (tmp);
	 mpz_set_ui (tmp, 0);

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
			   hash_add (H, old_p[j], old_i[j], NULL, 0, 0, NULL, 0, tmp);
	 free (old_p);
	 free (old_i);
	 mpz_clear (tmp);

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


/* first combination: 0, 1, 2, 3, \cdots */
void
first_comb ( unsigned long k,
			 unsigned long *r )
{
	 unsigned long i;
	 for (i = 0; i < k; ++i)
		  r[i] = i;
}


/* next combination */
unsigned long
next_comb ( unsigned long n,
			unsigned long k,
			unsigned long *r )
{
	 unsigned long j;

	 /* if the last combination */
	 if (r[0] == (n - k)) {
		  return k;
	 }

	 /* r[k-1] doesnot equal to the n-1, just increase it */
	 j = k - 1;
	 if (r[j] < (n-1)) {
		  r[j] ++;
		  return j;
	 }

	 /* find which one we should increase */
	 while ( (r[j] - r[j-1]) == 1)
		  j --;

	 unsigned long ret = j - 1;
	 unsigned long z = ++r[j-1];

	 while (j < k) {
		  r[j] = ++z;
		  j ++;
	 }
	 return ret;
}


/* debug */
void
print_comb ( unsigned long k,
			 unsigned long *r )
{
	 unsigned long i;
	 for (i = 0; i < k; i ++)
		  fprintf (stderr, "%lu ", r[i]);
	 fprintf (stderr, "\n");
}


/* return number of n choose k */
static inline unsigned long
binom ( unsigned long n,
		unsigned long k )
{
	 if (k > n)
		  return 0;
	 if (k == 0 || k == n)
		  return 1;
	 if (2*k > n)
		  k = n - k;

	 unsigned long tot = n - k + 1, f = tot, i;
	 for (i = 2; i <= k; i++) {
		  f ++;
		  tot *= f;
		  tot /= i;
	 }
	 return tot;
}


/* init the header struct */
static void
header_init ( header_t header,
			  mpz_t N,
			  unsigned long d,
			  unsigned long ad )
{
	 /* compute Ntilde, m0 */
	 mpz_init_set (header->N, N);
	 mpz_init (header->Ntilde);
	 mpz_init (header->m0);
	 header->d = d;
	 header->ad = ad;

	 /* compute Ntilde, ... from N, ... */
	 mpz_set_ui (header->Ntilde, header->ad);
	 mpz_mul_ui (header->Ntilde, header->Ntilde, header->d);
	 mpz_pow_ui (header->Ntilde, header->Ntilde, header->d - 1);
	 mpz_mul_ui (header->Ntilde, header->Ntilde, header->d);
	 mpz_mul (header->Ntilde, header->Ntilde, header->N); /* d^d * ad^(d-1) * N */
	 mpz_root (header->m0, header->Ntilde, header->d);
}


static void
header_clear ( header_t header )
{
	 mpz_clear (header->m0);
	 mpz_clear (header->Ntilde);
	 mpz_clear (header->N);
}


static void
collision_on_p ( header_t header,
				 roots_t R )
{
	 unsigned long i, j, nprimes, p, *rp, nrp;
	 long ppl = 0;
	 double pc1;
	 modulusul_t pp;
	 mpz_t *f, tmp;

	 mpz_init_set_ui (tmp, 0);
	 f = (mpz_t*) malloc ((header->d + 1) * sizeof (mpz_t));
	 if (f == NULL)
	 {
		  fprintf (stderr, "Error, cannot allocate memory in collision_on_p\n");
		  exit (1);
	 }
	 for (i = 0; i <= header->d; i++)
		  mpz_init (f[i]);
	 mpz_set_ui (f[header->d], 1);

	 rp = (unsigned long*) malloc (header->d * sizeof (unsigned long));
	 if (rp == NULL)
	 {
		  fprintf (stderr, "Error, cannot allocate memory in collision_on_p\n");
		  exit (1);
	 }

	 hash_t H;
	 hash_init (H);
	 for (nprimes = 0; (p = Primes[nprimes++]) != 0;)
	 {
		  if ((header->d * header->ad) % p == 0)
			   continue;
		  modul_initmod_ul (pp, p * p);
		  /* we want p^2 | N - (m0 + i)^d, thus
			 (m0 + i)^d = N (mod p^2) or m0 + i = N^(1/d) mod p^2 */
		  mpz_mod_ui (f[0], header->Ntilde, p);
		  mpz_neg (f[0], f[0]); /* f = x^d - N */
		  nrp = poly_roots_ulong (rp, f, header->d, p);
		  roots_lift (rp, header->Ntilde, header->d, header->m0, p, nrp);
		  roots_add (R, p, nrp, rp);
		  for (j = 0; j < nrp; j++)
		  {
			   /* only consider r[j] and r[j] - pp */
			   ppl = (long) modul_getmod_ul (pp);
			   hash_add (H, p, rp[j], header->m0, header->ad, header->d, header->N, 1, tmp);
			   hash_add (H, p, rp[j] - ppl, header->m0, header->ad, header->d, header->N, 1, tmp);
		  }
		  modul_clearmod (pp);
	 }

	 /* if the hash table contains n entries, each one smaller than (2P)^2,
		the number of potential collisions is about 1/2n^2/(2P)^2 */
	 fprintf (stderr, "# p hash_size: %lu for ad = %lu\n", H->size, header->ad);
	 pc1 = 0.5 * pow ((double) H->size / (double) Primes[nprimes - 2], 2.0);
	 hash_clear (H);
	 roots_realloc (R, R->size); /* free unused space */
	 for (i = 0; i <= header->d; i++)
		  mpz_clear (f[i]);
	 free (f);
	 free (rp);
	 mpz_clear (tmp);

	 pthread_mutex_lock (&lock);
	 potential_collisions += pc1;
	 pthread_mutex_unlock (&lock);
}


/* prepare special-q's roots */
static void
comp_sq_r ( header_t header,
			roots_t SQ_R )
{
	 unsigned long i, Q[] = SPECIAL_Q, q, *rq, nrq;
	 modulusul_t qq;
	 mpz_t *f;

	 f = (mpz_t*) malloc ((header->d + 1) * sizeof (mpz_t));
	 if (f == NULL)
	 {
		  fprintf (stderr, "Error, cannot allocate memory in comp_sq_r\n");
		  exit (1);
	 }
	 for (i = 0; i <= header->d; i++)
		  mpz_init (f[i]);
	 mpz_set_ui (f[header->d], 1);

	 rq = (unsigned long*) malloc (header->d * sizeof (unsigned long));
	 if (rq == NULL)
	 {
		  fprintf (stderr, "Error, cannot allocate memory in comp_sq_q\n");
		  exit (1);
	 }

	 /* prepare the special-q's */
	 for (i = 1; (q = Q[i]) < ULONG_MAX; i++)
	 {
		  if ((header->d * header->ad) % q == 0)
			   continue;

		  modul_initmod_ul (qq, q * q);
		  mpz_mod_ui (f[0], header->Ntilde, q);
		  mpz_neg (f[0], f[0]); /* f = x^d - N */
		  nrq = poly_roots_ulong (rq, f, header->d, q);
		  roots_lift (rq, header->Ntilde, header->d, header->m0, q, nrq);
		  roots_add (SQ_R, q, nrq, rq);
		  modul_clearmod (qq);
	 }

	 free(rq);
	 for (i = 0; i <= header->d; i++)
		  mpz_clear (f[i]);
	 free (f);
	 roots_realloc (SQ_R, SQ_R->size); /* free unused space */
}


/* crt, set r and qqz */
static inline void
crt_sq ( mpz_t qqz,
		 mpz_t r,
		 unsigned long *q,
		 unsigned long *rq )
{
	 mpz_t prod, pprod, mod, inv, sum;
	 int i;
	 unsigned long qq[lq];

	 mpz_init_set_ui (prod, 1);
	 mpz_init (pprod);
	 mpz_init (mod);
	 mpz_init (inv);
	 mpz_init_set_ui (sum, 0);

	 for (i = 0; i < lq; i ++) {
		  qq[i] = q[i] * q[i];
		  mpz_mul_ui (prod, prod, qq[i]);
	 }

	 for (i = 0; i < lq; i ++) {
		  mpz_divexact_ui (pprod, prod, qq[i]);
		  mpz_set_ui (mod, qq[i]);
		  mpz_invert (inv, pprod, mod);
		  mpz_mul_ui (inv, inv, rq[i]);
		  mpz_mul (inv, inv, pprod);
		  mpz_add (sum, sum, inv);
	 }

	 mpz_mod (sum, sum, prod);
	 mpz_set (r, sum);
	 mpz_set (qqz, prod);

	 mpz_clear (prod);
	 mpz_clear (pprod);
	 mpz_clear (mod);
	 mpz_clear (inv);
	 mpz_clear (sum);
}


/* collision on each special-q */
static void
collision_on_each_sq ( header_t header,
					   roots_t R,
					   unsigned long q,
					   mpz_t qqz,
					   mpz_t rqqz )
{
	 unsigned int nr, j;
	 unsigned long nprimes, p, qql, rqql, rp;
	 mpz_t tmp;
	 long ppl, u;
	 modulusul_t pp;
	 mpz_init (tmp);

#ifdef DEBUG
	 mpz_t tmp_debug, tmp_debug2;
	 mpz_init_set (tmp_debug, header->Ntilde);
	 mpz_mod (tmp_debug, tmp_debug, qqz);
	 mpz_init_set (tmp_debug2, header->m0);
	 mpz_add (tmp_debug2, tmp_debug2, rqqz);
	 mpz_pow_ui (tmp_debug2, tmp_debug2, header->d);
	 mpz_mod (tmp_debug2, tmp_debug2, qqz);
	 if (mpz_cmp (tmp_debug, tmp_debug2) != 0)
	 {
		  fprintf (stderr, "Error: crt root is wrong in collision_on_each_sq\n");
		  exit (1);
	 }
	 mpz_clear (tmp_debug);
	 mpz_clear (tmp_debug2);
#endif

	 hash_t H;
	 hash_init (H);

	 for (nprimes = 0; nprimes < R->size; nprimes ++)
	 {
		  residueul_t k, inv;

		  /* set p, p^2, ppl */
		  p = R->p[nprimes];
		  modul_initmod_ul (pp, p * p);
		  ppl = (long) modul_getmod_ul (pp);
		  nr = R->nr[nprimes];

		  /* set rq (mod p^2) */
		  mpz_mod_ui (tmp, rqqz, modul_getmod_ul (pp));
		  rqql = mpz_get_ui (tmp);

		  /* set q^2 (mod p^2) */
		  mpz_mod_ui (tmp, qqz, modul_getmod_ul (pp));
		  qql = mpz_get_ui (tmp);

		  for (j = 0; j < nr; j++)
		  {
			   rp = R->roots[nprimes][j];
			   /* we have p^2 | N - (m0 + rp)^d and
				  q^2 | N - (m0 + rq)^d
				  thus (p*q)^2 | N - (m0 + rq + k*q^2)^d
				  for rq[i] + k*q^2 = rp (mod p^2), i.e.,
				  k = (rp-rq[i])/q^2 mod p^2 */
			   modul_init (k, pp);
			   modul_init (inv, pp);
			   /* since q < p and rq[i] < q^2, we have rq[i] < p^2 */
			   u = (long) ( (rp >= rqql) ? rp - rqql : (rp + p * p) - rqql);
			   assert (u < ppl);
			   modul_set_ul_reduced (k, u, pp);
			   assert (qql < (unsigned long) ppl);
			   modul_set_ul_reduced (inv, qql, pp);
			   modul_inv (inv, inv, pp);
			   modul_mul (k, k, inv, pp);

			   /* only consider k and k - pp */
			   u = (long) modul_get_ul (k, pp);

			   //fprintf (stderr, "u: %ld\n", u);

#ifdef DEBUG
			   mpz_t Ntmp, m0tmp, tmp;
			   mpz_init_set (Ntmp, header->Ntilde);
			   mpz_mod_ui (Ntmp, Ntmp, ppl); // tildeN (mod p^2)

			   mpz_init_set (m0tmp, header->m0);
			   mpz_add (m0tmp, m0tmp, rqqz); // m0 + rq
			   mpz_init_set_ui (tmp, u);
			   mpz_mul (tmp, tmp, qqz); // i*q^2
			   mpz_add (m0tmp, m0tmp, tmp); // m0 + rq + i*q^2

			   mpz_pow_ui (m0tmp, m0tmp, header->d);
			   mpz_mod_ui (m0tmp, m0tmp, ppl);
			   if (mpz_cmp (m0tmp, Ntmp) != 0)
			   {
					fprintf (stderr, "Error: i computation is wrong in collision_on_each_sq\n");

					gmp_printf ("Details: (p=%lu, i(mod p^2)=%ld) for ad=%lu, q=%lu, rq=%Zd, rp=%lu\n",
								p, u, header->ad, q, rqqz, rp);
					gmp_printf ("m0=%Zd\n", header->m0);
					gmp_printf ("Ntilde=%Zd\n", header->Ntilde);

					exit (1);
			   }
			   mpz_clear (Ntmp);
			   mpz_clear (m0tmp);
			   mpz_clear (tmp);
#endif

			   hash_add (H, p, u, header->m0, header->ad, header->d, header->N, q, rqqz);
			   hash_add (H, p, u - ppl, header->m0, header->ad, header->d, header->N, q, rqqz);

			   modul_clear (k, pp);
			   modul_clear (inv, pp);
		  }  // next rp
		  modul_clearmod (pp);
	 } // next p

	 //fprintf (stderr, "# - q hash_size (q=%lu): %lu\n", q, H->size);
	 double pc2 = 0.5 * pow ((double) H->size / (double) R->p[nprimes-1],
							 2.0);

	 hash_clear (H);

	 pthread_mutex_lock (&lock);
	 potential_collisions += pc2;
	 pthread_mutex_unlock (&lock);

	 mpz_clear (tmp);
}

/* given individual q's, return crted rq */
static inline unsigned long
return_q_rq ( roots_t SQ_R,
			  unsigned long *idx_q,
			  unsigned long k,
			  mpz_t qqz,
			  mpz_t rqqz )
{
	 unsigned long i, j, idv_q[k], idv_rq[k], q = 1;

	 /* q and roots */
	 for (i = 0; i < k; i ++) {
		  idv_q[i] = SQ_R->p[idx_q[i]];
		  q = q * idv_q[i];
		  srand (time(NULL));
		  j = rand() % SQ_R->nr[idx_q[i]];
		  // j = 0;
		  idv_rq[i] = SQ_R->roots[idx_q[i]][j];
	 }

#if 0
	 for (i = 0; i < k; i ++) {
		  fprintf (stderr, "(%lu:%lu) ", idv_q[i], idv_rq[i]);
	 }
	 gmp_fprintf (stderr, "%Zd\n", rqqz);
#endif

	 /* crt roots */
	 crt_sq (qqz, rqqz, idv_q, idv_rq);

	 return q;
}


/* collision on special-q, call collisio_on_each_sq */
static void
collision_on_sq ( header_t header,
				  roots_t R )
{
	 mpz_t qqz, rqqz;
	 roots_t SQ_R;

	 mpz_init (qqz);
	 mpz_init (rqqz);

	 roots_init (SQ_R);
	 comp_sq_r (header, SQ_R);
	 unsigned long k = lq, n = SQ_R->size, q;
	 // roots_print (SQ_R);

	 /* less than lq special primes having roots for this ad */
	 if (n == 0 || n < k)
		  return;

	 unsigned long idx_q[n], tot, i;
	 /* if tot (n, k) < wanted, use tot as wanted */
	 tot =  binom (n, k);

     if (tot > (unsigned long) nq)
		  tot = (unsigned long) nq;
	 // fprintf (stderr, "n=%lu, k=%lu, (n,k)=%lu, nq:%d\n", n, k, tot, nq);

	 /* enumerate each combination */
	 first_comb (n, idx_q);
	 q = return_q_rq (SQ_R, idx_q, k, qqz, rqqz);
	 collision_on_each_sq ( header,
							R,
							q,
							qqz,
							rqqz );

	 //fprintf (stderr, "# - q ");
	 for (i = 1; i < tot; i ++) {
		  next_comb (n, k, idx_q); // print_comb (k, idx_q);

		  /* if (i % 100 == 0) */
		  /* 	   fprintf (stderr, "%lu ", i); */

		  q = return_q_rq (SQ_R, idx_q, k, qqz, rqqz);
		  collision_on_each_sq ( header,
								 R,
								 q,
								 qqz,
								 rqqz );
	 }

	 mpz_clear (qqz);
	 mpz_clear (rqqz);
	 roots_clear (SQ_R);
}

static void
newAlgo (mpz_t N, unsigned long d, unsigned long ad)
{
	 header_t header;
	 header_init (header, N, d, ad);

	 roots_t R;
	 roots_init (R);

	 collision_on_p (header, R);
	 collision_on_sq (header, R);

	 roots_clear (R);
	 header_clear (header);
}

void*
one_thread (void* args)
{
	 tab_t *tab = (tab_t*) args;
	 newAlgo (tab[0]->N, tab[0]->d, tab[0]->ad);
	 return NULL;
}

#define TARGET_TIME 10000000 /* print stats every TARGET_TIME milliseconds */

static void
usage (char *argv)
{
	 fprintf (stderr, "Usage: %s [options] P\n", argv);
	 fprintf (stderr, "Parameters and options:\n");
	 fprintf (stderr, "P            --- degree-1 coefficient of g(x) has\n");
	 fprintf (stderr, "                 two prime factors in [P,2P]\n");
	 fprintf (stderr, "-v           --- verbose mode\n");
	 fprintf (stderr, "-q           --- quiet mode\n");
	 fprintf (stderr, "-t nnn       --- use n threads (default 1)\n");
	 fprintf (stderr, "-admin nnn   --- start from ad=nnn (default 0)\n");
	 fprintf (stderr, "-admax nnn   --- stop at ad=nnn\n");
	 fprintf (stderr, "-incr nnn    --- forced factor of ad (default 60)\n");
	 fprintf (stderr, "-N nnn       --- input number\n");
	 fprintf (stderr, "-degree nnn  --- wanted polynomial degree\n");
	 fprintf (stderr, "-nq nnn      --- number of special-q's considered for each ad\n");
	 fprintf (stderr, "-lq nnn      --- number of factors in the special-q\n");
	 fprintf (stderr, "-kmax nnn    --- rotation bound (default %d)\n",
			  default_MAX_k);
	 fprintf (stderr, "-save xxx    --- save state in file xxx\n");
	 fprintf (stderr, "-resume xxx  --- resume state from file xxx\n");
	 fprintf (stderr, "-maxnorm xxx --- only optimize polynomials with norm <= xxx\n");
	 fprintf (stderr, "-maxtime xxx --- stop the search after xxx seconds\n");
	 fprintf (stderr, "-out xxx     --- for msieve-format output\n");
	 exit (1);
}

int
main (int argc, char *argv[])
{
	 int argc0 = argc;
	 char **argv0 = argv, *save = NULL, *resume = NULL;
	 double st0 = seconds (), maxtime = DBL_MAX;
	 mpz_t N;
	 unsigned int d = 0;
	 unsigned long P, admin = 0, admax = ULONG_MAX;
	 int tries = 0, i, nthreads = 1, st, target_time = TARGET_TIME;
	 tab_t *T;
	 FILE *fp;
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
	 default_MAX_k = MAX_k;
	 cado_poly_init (best_poly);
	 cado_poly_init (curr_poly);

	 if (argc == 1)
		  usage (argv0[0]);

	 while (argc >= 2 && argv[1][0] == '-')
	 {
		  if (strcmp (argv[1], "-t") == 0)
		  {
			   nthreads = atoi (argv[2]);
			   argv += 2;
			   argc -= 2;
		  }
		  else if (argc >= 3 && strcmp (argv[1], "-maxnorm") == 0)
		  {
			   max_norm = atof (argv[2]);
			   argv += 2;
			   argc -= 2;
		  }
		  else if (argc >= 3 && strcmp (argv[1], "-maxtime") == 0)
		  {
			   maxtime = atof (argv[2]);
			   argv += 2;
			   argc -= 2;
		  }
		  else if (argc >= 3 && strcmp (argv[1], "-admin") == 0)
		  {
			   double d;
			   d = strtod (argv[2], NULL);
			   if (d > (double) ULONG_MAX)
			   {
					fprintf (stderr, "Error, too large value of admin\n");
					exit (1);
			   }
			   admin = (unsigned long) d;
			   argv += 2;
			   argc -= 2;
		  }
		  else if (argc >= 3 && strcmp (argv[1], "-admax") == 0)
		  {
			   double d;
			   d = strtod (argv[2], NULL);
			   if (d > (double) ULONG_MAX)
			   {
					fprintf (stderr, "Error, too large value of admax\n");
					exit (1);
			   }
			   admax = (unsigned long) d;
			   argv += 2;
			   argc -= 2;
		  }
		  else if (argc >= 3 && strcmp (argv[1], "-incr") == 0)
		  {
			   incr = strtoul (argv[2], NULL, 10);
			   argv += 2;
			   argc -= 2;
		  }
		  else if (argc >= 3 && strcmp (argv[1], "-N") == 0)
		  {
			   mpz_set_str (N, argv[2], 10);
			   argv += 2;
			   argc -= 2;
		  }
		  else if (argc >= 3 && strcmp (argv[1], "-degree") == 0)
		  {
			   d = atoi (argv[2]);
			   argv += 2;
			   argc -= 2;
		  }
		  else if (argc >= 3 && strcmp (argv[1], "-nq") == 0)
		  {
			   nq = atoi (argv[2]);
			   argv += 2;
			   argc -= 2;
		  }
		  else if (argc >= 3 && strcmp (argv[1], "-lq") == 0)
		  {
			   lq = atoi (argv[2]);
			   argv += 2;
			   argc -= 2;
		  }
		  else if (argc >= 3 && strcmp (argv[1], "-save") == 0)
		  {
			   save = argv[2];
			   argv += 2;
			   argc -= 2;
		  }
		  else if (argc >= 3 && strcmp (argv[1], "-resume") == 0)
		  {
			   resume = argv[2];
			   argv += 2;
			   argc -= 2;
		  }
		  else if (argc >= 3 && strcmp (argv[1], "-out") == 0)
		  {
			   out = argv[2];
			   argv += 2;
			   argc -= 2;
		  }
		  else if (argc >= 3 && strcmp (argv[1], "-kmax") == 0)
		  {
			   MAX_k = atoi (argv[2]);
			   argv += 2;
			   argc -= 2;
		  }
		  else if (strcmp (argv[1], "-v") == 0)
		  {
			   verbose ++;
			   argv += 1;
			   argc -= 1;
		  }
		  else if (strcmp (argv[1], "-q") == 0)
		  {
			   verbose = -1;
			   argv += 1;
			   argc -= 1;
		  }
		  else
		  {
			   fprintf (stderr, "Invalid option: %s\n", argv[1]);
			   exit (1);
		  }
	 }

	 if (mpz_cmp_ui (N, 0) == 0)
	 {
		  int ret;

		  ret = gmp_scanf ("n: %Zd\n", N);
		  if (ret != 1)
		  {
			   fprintf (stderr, "Error, input number N cannot be read\n");
			   exit (1);
		  }
	 }

	 if (lq < 1 || nq < 1)
	 {
		  fprintf (stderr, "Error, number of factors in special-q should >= 1 and/or number of special-q's should >=1\n");
		  exit (1);
	 }
	 if (lq == 1)
	 {
		  nq = LEN_SPECIAL_Q;
	 }

	 mpz_set (best_poly->n, N);
	 mpz_set (curr_poly->n, N);

	 if (argc != 2)
		  usage (argv0[0]);

#ifdef MAX_THREADS
	 if (nthreads > MAX_THREADS)
	 {
		  fprintf (stderr, "Error, nthreads should be <= %d\n", MAX_THREADS);
		  exit (1);
	 }
#endif

	 if (mpz_cmp_ui (N, 0) <= 0)
	 {
		  fprintf (stderr, "Error, missing input number (-N option)\n");
		  exit (1);
	 }

	 if (d == 0)
	 {
		  fprintf (stderr, "Error, missing degree (-d option)\n");
		  exit (1);
	 }
	 best_poly->degree = d;
	 best_poly->degreeg = 1;
	 curr_poly->degree = d;
	 curr_poly->degreeg = 1;

	 if (resume != NULL)
	 {
		  fp = fopen (resume, "r");
		  if (fp == NULL)
		  {
			   fprintf (stderr, "Cannot open resume file %s\n", resume);
			   exit (1);
		  }
		  if (fscanf (fp, "%lu", &admin) != 1)
		  {
			   fprintf (stderr, "Cannot read ad value from resume file %s\n",
						resume);
			   exit (1);
		  }
		  fclose (fp);
	 }

	 for (i = 0; i <= 10; i++)
		  best_logmu[i] = 999.9;

	 P = atoi (argv[1]);
	 st = cputime ();
	 initPrimes (P);
	 printf ("# Initializing primes took %dms\n", cputime () - st);
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

	 if (incr <= 0)
	 {
		  fprintf (stderr, "Error, incr should be positive\n");
		  exit (1);
	 }

	 if (admin == 0)
		  admin = incr;

	 /* force admin to be divisible by incr */
	 admin = ((admin + incr - 1) / incr) * incr; /* incr * ceil (admin/incr) */

	 while (admin <= admax && seconds () - st0 <= maxtime)
	 {
		  for (i = 0; i < nthreads ; i++)
		  {
			   tries ++;
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

		  if (save != NULL)
		  {
			   fp = fopen (save, "w");
			   if (fp == NULL)
			   {
					fprintf (stderr, "Cannot open save file %s\n", save);
					exit (1);
			   }
			   if (fprintf (fp, "%lu\n", admin) < 0)
			   {
					fprintf (stderr, "Cannot print ad value to save file %s\n",
							 save);
					exit (1);
			   }
			   fclose (fp);
		  }

		  if (cputime () > target_time || verbose > 0)
		  {
			   printf ("# ad=%lu time=%dms exp. coll.=%1.1f, got %lu (%0.3f/s) with %lu good ones, av. lognorm=%1.2f\n",
					   admin,
					   cputime (),
					   potential_collisions,
					   collisions,
					   1000.0 * (double) collisions / cputime (),
					   collisions_good,
					   aver_lognorm / collisions );
			   fflush (stdout);
			   target_time += TARGET_TIME;
		  }
	 }

	 printf ("# Tried %d ad-value(s), found %d polynomial(s), %d below maxnorm\n",
			 tries, tot_found, found);
	 printf ("# potential collisions=%1.2e (%1.2e/s)\n",
			 potential_collisions, 1000.0 * potential_collisions
			 / (double) cputime ());
	 printf (" av. adm2:%1.0f", total_adminus2 / (double) collisions);

	 /* print best 10 values of logmu */
	 printf ("# best logmu:");
	 for (i = 0; i < 10; i++)
		  printf (" %1.2f", best_logmu[i]);
	 printf ("\n");

	 /* print total time (format for cpu_time.sh) */
	 printf ("# Total phase took %.2fs\n", seconds () - st0);

	 if (best_E == 0.0)
		  printf ("No polynomial found, please increase the ad range or decrease P\n");
	 else
		  print_poly (stdout, best_poly, argc0, argv0, st0, 1 /* raw */);

	 for (i = 0; i < nthreads ; i++)
		  mpz_clear (T[i]->N);
	 free (T);
	 mpz_clear (N);
	 clearPrimes ();
	 cado_poly_clear (best_poly);
	 cado_poly_clear (curr_poly);

	 return 0;
}
