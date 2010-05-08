/* this version implements Kleinjung's original algorithm, with a search in a
   hash table, and M = P^2 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <pthread.h>
#include "utils.h"
#include "auxiliary.h"

#define MAX_THREADS 16

typedef struct
{
  uint32_t *hash_p;         /* contains the primes */
  int64_t *hash_i;          /* contains the values of r such that p^2
                               divides N - (m0 + r)^2 */
  unsigned long hash_alloc; /* total allocated size */
  unsigned long hash_size;  /* number of entries in hash table */
} __hash_struct;
typedef __hash_struct hash_t[1];

/* read-only global variables */
int found = 0, verbose = 0;
unsigned int nr = 0; /* minimum number of real roots wanted */
double max_norm = DBL_MAX; /* maximal wanted norm (before rotation) */
pthread_mutex_t found_lock;
uint32_t *Primes;

static void
match (unsigned long p1, unsigned long p2, int64_t i, mpz_t m0,
       unsigned long ad, unsigned int d, mpz_t N)
{
  mpz_t l, mtilde, m, adm1, t, k, *f, g[2];
  unsigned int j, nroots;
  int cmp;
  double skew, logmu, alpha;

#ifdef DEBUG
  gmp_fprintf (stderr, "Found match: (%lu,%ld) (%lu,%ld)\n", p1, i, p2, i);
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
  for (j = 0; j <= d; j++)
    mpz_init (f[j]);
  /* we have l = p1*p2 */
  mpz_set_ui (l, p1);
  mpz_mul_ui (l, l, p2);
  /* mtilde = m0 + i */
  mpz_add_si (mtilde, m0, i);
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
  gmp_fprintf (stderr, "Raw polynomial:\n");
  gmp_fprintf (stderr, "Y1: %Zd\nY0: -%Zd\n", l, m);
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
    gmp_fprintf (stderr, "c%u: %Zd\n", j, f[j]);
  nroots = numberOfRealRoots (f, d, 0, 0);
  skew = SKEWNESS (f, d, SKEWNESS_DEFAULT_PREC);
  logmu = LOGNORM (f, d, skew);
  alpha = get_alpha (f, d, ALPHA_BOUND);
  fprintf (stderr, "# lognorm %1.2f, alpha %1.2f, E %1.2f, %u rroots\n",
           logmu, alpha, logmu + alpha, nroots);
  gmp_fprintf (stderr, "Optimized polynomial:\n");
#endif

  optimize (f, d, g, 0);
  nroots = numberOfRealRoots (f, d, 0, 0);
  skew = SKEWNESS (f, d, SKEWNESS_DEFAULT_PREC);
  logmu = LOGNORM (f, d, skew);
  if (nroots >= nr && logmu <= max_norm)
    {
      alpha = get_alpha (f, d, ALPHA_BOUND);
      gmp_fprintf (stderr, "n: %Zd\n", N);
      gmp_fprintf (stderr, "Y1: %Zd\nY0: %Zd\n", g[1], g[0]);
      for (j = d + 1; j -- != 0; )
        gmp_fprintf (stderr, "c%u: %Zd\n", j, f[j]);
      fprintf (stderr, "# lognorm %1.2f, alpha %1.2f, E %1.2f, %u rroots\n",
               logmu, alpha, logmu + alpha, nroots);
      fprintf (stderr, "\n");
      pthread_mutex_lock (&found_lock);
      found ++;
      pthread_mutex_unlock (&found_lock);
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

  H->hash_alloc = 1;
  H->hash_p = (uint32_t*) malloc (H->hash_alloc * sizeof (uint32_t));
  H->hash_i = (int64_t*) malloc (H->hash_alloc * sizeof (int64_t));
  for (j = 0; j < H->hash_alloc; j++)
    H->hash_p[j] = 0;
  H->hash_size = 0;
}

static void
hash_add (hash_t H, unsigned long p, int64_t i, mpz_t m0, unsigned long ad,
          unsigned int d, mpz_t N)
{
  unsigned long h;

  if (i >= 0)
    h = i % H->hash_alloc;
  else
    {
      h = H->hash_alloc - ((-i) % H->hash_alloc);
      if (h == H->hash_alloc)
        h = 0;
    }
  while (H->hash_p[h] != 0)
    {
      if (m0 != NULL && H->hash_i[h] == i && H->hash_p[h] != p)
        match (H->hash_p[h], p, i, m0, ad, d, N);
      if (++h == H->hash_alloc)
        h = 0;
    }
  H->hash_p[h] = p;
  H->hash_i[h] = i;
  H->hash_size ++;
}

static void
hash_clear (hash_t H)
{
  free (H->hash_p);
  free (H->hash_i);
}

static void
hash_grow (hash_t H)
{
  unsigned long j, old_hash_alloc;
  uint32_t *old_hash_p;
  int64_t *old_hash_i;

  old_hash_alloc = H->hash_alloc;
  old_hash_p = H->hash_p;
  old_hash_i = H->hash_i;

  H->hash_alloc = 2 * old_hash_alloc;
  H->hash_p = (uint32_t*) malloc (H->hash_alloc * sizeof (uint32_t));
  for (j = 0; j < H->hash_alloc; j++)
    H->hash_p[j] = 0;
  H->hash_i = (int64_t*) malloc (H->hash_alloc * sizeof (int64_t));
  H->hash_size = 0;
  for (j = 0; j < old_hash_alloc; j++)
    if (old_hash_p[j] != 0)
      hash_add (H, old_hash_p[j], old_hash_i[j], NULL, 0, 0, NULL);
  free (old_hash_p);
  free (old_hash_i);
}

static void
initPrimes (unsigned long P)
{
  unsigned long maxprimes = (unsigned long) 2.0 * (double) P /
    log (2.0 * (double) P) - (double) P / log ((double) P);
  unsigned long nprimes = 0;
  mpz_t p;

  Primes = (uint32_t*) malloc (maxprimes * sizeof (uint32_t));
  mpz_init (p);
  mpz_set_ui (p, P - 1);
  while (mpz_cmp_ui (p, 2 * P) <= 0)
    {
      mpz_nextprime (p, p);
      if (nprimes + 1 >= maxprimes)
        {
          maxprimes += maxprimes / 10;
          Primes = (uint32_t*) realloc (Primes, maxprimes * sizeof (uint32_t));
        }
      Primes[nprimes++] = mpz_get_ui (p);
    }
  Primes[nprimes++] = 0; /* end of list */
  Primes = (uint32_t*) realloc (Primes, nprimes * sizeof (uint32_t));
  mpz_clear (p);
}

static void
clearPrimes (void)
{
  free (Primes);
}

static void
newAlgo (mpz_t N, unsigned long d, unsigned long ad)
{
  mpz_t *f, lambda, tmp, m0, Ntilde, ump;
  unsigned long *r, i, j, p, nprimes = 0, nr;
  hash_t H;
  uint64_t pp;

  hash_init (H);
  mpz_init (Ntilde);
  mpz_init (m0);
  mpz_init (lambda);
  mpz_init (tmp);
  mpz_init (ump);
  r = (unsigned long*) malloc (d * sizeof (unsigned long));
  f = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
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

  while ((p = Primes[nprimes++]) != 0)
    {
      if ((d * ad) % p == 0)
        continue;
      pp = (uint64_t) p * (uint64_t) p;
      /* we want p^2 | N - (m0 + i)^d, thus
         (m0 + i)^d = N (mod p^2) or m0 + i = N^(1/d) mod p^2 */
      mpz_mod_ui (f[0], Ntilde, p);
      mpz_neg (f[0], f[0]); /* f = x^d - N */
      nr = poly_roots_ulong (r, f, d, p);
      for (j = 0; j < nr; j++)
        {
          /* we have r[j]^d = N (mod p), lift mod p^2:
             (r[j]+lambda*p)^d = N (mod p^2) implies
             r[j]^d + d*lambda*p*r[j]^(d-1) = N (mod p^2)
             lambda = (N - r[j]^d)/p/d/r[j]^(d-1) mod p */
          mpz_ui_pow_ui (tmp, r[j], d - 1);
          mpz_mul_ui (lambda, tmp, r[j]);
          mpz_sub (lambda, Ntilde, lambda);
#ifdef DEBUG
          if (mpz_divisible_ui_p (lambda, pui) == 0)
            {
              fprintf (stderr, "Error, N~ - r[j]^d not divisible by p\n");
              exit (1);
            }
#endif
          mpz_divexact_ui (lambda, lambda, p);
          mpz_mul_ui (tmp, tmp, d);
          mpz_set_ui (ump, p);
          mpz_invert (tmp, tmp, ump); /* FIXME: we can share the inversions */
          mpz_mul (lambda, lambda, tmp);
          mpz_mul_ui (lambda, lambda, p);
          mpz_add_ui (lambda, lambda, r[j]);
          mpz_sub (lambda, lambda, m0);
          mpz_mod_ui (lambda, lambda, pp);

          /* only consider lambda and lambda - pp */
          if (2 * H->hash_size + 1 >= H->hash_alloc)
            hash_grow (H);
          hash_add (H, p, mpz_get_si (lambda), m0, ad, d, N);
          mpz_sub_ui (lambda, lambda, pp);
          hash_add (H, p, mpz_get_si (lambda), m0, ad, d, N);
        }
    }
  for (i = 0; i <= d; i++)
    mpz_clear (f[i]);
  free (f);
  free (r);
  mpz_clear (lambda);
  mpz_clear (tmp);
  mpz_clear (ump);
  mpz_clear (m0);
  mpz_clear (Ntilde);
  hash_clear (H);
}

typedef struct
{
  mpz_t N;
  unsigned int d;
  unsigned long ad;
  int thread;
} __tab_struct;
typedef __tab_struct tab_t[1];

void*
one_thread (void* args)
{
  tab_t *tab = (tab_t*) args;

  //  printf ("thread %d deals with ad=%lu\n", tab[0]->thread, tab[0]->ad);
  fflush (stdout);
  newAlgo (tab[0]->N, tab[0]->d, tab[0]->ad);
  pthread_exit (NULL);
}

int
main (int argc, char *argv[])
{
  mpz_t N;
  unsigned int d;
  unsigned long P, admin = 60, admax = ULONG_MAX;
  int tries = 0, i, nthreads = 1, nr = 0;
  double exp_tries;
  tab_t *T;
#ifdef MAX_THREADS
  pthread_t tid[MAX_THREADS];
#endif

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
          admin = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
      else if (strcmp (argv[1], "-admax") == 0)
        {
          admax = atoi (argv[2]);
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
      fprintf (stderr, "Usage: %s [-t nthreads -nr nnn -admin nnn -admax nnn] P\n", argv[0]);
      exit (1);
    }

#ifdef MAX_THREADS
  if (nthreads > MAX_THREADS)
    {
      fprintf (stderr, "Error, nthreads should be <= %d\n", MAX_THREADS);
      exit (1);
    }
#endif

  mpz_init (N);
#if 1 /* RSA-768 */
  mpz_set_str (N, "1230186684530117755130494958384962720772853569595334792197322452151726400507263657518745202199786469389956474942774063845925192557326303453731548268507917026122142913461670429214311602221240479274737794080665351419597459856902143413", 10); /* RSA 768 */
  d = 6; /* degree */
  // admin = 1000000000020;
  // admin = 100000020;
#elif 0 /* RSA 180 */
  mpz_set_str (N, "191147927718986609689229466631454649812986246276667354864188503638807260703436799058776201365135161278134258296128109200046702912984568752800330221777752773957404540495707851421041", 10);
  d = 5;
  admin = 10000020;
#else /* aliq(564,3331) c166 */
  mpz_set_str (N, "2987193183856651631187095115639767024853796193560927033884323499962128575016083987608438436687647187187228725354979166326080044899692965798796577991358266793284662823", 10);
  d = 5;
#endif

  P = atoi (argv[1]);
  initPrimes (P);
  exp_tries = pow (log ((double) P), 2.0);
  printf ("P=%lu (%1.0f tries expected)\n", P, exp_tries);
  T = malloc (nthreads * sizeof (tab_t));
  for (i = 0; i < nthreads ; i++)
    {
      mpz_init_set (T[i]->N, N);
      T[i]->d = d;
      T[i]->thread = i;
    }
  while (admin <= admax)
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
          admin += 60;
        }
#ifdef MAX_THREADS
      for (i = 0 ; i < nthreads ; i++)
        pthread_join(tid[i], NULL);
#endif
    }
  fprintf (stderr, "Tried %d ad-values (log^2(P)=%1.0f)\n", tries, exp_tries);
  fprintf (stderr, "Found %d polynomials\n", found);
  for (i = 0; i < nthreads ; i++)
    mpz_clear (T[i]->N);
  free (T);
  mpz_clear (N);
  clearPrimes ();

  return 0;
}
