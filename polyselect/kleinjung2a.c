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
  unsigned long *hash_p;    /* contains the primes */
  mpz_t *hash_i;            /* contains the values of r such that p^2
                               divides N - (m0 + r)^2 */
  unsigned long hash_alloc; /* total allocated size */
  unsigned long hash_size;  /* number of entries in hash table */
} __hash_struct;
typedef __hash_struct hash_t[1];

int found = 0;
unsigned int nr = 0; /* minimum number of real roots wanted */
double max_norm = DBL_MAX; /* maximal wanted norm (before rotation) */
pthread_mutex_t found_lock;

static void
match (unsigned long p1, unsigned long p2, mpz_t i, mpz_t m0, unsigned long ad,
       unsigned int d, mpz_t N)
{
  mpz_t l, mtilde, m, adm1, t, k, *f, g[2];
  unsigned int j, nroots;
  int cmp;
  double skew, logmu, alpha;

  gmp_fprintf (stderr, "Found match: (%lu,%Zd) (%lu,%Zd)\n", p1, i, p2, i);

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
  mpz_add (mtilde, m0, i);
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
  gmp_fprintf (stderr, "Raw polynomial:\n");
  gmp_fprintf (stderr, "Y1: %Zd\nY0: -%Zd\n", l, m);
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
  for (j = d + 1; j -- != 0; )
    gmp_fprintf (stderr, "c%u: %Zd\n", j, f[j]);
  nroots = numberOfRealRoots (f, d, 0, 0);
  fprintf (stderr, "# %u real root(s)\n", nroots);

  skew = SKEWNESS (f, d, SKEWNESS_DEFAULT_PREC);
  logmu = LOGNORM (f, d, skew);
  alpha = get_alpha (f, d, ALPHA_BOUND);
  fprintf (stderr, "# lognorm: %1.2f, alpha: %1.2f E=%1.2f\n", logmu, alpha,
           logmu + alpha);

  optimize (f, d, g, 0);
  skew = SKEWNESS (f, d, SKEWNESS_DEFAULT_PREC);
  logmu = LOGNORM (f, d, skew);
  alpha = get_alpha (f, d, ALPHA_BOUND);
  gmp_fprintf (stderr, "Optimized polynomial:\n");
  gmp_fprintf (stderr, "n: %Zd\n", N);
  gmp_fprintf (stderr, "Y1: %Zd\nY0: %Zd\n", g[1], g[0]);
  for (j = d + 1; j -- != 0; )
    gmp_fprintf (stderr, "c%u: %Zd\n", j, f[j]);
  nroots = numberOfRealRoots (f, d, 0, 0);
  fprintf (stderr, "# %u real root(s)\n", nroots);
  fprintf (stderr, "# lognorm: %1.2f, alpha: %1.2f E=%1.2f\n", logmu, alpha,
           logmu + alpha);
  fprintf (stderr, "\n");

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

  if (nroots >= nr && logmu <= max_norm)
    {
      pthread_mutex_lock (&found_lock);
      found ++;
      pthread_mutex_unlock (&found_lock);
    }
}

static void
hash_init (hash_t H)
{
  unsigned long j;

  H->hash_alloc = 1;
  H->hash_p = (unsigned long*) malloc (H->hash_alloc * sizeof (unsigned long));
  H->hash_i = (mpz_t*) malloc (H->hash_alloc * sizeof (mpz_t));
  for (j = 0; j < H->hash_alloc; j++)
    {
      H->hash_p[j] = 0;
      mpz_init (H->hash_i[j]);
    }
  H->hash_size = 0;
}

static void
hash_add (hash_t H, unsigned long p, mpz_t i, mpz_t m0, unsigned long ad,
          unsigned int d, mpz_t N)
{
  unsigned long h = mpz_fdiv_ui (i, H->hash_alloc);

  while (H->hash_p[h] != 0)
    {
      if (m0 != NULL && mpz_cmp (H->hash_i[h], i) == 0 && H->hash_p[h] != p)
        match (H->hash_p[h], p, i, m0, ad, d, N);
      if (++h == H->hash_alloc)
        h = 0;
    }
  H->hash_p[h] = p;
  mpz_set (H->hash_i[h], i);
  H->hash_size ++;
}

static void
hash_clear (hash_t H)
{
  unsigned long j;

  for (j = 0; j < H->hash_alloc; j++)
    mpz_clear (H->hash_i[j]);
  free (H->hash_p);
  free (H->hash_i);
}

static void
hash_grow (hash_t H)
{
  unsigned long j, *old_hash_p, old_hash_alloc;
  mpz_t *old_hash_i;

  old_hash_alloc = H->hash_alloc;
  old_hash_p = H->hash_p;
  old_hash_i = H->hash_i;

  H->hash_alloc = 2 * old_hash_alloc;
  H->hash_p = (unsigned long*) malloc (H->hash_alloc *
                                       sizeof (unsigned long));
  for (j = 0; j < H->hash_alloc; j++)
    H->hash_p[j] = 0;
  H->hash_i = (mpz_t*) malloc (H->hash_alloc * sizeof (mpz_t));
  for (j = 0; j < H->hash_alloc; j++)
    mpz_init (H->hash_i[j]);
  H->hash_size = 0;
  for (j = 0; j < old_hash_alloc; j++)
    if (old_hash_p[j] != 0)
      hash_add (H, old_hash_p[j], old_hash_i[j], NULL, 0, 0, NULL);
  for (j = 0; j < old_hash_alloc; j++)
    mpz_clear (old_hash_i[j]);
  free (old_hash_p);
  free (old_hash_i);
}

static void
newAlgo (mpz_t N, unsigned long d, unsigned long P, unsigned long ad)
{
  mpz_t p, pp, pmax, *f, lambda, tmp, m0, Ntilde;
  unsigned long *r, i, j, pui, nprimes = 0, nr;
  hash_t H;

  hash_init (H);
  mpz_init (Ntilde);
  mpz_init (p);
  mpz_init (m0);
  mpz_init (pp);
  mpz_init (pmax);
  mpz_init (lambda);
  mpz_init (tmp);
  mpz_set_ui (p, P - 1);
  mpz_set_ui (pmax, P);
  mpz_mul_2exp (pmax, pmax, 1);
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

  while (mpz_cmp (p, pmax) <= 0)
    {
      mpz_nextprime (p, p);
      nprimes ++;
      mpz_mul (pp, p, p);
      pui = mpz_get_ui (p);
      /* we want p^2 | N - (m0 + i)^d, thus
         (m0 + i)^d = N (mod p^2) or m0 + i = N^(1/d) mod p^2 */
      mpz_mod (f[0], Ntilde, p);
      mpz_neg (f[0], f[0]); /* f = x^d - N */
      nr = poly_roots_ulong (r, f, d, pui);
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
          mpz_divexact_ui (lambda, lambda, pui);
          mpz_mul_ui (tmp, tmp, d);
          mpz_invert (tmp, tmp, p);
          mpz_mul (lambda, lambda, tmp);
          mpz_mul_ui (lambda, lambda, pui);
          mpz_add_ui (lambda, lambda, r[j]);
          mpz_sub (lambda, lambda, m0);
          mpz_mod (lambda, lambda, pp);

          /* only consider lambda and lambda - pp */
          if (2 * H->hash_size + 1 >= H->hash_alloc)
            hash_grow (H);
          hash_add (H, pui, lambda, m0, ad, d, N);
          mpz_sub (lambda, lambda, pp);
          hash_add (H, pui, lambda, m0, ad, d, N);
        }
    }
  for (i = 0; i <= d; i++)
    mpz_clear (f[i]);
  free (f);
  free (r);
  mpz_clear (p);
  mpz_clear (pp);
  mpz_clear (pmax);
  mpz_clear (lambda);
  mpz_clear (tmp);
  mpz_clear (m0);
  mpz_clear (Ntilde);
  hash_clear (H);
}

typedef struct
{
  mpz_t N;
  unsigned int d;
  unsigned long P;
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
  newAlgo (tab[0]->N, tab[0]->d, tab[0]->P, tab[0]->ad);
  pthread_exit (NULL);
}

int
main (int argc, char *argv[])
{
  mpz_t N;
  unsigned int d;
  unsigned long P, ad;
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
      else
        {
          fprintf (stderr, "Invalid option: %s\n", argv[1]);
          exit (1);
        }
    }

  if (argc != 2)
    {
      fprintf (stderr, "Usage: %s [-t nthreads -nr nnn] P\n", argv[0]);
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
  // ad = 1000000000020;
  ad = 100000020;
#elif 0 /* RSA 180 */
  mpz_set_str (N, "191147927718986609689229466631454649812986246276667354864188503638807260703436799058776201365135161278134258296128109200046702912984568752800330221777752773957404540495707851421041", 10);
  d = 5;
  ad = 10000020;
#else /* aliq(564,3331) c166 */
  mpz_set_str (N, "2987193183856651631187095115639767024853796193560927033884323499962128575016083987608438436687647187187228725354979166326080044899692965798796577991358266793284662823", 10);
  d = 5;
  ad = 1000020;
#endif

  /*    P        #ad tried    seconds   X4 (with ad=10000020):
     100000      66/133        13s       7.3e+24
     200000      31/149        11.7s     1.0e+26

     300000      569/159       310.948   8.2e+25
     400000      1688/166      1224s     7.2e+24
     500000      566/172       515s      2.3e+25
     1000000     1234/191      2232s     8.5e+23
     2000000     453/211       1637s     1.5e+24
     5000000     603/238       5374s     4.2e+23
  */

  //  P = 1000000;
  P = atoi (argv[1]);
  exp_tries = pow (log ((double) P), 2.0);
  printf ("P=%lu (%1.0f tries expected)\n", P, exp_tries);
  T = malloc (nthreads * sizeof (tab_t));
  for (i = 0; i < nthreads ; i++)
    {
      mpz_init_set (T[i]->N, N);
      T[i]->d = d;
      T[i]->P = P;
      T[i]->thread = i;
    }
  do
    {
      for (i = 0; i < nthreads ; i++)
        {
          tries ++;
          gmp_printf ("%d ad=%lu\r", tries, ad);
          fflush (stdout);
          T[i]->ad = ad;
#ifndef MAX_THREADS
          newAlgo (N, d, P, ad);
#else
          pthread_create (&tid[i], NULL, one_thread, (void *) (T+i));
#endif
          ad += 60;
        }
#ifdef MAX_THREADS
      for (i = 0 ; i < nthreads ; i++)
        pthread_join(tid[i], NULL);
#endif
    }
  while (found == 0);
  fprintf (stderr, "Tried %d ad-values (log^2(P)=%1.0f)\n", tries, exp_tries);
  for (i = 0; i < nthreads ; i++)
    mpz_clear (T[i]->N);
  free (T);
  mpz_clear (N);

  return 0;
}
