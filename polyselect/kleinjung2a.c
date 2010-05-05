/* this version implements Kleinjung's original algorithm, with a search in a
   hash table, and M = P^2 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <pthread.h>
#include "utils.h"

// #define MAX_THREADS 8

unsigned long *hash_p = NULL;
mpz_t *hash_i = NULL;
unsigned long hash_alloc = 0, hash_size = 0;

static void
match (unsigned long p1, unsigned long p2, mpz_t i, mpz_t m0, unsigned long ad,
       unsigned int d, mpz_t N)
{
  mpz_t l, mtilde, m, adm1, t, k;
  unsigned int j;
  int cmp;

  gmp_fprintf (stderr, "Found match: (%lu,%Zd) (%lu,%Zd)\n", p1, i, p2, i);

  mpz_init (l);
  mpz_init (m);
  mpz_init (t);
  mpz_init (k);
  mpz_init (adm1);
  mpz_init (mtilde);
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
  //  gmp_fprintf (stderr, "a_{d-1}=%Zd\n", adm1);
  mpz_mul (m, adm1, l);
  mpz_sub (m, mtilde, m);
  if (mpz_divisible_ui_p (m, d) == 0)
    {
      fprintf (stderr, "Error: m-a_{d-1}*l not divisible by d\n");
      exit (1);
    }
  mpz_divexact_ui (m, m, d);
  if (mpz_divisible_ui_p (m, ad) == 0)
    {
      fprintf (stderr, "Error: (m-a_{d-1}*l)/d not divisible by ad\n");
      exit (1);
    }
  mpz_divexact_ui (m, m, ad);
  gmp_fprintf (stderr, "Y1: %Zd\nY0: -%Zd\n", l, m);
  gmp_fprintf (stderr, "X%d: %lu\n", d, ad);
  mpz_pow_ui (t, m, d);
  mpz_mul_ui (t, t, ad);
  mpz_sub (t, N, t);
  gmp_fprintf (stderr, "X%d: %Zd\n", d-1, adm1);
  if (mpz_divisible_p (t, l) == 0)
    {
      fprintf (stderr, "Error: t not divisible by l\n");
      exit (1);
    }
  mpz_divexact (t, t, l);
  mpz_pow_ui (mtilde, m, d-1);
  mpz_mul (mtilde, mtilde, adm1);
  mpz_sub (t, t, mtilde);
  for (j = d - 2; j > 0; j--)
    {
      if (mpz_divisible_p (t, l) == 0)
        {
          fprintf (stderr, "Error: t not divisible by l\n");
          exit (1);
        }
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
      gmp_fprintf (stderr, "X%d: %Zd", j, adm1);
      if (j == d-2)
        fprintf (stderr, " (%.1e)", fabs (mpz_get_d (adm1)));
      fprintf (stderr, "\n");
      /* subtract adm1*m^j */
      mpz_submul (t, mtilde, adm1);
    }
  if (mpz_divisible_p (t, l) == 0)
    {
      fprintf (stderr, "Error: t not divisible by l\n");
      exit (1);
    }
  mpz_divexact (t, t, l);
  gmp_fprintf (stderr, "X0: %Zd\n", t);
  mpz_clear (l);
  mpz_clear (m);
  mpz_clear (t);
  mpz_clear (k);
  mpz_clear (adm1);
  mpz_clear (mtilde);
}

static void
hash_init (void)
{
  unsigned long j;

  hash_alloc = 1;
  hash_p = (unsigned long*) malloc (hash_alloc * sizeof (unsigned long));
  hash_i = (mpz_t*) malloc (hash_alloc * sizeof (mpz_t));
  for (j = 0; j < hash_alloc; j++)
    {
      hash_p[j] = 0;
      mpz_init (hash_i[j]);
    }
  hash_size = 0;
}

/* return 1 if a match is found, 0 otherwise */
static int
hash_add (unsigned long *hash_p, mpz_t *hash_i, unsigned long hash_alloc,
          unsigned long p, mpz_t i, mpz_t m0, unsigned long ad, unsigned int d,
          mpz_t N)
{
  unsigned long h = mpz_fdiv_ui (i, hash_alloc);
  int found = 0;

  while (hash_p[h] != 0)
    {
      if (m0 != NULL && mpz_cmp (hash_i[h], i) == 0 && hash_p[h] != p)
        {
          found = 1;
          match (hash_p[h], p, i, m0, ad, d, N);
        }
      if (++h == hash_alloc)
        h = 0;
    }
  hash_p[h] = p;
  mpz_set (hash_i[h], i);
  hash_size ++;
  return found;
}

static void
hash_clear (void)
{
  unsigned long j;

  for (j = 0; j < hash_alloc; j++)
    mpz_clear (hash_i[j]);
  free (hash_p);
  free (hash_i);
}

static void
hash_grow (void)
{
  unsigned long j, *new_hash_p, new_hash_alloc;
  mpz_t *new_hash_i;

  new_hash_alloc = 2 * hash_alloc;
  new_hash_p = (unsigned long*) malloc (new_hash_alloc *
                                        sizeof (unsigned long));
  for (j = 0; j < new_hash_alloc; j++)
    new_hash_p[j] = 0;
  new_hash_i = (mpz_t*) malloc (new_hash_alloc * sizeof (mpz_t));
  for (j = 0; j < new_hash_alloc; j++)
    mpz_init (new_hash_i[j]);
  hash_size = 0;
  for (j = 0; j < hash_alloc; j++)
    if (hash_p[j] != 0)
      hash_add (new_hash_p, new_hash_i, new_hash_alloc, hash_p[j], hash_i[j],
                NULL, 0, 0, NULL);
  hash_clear ();
  hash_p = new_hash_p;
  hash_i = new_hash_i;
  hash_alloc = new_hash_alloc;
}

/* return number of polynomials found */
static int
newAlgo (mpz_t N, unsigned long d, unsigned long P, unsigned long ad)
{
  mpz_t p, pp, pmax, *f, lambda, tmp, m0, Ntilde;
  unsigned long *r, nr, i, j, pui, nprimes = 0;
  int found = 0;

  hash_init ();
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
#if 0
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
          if (2 * hash_size + 1 >= hash_alloc)
            hash_grow ();
          found += hash_add (hash_p, hash_i, hash_alloc, pui, lambda, m0, ad, d, N);
          mpz_sub (lambda, lambda, pp);
          found += hash_add (hash_p, hash_i, hash_alloc, pui, lambda, m0, ad, d, N);
        }
    }
  // fprintf (stderr, "Number of primes: %lu\n", nprimes);
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
  hash_clear ();
  
  // fprintf (stderr, "Number of pairs: %lu (expected %1.0f)\n", hash_size, mpz_get_d (M) / (double) P / log ((double) P));

  return found;
}

typedef struct
{
  mpz_t N;
  unsigned int d;
  unsigned long P;
  unsigned long ad;
} __tab_struct;
typedef __tab_struct tab_t[1];

void*
one_thread (void* args)
{
  tab_t *tab = (tab_t*) args;
  newAlgo (tab[0]->N, tab[0]->d, tab[0]->P, tab[0]->ad);
  pthread_exit (NULL);
}

int
main (int argc, char *argv[])
{
  mpz_t N;
  unsigned int d;
  unsigned long P, ad;
  int found = 0, tries = 0, i, nthreads = 1;
  double exp_tries;
  tab_t *T;
#ifdef MAX_THREADS
  pthread_t tid[MAX_THREADS];
#endif

  if (argc >= 3 && strcmp (argv[1], "-t") == 0)
    {
      nthreads = atoi (argv[2]);
      argv += 2;
      argc -= 2;
    }

  if (argc != 2)
    {
      fprintf (stderr, "Usage: %s [-t nthreads] P\n", argv[0]);
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
  mpz_set_str (N, "1230186684530117755130494958384962720772853569595334792197322452151726400507263657518745202199786469389956474942774063845925192557326303453731548268507917026122142913461670429214311602221240479274737794080665351419597459856902143413", 10); /* RSA 768 */
  d = 6; /* degree */
  ad = 100000020;

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
          found += newAlgo (N, d, P, ad);
#else
          pthread_create (&tid[i], NULL, one_thread, (void *) (T+i));
#endif
          ad += 60;
        }
#ifdef MAX_THREADS
      for (i = 0 ; i < MAX_THREADS ; i++)
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
