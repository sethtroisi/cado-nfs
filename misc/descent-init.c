/* Experimental program to initialize the individual logarithm ("descent")
   in DLP.

   Input:
   ./descent-init [-seed s] p z
   s is the random seed used (if not given, getpid() is used)
   p is the prime defining the DLP group
   z is the target element

   Output:
   integers e, u0, v0 such that z^e = u0/v0 mod p and u0, v0 are smooth.

   Requires libecm.{a,so} in addition to -lgmp and -lm.

   Example: with the 180-digit prime and the target z=rsa1024 from
   http://caramel.loria.fr/p180.txt we get in a few minutes:

   e=1125822966 L=91 aver_gain[189]=131.337396 gain_u=131.337396 gain_v=207.443672 S=5223565.252471 Smax=751333.655143
   u0=272156503062481694016191511460637222570435895062537206732822181869087792349037223514408066
   v0=327142099764205173410624844875387710737486249289718351091450638348510670002468548292102115
   u=1989808957508784679468391539 (p91)
   v=1169415656325906683767999889 (p90)

   Note: even with the same random seed, the results might differ between two
   different runs since the random seed does not control the ECM runs, but only
   the choice of the exponents e.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/resource.h>
#include "gmp.h"
#include "ecm.h"

typedef struct {
  unsigned long e; /* exponent of z */
  mpz_t u0, v0;    /* z^e = u0/v0 mod p */
  mpz_t u, v;      /* unfactored parts of u, v */
  int uprime, vprime;
  double B1;       /* current B1 value */
  unsigned long l; /* potential size */
} cand_t;
typedef cand_t cand[1];

unsigned long
bsize (mpz_t u, int uprime)
{
  if (uprime)
    return mpz_sizeinbase (u, 2);
  else /* if u has n bits, we have 2^(n-1) <= u < 2^n: if u = p*q with p <= q,
          then q >= 2^((n-1)/2) and q has at least floor(((n-1)/2)+1 bits */
    return (mpz_sizeinbase (u, 2) + 1) / 2;
}

void
cand_init_set (cand c, unsigned long e, mpz_t u0, mpz_t v0, mpz_t u, mpz_t v,
               int uprime, int vprime, double B1, unsigned long l)
{
  c->e = e;
  mpz_init_set (c->u0, u0);
  mpz_init_set (c->v0, v0);
  mpz_init_set (c->u, u);
  mpz_init_set (c->v, v);
  c->uprime = uprime;
  c->vprime = vprime;
  c->B1 = B1;
  c->l = l;
}

void
cand_clear (cand c)
{
  mpz_clear (c->u0);
  mpz_clear (c->v0);
  mpz_clear (c->u);
  mpz_clear (c->v);
}

void
cand_scan (cand c, ecm_params params)
{
  mpz_t f;
  unsigned long lu, lv;

  mpz_init (f);
  if (c->uprime == 0)
    {
      mpz_set_d (params->sigma, c->B1);
      mpz_set_ui (f, 1);
      ecm_factor (f, c->u, c->B1, params);
      params->B1done = ECM_DEFAULT_B1_DONE; /* issue in 6.4.x */
      if (mpz_cmp_ui (f, 1) > 0)
        {
          mpz_divexact (c->u, c->u, f);
          c->uprime = mpz_probab_prime_p (c->u, 1);
          lu = bsize (c->u, c->uprime);
          lv = bsize (c->v, c->vprime);
          c->l = (lu > lv) ? lu : lv;
        }
    }
  if (c->vprime == 0)
    {
      mpz_set_d (params->sigma, c->B1);
      mpz_set_ui (f, 1);
      ecm_factor (f, c->v, c->B1, params);
      params->B1done = ECM_DEFAULT_B1_DONE; /* issue in 6.4.x */
      if (mpz_cmp_ui (f, 1) > 0)
        {
          mpz_divexact (c->v, c->v, f);
          c->vprime = mpz_probab_prime_p (c->v, 1);
          lu = bsize (c->u, c->uprime);
          lv = bsize (c->v, c->vprime);
          c->l = (lu > lv) ? lu : lv;
        }
    }
  c->B1 = c->B1 + sqrt (c->B1);
  mpz_clear (f);
}

void
cand_print (cand c)
{
  printf ("e=%lu\n", c->e);
  gmp_printf ("u0=%Zd\n", c->u0);
  gmp_printf ("v0=%Zd\n", c->v0);
  gmp_printf ("u=%Zd (%lu bits)\n", c->u, mpz_sizeinbase (c->u, 2));
  gmp_printf ("v=%Zd (%lu bits)\n", c->v, mpz_sizeinbase (c->v, 2));
  fflush (stdout);
}

void
cand_set (cand c, cand d)
{
  c->e = d->e;
  mpz_swap (c->u0, d->u0);
  mpz_swap (c->v0, d->v0);
  mpz_swap (c->u, d->u);
  mpz_swap (c->v, d->v);
  c->uprime = d->uprime;
  c->vprime = d->vprime;
  c->B1 = d->B1;
  c->l = d->l;
}

/*****************************************************************************/

typedef struct {
  cand *list;
  unsigned long n;
} pool_t;
typedef pool_t pool[1];

void
pool_init (pool p)
{
  p->list = NULL;
  p->n = 0;
}

void
pool_clear (pool p)
{
  unsigned long i;

  for (i = 0; i < p->n; i++)
    cand_clear (p->list[i]);
  free (p->list);
}

void
pool_add (pool p, unsigned long e, mpz_t u0, mpz_t v0, mpz_t u, mpz_t v,
          int uprime, int vprime, double B1, unsigned long l)
{
  p->list = realloc (p->list, (p->n + 1) * sizeof (cand));
  cand_init_set (p->list[p->n], e, u0, v0, u, v, uprime, vprime, B1, l);
  p->n += 1;
}

/* perform an ECM curve for each candidate in the pool */
unsigned long
pool_scan (pool p, ecm_params params, unsigned long L)
{
  unsigned long i, j;

  for (i = 0; i < p->n; i++)
    cand_scan (p->list[i], params);
  for (i = 0; i < p->n; i++)
    if (p->list[i]->uprime && p->list[i]->vprime)
      {
        if (p->list[i]->l < L) /* found better solution */
          {
            cand_print (p->list[i]);
            L = p->list[i]->l;
          }
      }
  for (i = j = 0; i < p->n; i++)
    if ((p->list[i]->uprime && p->list[i]->vprime) || p->list[i]->l >= L)
      {
      }
    else
      cand_set (p->list[j++], p->list[i]);
  if (j < p->n)
    {
      p->n = j;
      p->list = realloc (p->list, p->n * sizeof (cand));
      printf ("Remains %lu candidates in pool\n", p->n);
      fflush (stdout);
    }
  return L;
}

/*****************************************************************************/

int
cputime ()
{
  struct rusage rus;

  getrusage (0, &rus);
  return rus.ru_utime.tv_sec * 1000 + rus.ru_utime.tv_usec / 1000;
}

void
HalfGcd (mpz_t a, mpz_t b, mpz_t u)
{
  mpz_t v, w, x, q, r;
  unsigned long n;

  n = (mpz_cmpabs (a, b) > 0) ? mpz_sizeinbase (a, 2) : mpz_sizeinbase (b, 2);
  mpz_set_ui (u, 1);
  mpz_init_set_ui (w, 0);
  mpz_init_set_ui (v, 0);
  mpz_init_set_ui (x, 1);
  mpz_init (q);
  mpz_init (r);
  /* invariant: a = u*a0 + v*b0 */
  while (mpz_cmpabs (a, u) > 0)
    {
      mpz_tdiv_qr (q, r, a, b);
      mpz_swap (a, b);
      mpz_swap (b, r);
      mpz_submul (u, q, w);
      mpz_swap (u, w);
      mpz_submul (v, q, x);
      mpz_swap (v, x);
    }
  mpz_clear (v);
  mpz_clear (w);
  mpz_clear (x);
  mpz_clear (q);
  mpz_clear (r);
}

#define N 2048
double aver_gain[N], gain_u[N], gain_v[N];
unsigned long nb_test[N];

int
main (int argc, char *argv[])
{
  mpz_t p, z, ze, t, u, v, f, u0, v0;
  unsigned long e, L, l, lu, lv, i, seed = 0, target = 0;
  double B1, S, Smax = 1000.0;
  int uprime, vprime;
  ecm_params params;
  pool P;

  while (argc > 2 && argv[1][0] == '-')
    {
      if (strcmp (argv[1], "-seed") == 0)
        {
          seed = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else if (strcmp (argv[1], "-target") == 0)
        {
          target = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else
        {
          fprintf (stderr, "Unknown option %s\n", argv[1]);
          exit (1);
        }
    }

  assert (argc == 3);

  for (i = 0; i < N; i++)
    {
      aver_gain[i] = 0.0;
      nb_test[i] = 0;
    }

  pool_init (P);
  mpz_init_set_str (p, argv[1], 10);
  mpz_init_set_str (z, argv[2], 10);
  gmp_printf ("p=%Zd\n", p);
  gmp_printf ("z=%Zd\n", z);
  if (seed == 0)
    seed = getpid ();
  srand48 (seed);
  printf ("Using seed %lu\n", seed);
  mpz_init (ze);
  mpz_init (t);
  mpz_init (u);
  mpz_init (v);
  mpz_init (u0);
  mpz_init (v0);
  mpz_init (f);
  L = mpz_sizeinbase (p, 2);
  while (L > target)
    {
      e = lrand48 ();
      mpz_powm_ui (u, z, e, p);
      mpz_set (t, p);
      HalfGcd (u, t, v);
      mpz_set (u0, u);
      mpz_set (v0, v);
      /* u = v*z^e mod p thus z^e = u/v mod p */
      
      /* remove factors of 2 and 3 */
      gain_u[0] = gain_v[0] = 0;
      while (mpz_divisible_ui_p (u, 2))
        {
          mpz_divexact_ui (u, u, 2);
          gain_u[0] += 1.0;
        }
      while (mpz_divisible_ui_p (u, 3))
        {
          mpz_divexact_ui (u, u, 3);
          gain_u[0] += log2 (3.0);
        }
      while (mpz_divisible_ui_p (v, 2))
        {
          mpz_divexact_ui (v, v, 2);
          gain_v[0] += 1.0;
        }
      while (mpz_divisible_ui_p (v, 3))
        {
          mpz_divexact_ui (v, v, 3);
          gain_v[0] += log2 (3.0);
        }
      B1 = 100.0;
      S = 0.0;
      unsigned int i = 0;
      uprime = mpz_probab_prime_p (u, 1);
      if (uprime && mpz_sizeinbase (u, 2) >= L)
        continue;
      vprime = mpz_probab_prime_p (v, 1);
      if (vprime && mpz_sizeinbase (v, 2) >= L)
        continue;
      ecm_init (params);
      while (S < Smax && (uprime == 0 || vprime == 0))
        {
          i = i + 1;
	  assert (i < N);
          gain_u[i] = gain_u[i-1];
          gain_v[i] = gain_v[i-1];
          if (uprime == 0)
            {
              mpz_set_d (params->sigma, B1);
              mpz_set_ui (f, 1);
              ecm_factor (f, u, B1, params);
              params->B1done = ECM_DEFAULT_B1_DONE; /* fix issue in 6.4.x */
	      gain_u[i] += log2 (mpz_get_d (f));
	      aver_gain[i] = aver_gain[i] * nb_test[i] + gain_u[i];
	      nb_test[i] += 1;
	      aver_gain[i] /= nb_test[i];
              if (mpz_cmp_ui (f, 1) > 0)
                {
                  mpz_divexact (u, u, f);
                  uprime = mpz_probab_prime_p (u, 1);
                  if (uprime && mpz_sizeinbase (u, 2) >= L)
                    break;
                }
            }

          if (vprime == 0)
            {
              mpz_set_d (params->sigma, B1);
              mpz_set_ui (f, 1);
              ecm_factor (f, v, B1, params);
              params->B1done = ECM_DEFAULT_B1_DONE; /* fix issue in 6.4.x */
	      gain_v[i] += log2 (mpz_get_d (f));
	      aver_gain[i] = aver_gain[i] * nb_test[i] + gain_v[i];
	      nb_test[i] += 1;
	      aver_gain[i] /= nb_test[i];
              if (mpz_cmp_ui (f, 1) > 0)
                {
                  mpz_divexact (v, v, f);
                  vprime = mpz_probab_prime_p (v, 1);
                  if (vprime && mpz_sizeinbase (v, 2) >= L)
                    break;
                }
            }

          if (gain_u[i] < aver_gain[i] && gain_v[i] < aver_gain[i])
            break;
          B1 += sqrt (B1);
          S += B1;
        }
      lu = bsize (u, uprime);
      lv = bsize (v, vprime);
      l = (lu > lv) ? lu : lv;
      if (l < L) /* potential better relation */
        {
          pool_add (P, e, u0, v0, u, v, uprime, vprime, B1, l);
          L = pool_scan (P, params, L);
        }
      Smax += sqrt (Smax);
    }
  mpz_clear (p);
  mpz_clear (z);
  mpz_clear (ze);
  mpz_clear (t);
  mpz_clear (u);
  mpz_clear (v);
  mpz_clear (u0);
  mpz_clear (v0);
  mpz_clear (f);
  ecm_clear (params);
  pool_clear (P);
}
