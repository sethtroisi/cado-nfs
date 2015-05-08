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
  unsigned long e, L, l, lu, lv, i, seed = 0;
  double B1, S, Smax = 1000.0;
  int uprime, vprime;
  ecm_params params;

  if (argc > 2 && strcmp (argv[1], "-seed") == 0)
    {
      seed = strtoul (argv[2], NULL, 10);
      argc -= 2;
      argv += 2;
    }

  assert (argc == 3);

  for (i = 0; i < N; i++)
    {
      aver_gain[i] = 0.0;
      nb_test[i] = 0;
    }

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
  while (1)
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
      if (uprime && mpz_sizeinbase (u, 2) > L)
        continue;
      vprime = mpz_probab_prime_p (v, 1);
      if (vprime && mpz_sizeinbase (v, 2) > L)
        continue;
      while (S < Smax && (uprime == 0 || vprime == 0))
        {
          i = i + 1;
	  assert (i < N);
          gain_u[i] = gain_u[i-1];
          gain_v[i] = gain_v[i-1];
          if (uprime == 0)
            {
              ecm_init (params);
              mpz_set_d (params->sigma, B1);
              mpz_set_ui (f, 1);
              ecm_factor (f, u, B1, params);
              ecm_clear (params);
	      gain_u[i] += log2 (mpz_get_d (f));
	      aver_gain[i] = aver_gain[i] * nb_test[i] + gain_u[i];
	      nb_test[i] += 1;
	      aver_gain[i] /= nb_test[i];
              if (mpz_cmp_ui (f, 1) > 0)
                {
                  mpz_divexact (u, u, f);
                  uprime = mpz_probab_prime_p (u, 1);
                  if (uprime && mpz_sizeinbase (u, 2) > L)
                    break;
                }
            }

          if (vprime == 0)
            {
              ecm_init (params);
              mpz_set_d (params->sigma, B1);
              mpz_set_ui (f, 1);
              ecm_factor (f, v, B1, params);
              ecm_clear (params);
	      gain_v[i] += log2 (mpz_get_d (f));
	      aver_gain[i] = aver_gain[i] * nb_test[i] + gain_v[i];
	      nb_test[i] += 1;
	      aver_gain[i] /= nb_test[i];
              if (mpz_cmp_ui (f, 1) > 0)
                {
                  mpz_divexact (v, v, f);
                  vprime = mpz_probab_prime_p (v, 1);
                  if (vprime && mpz_sizeinbase (v, 2) > L)
                    break;
                }
            }

          if (gain_u[i] < aver_gain[i] && gain_v[i] < aver_gain[i])
            break;
          B1 += sqrt (B1);
          S += B1;
        }
      lu = (uprime) ? mpz_sizeinbase (u, 2) : mpz_sizeinbase (u, 2) / 2;
      lv = (vprime) ? mpz_sizeinbase (v, 2) : mpz_sizeinbase (v, 2) / 2;
      l = (lu > lv) ? lu : lv;
      if (l < L) /* potential better relation */
        {
          int s = cputime ();
          printf ("focus on u=%c%lu and v=%c%lu\n",
                  (uprime) ? 'p' : 'c', mpz_sizeinbase (u, 2),
                  (vprime) ? 'p' : 'c', mpz_sizeinbase (v, 2));
          while (S < 10 * Smax && (uprime == 0 || vprime == 0))
            {
              if (uprime == 0)
                {
                  ecm_init (params);
                  mpz_set_d (params->sigma, B1);
                  mpz_set_ui (f, 1);
                  ecm_factor (f, u, B1, params);
                  ecm_clear (params);
                  if (mpz_cmp_ui (f, 1) > 0)
                    {
                      gmp_printf ("found factor of u with B1=%.0f: %Zd\n",
                                  B1, f);
                      mpz_divexact (u, u, f);
                      uprime = mpz_probab_prime_p (u, 1);
                      if (uprime && mpz_sizeinbase (u, 2) > L)
                        break;
                    }
                }
              if (vprime == 0)
                {
                  ecm_init (params);
                  mpz_set_d (params->sigma, B1);
                  mpz_set_ui (f, 1);
                  ecm_factor (f, v, B1, params);
                  ecm_clear (params);
                  if (mpz_cmp_ui (f, 1) > 0)
                    {
                      gmp_printf ("found factor of v with B1=%.0f: %Zd\n",
                                  B1, f);
                      mpz_divexact (v, v, f);
                      vprime = mpz_probab_prime_p (v, 1);
                      if (vprime && mpz_sizeinbase (v, 2) > L)
                        break;
                    }
                }
              B1 += sqrt (B1);
              S += B1;
            }
          lu = (uprime) ? mpz_sizeinbase (u, 2) : mpz_sizeinbase (u, 2) / 2;
          lv = (vprime) ? mpz_sizeinbase (v, 2) : mpz_sizeinbase (v, 2) / 2;
          l = (lu > lv) ? lu : lv;
          printf ("focus took %dms: u=%c%lu and v=%c%lu\n",
                  cputime () - s,
                  (uprime) ? 'p' : 'c', mpz_sizeinbase (u, 2),
                  (vprime) ? 'p' : 'c', mpz_sizeinbase (v, 2));
          if (l < L)
            {
              L = l;
              printf ("e=%lu L=%lu aver_gain[%d]=%f gain_u=%f gain_v=%f S=%f Smax=%f\n",
                      e, l, i, aver_gain[i], gain_u[i], gain_v[i], S, Smax);
              gmp_printf ("u0=%Zd\n", u0);
              gmp_printf ("v0=%Zd\n", v0);
              gmp_printf ("u=%Zd (%c%lu)\n", u, (uprime) ? 'p' : 'c',
                          mpz_sizeinbase (u, 2));
              gmp_printf ("v=%Zd (%c%lu)\n", v, (vprime) ? 'p' : 'c',
                          mpz_sizeinbase (v, 2));
              fflush (stdout);
            }
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
}
