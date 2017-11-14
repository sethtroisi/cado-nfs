/* This implements Algorithm 1 from "Root Optimization of Polynomials in the
   Number Field Sieve" by Shi Bai, Richard P. Brent and Emmanuel Thomé,
   Mathematics of Computation, 2015.
   More precisely:
   * when the macro ORIGINAL is defined, it implements the original version
     described in the above reference;
   * when ORIGINAL is not defined (which is the default), it implements a
     variant which gives better results. Namely, when p^k is the largest
     power of p < B, the contribution is multiplied by p/(p-1) to take into
     account lifted roots mod p^(k+1), p^(k+2), ...
     On the RSA-768 polynomial, with B=V=W=200 and ORIGINAL defined, we get a
     maximal difference of 0.55 between the affine alpha-value and the
     computed estimation, and an average difference of 0.067, for all
     polynomials of the [-200,200]^2 grid.
     With ORIGINAL undefined, we get a maximal difference of 0.42 only,
     and an average of 0.0045 only.
*/

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "portability.h"
#include "utils.h"
#include "auxiliary.h" /* for common routines with polyselect.c */
#include "area.h"
#include "murphyE.h"
#include "size_optimization.h"

/* define ORIGINAL if you want the original algorithm from the paper */
// #define ORIGINAL

// #define TRACE_V 7
// #define TRACE_W 3

/* global variables */
int verbose = 0;                /* verbosity level */
long *Primes, nprimes;          /* primes less than B */
long *Q;                        /* largest p^k < B */
long bestu = 0, bestv = 0;      /* current best rotation */
double best_alpha = DBL_MAX;    /* alpha of best rotation */
double best_E = 0;              /* E of best rotation (with -E) */
mpz_t bestw;                    /* current best rotation in w */
int debug = 0;                  /* to print (max,avg) discrepancy with respect
                                   to real alpha */
double tot_pols = 0;            /* number of sieved polynomial */
double sieve_time = 0;
double max_discrepancy = 0;
double sum_discrepancy = 0;
double num_discrepancy = 0;
long u0 = 0, v0 = 0, w0 = 0;    /* initial translation */
int optimizeE = 0;              /* if not zero, optimize E instead of alpha */
double margin_alpha = 0.0;      /* margin with -E */
#define MARGIN_ALPHA 1.0

static void
usage_and_die (char *argv0)
{
  fprintf (stderr, "usage: %s [-area a] [-I n] [-Bf b] [-Bg c] [-margin x] [-v] [-debug] [-sopt] [-V vmax] [-W wmax] [-B bbb] poly\n", argv0);
  fprintf (stderr, "  poly: filename of polynomial\n");
  fprintf (stderr, "  -area a      area parameter for Murphy-E computation (default %.2e)\n", AREA);
  fprintf (stderr, "  -I nnn       I-value for Murphy-E computation (overrides area)\n");
  fprintf (stderr, "  -Bf b        Bf bound for Murphy-E computation (default %.2e)\n", BOUND_F);
  fprintf (stderr, "  -Bg c        Bg bound for Murphy-E computation (default %.2e)\n", BOUND_G);
  fprintf (stderr, "  -margin x    allows a lognorm increase of x (default %.2f)\n", NORM_MARGIN);
  fprintf (stderr, "  -v           verbose toggle\n");
  fprintf (stderr, "  -debug       prints avg/max discrepancy with true affine alpha\n");
  fprintf (stderr, "  -sopt        first size-optimize the given polynomial\n");
  fprintf (stderr, "  -B nnn       parameter for alpha computation (default %d)\n", ALPHA_BOUND);
  fprintf (stderr, "  -E           optimize E instead of alpha\n");
  exit (1);
}

static unsigned long
initPrimes (unsigned long B)
{
  unsigned long nprimes = 0, p, q, l;

  Primes = malloc (B * sizeof (unsigned long));
  ASSERT_ALWAYS(Primes != NULL);
  for (p = 2; p < B; p += 1 + (p > 2))
    if (ulong_isprime (p))
      Primes[nprimes++] = p;
  Primes = realloc (Primes, nprimes * sizeof (unsigned long));
  ASSERT_ALWAYS(Primes != NULL);

  /* compute prime powers */
  Q = malloc (nprimes * sizeof (unsigned long));
  ASSERT_ALWAYS(Q != NULL);
  for (l = 0; l < nprimes; l++)
    {
      p = Primes[l];
      for (q = p; q * p < B; q *= p);
      Q[l] = q;
    }

  return nprimes;
}

#if 0
/* put in roots[0], roots[1], ... the roots of f + g * w = 0 mod q,
   and return the number of roots. Naive implementation. */
static unsigned long
get_roots (unsigned long *roots, unsigned long f, unsigned long g,
           unsigned long q)
{
  unsigned long nroots = 0;

  for (unsigned long w = 0; w < q; w++)
    if (((f + g * w) % q) == 0)
      roots[nroots++] = w;
  return nroots;
}
#else
/* put in roots[0], roots[1], ... the roots of f + g * w = 0 mod q,
   and return the number of roots */
static unsigned long
get_roots (unsigned long *roots, unsigned long f, unsigned long g,
           unsigned long q)
{
  unsigned long nroots = 0;
  unsigned long h = gcd_ul (g, q);
  if (h == 1) /* only one root, namely -f/g mod q */
    {
      unsigned long invg = invert_ul (g, q);
      roots[0] = ((q - f) * invg) % q;
      nroots = 1;
    }
  else if ((f % h) != 0)
    nroots = 0;
  else
    {
      f /= h;
      g /= h;
      q /= h;
      unsigned long invg = invert_ul (g, q);
      roots[0] = ((q - f) * invg) % q;
      for (unsigned long j = 1; j < h; j++)
        roots[j] = roots[j-1] + q;
      nroots = h;
    }
  return nroots;
}
#endif

/* rotation for a fixed value of v */
static void
rotate_v (cado_poly_srcptr poly0, long v, long B,
          double maxlognorm, double Bf, double Bg, double area, long u)
{
  long w, wmin, wmax;
  cado_poly poly;
  long l;
  mpz_t g1, g0;
  mpz_t wminz, wmaxz;

  mpz_init (wminz);
  mpz_init (wmaxz);

  /* first make a local copy of the original polynomial */
  cado_poly_init (poly);
  cado_poly_set (poly, poly0);

  mpz_init_set (g1, poly->pols[RAT_SIDE]->coeff[1]);
  mpz_init_set (g0, poly->pols[RAT_SIDE]->coeff[0]);

  /* compute f + (v*x)*g */
  rotate_aux (poly->pols[ALG_SIDE]->coeff, g1, g0, 0, v, 1);

  double lognorm = L2_lognorm (poly->pols[ALG_SIDE], poly->skew);
  rotation_space r;
  expected_growth (&r, poly->pols[ALG_SIDE], poly->pols[RAT_SIDE], 0,
                   maxlognorm - lognorm);
  mpz_set_d (wminz, r.kmin);
  mpz_set_d (wmaxz, r.kmax);

  if (mpz_cmp (wminz, wmaxz) >= 0)
    goto end;

  tot_pols += mpz_get_d (wmaxz) - mpz_get_d (wminz);

  if (verbose)
    {
      gmp_printf ("v=%ld: wmin=%Zd wmax=%Zd\n", v, wminz, wmaxz);
      fflush (stdout);
    }

  /* compute the expected value sum(log(p)/(p-1), p < B) */
  double expected = 0.0;
  for (l = 0; l < nprimes; l++)
    {
      long p = Primes[l];
      expected += log ((double) p) / (double) (p - 1);
    }

  ASSERT_ALWAYS (mpz_fits_slong_p (wmaxz));
  ASSERT_ALWAYS (mpz_fits_slong_p (wminz));

  wmax = mpz_get_si (wmaxz);
  wmin = mpz_get_si (wminz);

  long len = wmax - wmin + 1;
  float *A = malloc (len * sizeof (float));
  /* A[w-wmin] corresponds to w */
  ASSERT_ALWAYS(A != NULL);
  for (long j = 0; j < len; j++)
    A[j] = expected;

#if defined(TRACE_V) && defined(TRACE_W)
  ASSERT_ALWAYS(wmin <= TRACE_W && TRACE_W <= wmax);
  if (v == TRACE_V)
    printf ("initialized A[%d] to %f\n", TRACE_W, A[TRACE_W - wmin]);
#endif

  mpz_t ump;
  mpz_init (ump);
  unsigned long *roots = malloc (B * sizeof (long));
  ASSERT_ALWAYS(roots != NULL);
  float *L;
  double nu;
  L = malloc (B * sizeof (float));
  for (l = 0; l < nprimes; l++)
    {
      long p = Primes[l], s, t, q;
      double logp = log ((double) p);
      memset (L, 0, B * sizeof (float));
      for (q = p; q <= Q[l]; q *= p)
        {
          /* the contribution is log(p)/p^(k-1)/(p+1) when the exponent k
             is not the largest one, and log(p)/p^(k-1)/(p+1)*p/(p-1) for the
             largest exponent k */
          nu = logp / (double) q * (double) p / (double) (p + 1);
#ifndef ORIGINAL
          if (q == Q[l])
            nu *= (double) p / (double) (p - 1);
#endif
          for (long x = 0; x < q; x++)
            {
              /* compute f(x) and g(x) mod p^k, where q = p^k */
              unsigned long fx, gx;
              mpz_poly_eval_ui (ump, poly->pols[ALG_SIDE], x);
              fx = mpz_fdiv_ui (ump, q);
              mpz_poly_eval_ui (ump, poly->pols[RAT_SIDE], x);
              gx = mpz_fdiv_ui (ump, q);
              /* search roots w of fx + w*gx = 0 mod q */
              unsigned long nroots = get_roots (roots, fx, gx, q);
              for (unsigned long i = 0; i < nroots; i++)
                {
                  long w = roots[i];
                  /* update for w+t*q */
                  for (t = 0; t < Q[l] / q; t++)
                    L[w + t * q] += nu;
                }
            }
        }

      /* now perform the sieve */
      sieve_time -= seconds ();
      q = Q[l];
      for (w = 0; w < q; w++)
        {
          nu = L[w];
          if (nu == 0.0)
            continue;
          /* compute s = w+t*q-wmin such that s - q < 0 <= s, i.e.,
             (t-1)*q < wmin-w <= t*q: t = ceil((wmin-w)/q) */
          if (wmin - w < 0)
            t = (wmin - w) / q;
          else
            t = (wmin - w + q - 1) / q;
          s = w + t * q - wmin;
          while (s < len)
            {
#if defined(TRACE_V) && defined(TRACE_W)
              if (v == TRACE_V && TRACE_W - wmin == s)
                printf ("q=%ld: update A[%d] from %f to %f\n",
                        q, TRACE_W, A[s], A[s] - nu);
#endif
              A[s] -= nu;
              s += q;
            }
        }
      sieve_time += seconds ();
    }
  free (L);
  free (roots);
  mpz_clear (ump);

#if defined(TRACE_V) && defined(TRACE_W)
  if (v == TRACE_V)
    printf ("A[%d] = %f\n", TRACE_W, A[TRACE_W - wmin]);
#endif


  /* check for the smallest A[s] */
  for (w = wmin; w <= wmax; w++)
    {
      /* print alpha and E of original polynomial */
      if (u == -u0 && v == -v0 && w == -w0)
        {
          rotate_aux (poly->pols[ALG_SIDE]->coeff, g1, g0, 0, w, 0);
          poly->skew = L2_skewness (poly->pols[ALG_SIDE],
                                    SKEWNESS_DEFAULT_PREC);
          double E = MurphyE (poly, Bf, Bg, area, MURPHY_K);
          gmp_printf ("u=%ld v=%ld w=%ld est_alpha_aff=%f E=%.2e [original]\n",
                      u, v, w, A[w - wmin], E);
          /* restore the original polynomial (w=0) */
          rotate_aux (poly->pols[ALG_SIDE]->coeff, g1, g0, w, 0, 0);
        }
      if (A[w - wmin] < best_alpha + margin_alpha)
        {
          /* compute E */
          rotate_aux (poly->pols[ALG_SIDE]->coeff, g1, g0, 0, w, 0);
          poly->skew = L2_skewness (poly->pols[ALG_SIDE],
                                    SKEWNESS_DEFAULT_PREC);
          double E = MurphyE (poly, Bf, Bg, area, MURPHY_K);
          /* restore the original polynomial (w=0) */
          rotate_aux (poly->pols[ALG_SIDE]->coeff, g1, g0, w, 0, 0);

          if (optimizeE == 0 || (optimizeE == 1 && E > best_E))
            {
              bestu = u;
              bestv = v;
              mpz_set_si (bestw, w);
              best_alpha = (double) A[w - wmin];
              best_E = E;
              gmp_printf ("u=%ld v=%ld w=%Zd est_alpha_aff=%f E=%.2e\n",
                          u, v, bestw, best_alpha, E);
            }
        }
    }

  if (debug)
    {
      long w0 = 0;
      for (w = wmin; w <= wmax; w++)
        {
          double alpha, disc;
          rotate_aux (poly->pols[ALG_SIDE]->coeff, g1, g0, w0, w, 0);
          w0 = w;
          alpha = get_biased_alpha_affine (poly->pols[ALG_SIDE], B);
          disc = fabs (alpha - A[w - wmin]);
          if (disc > max_discrepancy)
            {
              max_discrepancy = disc;
              // printf ("u=%ld v=%ld w=%ld alpha=%f est=%f error=%f\n", u, v, w, alpha, A[w - wmin], disc);
            }
          sum_discrepancy += disc;
          num_discrepancy += 1;
        }
    }

  free (A);

 end:
  cado_poly_clear (poly);
  mpz_clear (g1);
  mpz_clear (g0);
  mpz_clear (wminz);
  mpz_clear (wmaxz);
}

static void
rotate (cado_poly poly, long B, double maxlognorm, double Bf, double Bg,
        double area, long u)
{
  /* determine range [vmin,vmax] */
  double lognorm = L2_lognorm (poly->pols[ALG_SIDE], poly->skew);
  rotation_space r;
  expected_growth (&r, poly->pols[ALG_SIDE], poly->pols[RAT_SIDE], 1,
                   maxlognorm - lognorm);
  long vmin = (r.kmin < (double) LONG_MIN) ? LONG_MIN : r.kmin;
  long vmax = (r.kmax > (double) LONG_MAX) ? LONG_MAX : r.kmax;
  if (verbose)
    printf ("u=%ld: vmin=%ld vmax=%ld\n", u, vmin, vmax);

  for (long v = vmin; v <= vmax; v ++)
    rotate_v (poly, v, B, maxlognorm, Bf, Bg, area, u);
}

/* don't modify poly, which is the size-optimized polynomial
   (poly0 is the initial polynomial) */
static void
print_transformation (cado_poly_ptr poly0, cado_poly_srcptr poly)
{
  mpz_t k;
  int d = poly0->pols[ALG_SIDE]->deg;

  mpz_init (k);
  /* first compute the translation k: g(x+k) = g1*x + g1*k + g0 */
  mpz_sub (k, poly->pols[RAT_SIDE]->coeff[0], poly0->pols[RAT_SIDE]->coeff[0]);
  ASSERT_ALWAYS(mpz_divisible_p (k, poly0->pols[RAT_SIDE]->coeff[1]));
  mpz_divexact (k, k, poly0->pols[RAT_SIDE]->coeff[1]);
  gmp_printf ("translation %Zd, ", k);
  do_translate_z (poly0->pols[ALG_SIDE], poly0->pols[RAT_SIDE]->coeff, k);
  /* size_optimization might multiply f0 by some integer t */
  ASSERT_ALWAYS(mpz_divisible_p (poly->pols[ALG_SIDE]->coeff[d],
				 poly0->pols[ALG_SIDE]->coeff[d]));
  mpz_divexact (k, poly->pols[ALG_SIDE]->coeff[d],
		poly0->pols[ALG_SIDE]->coeff[d]);
  if (mpz_cmp_ui (k, 1) != 0)
    {
      gmp_printf ("multiplier %Zd, ", k);
      mpz_poly_mul_mpz (poly0->pols[ALG_SIDE], poly0->pols[ALG_SIDE], k);
    }
  /* now compute rotation by x^2 */
  mpz_sub (k, poly->pols[ALG_SIDE]->coeff[3], poly0->pols[ALG_SIDE]->coeff[3]);
  ASSERT_ALWAYS(mpz_divisible_p (k, poly0->pols[RAT_SIDE]->coeff[1]));
  mpz_divexact (k, k, poly0->pols[RAT_SIDE]->coeff[1]);
  gmp_printf ("rotation [%Zd,", k);
  assert (mpz_fits_slong_p (k));
  u0 = mpz_get_si (k);
  rotate_auxg_z (poly0->pols[ALG_SIDE]->coeff, poly0->pols[RAT_SIDE]->coeff[1],
                 poly0->pols[RAT_SIDE]->coeff[0], k, 2);
  mpz_sub (k, poly->pols[ALG_SIDE]->coeff[2], poly0->pols[ALG_SIDE]->coeff[2]);
  ASSERT_ALWAYS(mpz_divisible_p (k, poly0->pols[RAT_SIDE]->coeff[1]));
  mpz_divexact (k, k, poly0->pols[RAT_SIDE]->coeff[1]);
  gmp_printf ("%Zd,", k);
  assert (mpz_fits_slong_p (k));
  v0 = mpz_get_si (k);
  rotate_auxg_z (poly0->pols[ALG_SIDE]->coeff, poly0->pols[RAT_SIDE]->coeff[1],
                 poly0->pols[RAT_SIDE]->coeff[0], k, 1);
  mpz_sub (k, poly->pols[ALG_SIDE]->coeff[1], poly0->pols[ALG_SIDE]->coeff[1]);
  ASSERT_ALWAYS(mpz_divisible_p (k, poly0->pols[RAT_SIDE]->coeff[1]));
  mpz_divexact (k, k, poly0->pols[RAT_SIDE]->coeff[1]);
  gmp_printf ("%Zd]\n", k);
  assert (mpz_fits_slong_p (k));
  w0 = mpz_get_si (k);
  rotate_auxg_z (poly0->pols[ALG_SIDE]->coeff, poly0->pols[RAT_SIDE]->coeff[1],
                 poly0->pols[RAT_SIDE]->coeff[0], k, 0);
  ASSERT_ALWAYS(mpz_cmp (poly0->pols[ALG_SIDE]->coeff[0],
                         poly->pols[ALG_SIDE]->coeff[0]) == 0);
  mpz_clear (k);
}

int
main (int argc, char **argv)
{
    int argc0 = argc;
    char **argv0 = argv;
    cado_poly poly;
    int I = 0;
    double margin = NORM_MARGIN;
    long umin, umax;
    int sopt = 0;
    long B = ALPHA_BOUND;
    double time = seconds ();

    while (argc >= 2 && argv[1][0] == '-')
      {
        if (strcmp (argv[1], "-area") == 0)
          {
            area = atof (argv [2]);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-I") == 0)
          {
            I = atoi (argv [2]);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-Bf") == 0)
          {
            bound_f = atof (argv [2]);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-Bg") == 0)
          {
            bound_g = atof (argv [2]);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-margin") == 0)
          {
            margin = atof (argv [2]);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-v") == 0)
          {
            verbose ++;
            argv ++;
            argc --;
          }
        else if (strcmp (argv[1], "-E") == 0)
          {
            optimizeE = 1;
            margin_alpha = MARGIN_ALPHA;
            argv ++;
            argc --;
          }
        else if (strcmp (argv[1], "-debug") == 0)
          {
            debug = 1;
            argv ++;
            argc --;
          }
        else if (strcmp (argv[1], "-sopt") == 0)
          {
            sopt = 1;
            argv ++;
            argc --;
          }
        else if (strcmp (argv[1], "-B") == 0)
          {
            B = strtol (argv [2], NULL, 10);
            argv += 2;
            argc -= 2;
          }
        else
          break;
      }
    if (argc != 2)
        usage_and_die (argv[0]);

    if (I != 0)
      area = bound_f * pow (2.0, (double) (2 * I - 1));

    cado_poly_init (poly);
    if (!cado_poly_read (poly, argv[1]))
      {
        fprintf (stderr, "Problem when reading file %s\n", argv[1]);
        usage_and_die (argv[0]);
      }

    if (poly->skew == 0.0)
      poly->skew = L2_combined_skewness2 (poly->pols[0], poly->pols[1],
					  SKEWNESS_DEFAULT_PREC);

    /* if -sopt, size-optimize */
    if (sopt)
      {
        cado_poly c;
        cado_poly_init (c);
        cado_poly_set (c, poly);
        size_optimization (c->pols[ALG_SIDE], c->pols[RAT_SIDE],
                           poly->pols[ALG_SIDE], poly->pols[RAT_SIDE],
                           SOPT_DEFAULT_EFFORT, verbose);
        printf ("initial polynomial:\n");
        cado_poly_fprintf (stdout, poly, "");
        print_transformation (poly, c);
        cado_poly_set (poly, c);
        printf ("size-optimized polynomial:\n");
        cado_poly_fprintf (stdout, poly, "");
        cado_poly_clear (c);
      }

    nprimes = initPrimes (B);

    /* determine range [umin,umax] */
    rotation_space r;
    expected_growth (&r, poly->pols[ALG_SIDE], poly->pols[RAT_SIDE], 2,
                     margin);
    umin = (r.kmin < (double) LONG_MIN) ? LONG_MIN : r.kmin;
    umax = (r.kmax > (double) LONG_MAX) ? LONG_MAX : r.kmax;
    if (verbose)
      printf ("umin=%ld umax=%ld\n", umin, umax);

    double maxlognorm = L2_lognorm (poly->pols[ALG_SIDE], poly->skew) + margin;

    mpz_init (bestw);

    long u0 = 0;
    for (long u = umin; u <= umax; u++)
      {
        rotate_aux (poly->pols[ALG_SIDE]->coeff,
                    poly->pols[RAT_SIDE]->coeff[1],
                    poly->pols[RAT_SIDE]->coeff[0], u0, u, 2);
        u0 = u;

        rotate (poly, B, maxlognorm, bound_f, bound_g, area, u);
      }

    if (debug)
      printf ("discrepancy: avg=%f, max=%g\n",
              sum_discrepancy / num_discrepancy, max_discrepancy);
      
    /* restore original polynomial */
    rotate_aux (poly->pols[ALG_SIDE]->coeff, poly->pols[RAT_SIDE]->coeff[1],
                poly->pols[RAT_SIDE]->coeff[0], u0, 0, 2);

    /* perform the best rotation */
    gmp_printf ("best rotation: u=%ld v=%ld w=%Zd alpha=%1.2f\n",
            bestu, bestv, bestw, best_alpha);

    /* perform the best rotation */
    rotate_aux (poly->pols[ALG_SIDE]->coeff, poly->pols[RAT_SIDE]->coeff[1],
                poly->pols[RAT_SIDE]->coeff[0], 0, bestu, 2);
    rotate_aux (poly->pols[ALG_SIDE]->coeff, poly->pols[RAT_SIDE]->coeff[1],
                poly->pols[RAT_SIDE]->coeff[0], 0, bestv, 1);
    rotate_auxg_z (poly->pols[ALG_SIDE]->coeff, poly->pols[RAT_SIDE]->coeff[1],
                poly->pols[RAT_SIDE]->coeff[0], bestw, 0);

    /* recompute the skewness of the best polynomial */
    poly->skew = L2_combined_skewness2 (poly->pols[0], poly->pols[1],
                                        SKEWNESS_DEFAULT_PREC);

    print_cadopoly_extra (stdout, poly, argc0, argv0, 0);

    time = seconds () - time;
    printf ("Sieved %.2e polynomials in %.2f seconds (%.2es/p)\n",
            tot_pols, time, time / tot_pols);
    printf ("(sieve time %.2f = %.2f%%)\n", sieve_time,
            100.0 * sieve_time / time);

    free (Primes);
    free (Q);
    cado_poly_clear (poly);
    mpz_clear (bestw);

    return 0;
}
