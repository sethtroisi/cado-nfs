/*
   Test data
    P  = 31081938120519680804196101011964261019661412191103091971180537759
    (P - 1)/2 = 2 * Q where Q is a prime
    Q = 15540969060259840402098050505982130509830706095551545985590268879

    --run--
    ./dlpolyselect -df 3 -N 31081938120519680804196101011964261019661412191103091971180537759

    ./dlpolyselect -df 4 -dg 3 -N 191147927718986609689229466631454649812986246276667354864188503638807260703436799058776201365135161278134258296128109200046702912984568752800330221777752773957404540495707851421851

    ./dlpolyselect -df 3 -dg 2 -N 191147927718986609689229466631454649812986246276667354864188503638807260703436799058776201365135161278134258296128109200046702912984568752800330221777752773957404540495707851421851

*/

#include "cado.h"
#include "auxiliary.h"
#include "utils.h"
#include "portability.h"
#include "murphyE.h"
#include "ropt_param.h"
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <ctype.h>
#include <stdlib.h>
#include <time.h>

const double exp_rot[] = {0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 0};
#define SIZE(p) (((long *) (p))[1])
#define DATA(p) ((mp_limb_t *) (((long *) (p)) + 2))
#define LEN_QQ 11

const unsigned int QQ[LEN_QQ] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31};

mpz_poly *best_f = NULL;
double *best_alpha = NULL;
unsigned int keep = 100; /* number of polynomials f with best alpha to be kept */
unsigned int best_n = 0;

static void
save_f (mpz_t *f, unsigned int df)
{
  mpz_poly ff;
  double alpha;
  int i;

  ff->deg = df;
  ff->coeff = f;

  alpha = get_alpha (ff, ALPHA_BOUND);

  if (best_n == keep && alpha > best_alpha[best_n - 1])
    return;
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
  {
    for (i = best_n; i > 0 && alpha < best_alpha[i-1]; i--)
      {
        best_alpha[i] = best_alpha[i-1];
        mpz_poly_set (best_f[i], best_f[i-1]);
      }
    best_alpha[i] = alpha;
    mpz_poly_set (best_f[i], ff);
    if (i == 0)
      {
        printf ("best alpha %1.2f: ", alpha);
        mpz_poly_fprintf (stdout, ff);
      }
    best_n += (best_n < keep);
  }
}

/*
  Print two nonlinear poly info
*/
static void
print_nonlinear_poly_info ( mpz_t *f,
                            mpz_t *g,
                            unsigned int df,
                            unsigned int dg,
                            int format,
                            mpz_t n )
{
    unsigned int i;
    double skew, logmu[2], alpha[2], score;
    mpz_poly ff;
    ff->deg = df;
    ff->coeff = f;
    mpz_poly gg;
    gg->deg = dg;
    gg->coeff = g;
    static double best_score = DBL_MAX;
    /* the coefficients of g are O(n^(1/df)) */
    double target_score = log (mpz_get_d (n)) / (double) df;

    if (mpz_cmp_ui (g[dg], 0) == 0 || mpz_cmp_ui (g[0], 0) == 0)
      return;

    /* we use the skewness of polynomial g with large coefficients */
    skew = L2_skewness (gg, SKEWNESS_DEFAULT_PREC);
    logmu[1] = L2_lognorm (gg, skew);
    alpha[1] = get_alpha (gg, ALPHA_BOUND);
    logmu[0] = L2_lognorm (ff, skew);
    alpha[0] = get_alpha (ff, ALPHA_BOUND);

    score = logmu[1] + alpha[1] + logmu[0] + alpha[0];

    ASSERT_ALWAYS (0.5 * target_score < score && score < 1.5 * target_score);
    /* otherwise this might indicate that f is not irreducible */

    if (score >= best_score)
      return; /* only print record scores */

#ifdef HAVE_OPENMP
#pragma omp critical
#endif
    {
      best_score = score;

      if (format == 1)
	gmp_printf ("n: %Zd\n", n);
      else
	gmp_printf ("N %Zd\n", n);

      if (format == 1) {
        for (i = df + 1; i -- != 0; )
	  gmp_printf ("c%u: %Zd\n", i, f[i]);
      }
      else {
        for (i = df + 1; i -- != 0; )
	  gmp_printf ("X%u %Zd\n", i, f[i]);
      }
      if (format == 1) {
        for (i = dg + 1; i -- != 0; )
	  gmp_printf ("Y%u: %Zd\n", i, g[i]);
      }
      else {
        for (i = dg + 1; i -- != 0; )
	  gmp_printf ("Y%u %Zd\n", i, g[i]);
      }
      printf ("skew: %1.2f\n", skew);
      printf ("# f lognorm %1.2f, alpha %1.2f, E %1.2f\n",
	      logmu[0], alpha[0], logmu[0] + alpha[0]);
      printf ("# g lognorm %1.2f, alpha %1.2f, E %1.2f\n",
	      logmu[1], alpha[1], logmu[1] + alpha[1]);
      printf ("# f+g score %1.2f\n", score);
      printf ("\n");
    }
}

/* Check that f is irreducible.
   Let p be a prime such that f has a root r modulo p.
   Search by LLL a small linear combination between 1, r, ..., r^(d-1).
   If f is not irreducible, then p will be root of a factor of degree <= d-1,
   which will yield a linear dependency of the same order as the coefficients
   of f, otherwise if f is irreducible the linear dependency will be of order
   p^(1/d). */
static int
is_irreducible (mpz_poly f)
{
  mpz_t p;
  int d = f->deg;
  int dg = d - 1;
  int i, j, nr;
  mpz_t a, b, det, r, *roots;
  size_t normf;
  int ret;
  mat_Z g;

  mpz_init (p);
  mpz_poly_infinity_norm (p, f);
  normf = mpz_sizeinbase (p, 2);
#define MARGIN 16
  /* add some margin bits */
  mpz_mul_2exp (p, p, MARGIN);
  mpz_pow_ui (p, p, d);

  roots = malloc (d * sizeof (mpz_t));
  for (i = 0; i < d; i++)
    mpz_init (roots[i]);

  do {
    mpz_nextprime (p, p);
    nr = mpz_poly_roots_mpz (roots, f, p);
  } while (nr == 0);

  g.coeff = (mpz_t **) malloc ((dg + 2)*sizeof(mpz_t*));
  g.NumRows = g.NumCols = dg + 1;
  for (i = 0; i <= dg + 1; i ++) {
    g.coeff[i] = (mpz_t *) malloc ((dg + 2)*sizeof(mpz_t));
    for (j = 0; j <= dg + 1; j ++) {
      mpz_init (g.coeff[i][j]);
    }
  }

  mpz_init (det);
  mpz_init_set_ui (a, 3);
  mpz_init_set_ui (b, 4);
  mpz_init_set (r, roots[0]);
  for (i = 0; i <= dg + 1; i ++) {
    for (j = 0; j <= dg + 1; j ++) {
      mpz_set_ui (g.coeff[i][j], 0);
    }
  }

  for (i = 1;  i <= dg + 1; i++) {
    for (j = 1; j <= dg + 1; j++) {
      if (i == 1) {
        if (j == 1) {
          mpz_set (g.coeff[j][i], p);
        }
        else {
          mpz_neg (g.coeff[j][i], r);
          mpz_mul (r, r, roots[0]);
        }
      }
      else
        mpz_set_ui (g.coeff[j][i], i==j);
    }
  }

  LLL (det, g, NULL, a, b);

  for (j = 1; j <= dg + 1; j ++)
    {
      /* the coefficients of vector j are in g.coeff[j][i], 1 <= i <= dg + 1 */
      mpz_abs (a, g.coeff[j][1]);
      for (i = 2; i <= dg + 1; i++)
        if (mpz_cmpabs (g.coeff[j][i], a) > 0)
          mpz_abs (a, g.coeff[j][i]);
      /* now a = max (|g.coeff[j][i]|, 1 <= i <= dg+1) */
      if (j == 1 || mpz_cmpabs (a, b) < 0)
        mpz_set (b, a);
    }
  /* now b is the smallest infinity norm */
  if (mpz_sizeinbase (b, 2) < normf + MARGIN / 2)
    ret = 0;
  else
    ret = 1;

  for (i = 0; i <= dg + 1; i ++) {
    for (j = 0; j <= dg + 1; j ++)
      mpz_clear (g.coeff[i][j]);
    free(g.coeff[i]);
  }
  free (g.coeff);

  for (i = 0; i < d; i++)
    mpz_clear (roots[i]);
  free (roots);
  mpz_clear (det);
  mpz_clear (a);
  mpz_clear (b);
  mpz_clear (r);
  mpz_clear (p);

  return ret;
}

/*
  Generate polynomial f(x) of degree d with rank 'idx',
  with coefficients in [-bound, bound].
  The coefficient of degree d can be taken in [1, bound],
  due to the symmetry f(x) -> -f(x).
  The coefficient of degree d-1 can be taken in [0, bound],
  due to the symmetry f(x) -> f(-x).
  The coefficient of degree 0 should not be 0.
  Thus there are bound*(bound+1)*(2*bound+1)^(d-2)*(2*bound) possible values.
  Return 0 if the corresponding polynomial is not irreducible,
  or has no roots modulo n, otherwise return the number of roots
  modulo n.
*/
static int
polygen_JL_f ( mpz_t n,
               int d,
               unsigned int bound,
               mpz_t *f,
	       unsigned long idx )
{
    unsigned int i;
    unsigned long *rq;
    int nr = 0, *fint;
    mpz_t t;
    mpz_init (t);
    fint = (int *) malloc ((d + 1)*sizeof(int));
    rq = (unsigned long *) malloc ((d + 1)*sizeof(unsigned long));

    /* find irreducible polynomial f */
    {
        /* we take 1 <= f[d] <= bound */
        fint[d] = 1 + (idx % bound);
	idx = idx / bound;
        mpz_set_ui (f[d], fint[d]);
	/* we take 0 <= f[d-1] <= bound */
	fint[d-1] = idx % (bound + 1);
	idx = idx / (bound + 1);
        mpz_set_ui (f[d-1], fint[d-1]);
        for (i = d-2; i > 0; i --) {
	    fint[i] = (idx % (2 * bound + 1)) - bound;
	    idx = idx / (2 * bound + 1);
            mpz_set_si (f[i], fint[i]);
        }
	/* we take -bound <= f[0] < bound, f[0] <> 0,
	   which makes 2*bound possible values */
	ASSERT_ALWAYS(idx < 2 * bound);
	fint[0] = (idx < bound) ? idx - bound : idx - (bound - 1);
	mpz_set_si (f[0], fint[0]);

	/* since f and the reversed polynomial are equivalent, we can assume
	   |f[d]| < |f[0]| or (|f[d]| = |f[0]| and |f[d-1]| < |f[1]|) or ... */

	int ok = 1;
	for (int i = 0; 2 * i < d && ok; i++)
	  {
	    if (abs(fint[d-i]) > abs(fint[i]))
	      ok = 0;
	    else if (abs(fint[d-i]) < abs(fint[i]))
	      break;
	  }
	if (ok == 0)
	  goto end;

        /* content test (not necessary for now) */
        mpz_poly ff;
        ff->deg = d;
        ff->coeff = f;
        mpz_poly_content (t, ff);
        if (mpz_cmp_ui(t, 1) != 0 )
	  goto end; /* duplicate with f/t */

        /* irreducibility test */
        /*
        f_ZZ.SetLength(d+1);
        for (i = 0; i <= d; i ++)
            SetCoeff(f_ZZ, i, fint[i]);
        if (is_irreducible(f_ZZ) == 0)
            continue;
        */
        int test = 0;
        for (i = 0; i < LEN_QQ; i ++) {
          test = mpz_poly_roots_ulong (rq, ff, QQ[i]);
          ASSERT(test <= d);
          if (test == 0)
            break;
        }
        if (test != 0) /* f has roots for all primes in QQ[], thus is
			  probably not irreducible */
	  goto end;

        if (mpz_poly_squarefree_p (ff) == 0)
	  goto end;

        if (is_irreducible (ff) == 0)
          goto end;

	/* Note: this is not enough for degree 4 or more, since f might
	   have one factor of degree 2 and one factor of degree d-2.
	   When this is the case, those factors will pop up in LLL. */

        /* find roots mod n */
        nr = mpz_poly_roots_mpz (NULL, ff, n);
        ASSERT(nr <= d);
    }

 end:
    mpz_clear (t);
    free(fint);
    free(rq);
    return nr;
}

/*
  Generate polynomial g(x) of degree dg, given root
*/
static void
polygen_JL_g ( mpz_t N,
               int dg,
               mat_Z g,
               mpz_t root )
{
    int i, j;
    mpz_t a, b, det, r;
    mpz_init (det);
    mpz_init_set_ui(a, 3);
    mpz_init_set_ui(b, 4);
    mpz_init_set(r, root);
    for (i = 0; i <= dg + 1; i ++) {
        for (j = 0; j <= dg + 1; j ++) {
            mpz_set_ui (g.coeff[i][j], 0);
        }
    }

    for (i = 1;  i <= dg + 1; i++) {
        for (j = 1; j <= dg + 1; j++) {
            if (i == 1) {
                if (j == 1) {
                    mpz_set (g.coeff[j][i], N);
                }
                else {
                    mpz_neg (g.coeff[j][i], r);
                    mpz_mul (r, r, root);
                }
            }
            else {
                mpz_set_ui (g.coeff[j][i], i==j);
            }
        }
    }

    LLL (det, g, NULL, a, b);

    mpz_clear (det);
    mpz_clear (a);
    mpz_clear (b);
    mpz_clear (r);
}

/* JL method to generate d and d-1 polynomial.
   First pass: generate best 'keep' polynomials f with best alpha value. */
static void
polygen_JL1 ( mpz_t n,
             unsigned int df,
             unsigned int bound,
	     unsigned long idx )
{
    unsigned int i, nr;
    mpz_t *f;

    ASSERT_ALWAYS (df >= 3);
    f = (mpz_t *) malloc ((df + 1)*sizeof(mpz_t));
    for (i = 0; i <= df; i ++)
      mpz_init (f[i]);

    /* generate f of degree d with small coefficients */
    nr = polygen_JL_f (n, df, bound, f, idx);
    if (nr > 0)
      save_f (f, df);

    /* clear */
    for (i = 0; i <= df; i ++)
      mpz_clear (f[i]);
    free (f);
}

/* JL method to generate d and d-1 polynomial.
   Second pass: try best 'keep' polynomials f with best alpha value. */
static void
polygen_JL2 ( mpz_t n,
             unsigned int df,
             unsigned int dg,
	     unsigned long c )
{
    unsigned int i, j, nr, format = 1;
    mpz_t *rf;
    mat_Z g;

    ASSERT_ALWAYS (df >= 3);
    rf = (mpz_t *) malloc (df * sizeof(mpz_t));
    for (i = 0; i < df; i ++)
      mpz_init (rf[i]);
    g.coeff = (mpz_t **) malloc ((dg + 2)*sizeof(mpz_t*));
    g.NumRows = g.NumCols = dg + 1;
    for (i = 0; i <= dg + 1; i ++) {
        g.coeff[i] = (mpz_t *) malloc ((dg + 2)*sizeof(mpz_t));
        for (j = 0; j <= dg + 1; j ++) {
            mpz_init (g.coeff[i][j]);
        }
    }

    /* compute number of roots of the c-th best polynomial f */
    nr = mpz_poly_roots_mpz (rf, best_f[c], n);
    ASSERT(0 < nr && nr <= df);

    for (i = 0; i < nr; i ++) {
        /* generate g of degree dg */
        polygen_JL_g (n, dg, g, rf[i]);

        for (j = 1; j <= dg + 1; j ++) {
          print_nonlinear_poly_info (best_f[c]->coeff, &((g.coeff[j])[1]), df, dg,
                                     format, n);
        }
    }

    /* clear */
    for (i = 0; i < df; i ++)
      mpz_clear (rf[i]);
    for (i = 0; i <= dg + 1; i ++) {
        for (j = 0; j <= dg + 1; j ++)
            mpz_clear (g.coeff[i][j]);
        free(g.coeff[i]);
    }
    free (g.coeff);
    free (rf);
}

static void
usage ()
{
    fprintf (stderr, "./dlpolyselect -N xxx -df xxx -dg xxx -bound xxx [-t xxx]\n");
    exit (1);
}


int
main (int argc, char *argv[])
{
    int i;
    mpz_t N;
    unsigned int df = 0, dg = 0;
    int nthreads = 1;
    unsigned int bound = 1024; /* bound on the coefficients of f */
    unsigned long maxtries;
    mpz_init (N);

    /* printf command-line */
    printf ("#");
    for (i = 0; i < argc; i++)
        printf (" %s", argv[i]);
    printf ("\n");
    fflush (stdout);

    /* parsing */
    while (argc >= 3 && argv[1][0] == '-')
    {
        if (argc >= 3 && strcmp (argv[1], "-N") == 0) {
            mpz_set_str (N, argv[2], 10);
            argv += 2;
            argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-df") == 0) {
            df = atoi (argv[2]);
            argv += 2;
            argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-dg") == 0) {
            dg = atoi (argv[2]);
            argv += 2;
            argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-bound") == 0) {
            bound = atoi (argv[2]);
            argv += 2;
            argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-t") == 0) {
            nthreads = atoi (argv[2]);
            argv += 2;
            argc -= 2;
        }
        else {
            fprintf (stderr, "Invalid option: %s\n", argv[1]);
            usage();
            exit (1);
        }
    }

    if (mpz_cmp_ui (N, 0) <= 0) {
        fprintf (stderr, "Error, missing input number (-N option)\n");
        usage ();
    }

    if (df == 0) {
        fprintf (stderr, "Error, error degree (-df option)\n");
        usage ();
    }

    if (dg == 0 || dg >= df) {
        fprintf (stderr, "Error, missing or error degree (-dg option)\n");
        fprintf (stderr, "       only support dg < df.\n");
        usage ();
    }

    srand (time (NULL));

    double maxtries_double = (double) bound;
    maxtries_double *= (double) (bound + 1);
    maxtries_double *= pow ((double) (2 * bound + 1), (double) (df - 2));
    maxtries_double *= (double) (2 * bound);
    /* since each coefficient of f of degree 0 to df-1 is chosen randomly
       in [-bound, bound-1], we have (2*bound)^df possible values for f */
    if (maxtries_double >= (double) ULONG_MAX)
      maxtries = ULONG_MAX;
    else
      maxtries = (unsigned long) maxtries_double;

    best_f = malloc ((keep + 1) * sizeof (mpz_poly));
    best_alpha = malloc ((keep + 1) * sizeof (double));
    for (unsigned int i = 0; i <= keep; i++)
      mpz_poly_init (best_f[i], df);

#ifdef HAVE_OPENMP
    omp_set_num_threads (nthreads);
#pragma omp parallel for schedule(dynamic)
#endif
    for (unsigned long c = 0; c < maxtries; c++)
      polygen_JL1 (N, df, bound, c);

    printf ("best alpha %1.2f: ", best_alpha[0]);
    mpz_poly_fprintf (stdout, best_f[0]);
    printf ("worst alpha (%d) %1.2f: ", best_n, best_alpha[best_n-1]);
    mpz_poly_fprintf (stdout, best_f[best_n-1]);

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (unsigned long c = 0; c < best_n; c++)
      polygen_JL2 (N, df, dg, c);

    printf ("best alpha %1.2f: ", best_alpha[0]);
    mpz_poly_fprintf (stdout, best_f[0]);
    printf ("worst alpha (%d) %1.2f: ", best_n, best_alpha[best_n-1]);
    mpz_poly_fprintf (stdout, best_f[best_n-1]);

    mpz_clear (N);
    for (unsigned int i = 0; i <= keep; i++)
      mpz_poly_clear (best_f[i]);
    free (best_f);
    free (best_alpha);

    return 0;
}
