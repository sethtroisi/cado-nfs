/* dlpolyselect - discrete logarithm polynomial selection via Joux-Lercier's
   algorithm

Example with Apache 512-bit key (Table 1 of
https://weakdh.org/imperfect-forward-secrecy-ccs15.pdf):

# we first search a degree (3,2) pair
$ ./dlpolyselect -df 3 -dg 2 -N 8372421755538377327377912526045445423027732035562313241965800453667849685158691589507936013805295187219621475007123900107532269487803598942841993804845107 -bound 8
...
c3: 4
c2: 3
c1: -6
c0: -8
Y2: 430506467754036192750431535685941771886700493958222
Y1: 687078639517265884307374391507473735402330150842399
Y0: -533201855944366475689054404823391134825953676216635
skew: 1.23
# f lognorm 1.28, alpha 0.62, score 1.90
# g lognorm 116.59, alpha -2.18, score 114.41
# f+g score 116.31

# then we compare with a degree (4,3) pair
$ ./dlpolyselect -df 4 -dg 3 -N 8372421755538377327377912526045445423027732035562313241965800453667849685158691589507936013805295187219621475007123900107532269487803598942841993804845107 -bound 8
...
c4: 6
c3: 1
c2: -8
c1: -3
c0: 6
Y3: -64071264306884991611859009153886700616
Y2: -115884379190374676852348454783130883186
Y1: 382823080720299801084267253734861739469
Y0: -14529492984288436819691699253895818591
skew: 1.05
# f lognorm 1.15, alpha -0.04, score 1.12
# g lognorm 87.27, alpha -4.07, score 83.20
# f+g score 84.32

To recover the polynomial pair used for the DLP768 record
(https://listserv.nodak.edu/cgi-bin/wa.exe?A2=NMBRTHRY;a0c66b63.1606):

# -modr r -modm m enables one to only consider polynomials f
# with index r mod m in Stage 1 (default is r=0 and m=1).
# This is useful to distribute the work among several machines.
$ ./dlpolyselect -N 1219344858334286932696341909195796109526657386154251328029273656175766870980306505584577389125860826715201547225794072935883258868036433287217994721542199148182841505800433148410869683590659346847659519108393837414567892730579162319 -df 4 -dg 3 -bound 140 -keep 1000 -modr 66234274 -modm 100000000
...
c4: 55
c3: 5
c2: -86
c1: 34
c0: -140
Y3: 277260730400349522890422618473498148528706115003337935150
Y2: 217583293626947899787577441128333027617541095004734736415
Y1: -1937981312833038778565617469829395544065255938015920309679
Y0: -370863403886416141150505523919527677231932618184100095924
skew: 1.37
# f lognorm 3.93, alpha -1.96, score 1.97
# g lognorm 130.19, alpha -5.47, score 124.73
# f+g score 126.69
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

/* We assume a difference <= ALPHA_BOUND_GUARD between alpha computed
   with ALPHA_BOUND_SMALL and ALPHA_BOUND. In practice the largest value
   observed is 0.79. */
#define ALPHA_BOUND_GUARD 1.0

mpz_poly *best_f = NULL;
double *best_score = NULL, Best_score = DBL_MAX;
unsigned long keep = 0; /* target number of polynomials f with best alpha */
unsigned long best_n = 0; /* number of kept polynomials so far */
unsigned long f_irreducible = 0;
double max_guard = DBL_MIN;

#define PARENT(i) (((i)-1)/2) /* 1,2 -> 0, 3,4 -> 1, 5,6 -> 2, ... */
#define LEFT(i) (2*(i)+1)
#define RIGHT(i) (2*(i)+2)

/* Add (f, score) to heap. The top element of the heap has the largest score
   (i.e., it is the worst polynomial). */
static void MAYBE_UNUSED
heap_add (mpz_poly f, double score)
{
  unsigned long i, j;

  if (best_n < keep)
    {
      i = best_n;
      /* move towards top of heap to keep heap property */
      while (i > 0 && score > best_score[PARENT(i)])
        {
          j = PARENT(i);
          mpz_poly_swap (best_f[i], best_f[j]);
          best_score[i] = best_score[j];
          i = j;
        }
      mpz_poly_set (best_f[i], f);
      best_score[i] = score;
      best_n ++;
    }
  else /* heap is full, we replace the top element */
    {
      /* the new element should have a better (i.e., smaller) score than the
         top element of the heap */
      ASSERT (best_n == keep);
      i = 0;
      while (LEFT(i) < keep)
        {
          j = i; /* invariant: j is the element with the largest score
                    among i and its two sons, which should replace i */
          if (best_score[LEFT(i)] > score)
            j = LEFT(i);
          if (RIGHT(i) < keep && best_score[RIGHT(i)] > score &&
              best_score[RIGHT(i)] > best_score[LEFT(i)])
            j = RIGHT(i);
          if (j == i)
            break; /* both sons have a smaller score */
          mpz_poly_swap (best_f[i], best_f[j]);
          best_score[i] = best_score[j];
          i = j;
        }
      mpz_poly_set (best_f[i], f);
      best_score[i] = score;
    }
}

#if 0
static void
check_heap (void)
{
  unsigned long i;

  for (i = 0; LEFT(i) < best_n; i++)
    ASSERT_ALWAYS(best_score[i] >= best_score[LEFT(i)]);
  for (i = 0; RIGHT(i) < best_n; i++)
    ASSERT_ALWAYS(best_score[i] >= best_score[RIGHT(i)]);
}
#endif

static void
save_f (mpz_t *f, unsigned int df)
{
  mpz_poly ff;
  double skew, logmu, alpha, score;

  ff->deg = df;
  ff->coeff = f;

  skew = L2_skewness (ff, SKEWNESS_DEFAULT_PREC);
  logmu = L2_lognorm (ff, skew);
  alpha = get_alpha (ff, ALPHA_BOUND_SMALL);
  score = logmu + alpha;

  if (best_n == keep && score > best_score[0] + ALPHA_BOUND_GUARD)
    return;

  /* refine alpha */
  alpha = get_alpha (ff, ALPHA_BOUND);
  score = logmu + alpha;

  if (best_n == keep && score >= best_score[0])
    return;

#ifdef HAVE_OPENMP
#pragma omp critical
#endif
  {
    heap_add (ff, score);
    // check_heap ();
    if (score < Best_score)
      {
        Best_score = score;
        printf ("lognorm %1.2f, alpha %1.2f, score %1.2f: ",
                logmu, alpha, score);
        mpz_poly_fprintf (stdout, ff);
        fflush (stdout);
      }
  }
}

/*
  Print two nonlinear poly info. Return non-zero for record polynomials.
*/
static int
print_nonlinear_poly_info (mpz_poly ff, double alpha_f, mpz_poly gg,
                           int format,  mpz_t n)
{
    unsigned int i;
    double skew, logmu[2], alpha_g_approx, alpha_g, score, score_approx;
    int df = ff->deg;
    mpz_t *f = ff->coeff;
    int dg = gg->deg;
    mpz_t *g = gg->coeff;
    static double best_score = DBL_MAX;
    /* the coefficients of g are O(n^(1/df)) */

    /* we use the skewness minimizing the sum of lognorms */
    skew = L2_combined_skewness2 (ff, gg, SKEWNESS_DEFAULT_PREC);
    logmu[1] = L2_lognorm (gg, skew);
    logmu[0] = L2_lognorm (ff, skew);
    /* first estimate alpha with a small bound */
    alpha_g_approx = get_alpha (gg, ALPHA_BOUND_SMALL);

    score_approx = logmu[1] + alpha_g_approx + logmu[0] + alpha_f;

    if (score_approx >= best_score + ALPHA_BOUND_GUARD)
      return 0;

    /* now get a more precise alpha value */
    alpha_g = get_alpha (gg, ALPHA_BOUND);

    score = logmu[1] + alpha_g + logmu[0] + alpha_f;
    if (score_approx - score > max_guard)
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
      max_guard = score_approx - score;

    if (score >= best_score)
      return 0; /* only print record scores */

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
      printf ("# f lognorm %1.2f, alpha %1.2f, score %1.2f\n",
	      logmu[0], alpha_f, logmu[0] + alpha_f);
      printf ("# g lognorm %1.2f, alpha %1.2f, score %1.2f\n",
	      logmu[1], alpha_g, logmu[1] + alpha_g);
      printf ("# f+g score %1.2f\n", score);
      printf ("\n");
      fflush (stdout);
    }
    return 1;
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
  mpz_init_set_ui (a, 1);
  mpz_init_set_ui (b, 1);
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

        /* content test */
        mpz_poly ff;
        ff->deg = d;
        ff->coeff = f;
        mpz_poly_content (t, ff);
        if (mpz_cmp_ui (t, 1) != 0 )
	  goto end; /* duplicate with f/t */

        /* irreducibility test */

        if (mpz_poly_squarefree_p (ff) == 0)
	  goto end;

        if (is_irreducible (ff) == 0)
          goto end;

#ifdef HAVE_OPENMP
#pragma omp critical
        f_irreducible ++;
#endif

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

/* Generate polynomial g(x) of degree dg, given root 'root' of f.
   It might be better to take into account the skewness of f in the LLL
   lattice, but experimentally this does not give better results (probably
   because LLL is not very sensible to a small change of the skewness). */
static void
polygen_JL_g (mpz_t N, int dg, mat_Z g, mpz_t root)
{
    int i, j;
    mpz_t a, b, det, r;

    mpz_init (det);
    mpz_init_set_ui (a, 1);
    mpz_init_set_ui (b, 1);
    mpz_init_set (r, root);
    for (i = 0; i <= dg + 1; i ++) {
        for (j = 0; j <= dg + 1; j ++) {
            mpz_set_ui (g.coeff[i][j], 0);
        }
    }

    for (i = 1;  i <= dg + 1; i++) {
        for (j = 1; j <= dg + 1; j++) {
            if (i == 1)
              {
                if (j == 1)
                    mpz_set (g.coeff[j][i], N);
                else
                  {
                    mpz_neg (g.coeff[j][i], r);
                    mpz_mul (r, r, root);
                  }
              }
            else
              mpz_set_ui (g.coeff[j][i], i==j);
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
polygen_JL2 (mpz_t n, unsigned int df, unsigned int dg, long bound,
	     unsigned long nb_comb, unsigned long c)
{
    unsigned int i, j, nr, format = 1;
    mpz_t *rf;
    mat_Z g;
    mpz_poly *v, u;
    long *a;
    double alpha_f;

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
    v = malloc ((dg + 1) * sizeof (mpz_poly));
    for (j = 0; j <= dg; j++)
      {
        v[j]->deg = dg;
        v[j]->coeff = g.coeff[j + 1] + 1;
      }
    mpz_poly_init (u, dg);
    u->deg = dg;
    a = malloc ((dg + 1) * sizeof (long));

    /* compute number of roots of the c-th best polynomial f */
    nr = mpz_poly_roots_mpz (rf, best_f[c], n);
    ASSERT(0 < nr && nr <= df);

    alpha_f = get_alpha (best_f[c], ALPHA_BOUND);

    for (i = 0; i < nr; i ++) {
        /* generate g of degree dg */
        polygen_JL_g (n, dg, g, rf[i]);

        /* we skip idx = 0 which should correspond to c[0] = ... = c[dg] = 0 */
        for (unsigned long idx = 1; idx < nb_comb; idx ++)
          {
            unsigned long k = idx;
	    int64_t g;

            /* compute first index */
            a[0] = k % (bound + 1);
            k = k / (bound + 1);
	    g = a[0];
            for (j = 1; j <= dg; j++)
              {
                a[j] = k % (2 * bound + 1);
                k = k / (2 * bound + 1);
                a[j] = (a[j] <= bound) ? a[j] : a[j] - (2 * bound + 1);
		g = gcd_int64 (g, a[j]);
              }
            ASSERT_ALWAYS(k == 0);
	    if (g != 1)
	      continue;
            mpz_poly_mul_si (u, v[0], a[0]);
	    for (j = 1; j <= dg; j++)
	      mpz_poly_addmul_si (u, v[j], a[j]);
            if (print_nonlinear_poly_info (best_f[c], alpha_f, u, format, n))
              {
#if 0 /* print coefficients of record combination */
                for (j = 0; j <= dg; j++)
                  printf ("%ld ", a[j]);
                printf ("\n");
#endif
              }
        }
    }

    /* clear */
    free (a);
    mpz_poly_clear (u);
    free (v);
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
    fprintf (stderr, "./dlpolyselect -N xxx -df xxx -dg xxx -bound xxx [-keep xxx] [-t xxx]\n");
    exit (1);
}


int
main (int argc, char *argv[])
{
    int i;
    mpz_t N;
    unsigned int df = 0, dg = 0;
    int nthreads = 1;
    unsigned int bound = 4; /* bound on the coefficients of f */
    unsigned long maxtries;
    double t1, t2;
    unsigned long modr = 0, modm = 1;

    t1 = seconds ();
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
        else if (argc >= 3 && strcmp (argv[1], "-keep") == 0) {
            keep = atoi (argv[2]);
            argv += 2;
            argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-modm") == 0) {
	    modm = strtoul (argv[2], NULL, 10);
            argv += 2;
            argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-modr") == 0) {
	    modr = strtoul (argv[2], NULL, 10);
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

    ASSERT_ALWAYS (bound >= 1);

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

    unsigned int bound2 = 1; /* bound on the coefficients of linear
                                combinations from the LLL short vectors */
    double nb_comb_f;
    unsigned long nb_comb;
    /* we try all combinations u = c[0]*v[0] + ... + c[dg]*v[dg] with
       -bound2 <= c[j] <= bound2, except 0 <= c[0] <= bound2 since u and -u
       are equivalent, and except c[0] = ... = c[dg] = 0, this gives a total
       of (bound2+1)*(2*bound2+1)^dg-1 */
    nb_comb_f = (double) (bound2 + 1) * pow ((double) (2 * bound2 + 1),
                                            (double) dg);
    nb_comb = (nb_comb_f > (double) ULONG_MAX) ? ULONG_MAX
      : (unsigned long) nb_comb_f;

    if (keep == 0)
      keep = (10 * (maxtries / modm)) / nb_comb;

    printf ("tries %lu, nb_comb %lu, keep %lu\n",
            maxtries / modm, nb_comb, keep);

    best_f = malloc ((keep + 1) * sizeof (mpz_poly));
    best_score = malloc ((keep + 1) * sizeof (double));
    for (unsigned int i = 0; i <= keep; i++)
      mpz_poly_init (best_f[i], df);

#ifdef HAVE_OPENMP
    omp_set_num_threads (nthreads);
#pragma omp parallel for schedule(dynamic)
#endif
    for (unsigned long c = modr; c < maxtries; c += modm)
      polygen_JL1 (N, df, bound, c);

    t1 = seconds () - t1;

    t2 = seconds ();

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (unsigned long c = 0; c < best_n; c++)
      polygen_JL2 (N, df, dg, bound2, nb_comb, c);

    t2 = seconds () - t2;

    printf ("found %lu irreducible f-polynomials, ", f_irreducible);
    printf ("best score %1.2f, worst %1.2f\n",
            Best_score, best_score[0]);
    printf ("tries %lu, nb_comb %lu, keep %lu, ",
            maxtries / modm, nb_comb, keep);
    printf ("max guard %1.2f\n", max_guard);
    printf ("Stage 1: %.0fs, Stage 2: %.0fs\n", t1, t2);

    mpz_clear (N);
    for (unsigned int i = 0; i <= keep; i++)
      mpz_poly_clear (best_f[i]);
    free (best_f);
    free (best_score);

    return 0;
}
