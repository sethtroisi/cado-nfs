/* polynomial selection with Kleinjung's algorithm

   Reference:
   "On polynomial selection for the general number field sieve",
   Thorsten Kleinjung, Mathematics of Computation 75 (2006), p. 2037-2047.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <values.h> /* for DBL_MAX */
#include "cado.h"
#include "utils/utils.h"

#include "aux.c"

/* Compute the sup-norm of f(x) = A[d]*x^d + ... + A[0]
   as in Definition 3.1 of Kleinjung's paper
*/
double
sup_norm (mpz_t *A, int d)
{
  double norm, min_norm, *a, s, t;
  int i, j, k;

  min_norm = DBL_MAX;
  a = (double*) malloc ((d + 1) * sizeof (double));
  for (i = 0; i <= d; i++)
    a[i] = fabs (mpz_get_d (A[i]));
  for (i = 0; i <= d; i++)
    {
      if (a[i] == 0.0)
        continue;
      for (j = i + 1; j <= d; j++)
        {
          s = pow (a[j] / a[i], 1.0 / (double) (i - j));
          norm = a[i] * pow (s, (double) i - (double) d / 2.0);
          for (k = 0; norm < min_norm && k <= d; k++)
            if (k != i && k != j)
              {
                t = a[k] * pow (s, (double) k - (double) d / 2.0);
                if (t > norm)
                  norm = min_norm; /* will exit the loop */
              }
          if (norm < min_norm)
            min_norm = norm;
        }
    }
  free (a);
  return min_norm; 
}

/* Implements Lemma 2.1 from Kleinjung's paper, assumes a[d] is already set. */
void
Lemma21 (mpz_t *a, mpz_t N, int d, mpz_t p, mpz_t m)
{
  mpz_t r, mi;
  int i;

  mpz_init (r);
  mpz_init (mi);
  mpz_set (r, N);
  mpz_pow_ui (mi, m, d);
  for (i = d - 1; i >= 0; i--)
    {
      /* invariant: mi = m^(i+1) */
      mpz_mul (a[i], a[i+1], mi);
      mpz_sub (r, r, a[i]);
      ASSERT (mpz_divisible_p (r, p));
      mpz_divexact (r, r, p);
      mpz_divexact (mi, mi, m); /* now mi = m^i */
      mpz_invert (a[i], p, mi); /* 1/p mod m^i */
      mpz_mul (a[i], a[i], r);
      mpz_neg (a[i], a[i]);
      mpz_mod (a[i], a[i], mi); /* -r/p mod m^i */
      mpz_mul_2exp (a[i], a[i], 1);
      if (mpz_cmp (a[i], mi) >= 0)
        {
          mpz_div_2exp (a[i], a[i], 1);
          mpz_sub (a[i], a[i], mi);
        }
      else
        mpz_div_2exp (a[i], a[i], 1);
      mpz_mul (a[i], a[i], p);
      mpz_add (a[i], a[i], r);
      ASSERT (mpz_divisible_p (a[i], mi));
      mpz_divexact (a[i], a[i], mi);
    }
  mpz_clear (r);
  mpz_clear (mi);
}

/* Outputs all (mu[0], ..., mu[l-1]), 0 <= mu_i < d, such that S is at distance
   less than eps from an integer, with
   S = f0 + f[0][mu[0]] + ... + f[l-1][mu[l-1]].
   Assumes a[d] is set to the current search value.
*/
void
naive_search (double f0, double **f, int l, int d, double eps, mpz_t *a,
	      mpz_t P, mpz_t N, double M, mpz_t **x, mpz_t m0)
{
  int *mu, i;
  double *s, fr, norm;
  mpz_t t;

  mu = (int*) malloc (l * sizeof (int));
  s = (double*) malloc ((l + 1) * sizeof (double));
  /* s[i] contains f0 + f[0][mu[0]] + ... + f[i-1][mu[i-1]] */

  /* initializes mu[] to (0,...,0) */
  for (i = 0; i < l; i++)
    mu[i] = 0;
  s[0] = f0;
  for (i = 1; i <= l; i++)
    s[i] = s[i-1] + f[i-1][mu[i-1]];
  mpz_init (t);

  while (1)
    {
      /* check current sum */
      fr = fabs (s[l] - round (s[l]));
      if (fr <= eps)
	{
	  //mpz_set (t, m[0][mu[0]]);
          mpz_add (t, m0, x[0][mu[0]]);
	  for (i = 1; i < l; i++)
	    // mpz_add (t, t, m[i][mu[i]]);
            mpz_add (t, t, x[i][mu[i]]);
          Lemma21 (a, N, d, P, t);
          norm = sup_norm (a, d);
          if (norm <= M)
            {
              gmp_printf ("p=%Zd m=%Zd norm=%1.2e\n", P, t, norm);
              fprint_polynomial (stdout, a, d);
              printf ("\n");
              fflush (stdout);
            }
	}
      
      /* go to next combination */
      for (i = l - 1; i >= 0 && mu[i] == d - 1; i--);
      /* now either i < 0 and we are done, or mu[i] < d-1 */
      if (i < 0)
	break;
      mu[i] ++;
      s[i+1] = s[i] + f[i][mu[i]];
      while (++i < l)
	{
	  mu[i] = 0;
	  s[i+1] = s[i] + f[i][mu[i]];
	}
    }

  mpz_clear (t);
  free (mu);
  free (s);
}

/* Enumerates all subsets of kk elements of Q, where Q has lQ elements,
   such that the product does not exceed max_adm1.
   Assumes a[d] is set to the current search value.
*/
void
enumerate (unsigned int *Q, int lQ, int l, double max_adm1, double max_adm2,
           mpz_t *a, mpz_t N, int d, mpz_t *g, mpz_t mtilde, double M)
{
  int *p, k, i, j;
  mpz_t **x, **e, t, u, P, P_over_pi, m0, invN, M0;
  unsigned int pi;
  LONG *roots;
  double eps, f0, **f, one_over_P2;

  p = (int*) malloc (l * sizeof (int));
  for (k = 0; k < l; k++)
    p[k] = k;
  x = (mpz_t**) malloc (l * sizeof(mpz_t*));
  e = (mpz_t**) malloc (l * sizeof(mpz_t*));
  f = (double**) malloc (l * sizeof(double*));
  for (i = 0; i < l; i++)
    {
      x[i] = (mpz_t*) malloc (d * sizeof(mpz_t));
      e[i] = (mpz_t*) malloc (d * sizeof(mpz_t));
      f[i] = (double*) malloc (d * sizeof(double));
      for (j = 0; j < d; j++)
        {
          mpz_init (x[i][j]);
          mpz_init (e[i][j]);
        }
    }
  mpz_init (t);
  mpz_init (u);
  mpz_init (P_over_pi);
  mpz_init (P);
  mpz_init (m0);
  mpz_init (M0);
  mpz_init (invN);
  roots = (LONG*) malloc (d * sizeof(LONG));
  do
    {
      /* compute product of current subset */
      mpz_set_ui (P, 1);
      for (k = 0; k < l; k++)
        mpz_mul_ui (P, P, Q[p[k]]);

      if (mpz_get_d (P) <= max_adm1)
        {
#if 0
          /* print current subset */
          gmp_printf ("%Zd", a[d]);
          for (k = 0; k < l; k++)
            printf (" %u", Q[p[k]]);
          printf (": %e\n", mpz_get_d (P));
#endif

          /* compute 1/N mod P */
          mpz_invert (invN, N, P);

          /* m0 is the smallest integer bigger than mtilde and divisible by P:
             m0 = s*P, where s = ceil(mtilde/P) = floor((mtilde + P - 1)/P) */
          mpz_add (t, mtilde, P);
          mpz_sub_ui (t, t, 1);
          mpz_tdiv_q (t, t, P);
          mpz_mul (m0, t, P);
          eps = max_adm2 / mpz_get_d (m0);

          /* compute f0 */
          mpz_pow_ui (t, m0, d);
          mpz_mul (t, t, a[d]);
          mpz_sub (t, N, t); /* N - a[d]*m0^d */
          f0 = mpz_get_d (t);
          mpz_pow_ui (t, m0, d - 1);
          mpz_mul (t, t, P);
          mpz_mul (t, t, P);
          f0 = f0 / mpz_get_d (t);
          
          /* compute the x[i][j] from (3.2) */
          for (i = 0; i < l; i++)
            {
              /* put in x[i][0..d-1] the d roots of x^d = N/a[d] mod Q[p[i]] */
              pi = Q[p[i]];
              mpz_set_ui (u, pi);
              mpz_invert (t, a[d], u);
              mpz_mul (t, t, N);
              mpz_mod_ui (t, t, pi);
              mpz_ui_sub (g[0], pi, t);
              if (roots_mod_long (roots, g, d, pi) != d)
                {
                  fprintf (stderr, "Error, d roots expected\n");
                  exit (1);
                }
              mpz_divexact_ui (P_over_pi, P, pi);
              mpz_invert (t, P_over_pi, u);       /* 1 / (P/pi) mod pi */
              for (j = 0; j < d; j++)
                {
                  /* we want x[i][j] = c*(P/pi) and x[i][j] = roots[j] mod pi,
                     thus c = roots[j] / (P/pi) mod pi */
                  mpz_mul_si (x[i][j], t, roots[j]);
                  mpz_mod_ui (x[i][j], x[i][j], pi);
                  mpz_mul (x[i][j], x[i][j], P_over_pi);
                }
            }

          /* The m[i][j], cf (3.3) are not needed, since m[i][j] = x[i][j]
             for i >= 1, and m[0][j] = m0 + x[0][j] */

          /* compute the e[i][j] from (3.6) */
          /* first compute e[0][j] = a_{d-1, (j,...,1)} */
          mpz_set (M0, m0);
          for (i = 0; i < l; i++)
            mpz_add (M0, M0, x[i][0]);
          mpz_set (t, M0);
          /* t = m0 + x_{(1,...,1)} = m_{(1,...,1)} */
          for (j = 0; j < d; j++)
            {
              if (j > 0)
                {
                  /* m_{(j,1,...,1)} = m_{(j-1,1,...,1)} + x_{1,j}-x_{1,j-1} */
                  mpz_sub (t, t, x[0][j - 1]);
                  mpz_add (t, t, x[0][j]);
                }
              mpz_pow_ui (u, t, d);
              mpz_mul (u, u, a[d]);
              mpz_sub (u, N, u);
              ASSERT (mpz_divisible_p (u, P));
              mpz_divexact (u, u, P);
              mpz_mul (u, u, invN);
              mpz_mul (u, u, a[d]);
              mpz_mul (u, u, t);
              mpz_mod (e[0][j], u, P);
            }
          /* now compute e[i][j] for i > 0 */
          for (i = 1; i < l; i++)
            {
              mpz_set_ui (e[i][0], 0); /* e_{i,1} = 0 */
              mpz_set (t, M0);         /* m_{(1,...,1)} */
              for (j = 1; j < d; j++)
                {
                  /* since i >= 1, m[i][.] = x[i][.] */
                  mpz_sub (t, t, x[i][j - 1]);
                  mpz_add (t, t, x[i][j]);
                  /* now t = m_{(1,...,j,...,1)} where the j is at the
                     ith place */
                  mpz_pow_ui (u, t, d);
                  mpz_mul (u, u, a[d]);
                  mpz_sub (u, N, u);
                  ASSERT (mpz_divisible_p (u, P));
                  mpz_divexact (u, u, P);
                  mpz_mul (u, u, invN);
                  mpz_mul (u, u, a[d]);
                  mpz_mul (e[i][j], u, t);
                  mpz_sub (e[i][j], e[i][j], e[0][0]);
                  mpz_mod (e[i][j], e[i][j], P);
                }
            }

          /* finally compute the f[i][j] */
          mpz_mul_ui (t, a[d], d);
          one_over_P2 = mpz_get_d (P);
          one_over_P2 = -1.0 / (one_over_P2 * one_over_P2);
          for (i = 0; i < l; i++)
            for (j = 0; j < d; j++)
              {
                mpz_mul (u, t, x[i][j]);
                mpz_addmul (u, e[i][j], P);
                f[i][j] = mpz_get_d (u) * one_over_P2;
              }

          /* now search for a small combination */
          naive_search (f0, f, l, d, eps, a, P, N, M, x, m0);
        }
      
      /* go to next subset */
      for (k = l - 1; k >= 0 && p[k] == lQ - l + k; k--);
      if (k < 0)
        break;
      p[k] ++;
      while (++k < l)
        p[k] = p[k - 1] + 1;
    }
  while (1);
  for (i = 0; i < l; i++)
    {
      for (j = 0; j < d; j++)
        {
          mpz_clear (x[i][j]);
          mpz_clear (e[i][j]);
        }
      free (x[i]);
      free (e[i]);
      free (f[i]);
    }
  free (x);
  free (e);
  free (f);
  mpz_clear (t);
  mpz_clear (u);
  mpz_clear (P_over_pi);
  mpz_clear (P);
  mpz_clear (m0);
  mpz_clear (M0);
  mpz_clear (invN);
  free (roots);
  free (p);
}

/* N is the number to factor
   d is the wanted degree
   M is the sup-norm bound
   l is the number of primes = 1 mod d in p
   pb is the prime bound for those primes
   incr is the increment for a[d]
*/
void
Algo36 (mpz_t N, unsigned int d, double M, unsigned int l, unsigned int pb,
        int incr)
{
  unsigned int *P = NULL, lP = 0;
  unsigned int *Q = NULL, lQ = 0;
  unsigned int r, i;
  mpz_t *a, *g, t, mtilde;
  double Nd, max_ad, max_adm1, max_adm2;

  ASSERT_ALWAYS (d >= 4);

  /* step 1 */
  for (r = 1; r < pb; r += d)
    if (isprime (r))
      {
        if (mpz_divisible_ui_p (N, r))
          fprintf (stderr, "Warning, N is divisible by %u\n", r);
        else /* add r to P */
          {
            lP ++;
            P = (unsigned int*) realloc (P, lP * sizeof (unsigned int));
            P[lP - 1] = r;
          }
      }
  //  fprintf (stderr, "P has %u primes\n", lP);

  a = alloc_mpz_array (d + 1);
  g = alloc_mpz_array (d + 1);
  /* g will store the polynomial x^d - t */
  mpz_set_ui (g[d], 1);
  for (i = 1; i < d; i++)
    mpz_set_ui (g[i], 0);

  Nd = mpz_get_d (N);

  max_ad = pow(pow (M, (double) (2 * d - 2)) / Nd, 1.0 / (double) (d - 3));
  fprintf (stderr, "max ad=%1.2e\n", max_ad);
  fflush (stderr);

  mpz_set_ui (a[d], incr);

  Q = (unsigned int*) malloc (lP * sizeof (unsigned int));
  mpz_init (t);
  mpz_init (mtilde);

  /* step 2 */
  do
    {
      for (i = lQ = 0; i < lP; i++)
        {
          /* add r = P[i] to Q if a[d]/N <> 0 and is a dth power modulo r */
          r = P[i];
          mpz_set_ui (t, r);
          mpz_invert (t, N, t);
          mpz_mul (t, t, a[d]);
          mpz_mod_ui (t, t, r);
          if (mpz_cmp_ui (t, 0) != 0)
            {
              mpz_set_si (g[0], - mpz_get_si (t));
              if (roots_mod_long (NULL, g, d, r) > 0)
                Q[lQ++] = r;
            }
        }
      if (lQ < l)
        goto next_ad;

      mpz_tdiv_q (mtilde, N, a[d]);
      mpz_root (mtilde, mtilde, d);
      max_adm1 = M * M / mpz_get_d (mtilde);
      max_adm2 = pow (pow (M, (double) (2 * d - 6)) /
                      pow (mpz_get_d (mtilde), (double) (d - 4)),
                      1.0 / (double) (d - 2));

      /* enumerate all subsets Pprime of at least l elements of Q such that
         prod(r, r in Pprime) <= max_adm1 */
      //      fprintf (stderr, "max_adm1=%e\n", max_adm1);
      for (i = l; i <= lQ; i++)
        enumerate (Q, lQ, i, max_adm1, max_adm2, a, N, d, g, mtilde, M);

    next_ad:
      mpz_add_ui (a[d], a[d], incr);
    }
  while (mpz_cmp_d (a[d], max_ad) <= 0);

  free (P);
  free (Q);
  clear_mpz_array (a, d + 1);
  clear_mpz_array (g, d + 1);
  mpz_clear (t);
  mpz_clear (mtilde);
}

static void
usage ()
{
  fprintf (stderr, "Usage: kleinjung [-v] [-degree d] [-e e] [-incr i] [-l l] [-M M] [-pb p] < in\n\n");
  fprintf (stderr, "       -v        - verbose\n");
  fprintf (stderr, "       -degree d - use algebraic polynomial of degree d (default 5)\n");
  fprintf (stderr, "       -e e      - use effort e (default 1e6)\n");
  fprintf (stderr, "       -incr i   - ad is incremented by i (default 1)\n");
  fprintf (stderr, "       -l l      - leading coefficient of g(x) has l prime factors (default 7)\n");
  fprintf (stderr, "       -M M      - keep polynomials with sup-norm <= M (default 1e25)\n");
  fprintf (stderr, "       -pb p   - prime factors are bounded by p (default 256)\n");
  fprintf (stderr, "       in        - input file (n:...)\n");
  exit (1);
}

int
main (int argc, char *argv[])
{
  cado_input in;
  int degree = 5;
  int verbose = 0;
  double effort = 1e6;
  double M = 1e25;
  int l = 7;
  int pb = 256;
  int incr = 1; /* Implements remark (1) following Algorithm 3.6:
                   try a[d]=incr, then 2*incr, 3*incr, ... */
  int i;

  /* print command line */
  fprintf (stderr, "# %s.r%s", argv[0], REV);
  for (i = 1; i < argc; i++)
    fprintf (stderr, " %s", argv[i]);
  fprintf (stderr, "\n");
  init (in);

  /* parse options */
  while (argc > 1 && argv[1][0] == '-')
    {
      if (strcmp (argv[1], "-v") == 0)
	{
	  verbose ++;
	  argv ++;
	  argc --;
	}
      else if (argc > 2 && strcmp (argv[1], "-e") == 0)
	{
	  effort = atof (argv[2]);
	  argv += 2;
	  argc -= 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-degree") == 0)
	{
          degree = atoi (argv[2]);
	  argv += 2;
	  argc -= 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-incr") == 0)
	{
          incr = atoi (argv[2]);
	  argv += 2;
	  argc -= 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-l") == 0)
	{
          l = atoi (argv[2]);
	  argv += 2;
	  argc -= 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-M") == 0)
	{
          M = atof (argv[2]);
	  argv += 2;
	  argc -= 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-pb") == 0)
	{
          pb = atoi (argv[2]);
	  argv += 2;
	  argc -= 2;
	}
      else
	{
	  fprintf (stderr, "Unknown option: %s\n", argv[1]);
	  usage ();
	}
    }

  if (argc != 1)
    usage ();

  parse_input (in);

  Algo36 (in->n, degree, M, l, pb, incr);
  
  clear (in);
  return 0;
}
