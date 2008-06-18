/* polynomial selection with Kleinjung's algorithm

   Reference:
   "On polynomial selection for the general number field sieve",
   Thorsten Kleinjung, Mathematics of Computation 75 (2006), p. 2037-2047.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cado.h"
#include "utils/utils.h"

#include "aux.c"

/* Outputs all (mu_1, ..., mu_l), 0 <= mu_i < d, such that S is at distance
   less than eps from an integer, with
   S = f0 + f[1][mu_1] + ... + f[l][mu_l].
*/
void
naive_search (double f0, double **f, int l, int d, double eps)
{
}

/* enumerates all subsets of kk elements of Q, where Q has lQ elements,
   such that the product does not exceed max_adm1 */
void
enumerate (unsigned int *Q, int lQ, int l, double max_adm1, double max_adm2,
           mpz_t ad, mpz_t N, int d, mpz_t *g, double mtilde)
{
  int *p, k, i, j;
  mpz_t **x, **m, **e, t, u, P, P_over_pi, m0, invN, M0;
  unsigned int pi;
  LONG *roots;
  double eps, f0, **f, one_over_P2;

  p = (int*) malloc (l * sizeof (int));
  for (k = 0; k < l; k++)
    p[k] = k;
  x = (mpz_t**) malloc (l * sizeof(mpz_t*));
  m = (mpz_t**) malloc (l * sizeof(mpz_t*));
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
  /* m[0][j] = m0 + x[0][j], and m[i][j] = x[i][j] for i > 0 */
  m[0] = (mpz_t*) malloc (d * sizeof(mpz_t));
  for (j = 0; j < d; j++)
    mpz_init (m[0][j]);
  for (i = 1; i < l; i++)
    m[i] = x[i];
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
          /* print current subset */
          gmp_printf ("%Zd", ad);
          for (k = 0; k < l; k++)
            printf (" %u", Q[p[k]]);
          printf (": %e\n", mpz_get_d (P));

          /* compute 1/N mod P */
          mpz_invert (invN, N, P);

          /* m0 is the smallest integer bigger than mtilde and divisible by P:
             m0 = s*P, where s = ceil(mtilde/P) = floor((mtilde + P - 1)/P) */
          mpz_set_d (t, mtilde);
          mpz_add (t, t, P);
          mpz_sub_ui (t, t, 1);
          mpz_tdiv_q (t, t, P);
          mpz_mul (m0, t, P);
          eps = max_adm2 / mpz_get_d (m0);

          /* compute f0 */
          mpz_pow_ui (t, m0, d);
          mpz_mul (t, t, ad);
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
              mpz_invert (t, ad, u);
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

          /* compute the m[i][j] from (3.3): we only need to compute m[0][j],
             since m[i][j] = x[i][j] for i > 0 */
          for (j = 0; j < d; j++)
            mpz_add (m[0][j], m0, x[0][j]);

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
              mpz_mul (u, u, ad);
              mpz_sub (u, N, u);
              ASSERT (mpz_divisible_p (u, P));
              mpz_divexact (u, u, P);
              mpz_mul (u, u, invN);
              mpz_mul (u, u, ad);
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
                  mpz_sub (t, t, m[i][j - 1]);
                  mpz_add (t, t, m[i][j]);
                  /* now t = m_{(1,...,j,...,1)} where the j is at the
                     ith place */
                  mpz_pow_ui (u, t, d);
                  mpz_mul (u, u, ad);
                  mpz_sub (u, N, u);
                  ASSERT (mpz_divisible_p (u, P));
                  mpz_divexact (u, u, P);
                  mpz_mul (u, u, invN);
                  mpz_mul (u, u, ad);
                  mpz_mul (e[i][j], u, t);
                  mpz_sub (e[i][j], e[i][j], e[0][0]);
                  mpz_mod (e[i][j], e[i][j], P);
                }
            }

          /* finally compute the f[i][j] */
          mpz_mul_ui (t, ad, d);
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
          naive_search (f0, f, l, d, eps);
        }
      
      /* go to next subset */
      for (k = l - 1; p[k] == lQ - l + k; k--);
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
  for (j = 0; j < d; j++)
    mpz_clear (m[0][j]);
  free (x);
  free (m);
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
}

void
Algo36 (mpz_t N, unsigned int d, double M, unsigned int l, unsigned int pb)
{
  unsigned int *P = NULL, lP = 0;
  unsigned int *Q = NULL, lQ = 0;
  unsigned int r, i;
  mpz_t *a, *g, t;
  double Nd, max_ad, mtilde, max_adm1, max_adm2;

  ASSERT(d >= 4);

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
  fprintf (stderr, "P has %u primes\n", lP);

  a = alloc_mpz_array (d + 1);
  g = alloc_mpz_array (d + 1);
  /* g will store the polynomial x^d - t */
  mpz_set_ui (g[d], 1);
  for (i = 1; i < d; i++)
    mpz_set_ui (g[i], 0);

  Nd = mpz_get_d (N);

  max_ad = pow(pow (M, (double) (2 * d - 2)) / Nd, 1.0 / (double) (d - 3));
  fprintf (stderr, "max a[d]=%e\n", max_ad);

  mpz_set_ui (a[d], 1);

  Q = (unsigned int*) malloc (lP * sizeof (unsigned int));
  mpz_init (t);

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

      mtilde = pow (Nd / mpz_get_d (a[d]), 1.0 / (double) d);
      max_adm1 = M * M / mtilde;
      max_adm2 = pow (pow (M, (double) (2 * d - 6)) /
                      pow (mtilde, (double) (d - 4)), 1.0 / (double) (d - 2));

      /* enumerate all subsets Pprime of at least l elements of Q such that
         prod(r, r in Pprime) <= max_adm1 */
      //      fprintf (stderr, "max_adm1=%e\n", max_adm1);
      for (i = l; i <= lQ; i++)
        enumerate (Q, lQ, i, max_adm1, max_adm2, a[d], N, d, g, mtilde);

    next_ad:
      mpz_add_ui (a[d], a[d], 1);
    }
  while (mpz_cmp_d (a[d], max_ad) <= 0);

  free (P);
  free (Q);
  clear_mpz_array (a, d + 1);
  clear_mpz_array (g, d + 1);
  mpz_clear (t);
}

int
main (int argc, char *argv[])
{
  cado_input in;

  init (in);

  parse_input (in);
  gmp_printf ("N=%Zd\n", in->n);

  /* for Elie: d=5, M=1e21, l=8, pb=256 */
  Algo36 (in->n, 5, 1e21, 7, 256);
  
  clear (in);
  return 0;
}
