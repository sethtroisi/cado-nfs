#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include "portability.h"
#include "utils.h"
#include "auxiliary.h"
#include "size_optimization.h"


/******************************************************************************/
/************************ Internal functions **********************************/
/******************************************************************************/

/***************** List of mpz_t, to handle list of translations **************/
typedef struct
{
  mpz_t *tab;
  uint64_t len;
  uint64_t alloc;
} list_mpz_s;
typedef list_mpz_s list_mpz_t[1];
typedef list_mpz_s * list_mpz_ptr;
typedef const list_mpz_s * list_mpz_srcptr;

static inline void
list_mpz_realloc (list_mpz_ptr l, uint64_t newalloc)
{
#if DEBUG >= 2
  fprintf (stderr, "debug: %s in %s: l->alloc = %" PRIu64 ", l->len = "
                   "%" PRIu64 ", will be reallocated to newalloc = %" PRIu64
                   "\n", __FILE__, __func__, l->alloc, l->len, newalloc);
#endif
  l->tab = (mpz_t *) realloc (l->tab, newalloc * sizeof (mpz_t));
  ASSERT_ALWAYS (l->tab != NULL);
  for (uint64_t i = l->alloc; i < newalloc; i++)
    mpz_init (l->tab[i]);
  l->alloc = newalloc;
}

static inline void
list_mpz_init (list_mpz_ptr l, uint64_t init_alloc)
{
  l->len = 0;
  l->alloc = 0;
  l->tab = NULL;
  list_mpz_realloc (l, init_alloc);
}

static inline void
list_mpz_clear (list_mpz_t l)
{
  for (uint64_t i = 0; i < l->alloc; i++)
    mpz_clear (l->tab[i]);
  free (l->tab);
}

static inline void
list_mpz_append (list_mpz_ptr l, mpz_t e)
{
  if (l->len == l->alloc)
    list_mpz_realloc (l, 2*l->alloc);
  mpz_set (l->tab[l->len], e);
  l->len++;
}

static inline void
list_mpz_append_from_rounded_double (list_mpz_ptr l, double e)
{
  if (l->len == l->alloc)
    list_mpz_realloc (l, 2*l->alloc);
  mpz_set_d (l->tab[l->len], e > 0 ? e + 0.5 : e - 0.5);
  l->len++;
}

static inline void
list_mpz_sort_and_remove_dup (list_mpz_ptr l, const int verbose)
{
  /* Sort list_k by increasing order */
  qsort (l->tab, l->len, sizeof (mpz_t), (int(*)(const void*,const void*)) mpz_cmp);

  /* Remove duplicates */
  uint64_t len = l->len;
  l->len = 1;
  for (unsigned int i = 1; i < len; i++)
    if (mpz_cmp (l->tab[i], l->tab[l->len-1]) != 0)
      list_mpz_append (l, l->tab[i]);
    else if (verbose)
      gmp_fprintf (stderr, "# sopt: Remove duplicate %Zd\n", l->tab[i]);
}


/****************************** Hash table for uint64_t ***********************/
typedef struct {
  uint64_t *tab;
  unsigned long size;
  unsigned long used;
} hash_table_uint64_s;

typedef hash_table_uint64_s hash_table_uint64_t[1];
typedef hash_table_uint64_s * hash_table_uint64_ptr;
typedef const hash_table_uint64_s * hash_table_uint64_srcptr;

/* Internal function. Do not use it directly, use hash_uint64_insert instead. */
/* return 1 if new, 0 if already present in table */
static inline int
hash_uint64_insert_internal (uint64_t *H, unsigned long size, uint64_t h)
{
  uint64_t i = h % size;

  while (H[i] != 0 && H[i] != h)
    if (++i == size)
      i = 0;
  if (H[i] == h) /* already present */
    return 0;
  else
  {
    ASSERT_ALWAYS(H[i] == 0);
    H[i] = h;
    return 1;
  }
}

static inline void
hash_uint64_realloc (hash_table_uint64_ptr H, unsigned long wanted_size)
{
  uint64_t *new;
  unsigned long new_used = 0, new_size;

  new_size = ulong_nextprime (wanted_size);
  new = (uint64_t *) malloc (new_size * sizeof (uint64_t));
  memset (new, 0, new_size * sizeof (uint64_t));
  for (unsigned long i = 0; i < H->size; i++)
    if (H->tab[i])
      new_used += hash_uint64_insert_internal (new, new_size, H->tab[i]);
  free (H->tab);
  H->tab = new;
  H->size = new_size;
  ASSERT_ALWAYS(new_used == H->used);
}

static inline void
hash_uint64_init (hash_table_uint64_ptr H, unsigned long init_size)
{
  H->used  = 0;
  H->size = 0;
  H->tab = NULL;
  hash_uint64_realloc (H, init_size);
}

static inline void
hash_uint64_clear (hash_table_uint64_ptr H)
{
  free (H->tab);
  H->tab = NULL;
  H->size = 0;
  H->used = 0;
}

static inline void
hash_uint64_reset (hash_table_uint64_ptr H)
{
  memset (H->tab, 0, H->size * sizeof (uint64_t));
  H->used = 0;
}

/* return 1 if new, 0 if already present in table */
static inline int
hash_uint64_insert (hash_table_uint64_ptr H, uint64_t h)
{
  if (2 * H->used >= H->size) /* realloc */
    hash_uint64_realloc (H, 2*H->size + 1);

  if (hash_uint64_insert_internal (H->tab, H->size, h))
  {
    H->used++;
    return 1;
  }
  else
    return 0;
}


/******************************************************************************/

/* store in Q[0..nb_approx-1] the nb_approx best rational approximations of q
   with denominator <= bound. nb_approx should be greater than 0
   Return the number of found approximations. */
static inline unsigned int
compute_rational_approximation (double *Q, double q, unsigned int nb_approx,
                                double bound)
{
  unsigned int n = 0;
  double *E = (double *) malloc (nb_approx * sizeof(double));
  ASSERT_ALWAYS (E != NULL);

  for (double den = 2.0; den <= bound; den += 1.0)
  {
    double num = floor (den * q + 0.5);
    double e = fabs (q - num / den);
    unsigned int j;
    /* search for duplicate before inserting */
    for (j = 0; j < n && e > E[j]; j++);
    if (j < n && e == E[j])
      continue;
    /* If we arrive here, it means that num/den and e should be inserted in
       position j in Q and E respectively */
    if (n < nb_approx - 1) /* We are not at maximum size (we need one space left
    for approximation with denominator 1) */
      n++;

    for (unsigned int l = n; l > j; l--)
    {
      Q[l] = Q[l-1];
      E[l] = E[l-1];
    }
    Q[j] = num / den;
    E[j] = e;
  }
  /* force approximation with denominator 1 at the end */
  Q[n++] = floor (q + 0.5);
  return n;
}

/* Compute the skewness value that are going to be used for LLL */
static inline void
sopt_list_skewness_values (mpz_t *skew, mpz_poly_srcptr f, mpz_poly_srcptr g,
                           const int verbose)
{
  const int d = f->deg;
  unsigned int i = 0;
  /* skew0 = (|g0/ad|)^(1/d) */
  double skew0 = pow (fabs (mpz_get_d (g->coeff[0]) / mpz_get_d (f->coeff[d])),
              1.0 / (double) d);
  if (verbose)
    fprintf (stderr, "# sopt: skew0 = %f\n", skew0);
  for (int e = SOPT_NSKEW; e <= 3 * SOPT_NSKEW; e++, i++)
  {
    double s = pow (skew0, (double) e / (double) (2*SOPT_NSKEW));
    mpz_set_d (skew[i], s + 0.5); /* round to nearest */
    if (verbose)
      gmp_fprintf (stderr, "# sopt: skew[%u] = %Zd\n", i, skew[i]);
  }
}

/* Computes the roots of R below bound_on_roots_and_extrema, and the extrema
   of R which are below the same bound and closer to 0 than their neighbourhood.
   Stores those roots and extrema in roots_and_extrema, and returns their
   number.
   Note: roots_and_extrema must have enough storage for all roots and ALL
   extrema, whether close to 0 or not; i.e., 2*deg - 1 entries is enough.
*/
static inline unsigned int
sopt_compute_roots_and_extrema_close_to_0 (double *roots_and_extrema,
                                           double_poly_srcptr R,
                                           double bound_on_roots_and_extrema,
                                           int verbose)
{
  unsigned int nr, ne;
  nr = double_poly_compute_all_roots_with_bound (roots_and_extrema, R,
                                                 bound_on_roots_and_extrema);

  if (R->deg < 2) /* No extremum of f' in this case */
    return nr;
  else
  {
    /* add roots of the derivative */
    double_poly_t dR;
    double_poly_init (dR, R->deg - 1);
    double_poly_derivative (dR, R);
    ne = double_poly_compute_all_roots_with_bound (roots_and_extrema + nr, dR,
                                                   bound_on_roots_and_extrema);

    /* Keep only those extrema which are closer to 0 than points in their
       neighbourhood. They are exactly those for which f(x) * f''(x) > 0. */

    double_poly_derivative (dR, dR); /* dR is now the second derivative of R */

    unsigned int kept = nr;
    for (unsigned int i = nr; i < nr + ne; i++)
    {
      const double x = roots_and_extrema[i];
      const double Rx = double_poly_eval(R, x);
      const double ddRx = double_poly_eval(dR, x);
      if (Rx * ddRx > 0.)
      {
        roots_and_extrema[kept++] = x;
        if (verbose)
          printf ("Keeping x = %f, f(x) = %f, f''(x) = %f\n", x, Rx, ddRx);
      }
      else if (verbose)
        printf ("Not keeping x = %f, f(x) = %f, f''(x) = %f\n", x, Rx, ddRx);
    }

    double_poly_clear (dR);
    return kept;
  }
}

/* For deg(f) = 6 and deg(g) = 1 */
/* Return in list_k good values of k such that f(x+k) has small coefficients of
   degree 4 and 3. */
static void
sopt_find_translations_deg6 (list_mpz_t list_k, mpz_poly_srcptr f,
                             mpz_poly_srcptr g, const int verbose)
{
  ASSERT_ALWAYS (f->deg == 6);
  ASSERT_ALWAYS (g->deg == 1);

  unsigned int i, nb_q3roots;

  /* Let f(x) = a6 x^6 + ... + a0 and g(x) = g1 x + g0. */
  double a6, a5, a4, a3, a2, g1, g0;
  a6 = mpz_get_d (f->coeff[6]);
  a5 = mpz_get_d (f->coeff[5]);
  a4 = mpz_get_d (f->coeff[4]);
  a3 = mpz_get_d (f->coeff[3]);
  a2 = mpz_get_d (f->coeff[2]);
  g1 = mpz_get_d (g->coeff[1]);
  g0 = mpz_get_d (g->coeff[0]);

  /* Let f~(x,k,q3) = f(x+k)+q3*x^3*g(x+k).
     Let c_i(k,q3) be the coefficient of x^i in f~. Then
        c_3(k,q3) = 20*a6*k^3 + 10*a5*k^2 + 4*a4*k + q3*g1*k + a3 + q3*g0
        c_4(k,q3) = 15*a6*k^2 + 5*a5*k + q3*g1 + a4

     Let res = Res (c_3(k,q3), c_4(k,q3)), where the resultant is taken to
     eliminate the variable k.
     We compute the roots of res and of the derivative of res.

     #Sage code to compute Res(c3(k), c4(k)):
        R.<x,k,q3,g1,g0,a6,a5,a4,a3,a2,a1,a0> = ZZ[]
        f = a6*x^6 + a5*x^5 + a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0
        g = g1*x + g0
        ff=f(x=x+k)+q3*x^3*g(x=x+k)
        c3=ff.coefficient({x:3})
        c4=ff.coefficient({x:4})
        r = c3.resultant(c4,k)//(25*a6)
        SR(r.coefficient({q3:3})).factor()
        SR(r.coefficient({q3:2})).factor()
        SR(r.coefficient({q3:1})).factor()
        SR(r.coefficient({q3:0})).factor()
  */

  double_poly_t res;
  double roots_q3[5]; /* 5 = 3 roots for res + 2 roots for the derivative */
  double_poly_init (res, 3);

  /* degree 3: a6*g1^3 */
  res->coeff[3] = a6 * g1 * g1 * g1;

 /* degree 2: 135*a6^2*g0^2 - 45*a5*a6*g0*g1 + 10*a5^2*g1^2 - 15*a4*a6*g1^2 */
  res->coeff[2] = 135.0 * a6 * a6 * g0 * g0 - 45.0 * a5 * a6 * g0 * g1
                  + 10.0 * a5 * a5 * g1 * g1 - 15.0 * a4 * a6 * g1 * g1;

  /* degree 1: 50*a5^3*g0 - 180*a4*a5*a6*g0 + 270*a3*a6^2*g0 - 10*a4*a5^2*g1
     + 48*a4^2*a6*g1 - 45*a3*a5*a6*g1 */
  res->coeff[1] = 50.0 * a5 * a5 * a5 * g0 - 180.0 * a4 * a5 * a6 * g0
                  + 270.0 * a3 * a6 * a6 * g0 - 10.0 * a4 * a5 * a5 * g1
                  + 48.0 * a4 * a4 * a6 * g1 - 45.0 * a3 * a5 * a6 * g1;

  /* degree 0: -20*a4^2*a5^2 + 50*a3*a5^3 + 64*a4^3*a6 - 180*a3*a4*a5*a6
     + 135*a3^2*a6^2 */
  res->coeff[0] = -20.0 * a4 * a4 * a5 * a5 + 50.0 * a3 * a5 * a5 * a5
                  + 64.0 * a4 * a4 * a4 * a6 - 180.0 * a3 * a4 * a5 * a6
                  + 135.0 * a3 * a3 * a6 * a6;

  /* compute roots of res and some roots of the derivative */
  nb_q3roots = sopt_compute_roots_and_extrema_close_to_0 (roots_q3, res,
                                         SOPT_MAX_VALUE_FOR_Q_ROOTS, verbose);
  if (verbose)
  {
    fprintf (stderr, "# sopt: q-roots of Res(c3,c4) or Res(c3,c4)' = {");
    for (i = 0; i < nb_q3roots; i++)
      fprintf (stderr, " %f%c", roots_q3[i], (i+1==nb_q3roots)?' ':',');
    fprintf (stderr, "}\n");
  }

  double_poly_t C;
  double_poly_init (C, 4);
  for (i = 0; i < nb_q3roots; i++)
  {
    unsigned int t;
    double q3_rat_approx[SOPT_NB_RAT_APPROX_OF_Q_ROOTS];
    t = compute_rational_approximation (q3_rat_approx, roots_q3[i],
                                        SOPT_NB_RAT_APPROX_OF_Q_ROOTS,
                                        SOPT_MAX_DEN_IN_RAT_APPROX_OF_Q_ROOTS);
    for (unsigned int j = 0; j < t; j++)
    {
      double q3_rat = q3_rat_approx[j];
      double double_roots_k[5];
      unsigned int nb_k_roots;

      if (verbose)
        fprintf (stderr, "# sopt: process rational approximation %f of %f\n",
                         q3_rat, roots_q3[i]);

      /* find roots k of c4(k, q3=q3_rat) */
      C->deg = 2;
      C->coeff[2] = 15.0 * a6;
      C->coeff[1] = 5.0 * a5;
      C->coeff[0] = g1 * q3_rat + a4;

      nb_k_roots = double_poly_compute_all_roots (double_roots_k, C);
      for (unsigned int l = 0; l < nb_k_roots; l++)
      {
        list_mpz_append_from_rounded_double (list_k, double_roots_k[l]);
        if (verbose)
          gmp_fprintf (stderr, "# sopt:   Adding k = %Zd (roots of c4)\n",
                       list_k->tab[list_k->len-1]);
      }

      /* find roots k of c3(k, q3=q3_rat) */
      C->deg = 3;
      C->coeff[3] = 20.0 * a6;
      C->coeff[2] = 10.0 * a5;
      C->coeff[1] = g1 * q3_rat + 4.0 * a4;
      C->coeff[0] = q3_rat * g0 + a3;

      nb_k_roots = double_poly_compute_all_roots (double_roots_k, C);
      for (unsigned int l = 0; l < nb_k_roots; l++)
      {
        list_mpz_append_from_rounded_double (list_k, double_roots_k[l]);
        if (verbose)
          gmp_fprintf (stderr, "# sopt:   Adding k = %Zd (roots of c3)\n",
                       list_k->tab[list_k->len-1]);
      }

      /* Let f~~(x,k,q3,q2) f(x+k) + (q3*x^3 + q2*x^2)*g(x,k)
         Add roots of Res(c3(k, q3=q3_rat, q2), c2(k, q3_q3rat, q2)) and its
         derivative, where the resultant is taken to eliminate the variable q2.

         #Sage code to compute Res(c3(k,q3,q2), c2(k,q3,q2)):
            R.<a6,k,a5,g1,q,q2,g0,a3,a4,a2,x> = ZZ[]
            f = a6*x^6+a5*x^5+a4*x^4+a3*x^3+a2*x^2
            g = g1*x+g0
            ff = f(x=x+k)+q*x^3*g(x=x+k)+q2*x^2*g(x=x+k)
            c3 = ff.coefficient({x:3})
            c2 = ff.coefficient({x:2})
            c3.resultant(c2,q2)
      */
      C->deg = 4;
      C->coeff[4] = -5.0 * a6 * g1;
      C->coeff[3] = -20.0 * a6 * g0;
      C->coeff[2] = -g1 * g1 * q3_rat - 10.0 * a5 * g0 + 2.0 * g1 * a4;
      C->coeff[1] = 2.0 * g1 * a3 - 2.0 * g1 * q3_rat * g0 - 4.0 * g0 * a4;
      C->coeff[0] = -q3_rat * g0 * g0 - g0 * a3 + g1 * a2;

      nb_k_roots = double_poly_compute_all_roots (double_roots_k, C);
      for (unsigned int l = 0; l < nb_k_roots; l++)
      {
        list_mpz_append_from_rounded_double (list_k, double_roots_k[l]);
        if (verbose)
          gmp_fprintf (stderr, "# sopt:   Adding k = %Zd (roots of "
                               "Res(c3,c2))\n", list_k->tab[list_k->len-1]);
      }

      double_poly_derivative (C, C);
      nb_k_roots = double_poly_compute_all_roots (double_roots_k, C);
      for (unsigned int l = 0; l < nb_k_roots; l++)
      {
        list_mpz_append_from_rounded_double (list_k, double_roots_k[l]);
        if (verbose)
          gmp_fprintf (stderr, "# sopt:   Adding k = %Zd (roots of derivative "
                               "of Res(c3,c2)\n", list_k->tab[list_k->len-1]);
      }

    }
#if 0
    double_poly_init (C, 4);
    int n, i, j, t, m, ret = 0;
    //double_poly_t R, C;
#define MAX_ROOTSK 40
    //double roots_k[MAX_ROOTSK];
    //double *r, roots_q[5];
    //R->deg = 3;
    //r = R->coeff;
      while (t) {
      roots_q[i] = Q[--t];


      /* add roots of Res(c3(q3=roots_q[i](k),c2(k)) (eliminate k) */
      double m0 = -g0;
      double roots_q2[7];
      int ii, nb_roots_q2 = 0;

      C->deg = 4;
      C->coeff[4] = -125*a6*a6*a6*g1*g1*g1*g1;
      C->coeff[3] = (3000*a6*a6*a6*m0*g1*g1*g1 + 500*a5*a6*a6*g1*g1*g1*g1)*roots_q[i] - 160000*a6*a6*a6*a6*m0*m0*m0 - 80000*a5*a6*a6*a6*m0*m0*g1 - 5000*a5*a5*a6*a6*m0*g1*g1 - 20000*a4*a6*a6*a6*m0*g1*g1 - 1000*a4*a5*a6*a6*g1*g1*g1 - 3500*a3*a6*a6*a6*g1*g1*g1;
      C->coeff[2] = -225*a6*a6*g1*g1*g1*g1*g1*roots_q[i]*roots_q[i]*roots_q[i] + (15750*a6*a6*a6*m0*m0*g1*g1 + 5250*a5*a6*a6*m0*g1*g1*g1 - 500*a5*a5*a6*g1*g1*g1*g1 + 2250*a4*a6*a6*g1*g1*g1*g1)*roots_q[i]*roots_q[i] + (20000*a5*a5*a6*a6*m0*m0*g1 - 48000*a4*a6*a6*a6*m0*m0*g1 + 10000*a5*a5*a5*a6*m0*g1*g1 - 28000*a4*a5*a6*a6*m0*g1*g1 + 18000*a3*a6*a6*a6*m0*g1*g1 + 2000*a4*a5*a5*a6*g1*g1*g1 - 7200*a4*a4*a6*a6*g1*g1*g1 + 5000*a3*a5*a6*a6*g1*g1*g1 - 4000*a2*a6*a6*a6*g1*g1*g1)*roots_q[i] - 50000*a5*a5*a5*a5*a6*m0*m0 + 240000*a4*a5*a5*a6*a6*m0*m0 - 192000*a4*a4*a6*a6*a6*m0*m0 - 240000*a3*a5*a6*a6*a6*m0*m0 + 480000*a2*a6*a6*a6*a6*m0*m0 - 20000*a4*a5*a5*a5*a6*m0*g1 + 80000*a4*a4*a5*a6*a6*m0*g1 + 10000*a3*a5*a5*a6*a6*m0*g1 - 216000*a3*a4*a6*a6*a6*m0*g1 + 160000*a2*a5*a6*a6*a6*m0*g1 - 2000*a4*a4*a5*a5*a6*g1*g1 + 7200*a4*a4*a4*a6*a6*g1*g1 - 1000*a3*a4*a5*a6*a6*g1*g1 + 5000*a2*a5*a5*a6*a6*g1*g1 - 33750*a3*a3*a6*a6*a6*g1*g1 + 20000*a2*a4*a6*a6*a6*g1*g1;
      C->coeff[1] = (-18000*a6*a6*a6*m0*m0*m0*g1 - 9000*a5*a6*a6*m0*m0*g1*g1 - 3600*a4*a6*a6*m0*g1*g1*g1 - 900*a3*a6*a6*g1*g1*g1*g1)*roots_q[i]*roots_q[i]*roots_q[i] + (-15000*a5*a5*a6*a6*m0*m0*m0 + 36000*a4*a6*a6*a6*m0*m0*m0 - 20000*a5*a5*a5*a6*m0*m0*g1 + 63000*a4*a5*a6*a6*m0*m0*g1 - 67500*a3*a6*a6*a6*m0*m0*g1 - 8000*a4*a5*a5*a6*m0*g1*g1 + 25200*a4*a4*a6*a6*m0*g1*g1 - 7500*a3*a5*a6*a6*m0*g1*g1 - 30000*a2*a6*a6*a6*m0*g1*g1 - 2000*a3*a5*a5*a6*g1*g1*g1 + 6300*a3*a4*a6*a6*g1*g1*g1 - 5000*a2*a5*a6*a6*g1*g1*g1)*roots_q[i]*roots_q[i] + (8000*a4*a4*a5*a5*a6*m0*g1 - 20000*a3*a5*a5*a5*a6*m0*g1 - 28800*a4*a4*a4*a6*a6*m0*g1 + 84000*a3*a4*a5*a6*a6*m0*g1 - 20000*a2*a5*a5*a6*a6*m0*g1 - 81000*a3*a3*a6*a6*a6*m0*g1 + 48000*a2*a4*a6*a6*a6*m0*g1 + 2000*a3*a4*a5*a5*a6*g1*g1 - 10000*a2*a5*a5*a5*a6*g1*g1 - 7200*a3*a4*a4*a6*a6*g1*g1 + 4500*a3*a3*a5*a6*a6*g1*g1 + 32000*a2*a4*a5*a6*a6*g1*g1 - 36000*a2*a3*a6*a6*a6*g1*g1)*roots_q[i] + 16000*a4*a4*a4*a5*a5*a6*m0- 60000*a3*a4*a5*a5*a5*a6*m0+ 100000*a2*a5*a5*a5*a5*a6*m0- 57600*a4*a4*a4*a4*a6*a6*m0+ 240000*a3*a4*a4*a5*a6*a6*m0+ 15000*a3*a3*a5*a5*a6*a6*m0- 480000*a2*a4*a5*a5*a6*a6*m0- 324000*a3*a3*a4*a6*a6*a6*m0+ 384000*a2*a4*a4*a6*a6*a6*m0+ 480000*a2*a3*a5*a6*a6*a6*m0- 480000*a2*a2*a6*a6*a6*a6*m0+ 4000*a3*a4*a4*a5*a5*a6*g1 - 20000*a3*a3*a5*a5*a5*a6*g1 + 20000*a2*a4*a5*a5*a5*a6*g1 - 14400*a3*a4*a4*a4*a6*a6*g1 + 81000*a3*a3*a4*a5*a6*a6*g1 - 80000*a2*a4*a4*a5*a6*a6*g1 - 10000*a2*a3*a5*a5*a6*a6*g1 - 121500*a3*a3*a3*a6*a6*a6*g1 + 216000*a2*a3*a4*a6*a6*a6*g1 - 80000*a2*a2*a5*a6*a6*a6*g1;
      C->coeff[0] = (3375*a6*a6*a6*m0*m0*m0*m0 + 2250*a5*a6*a6*m0*m0*m0*g1 + 1350*a4*a6*a6*m0*m0*g1*g1 + 675*a3*a6*a6*m0*g1*g1*g1 + 225*a2*a6*a6*g1*g1*g1*g1)*roots_q[i]*roots_q[i]*roots_q[i]*roots_q[i] + (5000*a5*a5*a5*a6*m0*m0*m0 - 18000*a4*a5*a6*a6*m0*m0*m0 + 27000*a3*a6*a6*a6*m0*m0*m0 + 3000*a4*a5*a5*a6*m0*m0*g1 - 10800*a4*a4*a6*a6*m0*m0*g1 + 4500*a3*a5*a6*a6*m0*m0*g1 + 18000*a2*a6*a6*a6*m0*m0*g1 + 1500*a3*a5*a5*a6*m0*g1*g1 - 5400*a3*a4*a6*a6*m0*g1*g1 + 6000*a2*a5*a6*a6*m0*g1*g1 + 500*a2*a5*a5*a6*g1*g1*g1 - 675*a3*a3*a6*a6*g1*g1*g1)*roots_q[i]*roots_q[i]*roots_q[i] + (-6000*a4*a4*a5*a5*a6*m0*m0 + 15000*a3*a5*a5*a5*a6*m0*m0 + 21600*a4*a4*a4*a6*a6*m0*m0 - 63000*a3*a4*a5*a6*a6*m0*m0 + 15000*a2*a5*a5*a6*a6*m0*m0 + 60750*a3*a3*a6*a6*a6*m0*m0 - 36000*a2*a4*a6*a6*a6*m0*m0 - 3000*a3*a4*a5*a5*a6*m0*g1 + 15000*a2*a5*a5*a5*a6*m0*g1 + 10800*a3*a4*a4*a6*a6*m0*g1 - 6750*a3*a3*a5*a6*a6*m0*g1 - 48000*a2*a4*a5*a6*a6*m0*g1 + 54000*a2*a3*a6*a6*a6*m0*g1 - 1500*a3*a3*a5*a5*a6*g1*g1 + 3000*a2*a4*a5*a5*a6*g1*g1 + 4050*a3*a3*a4*a6*a6*g1*g1 - 7200*a2*a4*a4*a6*a6*g1*g1 - 3000*a2*a3*a5*a6*a6*g1*g1 + 12000*a2*a2*a6*a6*a6*g1*g1)*roots_q[i]*roots_q[i] + 6000*a3*a3*a4*a4*a5*a5*a6 - 16000*a2*a4*a4*a4*a5*a5*a6 - 20000*a3*a3*a3*a5*a5*a5*a6 + 60000*a2*a3*a4*a5*a5*a5*a6 - 50000*a2*a2*a5*a5*a5*a5*a6 - 21600*a3*a3*a4*a4*a4*a6*a6 + 57600*a2*a4*a4*a4*a4*a6*a6 + 81000*a3*a3*a3*a4*a5*a6*a6 - 240000*a2*a3*a4*a4*a5*a6*a6 - 15000*a2*a3*a3*a5*a5*a6*a6 + 240000*a2*a2*a4*a5*a5*a6*a6 - 91125*a3*a3*a3*a3*a6*a6*a6 + 324000*a2*a3*a3*a4*a6*a6*a6 - 192000*a2*a2*a4*a4*a6*a6*a6 - 240000*a2*a2*a3*a5*a6*a6*a6 + 160000*a2*a2*a2*a6*a6*a6*a6;


      nb_roots_q2 += double_poly_compute_all_roots_with_bound (roots_q2, C,
                                                               SOPT_MAX_VALUE_FOR_Q_ROOTS);
      ASSERT_ALWAYS(nb_roots_q2 <= 4);

      /* roots of the derivative of Res(c3,c2) */
      C->deg = 3;
      C->coeff[0] = C->coeff[1];
      C->coeff[1] = C->coeff[2] * 2.0;
      C->coeff[2] = C->coeff[3] * 3.0;
      C->coeff[3] = C->coeff[4] * 4.0;
      nb_roots_q2 += double_poly_compute_all_roots_with_bound (roots_q2
                                                  + nb_roots_q2, C, SOPT_MAX_VALUE_FOR_Q_ROOTS);
      ASSERT_ALWAYS(nb_roots_q2 <= 7);

      for (ii = 0; ii < nb_roots_q2; ii++)
      {
        /* Compute roots of c2(q2=roots_q2[ii]) */
        C->deg = 4;
        C->coeff[4] = 15.0 * a6;
        C->coeff[3] = 10.0 * a5;
        C->coeff[2] = 6.0 * a4;
        C->coeff[1] = g1 * roots_q2[ii] + 3.0 * a3;
        C->coeff[0] = a2 + g0 * roots_q2[ii];

        m += double_poly_compute_all_roots (roots_k + m, C);
        ASSERT_ALWAYS(nb_roots_q2 <= 12+4*(ii+1));
      }
#endif
  }
  double_poly_clear (C);
  double_poly_clear (res);
}

/* For deg(f) = 5 and deg(g) = 1 */
/* Return in list_k good values of k such that f(x+k) has small coefficients of
   degree 3 and 2. */
static void
sopt_find_translations_deg5 (list_mpz_t list_k, mpz_poly_srcptr f,
                             mpz_poly_srcptr g, const int verbose)
{
  ASSERT_ALWAYS (f->deg == 5);
  ASSERT_ALWAYS (g->deg == 1);

  unsigned int i, nb_q2roots;

  /* Let f(x) = a5 x^5 + ... + a0 and g(x) = g1 x + g0. */
  double a5, a4, a3, a2, g1, g0;
  a5 = mpz_get_d (f->coeff[5]);
  a4 = mpz_get_d (f->coeff[4]);
  a3 = mpz_get_d (f->coeff[3]);
  a2 = mpz_get_d (f->coeff[2]);
  g1 = mpz_get_d (g->coeff[1]);
  g0 = mpz_get_d (g->coeff[0]);

  /* Let f~(x,k,q2) = f(x+k)+q2*x^2*g(x+k).
     Let c_i(k,q2) be the coefficient of x^i in f~. Then
        c_2(k,q2) = 10*a5*k^3 + 6*a4*k^2 + 3*a3*k + q2*g1*k + a2 + q2*g0
        c_3(k,q2) = 10*a5*k^2 + 4*a4*k + q2*g1 + a3

     Let res = Res (c_2(k,q3), c_3(k,q3)), where the resultant is taken to
     eliminate the variable k.
     We compute the roots of res and of the derivative of res.
     Degree 5 is the only case where res is of degree 2 in q2 (instead of 3)

     #Sage code to compute Res(c2(k), c3(k)):
        R.<x,k,q2,g1,g0,a5,a4,a3,a2,a1,a0> = ZZ[]
        f = a5*x^5 + a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0
        g = g1*x + g0
        ff=f(x=x+k)+q2*x^2*g(x=x+k)
        c2=ff.coefficient({x:2})
        c3=ff.coefficient({x:3})
        r = c2.resultant(c3,k)//(40*a5)
        SR(r.coefficient({q2:2})).factor()
        SR(r.coefficient({q2:1})).factor()
        SR(r.coefficient({q2:0})).factor()
  */

  double_poly_t res;
  double roots_q2[3]; /* 3 = 2 roots for res + 1 roots for the derivative */
  double_poly_init (res, 2);

  /* degree 2: (5*a5*g0 - a4*g1)^2 */
  res->coeff[2] =  (5.0*a5*g0 - a4*g1)*(5.0*a5*g0 - a4*g1);

  /* degree 1: 8*a4^3*g0 - 30*a3*a4*a5*g0 + 50*a2*a5^2*g0 - 2*a3*a4^2*g1
               + 10*a3^2*a5*g1 - 10*a2*a4*a5*g1 */
  res->coeff[1] = 8.0*a4*a4*a4*g0 - 30.0*a3*a4*a5*g0 + 50.0*a2*a5*a5*g0
                  - 2.0*a3*a4*a4*g1 + 10.0*a3*a3*a5*g1 - 10.0*a2*a4*a5*g1;

  /* degree 0: -3*a3^2*a4^2 + 8*a2*a4^3 + 10*a3^3*a5 - 30*a2*a3*a4*a5
               + 25*a2^2*a5^2 */
  res->coeff[0] = -3.0*a3*a3*a4*a4 + 8.0*a2*a4*a4*a4 + 10.0*a3*a3*a3*a5
                  - 30.0*a2*a3*a4*a5 + 25.0*a2*a2*a5*a5;

  /* compute roots of res and some roots of the derivative */
  nb_q2roots = sopt_compute_roots_and_extrema_close_to_0 (roots_q2, res,
                                         SOPT_MAX_VALUE_FOR_Q_ROOTS, verbose);
  if (verbose)
  {
    fprintf (stderr, "# sopt: q-roots of Res(c2,c3) or Res(c2,c3)' = {");
    for (i = 0; i < nb_q2roots; i++)
      fprintf (stderr, " %f%c", roots_q2[i], (i+1==nb_q2roots)?' ':',');
    fprintf (stderr, "}\n");
  }

  double_poly_t C;
  double_poly_init (C, 3);
  for (i = 0; i < nb_q2roots; i++)
  {
    unsigned int t;
    double q2_rat_approx[SOPT_NB_RAT_APPROX_OF_Q_ROOTS];
    t = compute_rational_approximation (q2_rat_approx, roots_q2[i],
                                        SOPT_NB_RAT_APPROX_OF_Q_ROOTS,
                                        SOPT_MAX_DEN_IN_RAT_APPROX_OF_Q_ROOTS);
    for (unsigned int j = 0; j < t; j++)
    {
      double q2_rat = q2_rat_approx[j];
      double double_roots_k[5];
      unsigned int nb_k_roots;

      if (verbose)
        fprintf (stderr, "# sopt: process rational approximation %f of %f\n",
                         q2_rat, roots_q2[i]);

      /* find roots k of c3(k, q2=q2_rat) */
      C->deg = 2;
      C->coeff[2] = 10.0 * a5;
      C->coeff[1] = 4.0 * a4;
      C->coeff[0] = g1 * q2_rat + a3;

      nb_k_roots = double_poly_compute_all_roots (double_roots_k, C);
      for (unsigned int l = 0; l < nb_k_roots; l++)
      {
        list_mpz_append_from_rounded_double (list_k, double_roots_k[l]);
        if (verbose)
          gmp_fprintf (stderr, "# sopt:   Adding k = %Zd (roots of c3)\n",
                       list_k->tab[list_k->len-1]);
      }

      /* find roots k of c2(k, q2=q2_rat) */
      C->deg = 3;
      C->coeff[3] = 10.0 * a5;
      C->coeff[2] = 6.0 * a4;
      C->coeff[1] = g1 * q2_rat + 3.0 * a3;
      C->coeff[0] = g0 * q2_rat + a2;

      nb_k_roots = double_poly_compute_all_roots (double_roots_k, C);
      for (unsigned int l = 0; l < nb_k_roots; l++)
      {
        list_mpz_append_from_rounded_double (list_k, double_roots_k[l]);
        if (verbose)
          gmp_fprintf (stderr, "# sopt:   Adding k = %Zd (roots of c2)\n",
                       list_k->tab[list_k->len-1]);
      }
    }
  }
  double_poly_clear (C);
  double_poly_clear (res);
}

/* If sopt_effort is not 0, add translations to the list of translations that
   are going to be tried in LLL */
static inline void
improve_list_k (list_mpz_ptr list_k, const int sopt_effort, const int verbose)
{
  mpz_t tmp, tmp2;
  mpz_init (tmp);
  mpz_init (tmp2);

  uint64_t len;

  int effort_sqrt = (int) sqrt ((double) sopt_effort);
  int effort_sqrt_max = effort_sqrt/2;
  int effort_sqrt_min = -effort_sqrt_max;
  int effort_max = (((int) sopt_effort) + 1)/2;
  int effort_min = -effort_max;

  unsigned int emax; /* emax is the number of digits of max(|k|) */
  mpz_abs (tmp, list_k->tab[0]);
  mpz_abs (tmp2, list_k->tab[list_k->len-1]);
  if  (mpz_cmp (tmp, tmp2) >= 0)
    emax = mpz_sizeinbase (tmp, 10);
  else
    emax = mpz_sizeinbase (tmp2, 10);

  /* For each value of k computed earlier, add k + delta_k, where
      delta_k is in [ effort_sqrt_min..effort_sqrt_max ] */
  /* Then add k = j*10^e for e in [0..emax], where emax = log_10 (max |k|) + 1,
                         and j in [ effort_min..effort_max ] */
  if (verbose)
    fprintf (stderr, "# sopt: sopt_effort_sqrt = "
                     "%d\n# sopt: adding delta_k to the %" PRIu64 " values "
                     "of k in the table, for delta_k in [%d..%d]\n",
                     effort_sqrt, list_k->len, effort_sqrt_min,
                     effort_sqrt_max);

  len = list_k->len;
  for (unsigned int i = 0; i < len; i++)
  {
    mpz_set (tmp, list_k->tab[i]);
    for (int delta_k = effort_sqrt_min; delta_k <= effort_sqrt_max; delta_k++)
    {
      mpz_add_si (tmp2, tmp, delta_k);
      list_mpz_append (list_k, tmp2);
    }
  }

  if (verbose)
  {
    fprintf (stderr, "# sopt: %" PRIu64 " values for the translations were "
                     "added\n# sopt: add k = j*10^e for e in [0..%d] and j "
                     "in [%d..%d]\n", list_k->len-len, emax, effort_min,
                     effort_max);
    len = list_k->len;
  }

  for (unsigned int e = 0; e <= emax; e++)
  {
    mpz_ui_pow_ui (tmp, 10, e);
    for (int j = effort_min; j <= effort_max; j++)
    {
      mpz_mul_si (tmp2, tmp, j);
      list_mpz_append (list_k, tmp2);
    }
  }
  if (verbose)
    fprintf (stderr, "# sopt: %" PRIu64 " values for the translations were "
                     "added\n", list_k->len-len);

  /* Sort list_k by increasing order and remove duplicates */
  list_mpz_sort_and_remove_dup (list_k, verbose);
  if (verbose)
    fprintf (stderr, "# sopt: It remains %" PRIu64 " values after sorting "
                     " and removing duplicates\n", list_k->len);

  mpz_clear (tmp);
  mpz_clear (tmp2);
}

/* Construct the LLL matrix from the polynomial pair. */
static inline void
LLL_set_matrix_from_polys (mat_Z *m, mpz_poly_srcptr ft, mpz_poly_srcptr gt,
                           mpz_srcptr skew, mpz_ptr tmp)
{
  for (int k = 1; k <= m->NumRows; k++)
    for (int l = 1; l <= m->NumCols; l++)
      mpz_set_ui (m->coeff[k][l], 0);
  mpz_set_ui (tmp, 1);
  for (int k = 0; k < m->NumCols; k++)
  {
    if (k > 0)
      mpz_mul (tmp, tmp, skew); /* skew^k */
    if (k <= m->NumRows - 2)
      mpz_set (m->coeff[k+2][k+1], tmp);
    mpz_mul (m->coeff[1][k+1], tmp, ft->coeff[k]);
  }
  for (int k = 0; k <= m->NumRows - 2; k++)
  {
    /* m.coeff[k+2][k+1] is already skew^k */
    mpz_mul (m->coeff[k+2][k+2], m->coeff[k+2][k+1], skew); /* skew^(k+1) */
    mpz_mul (m->coeff[k+2][k+1], m->coeff[k+2][k+1], gt->coeff[0]);
    mpz_mul (m->coeff[k+2][k+2], m->coeff[k+2][k+2], gt->coeff[1]);
  }
}

/* Construct poly from i-th row vector of length d of m. */
static inline void
LLL_set_poly_from_vector (mpz_poly_ptr flll, mat_Z U, int k, mpz_poly_srcptr ft,
                          mpz_poly_srcptr gt)
{
  int d = ft->deg;
  mpz_poly_realloc (flll, d+1);
  for (int i = 0; i <= d; i++)
  {
    mpz_mul (flll->coeff[i], ft->coeff[i], U.coeff[k][1]);
    if (i <= d-2)
      mpz_addmul (flll->coeff[i], U.coeff[k][i+2], gt->coeff[0]);
    if (i >= 1 && i <= d-1)
      mpz_addmul (flll->coeff[i], U.coeff[k][i+1], gt->coeff[1]);
  }
  flll->deg = d;
  /* force leading coefficient of f to be positive */
  if (mpz_sgn (flll->coeff[d]) < 0)
    mpz_poly_neg (flll, flll);
}

static inline uint64_t
get_LLL_uid (mat_Z U, int k, int d)
{
  uint64_t v1 = ((uint64_t) mpz_get_ui (U.coeff[k][1])) & 0x00000000FFFFFFFF;
  uint64_t v2 = ((uint64_t) mpz_get_ui (U.coeff[k][d])) & 0x00000000FFFFFFFF;
  return  (v1 << 32) | v2;
}

/* Print only the two coefficients of higher degree of a polynomial */
static inline void
mpz_poly_fprintf_short (FILE *out, mpz_poly_srcptr f)
{
  int d = f->deg;
  gmp_fprintf (out, "%Zd*x^%d + %Zd*x^%d%s\n", f->coeff[d], d,
                    f->coeff[d-1], d-1, (d > 1) ? " + ...":"");
}

/* Print poly: use mpz_poly_fprintf or mpz_poly_fprintf_short depending on
   verbose */
static inline void
mpz_poly_fprintf_verbose (FILE *out, mpz_poly_srcptr f, int verbose)
{
  if (verbose > 1)
    mpz_poly_fprintf (out, f);
  else
    mpz_poly_fprintf_short (out, f);
}

/******************************************************************************/
/******************************************************************************/
/********************** End Internal functions ********************************/
/******************************************************************************/

/* Use rotation and translation to find a polynomial with smaller norm.
   Use local descent algorithm to find a local minimum.
   Maximum degree for the rotation is given by deg_rotation (deg_rotation < 0
   means no rotation).
   Translations are used only if use_translation is non-zero.
   Will look for polynomials of the form
   [f_raw + (k[deg_rotation]*x^deg_rotation + ... + k[1]*x + k[0])*g_raw](x+kt)
   and g(x+kt)
   Optimized polynomial are return in f_opt and g_opt, and the function returns
   the skew lognorm of f_opt.

   To use only translation: use_translation = 1 and deg_rotation = -1
   To use only rotation   : use_transaltion = 0 and deg_rotation >= 0
   To use both            : use_translation = 1 and deg_rotation >= 0
   TODO: _mp version like old optimize_aux_mp in auxiliary.c
*/
double
sopt_local_descent (mpz_poly_ptr f_opt, mpz_poly_ptr g_opt,
                    mpz_poly_srcptr f_raw, mpz_poly_srcptr g_raw,
                    int use_translation, int deg_rotation,
                    unsigned int max_iter, int verbose)
{
  mpz_t kt, k[SOPT_MAX_DEGREE_ROTATION+1]; /* current offset */
  mpz_t tmp;
  int changedt, changed[SOPT_MAX_DEGREE_ROTATION+1];
  double logmu_opt, logmu;
  unsigned int iter = 0;
  mpz_poly_t ftmp;

  ASSERT_ALWAYS(deg_rotation <= SOPT_MAX_DEGREE_ROTATION);

  logmu_opt = L2_skew_lognorm ((mpz_poly_ptr) f_raw, SKEWNESS_DEFAULT_PREC);
  mpz_init (tmp);
  mpz_poly_init (ftmp, f_raw->deg);

  if (verbose)
  {
    fprintf (stderr, "# sopt:       starting local descent with polynomials:\n"
                     "# sopt:       f_raw = ");
    mpz_poly_fprintf_verbose (stderr, f_raw, verbose);
    fprintf (stderr, "# sopt:       g_raw = ");
    mpz_poly_fprintf_verbose (stderr, g_raw, verbose);
    fprintf (stderr, "# sopt:       with lognorm = %f\n", logmu_opt);
  }

  /* initialize k[i] to the smallest power of two that increases the lognorm
     by SOPT_LOCAL_DESCENT_GUARD with respect to the initial polynomial, to
     avoid being stuck in a very flat region */
  for (int i = 0; i <= deg_rotation; i++)
  {
    mpz_init_set_ui (k[i], 1);
    while (1)
    {
      mpz_poly_rotation (ftmp, f_raw, g_raw, k[i], i);
      logmu = L2_skew_lognorm (ftmp, SKEWNESS_DEFAULT_PREC);
      if (logmu > logmu_opt + SOPT_LOCAL_DESCENT_GUARD)
      {
        mpz_neg (tmp, k[i]);
        mpz_poly_rotation (ftmp, f_raw, g_raw, tmp, i);
        logmu = L2_skew_lognorm (ftmp, SKEWNESS_DEFAULT_PREC);
        if (logmu > logmu_opt + SOPT_LOCAL_DESCENT_GUARD)
          break;
      }
    mpz_mul_2exp (k[i], k[i], 1);
    }
  }

  /* initialize kt likewise */
  mpz_init_set_ui (kt, 1);
  while (1)
  {
    mpz_poly_translation (ftmp, f_raw, kt);
    logmu = L2_skew_lognorm (ftmp, SKEWNESS_DEFAULT_PREC);
    if (logmu > logmu_opt + SOPT_LOCAL_DESCENT_GUARD)
    {
      mpz_neg (tmp, kt);
      mpz_poly_translation (ftmp, f_raw, tmp);
      logmu = L2_skew_lognorm (ftmp, SKEWNESS_DEFAULT_PREC);
      if (logmu > logmu_opt + SOPT_LOCAL_DESCENT_GUARD)
        break;
    }
    mpz_mul_2exp (kt, kt, 1);
  }

  mpz_poly_set (f_opt, f_raw);
  mpz_poly_set (g_opt, g_raw);

  /* abort in cases where the descent procedure takes too long to converge;
   * 300 iterations seems largely enough in most cases */
  while (iter < max_iter)
  {
    iter++;
    changedt = 0;
    if (deg_rotation >= 0)
      memset (changed, 0, (deg_rotation+1) * sizeof (int));

    if (use_translation != 0)
    {
      /* first try translation by kt */
      mpz_poly_translation (ftmp, f_opt, kt); /* f(x+kt) */
      logmu = L2_skew_lognorm (ftmp, SKEWNESS_DEFAULT_PREC);
      if (logmu < logmu_opt)
      {
        changedt = 1;
        logmu_opt = logmu;
        mpz_poly_swap (f_opt, ftmp);
        mpz_poly_translation (g_opt, g_opt, kt);
      }
      else
      {
        mpz_neg (tmp, kt);
        mpz_poly_translation (ftmp, f_opt, tmp); /* f(x-kt) */
        logmu = L2_skew_lognorm (ftmp, SKEWNESS_DEFAULT_PREC);
        if (logmu < logmu_opt)
        {
          changedt = 1;
          logmu_opt = logmu;
          mpz_poly_swap (f_opt, ftmp);
          mpz_poly_translation (g_opt, g_opt, tmp);
        }
      }
    }

    for (int i = 0; i <= deg_rotation; i++)
    {
      mpz_poly_rotation (ftmp, f_opt, g_opt, k[i], i); /* f + k[i]*x^i*g */
      logmu = L2_skew_lognorm (ftmp, SKEWNESS_DEFAULT_PREC);
      if (logmu < logmu_opt)
      {
        changed[i] = 1;
        logmu_opt = logmu;
        mpz_poly_swap (f_opt, ftmp);
      }
      else
      {
        mpz_neg (tmp, k[i]);
        mpz_poly_rotation (ftmp, f_opt, g_opt, tmp, i); /* f - k[i]*x^i*g */
        logmu = L2_skew_lognorm (ftmp, SKEWNESS_DEFAULT_PREC);
        if (logmu < logmu_opt)
        {
          changed[i] = 1;
          logmu_opt = logmu;
          mpz_poly_swap (f_opt, ftmp);
        }
      }
    }

    int is_finished = 1;
    if (use_translation)
      is_finished = (changedt == 0 && mpz_cmp_ui (kt, 1) == 0);
    for (int i = 0; (i <= deg_rotation) && is_finished; i++)
      is_finished = (changed[i] == 0 && mpz_cmp_ui (k[i], 1) == 0);

    if (is_finished)
      break;

    if (changedt == 1)
      mpz_mul_2exp (kt, kt, 1);       /* kt <- 2*kt */
    else if (mpz_cmp_ui (kt, 1) > 0)
      mpz_div_2exp (kt, kt, 1);       /* kt <- kt/2 */

    for (int i = 0; i <= deg_rotation; i++)
    {
      if (changed[i] == 1)
        mpz_mul_2exp (k[i], k[i], 1);       /* k[i] <- 2*k[i] */
      else if (mpz_cmp_ui (k[i], 1) > 0)
        mpz_div_2exp (k[i], k[i], 1);       /* k[i] <- k[i]/2 */
    }
  }

  if (verbose)
  {
    fprintf (stderr, "# sopt:       end of local descent after %u iterations "
                     "with polynomials:\n" "# sopt:       f_opt = ", iter);
    mpz_poly_fprintf_verbose (stderr, f_opt, verbose);
    fprintf (stderr, "# sopt:       g_opt = ");
    mpz_poly_fprintf_verbose (stderr, g_opt, verbose);
    fprintf (stderr, "# sopt:       with lognorm = %f\n", logmu_opt);
  }

  for (int i = 0; i <= deg_rotation; i++)
    mpz_clear (k[i]);
  mpz_clear (kt);
  mpz_clear (tmp);
  mpz_poly_clear (ftmp);

  return logmu_opt;
}


#define SOPT_INIT_SIZE_ALLOCATED_TRANSLATIONS 1024



/* Size optimize the polynomial pair (f_raw, g_raw) with rotations and
   translations.
   Assume that deg(g) = 1.
   Return the size-optimized pair (f_opt, g_opt) and the skew lognorm of f_opt.
   The sopt_effort parameter is used to increase the number of translations that
   are considered. sopt_effort = 0 means that we consider only translations
   found by sopt_find_translation_deg* and the translation k = 0.
   TODO: LLL on gram matrix to be faster
   TODO: precompute skew^i for i in [0..d]
*/
double
size_optimization (mpz_poly_ptr f_opt, mpz_poly_ptr g_opt,
                   mpz_poly_srcptr f_raw, mpz_poly_srcptr g_raw,
                   const unsigned int sopt_effort, const int verbose)
{
  ASSERT_ALWAYS (g_raw->deg == 1);
  ASSERT_ALWAYS (sopt_effort <= SOPT_MAX_EFFORT);
  const int d = f_raw->deg;
  double lognorm_opt =
      L2_skew_lognorm ((mpz_poly_ptr) f_raw, SKEWNESS_DEFAULT_PREC);


  if (verbose)
  {
    fprintf (stderr, "\n# sopt: start size-optimization with polynomials:\n"
                     "# sopt: f_raw = ");
    mpz_poly_fprintf_verbose (stderr, f_raw, verbose);
    fprintf (stderr, "# sopt: g_raw = ");
    mpz_poly_fprintf_verbose (stderr, g_raw, verbose);
    fprintf (stderr, "# sopt: with lognorm = %f\n", lognorm_opt);
  }

  /****************************** init *******************************/
  mpz_t tmp, tmp2, list_skew[SOPT_NB_OF_SKEWNESS_VALUES], a, b;
  mpz_poly_t ft, gt, flll, fld, gld;
  list_mpz_t list_k;
  mat_Z m, U;
  hash_table_uint64_t H;

  mpz_init (tmp);
  mpz_init (tmp2);
  for (unsigned int i = 0; i < SOPT_NB_OF_SKEWNESS_VALUES; i++)
    mpz_init (list_skew[i]);
  /* 1/4 < a/b < 1: the closer a/b is from 1, the better the reduction is.
     Some experiments suggest that increasing a/b does not yield better
     polynomials, thus we take the smallest possible value. */
  mpz_init_set_ui (a, 1);
  mpz_init_set_ui (b, 4);

  mpz_poly_init (ft, d);
  mpz_poly_init (gt, 1);
  mpz_poly_init (flll, d);
  mpz_poly_init (fld, d);
  mpz_poly_init (gld, 1);

  list_mpz_init (list_k, SOPT_INIT_SIZE_ALLOCATED_TRANSLATIONS);

  /*************** generate list of translations to try *****************/
  if (d == 6 || d == 5)
  {
    if (d == 6)
      sopt_find_translations_deg6 (list_k, f_raw, g_raw, verbose);
    else if (d == 5)
      sopt_find_translations_deg5 (list_k, f_raw, g_raw, verbose);
    if (verbose)
      fprintf (stderr, "# sopt: %" PRIu64 " values for the translations were "
                       "computed and added in list_k\n", list_k->len);

    /* Add 0 to the list of translations to try */
    mpz_set_ui (tmp, 0);
    list_mpz_append (list_k, tmp);
    if (verbose)
      fprintf (stderr, "# sopt: k = 0 was added list_k\n");

    /* Sort list_k by increasing order and remove duplicates */
    list_mpz_sort_and_remove_dup (list_k, verbose);
    if (verbose)
      fprintf (stderr, "# sopt: It remains %" PRIu64 " values after sorting and"
                       " removing duplicates\n", list_k->len);
  }
  else
  {
    /* Start with list_k = {0} */
    mpz_set_ui (tmp, 0);
    list_mpz_append (list_k, tmp);
    if (verbose)
      fprintf (stderr, "# sopt: Start with list_k = { 0 }\n");
  }

  /******* Improve list of translation (depending on sopt_effort) *******/
  if (sopt_effort > 0)
  {
    if (verbose)
      fprintf (stderr, "# sopt: sopt_effort = %u > 0. Try to improve list_k\n",
                        sopt_effort);
    improve_list_k (list_k, sopt_effort, verbose);
  }
  else if (verbose)
    fprintf (stderr, "# sopt: sopt_effort = 0. list_k is not modify.\n");


  /****************************** Main loop *******************************/
  /* For each value of k, call LLL on translated polynomial, and for each
     polynomial returned by LLL, call size_optimize_local_descent */

  /* Init hash table */
  hash_uint64_init (H, 2*SOPT_NB_OF_SKEWNESS_VALUES*SOPT_MAX_LLL_POLY_PROCESS);

  /* Compute values of skewness that are going to be used in LLL */
  sopt_list_skewness_values (list_skew, f_raw, g_raw, verbose);

  /* Init matrix for LLL. We consider rotation up to x^(d-2)*g(x), so m has d
     row vectors (d-1 for the rotation + 1 for the polynomial f) of d+1
     coefficients each. */
  LLL_init (&m, d, d+1);
  /* The transformation matrix is a square matrix of size d (because m has d
     row vectors). */
  LLL_init (&U, d, d);

  if (verbose)
    fprintf (stderr, "# sopt: Start processing all k in list_k of length "
                     "%" PRIu64 "\n", list_k->len);

  for (unsigned int i = 0; i < list_k->len; i++)
  {
    mpz_poly_translation (ft, f_raw, list_k->tab[i]);
    mpz_poly_translation (gt, g_raw, list_k->tab[i]);
    hash_uint64_reset (H);
    if (verbose)
    {
      gmp_fprintf (stderr, "# sopt: process k = %Zd\n# sopt:   translated "
                           "polynomials are:\n# sopt:   ft = ", list_k->tab[i]);
      mpz_poly_fprintf_verbose (stderr, ft, verbose);
      fprintf (stderr, "# sopt:   gt = ");
      mpz_poly_fprintf_verbose (stderr, gt, verbose);
    }

    for (unsigned int j = 0; j < SOPT_NB_OF_SKEWNESS_VALUES; j++)
    {
      LLL_set_matrix_from_polys (&m, ft, gt, list_skew[j], tmp);
      if (verbose)
        gmp_fprintf (stderr, "# sopt:   calling LLL with skew = %Zd\n",
                             list_skew[j]);

      LLL (tmp, m, &U, a, b);

      int npoly = 0;
      for (int k = 1; k <= m.NumRows && npoly < SOPT_MAX_LLL_POLY_PROCESS; k++)
      {
        /* We want a polynomial that is not a multiple of gt. */
        if (mpz_sgn (U.coeff[k][1]) != 0)
        {
          npoly++;
          LLL_set_poly_from_vector (flll, U, k, ft, gt);
          if (verbose)
          {
            fprintf (stderr, "# sopt:     LLL return the following polynomial "
                             "of degree %d:\n# sopt:       flll = ", d);
            mpz_poly_fprintf_verbose (stderr, flll, verbose);
          }

          if (hash_uint64_insert (H, get_LLL_uid (U, k, d)))
          {
            double lognorm = sopt_local_descent (fld, gld, flll, gt, 1, d-2,
                                             SOPT_DEFAULT_MAX_STEPS, verbose);
            if (lognorm < lognorm_opt)
            {
              if (verbose)
              {
                fprintf (stderr, "# sopt:       better lognorm = %f (previous "
                                 "was %f)\n", lognorm, lognorm_opt);
              }
              lognorm_opt = lognorm;
              mpz_poly_set (f_opt, fld);
              mpz_poly_set (g_opt, gld);
            }
            else if (verbose)
            {
              fprintf (stderr, "# sopt:       lognorm = %f (better is %f)\n",
                               lognorm, lognorm_opt);
            }
          }
          else if (verbose)
            fprintf (stderr, "# sopt:       local descent already run on a "
                             "similar polynomial\n");
        }
      }
    }
  }

  /************************** Clear everything ***************************/
  LLL_clear (&m);
  LLL_clear (&U);
  for (unsigned int i = 0; i < SOPT_NB_OF_SKEWNESS_VALUES; i++)
    mpz_clear (list_skew[i]);
  hash_uint64_clear (H);
  list_mpz_clear (list_k);
  mpz_poly_clear (flll);
  mpz_poly_clear (fld);
  mpz_poly_clear (gld);
  mpz_poly_clear (ft);
  mpz_poly_clear (gt);
  mpz_clear (a);
  mpz_clear (b);
  mpz_clear (tmp);
  mpz_clear (tmp2);

  if (verbose)
  {
    fprintf (stderr, "# sopt: end of size-optimization, best polynomial are\n");
    fprintf (stderr, "# sopt: f_opt = ");
    mpz_poly_fprintf_verbose (stderr, f_opt, verbose);
    fprintf (stderr, "# sopt: g_opt = ");
    mpz_poly_fprintf_verbose (stderr, g_opt, verbose);
    fprintf (stderr, "# sopt: with lognorm = %f\n\n", lognorm_opt);
  }

  return lognorm_opt;
}
