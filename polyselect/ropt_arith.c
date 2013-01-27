/**
 * @file ropt_arith.c
 * Some arithmetics used in ropt.
 */


#include "cado.h"
#include "ropt_arith.h"


/**
 * Solve x in a + b*x = c (mod p)
 */
unsigned long
solve_lineq ( unsigned long a,
              unsigned long b,
              unsigned long c,
              unsigned long p )
{
  /* in general, we should know that gcd(b, p) = 1 */
  if (b % p == 0) {
    fprintf (stderr, "Error, impossible inverse in solve_lineq().\n");
    exit (1);
  }

  unsigned long tmp;
  modulusul_t mod;
  residueul_t tmpr, ar, cr;
  modul_initmod_ul (mod, p);
  modul_init (tmpr, mod);
  modul_init (cr, mod);
  modul_init (ar, mod);
  modul_set_ul (cr, c, mod);
  modul_set_ul (ar, a, mod);
  modul_sub (cr, cr, ar, mod);
  modul_set_ul (tmpr, b, mod);
  modul_inv (tmpr, tmpr, mod);
  modul_mul (tmpr, tmpr, cr, mod);
  tmp = modul_get_ul(tmpr, mod);
  modul_clear (tmpr, mod);
  modul_clear (cr, mod);
  modul_clear (ar, mod);
  modul_clearmod (mod);
  return tmp;
}


/**
 * Change coordinate from (a, b) to (u, v),
 * where A + MOD*a = u.
 */
void
ab2uv ( mpz_t A,
        mpz_t MOD,
        long a,
        mpz_t u )
{
  mpz_mul_si (u, MOD, a);
  mpz_add (u, u, A);
}


/**
 * Change coordinate from (a, b) to the index of
 * the sieving array, where index = a - Amin,
 * where Amin is negative.
 */
long
ab2ij ( long Amin,
        long a )
{
  return ( a - Amin );
}


/**
 * Change coordinate from (i, j) to (a, b).
 */
long
ij2ab ( long Amin,
        long i )
{
  return ( i + Amin );
}


/**
 * Change coordinate from (i, j) to (u, v).
 */
void
ij2uv ( mpz_t A,
        mpz_t MOD,
        long Amin,
        long i,
        mpz_t u )
{
  ab2uv(A, MOD, ij2ab(Amin, i), u);
}


/**
 * Find coordinate a such that
 * A + MOD*a = u (mod p).
 */
unsigned int
uv2ab_mod ( mpz_t A,
            mpz_t MOD,
            unsigned int U,
            unsigned int p )
{
  unsigned long a = mpz_fdiv_ui (A, p);
  unsigned long mod = mpz_fdiv_ui (MOD, p);
  unsigned long u = U % p;
  /* compute the A + MOD * a = u (mod p) */
  return (unsigned int) solve_lineq(a, mod, u, p);
}


/**
 * Same as above, but return the
 * position of a in the array.
 */
long
uv2ij_mod ( mpz_t A,
            long Amin,
            mpz_t MOD,
            unsigned int U,
            unsigned int p )
{
  long i = (long) uv2ab_mod (A, MOD, U, p);

  /* smallest k + p*i such that A0 < k + p*i, where A0 < 0,
     hence i = ceil((A0-tmp)/p). Note, this should be negative. */
  i = (long) ceil (((double) Amin - (double) i) / (double) p)
    * (long) p + i;

  /* compute the position of this (u, v) in the array. */
  i = ab2ij (Amin, i);

  return i;
}


/**
 * Given numerator = f(x)*g'(x) - f'(x)*g(x) and g(x),
 * compute special_u = numerator/g^2(x) (mod p).
 */
static inline unsigned long
compute_special_u ( unsigned long numerator,
                    unsigned long gx,
                    unsigned long p )
{
  modulusul_t mod;
  residueul_t tmp, tmp1;
  unsigned long special_u;
  modul_initmod_ul (mod, p);
  modul_init (tmp, mod);

  modul_set_ul (tmp1, numerator, mod);
  gx = gx * gx; /* This should be ok since gx < p < 2000 */
  modul_set_ul (tmp, gx, mod);

  modul_inv (tmp, tmp, mod);
  modul_mul (tmp, tmp, tmp1, mod);
  special_u = modul_get_ul (tmp, mod);

  modul_clear (tmp, mod);
  modul_clear (tmp1, mod);
  modul_clearmod (mod);

  return special_u;
}


/**
 * Similar to above, but test whether r is a multiple root.
 * Given numerator = f(x)*g'(x) - f'(x)*g(x) and g(x),
 * test  u*g^2(x) == numerator (mod p).
 */
static inline unsigned long
test_special_u ( unsigned long numerator,
                 unsigned long gx,
                 unsigned long u,
                 unsigned long p )
{
  modulusul_t mod;
  residueul_t tmp, tmp1;
  unsigned long re;
  modul_initmod_ul (mod, p);
  modul_init (tmp, mod);
  modul_init (tmp1, mod);

  gx = gx * gx; /* This should be ok since gx < p < 2000 */
  modul_set_ul (tmp, gx, mod);
  modul_set_ul (tmp1, u, mod);
  modul_mul (tmp, tmp, tmp1, mod);
  re = (modul_get_ul (tmp, mod) == numerator);

  modul_clear (tmp, mod);
  modul_clear (tmp1, mod);
  modul_clearmod (mod);
  return re;
}


/**
 * Compute fuv = f+(u*x+v)*g,
 * f(r) + u*r*g(r) + v*g(r) = 0
 * The inputs for f and g are mpz.
 */
void
compute_fuv_mp ( mpz_t *fuv,
                 mpz_t *f,
                 mpz_t *g,
                 int d,
                 mpz_t u,
                 mpz_t v )
{
  mpz_t tmp, tmp1;
  mpz_init (tmp);
  mpz_init (tmp1);
  int i = 0;

  for (i = 3; i <= d; i ++)
    mpz_set (fuv[i], f[i]);

  /* f + u*g1*x^2
     + (g0*u* + v*g1)*x
     + v*g0 */

  /* Note, u, v are signed long! */
  /* u*g1*x^2 */
  mpz_mul (tmp, g[1], u);
  mpz_add (fuv[2], f[2], tmp);

  /* (g0*u* + v*g1)*x */
  mpz_mul (tmp, g[0], u);
  mpz_mul (tmp1, g[1], v);
  mpz_add (tmp, tmp, tmp1);
  mpz_add (fuv[1], f[1], tmp);

  /* v*g0 */
  mpz_mul (tmp, g[0], v);
  mpz_add (fuv[0], f[0], tmp);

  mpz_clear (tmp);
  mpz_clear (tmp1);
}


/**
 * Compute fuv = f+(u*x+v)*g,
 * f(r) + u*r*g(r) + v*g(r) = 0
 * The inputs for f and g are mpz.
 */
static inline void
compute_fuv ( mpz_t *fuv,
              mpz_t *f,
              mpz_t *g,
              int d,
              long u,
              long v)
{
  mpz_t tmp, tmp1;
  mpz_init (tmp);
  mpz_init (tmp1);
  int i = 0;

  for (i = 3; i <= d; i ++)
    mpz_set (fuv[i], f[i]);

  /* f + u*g1*x^2
     + (g0*u* + v*g1)*x
     + v*g0 */

  /* Note, u, v are signed long! */
  /* u*g1*x^2 */
  mpz_mul_si (tmp, g[1], u);
  mpz_add (fuv[2], f[2], tmp);

  /* (g0*u* + v*g1)*x */
  mpz_mul_si (tmp, g[0], u);
  mpz_mul_si (tmp1, g[1], v);
  mpz_add (tmp, tmp, tmp1);
  mpz_add (fuv[1], f[1], tmp);

  /* v*g0 */
  mpz_mul_si (tmp, g[0], v);
  mpz_add (fuv[0], f[0], tmp);

  mpz_clear (tmp);
  mpz_clear (tmp1);
}


/**
 * Compute fuv = f+(u*x+v)*g,
 * The inputs for f and g are unsigne long.
 * Note, u, v are unsigned int.
 * So they should be reduce (mod p) if necessary.
 */
void
compute_fuv_ui ( unsigned int *fuv_ui,
                 unsigned int *f_ui,
                 unsigned int *g_ui,
                 int d,
                 unsigned int u,
                 unsigned int v,
                 unsigned int p )
{
  int i;
  modulusul_t mod;
  residueul_t tmp, tmp1, tmp2;
  modul_initmod_ul (mod, p);
  modul_init (tmp, mod);
  modul_init (tmp1, mod);
  modul_init (tmp2, mod);

  for (i = 3; i <= d; i ++)
    fuv_ui[i] = f_ui[i];

  /* f + u*g1*x^2
     + (g0*u* + v*g1)*x
     + v*g0 */

  /* u*g1*x^2 */
  modul_set_ul (tmp, g_ui[1], mod);
  modul_set_ul (tmp2, u, mod);
  modul_mul (tmp, tmp, tmp2, mod);
  modul_set_ul (tmp1, f_ui[2], mod);
  modul_add (tmp, tmp, tmp1, mod);
  fuv_ui[2] = (unsigned int) modul_get_ul(tmp, mod);

  /* (g0*u* + v*g1)*x */
  modul_set_ul (tmp, g_ui[1], mod);
  modul_set_ul (tmp1, v, mod);
  modul_mul (tmp, tmp, tmp1, mod);
  modul_set_ul (tmp1, g_ui[0], mod);
  // tmp2 = u as set above.
  modul_mul (tmp1, tmp1, tmp2, mod);
  modul_add (tmp, tmp, tmp1, mod);
  modul_set_ul (tmp1, f_ui[1], mod);
  modul_add (tmp, tmp, tmp1, mod);
  fuv_ui[1] = (unsigned int) modul_get_ul(tmp, mod);

  /* v*g0 */
  modul_set_ul (tmp1, v, mod);
  modul_set_ul (tmp2, g_ui[0], mod);
  modul_mul (tmp1, tmp1, tmp2, mod);
  modul_set_ul (tmp, f_ui[0], mod);
  modul_add (tmp, tmp, tmp1, mod);
  fuv_ui[0] = (unsigned int) modul_get_ul(tmp, mod);
}


/**
 * Compute v (mod p) by
 * f(r) + u*r*g(r) + v*g(r) = 0 (mod p).
 * The inputs for f(r) and g(r) are unsigned int.
 */
unsigned int
compute_v_ui ( unsigned int fx,
               unsigned int gx,
               unsigned int r,
               unsigned int u,
               unsigned int p)
{
  modulusul_t mod;
  residueul_t tmp, tmp1;
  unsigned long v;
  modul_initmod_ul (mod, p);
  modul_init (tmp, mod);
  modul_init (tmp1, mod);

  /* g(r)*r*u + f(r) */
  modul_set_ul (tmp, gx, mod);
  modul_set_ul (tmp1, r, mod);
  modul_mul (tmp, tmp, tmp1, mod);
  modul_set_ul (tmp1, u, mod);
  modul_mul (tmp, tmp, tmp1, mod);
  modul_set_ul (tmp1, fx, mod);
  modul_add (tmp, tmp, tmp1, mod);
  v = modul_get_ul(tmp, mod);

  /* solve v in tmp2 + v*g(r) = 0 (mod p) */
  v = solve_lineq (v, gx, 0, p);

  modul_clear (tmp, mod);
  modul_clear (tmp1, mod);
  modul_clearmod (mod);
  return (unsigned int) v;
}


/**
 * Compute v = f(r) (mod pe), where f is of degree d.
 * The input f should be unsigned int.
 */
unsigned int
eval_poly_ui_mod ( unsigned int *f,
                   int d,
                   unsigned int r,
                   unsigned int pe )
{
  int i;
  modulusul_t mod;
  residueul_t vtmp, rtmp, tmp;
  unsigned int v;

  modul_initmod_ul (mod, pe);
  modul_init (vtmp, mod);
  modul_init (rtmp, mod);
  modul_init (tmp, mod);

  /* set vtmp = f[d] (mod p) and rtmp = r (mod p) */
  modul_set_ul (vtmp, f[d], mod);
  modul_set_ul (rtmp, r, mod);

  for (i = d - 1; i >= 0; i--) {
    modul_mul (vtmp, vtmp, rtmp, mod);
    modul_set_ul (tmp, f[i], mod);
    modul_add (vtmp, tmp, vtmp, mod);
  }

  v = (unsigned int) modul_get_ul (vtmp, mod);
  modul_clear (vtmp, mod);
  modul_clear (rtmp, mod);
  modul_clear (tmp, mod);
  modul_clearmod (mod);

  return v;
}


/**
 * Compute v = f'(r) (mod pe), where f is of degree d.
 * The input f should be unsigned long.
 */
static inline unsigned long
eval_poly_diff_ui_mod ( unsigned long *f,
                        int d,
                        unsigned long r,
                        unsigned long pe )
{
  int i;
  modulusul_t mod;
  residueul_t vtmp, rtmp, tmp, itmp;
  unsigned long v;

  modul_initmod_ul (mod, pe);
  modul_init (vtmp, mod);
  modul_init (rtmp, mod);
  modul_init (itmp, mod);
  modul_init (tmp, mod);

  /* set vtmp = f[d] (mod p) and rtmp = r (mod p) */
  modul_set_ul (vtmp, f[d], mod);
  modul_set_ul (tmp, (unsigned long) d, mod);
  modul_mul (vtmp, vtmp, tmp, mod);
  modul_set_ul (rtmp, r, mod);

  for (i = d - 1; i >= 1; i--) {
    modul_mul (vtmp, vtmp, rtmp, mod);
    /* vtmp <- vtmp + i*f[i] */
    modul_set_ul (tmp, f[i], mod);
    modul_set_ul (itmp, i, mod);
    modul_mul (tmp, itmp, tmp, mod);
    modul_add (vtmp, tmp, vtmp, mod);
  }

  v = modul_get_ul (vtmp, mod);
  modul_clear (vtmp, mod);
  modul_clear (itmp, mod);
  modul_clear (rtmp, mod);
  modul_clear (tmp, mod);
  modul_clearmod (mod);

  return v;
}


/**
 * Test whether r is a root for f(r) + g(r)*(u*x+v) = 0 (mod p^e);
 * Note that, the condition to call this function is that, r is a root
 * for f(r) + g(r)*(u*x+v) = 0 (mod p^(e-1)).
 *
 * (1) If it can be lifted, then it could be single or multiple root
 * -- single root return 1, -- note, the lifted root is not r in general.
 * Hence, we need to compute the lifted root and save it in r_lifted.
 * -- multiple root return 2.
 * (2) If it is can not be lifted, return 0;
 */
static inline int
isroot_fuv ( mpz_t *f,
             mpz_t *g,
             int d,
             unsigned long u,
             unsigned long v,
             unsigned long r,
             unsigned long pe,
             unsigned long p,
             unsigned long *r_lifted )
{
  mpz_t tmp, tmp2, fr, gr;
  int type;
  mpz_init (fr);
  mpz_init (gr);
  mpz_init (tmp);
  mpz_init (tmp2);

  /* r is a root over p^(e-1), we need to see whether it is a
     single or multiple;
     - If it is a single, computer r_lifted
     - If it is a multiple, see whether it can be lifted
     -- consider f(r) (mod p^e)	 */

  /* f'(x)
     + u*g(x)
     + g'(x)*(ux + v) */
  mpz_set_ui (tmp, r);
  mpz_mul_ui (tmp, tmp, u);
  mpz_add_ui (tmp, tmp, v);
  mpz_mul (tmp, tmp, g[1]);
  eval_poly_ui (gr, g, 1, r);
  mpz_mul_ui (gr, gr, u);
  mpz_add (gr, gr, tmp);
  /* keep fr unchanged for further reference */
  eval_poly_diff_ui (fr, f, d, r);
  mpz_add (tmp, fr, gr); // tmp is f'_uv(r)

  /* cmpute f(x) + g(x)*(ux + v) */
  mpz_set_ui (tmp2, r);
  mpz_mul_ui (tmp2, tmp2, u);
  mpz_add_ui (tmp2, tmp2, v);
  eval_poly_ui (gr, g, 1, r);
  mpz_mul (gr, gr, tmp2);
  eval_poly_ui (fr, f, d, r);
  mpz_add (fr, fr, gr);  // fr is f_uv(r)

  /* if f_uv'(r) = 0 (mod p), then it is a double root */
  if ( mpz_fdiv_ui (tmp, p) == 0) {
    if (mpz_fdiv_ui (fr, pe) == 0)
      type = 2;
    else
      type = 0;
  }
  /* r is a single root, computer r_lifted. */
  else {
    unsigned long tmp3;
    modulusul_t mod;
    residueul_t res1, res2;

    /* 1/f_uv'(r) (mod p) */
    tmp3 = mpz_fdiv_r_ui (tmp, tmp, p);
    modul_initmod_ul (mod, p);
    modul_init (res1, mod);
    modul_set_ul (res1, tmp3, mod);
    modul_inv (res1, res1, mod);
    /* -f_uv(r)/p^(e-1) (mod p) */
    mpz_fdiv_q_ui (tmp2, fr, pe/p);
    mpz_neg (tmp2, tmp2);
    tmp3 = mpz_fdiv_r_ui (tmp2, tmp2, p);
    modul_init (res2, mod);
    modul_set_ul (res2, tmp3, mod);
    modul_mul (res1, res1, res2, mod);

    (*r_lifted) = modul_get_ul(res1, mod);
    (*r_lifted) = (*r_lifted) * (pe/p) + r;
    //printf ("lifted_r: %lu\n", *r_lifted);
    type = 1;

    modul_clear (res1, mod);
    modul_clear (res2, mod);
    modul_clearmod (mod);
  }

  mpz_clear (fr);
  mpz_clear (gr);
  mpz_clear (tmp);
  mpz_clear (tmp2);
  return type;
}


/**
 * Given f (mod pe) and (u, v) (mod pe) pair.
 * Lift a single root r (mod pe/p) to r (mod p)
 * Please ensure that f_{u, v}(r) = 0 (mod pe/p)
 * when call this function.
 */
static inline void
liftroot_fuv_ul ( unsigned long *f,
                  int d,
                  unsigned long r,
                  unsigned long p,
                  unsigned long pe,
                  unsigned long *r_lifted )
{
  modulusul_t mod;
  residueul_t tmp, tmp1;
  unsigned long fr;
  modul_initmod_ul (mod, p);
  modul_init (tmp, mod);
  modul_init (tmp1, mod);

  /* f'(x) */
  fr = eval_poly_diff_ui_mod (f, d, r, p);
  modul_set_ul (tmp1, fr, mod);
  /* tmp should be kept unchanged from now on */
  modul_add (tmp, tmp, tmp1, mod); 

  /* 1/f'(r) (mod p) */
  modul_inv (tmp, tmp, mod);

  /* -f(r)/p^(e-1) (mod p) */
  fr = fr / (pe/p);
  modul_set_ul (tmp1, fr, mod);
  modul_neg (tmp1, tmp1, mod);
  modul_mul (tmp, tmp, tmp1, mod);

  (*r_lifted) = modul_get_ul(tmp, mod);
  (*r_lifted) = (*r_lifted) * (pe/p) + r;

  modul_clear (tmp1, mod);
  modul_clear (tmp, mod);
  modul_clearmod (mod);
}


/**
 * Reduce mpz_t *f to unsigned int *f_mod;
 * Given modulus pe, return f (mod pe).
 */
inline void
reduce_poly_ul ( unsigned int *f_ui,
                 mpz_t *f,
                 int d,
                 unsigned int pe )
{
  int i;
  for (i = 0; i <= d; i ++) {
    f_ui[i] = (unsigned int) mpz_fdiv_ui (f[i], pe);
  }
}


/**
 * From polyselect.c
 * Implements Lemma 2.1 from Kleinjung's paper.
 * If a[d] is non-zero, it is assumed it is already set, otherwise it is
 * determined as a[d] = N/m^d (mod p).
 */
void
Lemma21 ( mpz_t *a,
          mpz_t N,
          int d,
          mpz_t p,
          mpz_t m )
{
  mpz_t r, mi, invp;
  int i;

  mpz_init (r);
  mpz_init (mi);
  mpz_init (invp);
  mpz_set (r, N);
  mpz_pow_ui (mi, m, d);
  if (mpz_cmp_ui (a[d], 0) == 0)
  {
    mpz_invert (a[d], mi, p); /* 1/m^d mod p */
    mpz_mul (a[d], a[d], N);
    mpz_mod (a[d], a[d], p);
  }
  for (i = d - 1; i >= 0; i--)
  {
    /* invariant: mi = m^(i+1) */
    mpz_mul (a[i], a[i+1], mi);
    mpz_sub (r, r, a[i]);
    ASSERT (mpz_divisible_p (r, p));
    mpz_divexact (r, r, p);
    mpz_divexact (mi, mi, m); /* now mi = m^i */
    if (i == d - 1)
    {
      mpz_invert (invp, p, mi); /* 1/p mod m^i */
      mpz_sub (invp, mi, invp); /* -1/p mod m^i */
    }
    else
      mpz_mod (invp, invp, mi);
    mpz_mul (a[i], invp, r);
    mpz_mod (a[i], a[i], mi); /* -r/p mod m^i */
    /* round to nearest in [-m^i/2, m^i/2] */
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
  mpz_clear (invp);
}
