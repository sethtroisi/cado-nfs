#ifndef _EC_ARITH_WEIERSTRASS_H_
#define _EC_ARITH_WEIERSTRASS_H_

#include "ec_arith_common.h"

/********************* Short Weierstrass elliptic curves **********************/

/* Equation:
 *    y^2 = x^3 + a*x + b                for affine
 *    Y^2*Z = X^3 + a*X*Z^Z + b*Z^3      for projective
 *
 * Curve coefficient needed in computation: a
 *
 * Implemented functions:
 *  - weierstrass_aff_add (R:SHORT_WEIERSTRASS_aff, P:SHORT_WEIERSTRASS_aff,
 *                                                      Q:SHORT_WEIERSTRASS_aff)
 *      R <- P+Q
 *
 *  - weierstrass_aff_dbl (R:SHORT_WEIERSTRASS_aff, P:SHORT_WEIERSTRASS_aff)
 *      R <- 2*P
 *
 *  - weierstrass_proj_add (R:SHORT_WEIERSTRASS_proj, P:SHORT_WEIERSTRASS_proj,
 *                                                    Q:SHORT_WEIERSTRASS_proj)
 *      R <- P+Q
 *
 *  - weierstrass_proj_dbl (R:SHORT_WEIERSTRASS_proj, P:SHORT_WEIERSTRASS_proj)
 *      R <- 2*P
 */


/* Computes R=2P, with 1 inv, 4 muls (2 muls and 2 squares) and 8 add/sub.
 *    - m : modulus
 *    - a : curve coefficient
 *
 * Return 1 if it worked, 0 if a modular inverse failed.
 * If modular inverse failed, return non-invertible value in R->x.
 *
 * It is permissible to let P and Q use the same memory.
 */
static int
weierstrass_aff_dbl (ec_point_t R, const ec_point_t P, const residue_t a,
                     const modulus_t m)
{
  int ret = 1;
  residue_t lambda, u, v;

  mod_init_noset0 (lambda, m);
  mod_init_noset0 (u, m);
  mod_init_noset0 (v, m);

  mod_sqr (u, P->x, m);
  mod_add (v, u, u, m);
  mod_add (v, v, u, m);
  mod_add (v, v, a, m); /* 3x^2 + a */
  mod_add (u, P->y, P->y, m);
  ret = mod_inv (lambda, u, m);    /* 1/(2*y) */
  if (ret != 0)
  {
    mod_mul (lambda, lambda, v, m);
    mod_sqr (u, lambda, m);
    mod_sub (u, u, P->x, m);
    mod_sub (u, u, P->x, m);    /* x3 = u = lambda^2 - 2*x */
    mod_sub (v, P->x, u, m);
    mod_mul (v, v, lambda, m);
    mod_sub (R->y, v, P->y, m);
    mod_set (R->x, u, m);
  }
  else
    mod_set (R->x, u, m);

  mod_clear (v, m);
  mod_clear (u, m);
  mod_clear (lambda, m);
  return ret;
}

/* Computes R=P+Q, with 1 inv, 3 muls (2 muls and 1 square) and 6 add/sub.
 *    - m : modulus
 *    - a : curve coefficient
 *
 * Return 1 if it worked, 0 if a modular inverse failed.
 * If modular inverse failed, return non-invertible value in R->x.
 *
 * It is permissible to let R and P (or R and Q) use the same memory.
 */
static int
weierstrass_aff_add (ec_point_t R, const ec_point_t P, const ec_point_t Q,
                     const residue_t a, const modulus_t m)
{
  residue_t lambda, u, v;
  int ret;

  mod_init_noset0 (u, m);
  mod_init_noset0 (v, m);
  mod_init_noset0 (lambda, m);

  mod_sub (u, Q->x, P->x, m);
  ret = mod_inv (v, u, m);
  if (ret != 0)
  {
    mod_sub (u, Q->y, P->y, m);
    mod_mul (lambda, u, v, m);
    mod_sqr (u, lambda, m);
    mod_sub (u, u, P->x, m);
    mod_sub (u, u, Q->x, m);    /* x3 = u = lambda^2 - P->x - Q->x */
    mod_sub (v, P->x, u, m);
    mod_mul (v, v, lambda, m);
    mod_sub (R->y, v, P->y, m);
    mod_set (R->x, u, m);
  }
  else if (mod_equal (P->x, Q->x, m) && mod_equal (P->y, Q->y, m))
    ret = weierstrass_aff_dbl (R, P, a, m);
  else
    mod_set (R->x, u, m);

  mod_clear (lambda, m);
  mod_clear (v, m);
  mod_clear (u, m);
  return ret;
}


/* Convert the Montgomery curve B*Y^2*Z = X^3 + A*X^2*Z + X*Z^2 with a valid
 * point Pm into a affine Weierstrass curve y^2 = x^3 + a*x + b with a
 * valid affine point Pw.
 *
 * Return 1 if it worked, 0 if a modular inverse failed.
 * If modular inverse failed, return non-invertible value in Pw->x.
 *
 * The curve coefficient b of the short Weierstrass curve will not be computed.
 * a and A can be the same variable.
 * Pm and Pw can be the same variable.
 */
static int
weierstrass_aff_from_montgomery (residue_t a, ec_point_t Pw, const residue_t A,
                                 ec_point_t Pm, const modulus_t m)
{
  residue_t B, one, t, x;
  int ret;

  mod_init_noset0 (B, m);
  mod_init_noset0 (one, m);
  mod_init_noset0 (t, m);
  mod_init_noset0 (x, m);
  mod_set1 (one, m);

  ret = mod_inv (t, Pm->z, m);
  if (ret == 0)
  {
    fprintf (stderr, "%s: could not invert Z\n", __func__);
    mod_set (Pw->x, Pm->z, m);
  }
  else
  {
    mod_mul (x, Pm->x, t, m); /* x = X/Z */
    mod_add (B, x, A, m);
    mod_mul (B, B, x, m);
    mod_add (B, B, one, m);
    mod_mul (B, B, x, m); /* B = x^3 + A*x^2 + x */

    /* Now (x,1) is on the curve B*y^2 = x^3 + A*x^2 + x. */
    ret = mod_inv (Pw->y, B, m);    /* y = 1/B */
    if (ret == 0)
    {
      fprintf (stderr, "%s: could not invert B\n", __func__);
      mod_set (Pw->x, B, m);
    }
    else
    {
      mod_div3 (t, A, m);
      mod_add (Pw->x, x, t, m);
      mod_mul (Pw->x, Pw->x, Pw->y, m); /* x = (X + A/3)/B */
      mod_mul (a, t, A, m);
      mod_sub (a, one, a, m);
      mod_mul (a, a, Pw->y, m);
      mod_mul (a, a, Pw->y, m);         /* a = (1 - (A^2)/3)/B^2 */
    }
  }

  mod_clear (one, m);
  mod_clear (B, m);
  mod_clear (t, m);
  mod_clear (x, m);

  return ret;
}


/* Computes P<-eP, with double-and-add algorithm.
 *    - m : modulus
 *    - a : curve coefficient
 *
 * Return 0 if e*P is the point at infinity, else return nonzero.
 * If the point at infinity is due to a failed inversion, the non-invertible
 * value is returned in P->x.
 */
static int
weierstrass_aff_smul_ui (ec_point_t P, const unsigned long e, const residue_t a,
                         const modulus_t m)
{
  unsigned long i;
  ec_point_t T;
  int tfinite; /* Nonzero iff T is NOT point at infinity */

  if (e == 0)
    return 0; /* signal point at infinity */

  ec_point_init (T, m);

  i = ~(0UL);
  i -= i/2;   /* Now the most significant bit of i is set */
  while ((i & e) == 0)
    i >>= 1;

  ec_point_set (T, P, m, SHORT_WEIERSTRASS_aff);
  tfinite = 1;
  i >>= 1;

  while (i > 0)
  {
    if (tfinite)
      tfinite = weierstrass_aff_dbl (T, T, a, m);
    if (e & i)
    {
      if (tfinite)
        tfinite = weierstrass_aff_add (T, T, P, a, m);
      else
      {
        ec_point_set (T, P, m, SHORT_WEIERSTRASS_aff);
        tfinite = 1;
      }
    }
    i >>= 1;
  }

  if (tfinite)
    ec_point_set (P, T, m, SHORT_WEIERSTRASS_aff);
  else
    mod_set (P->x, T->x, m);

  ec_point_clear (T, m);

  return tfinite;
}

static inline void
weierstrass_proj_point_set_zero (ec_point_t P, const modulus_t m)
{
  mod_set0 (P->x, m);
  mod_set1 (P->y, m);
  mod_set0 (P->z, m);
}
/* Computes R=2P, with ? muls (? muls and ? squares) and ? add/sub.
 *    - m : modulus
 *    - a : curve coefficient
 *
 * It is permissible to let P and Q use the same memory.
 */
static void
weierstrass_proj_dbl (ec_point_t R, const ec_point_t P, const residue_t a,
                      const modulus_t m)
{
  residue_t xx, zz, w, u, s, ss, sss, r, rr, B, t, h;
  /* TODO reduce number of var */

  mod_init_noset0 (xx, m);
  mod_init_noset0 (zz, m);
  mod_init_noset0 (w, m);
  mod_init_noset0 (u, m);
  mod_init_noset0 (s, m);
  mod_init_noset0 (ss, m);
  mod_init_noset0 (sss, m);
  mod_init_noset0 (r, m);
  mod_init_noset0 (rr, m);
  mod_init_noset0 (B, m);
  mod_init_noset0 (t, m);
  mod_init_noset0 (h, m);

  mod_sqr (xx, P->x, m);
  mod_sqr (zz, P->z, m);

  mod_mul (w, a, zz, m);
  mod_add (t, xx, xx, m);
  mod_add (t, t, xx, m);
  mod_add (w, w, t, m);

  mod_mul (s, P->y, P->z, m);
  mod_add (s, s, s, m);
  mod_sqr (ss, s, m);
  mod_mul (sss, ss, s, m);
  mod_mul (r, P->y, s, m);
  mod_sqr (rr, r, m);

  mod_add (B, P->x, r, m);
  mod_sqr (B, B, m);
  mod_sub (B, B, xx, m);
  mod_sub (B, B, rr, m);
  mod_sqr (h, w, m);
  mod_sub (h, h, B, m);
  mod_sub (h, h, B, m);

  mod_mul (R->x, h, s, m);
  mod_sub (R->y, B, h, m);
  mod_mul (R->y, R->y, w, m);
  mod_sub (R->y, R->y, rr, m);
  mod_sub (R->y, R->y, rr, m);
  mod_set (R->z, sss, m);

  mod_clear (xx, m);
  mod_clear (zz, m);
  mod_clear (w, m);
  mod_clear (u, m);
  mod_clear (s, m);
  mod_clear (ss, m);
  mod_clear (sss, m);
  mod_clear (r, m);
  mod_clear (rr, m);
  mod_clear (B, m);
  mod_clear (t, m);
  mod_clear (h, m);
}

/* Computes R=P+Q, with 14 muls (12 muls and 2 squares) and 7 add/sub.
 *    - m : modulus
 *    - a : curve coefficient
 *
 * It is permissible to let R and P (or R and Q) use the same memory.
 */
static void
weierstrass_proj_add (ec_point_t R, const ec_point_t P, const ec_point_t Q,
                      const modulus_t m)
{
  residue_t t1, t2, t3, u, uu, v, vv, vvv, r, A; /* TODO reduce number of var */

  mod_init_noset0 (t1, m);
  mod_init_noset0 (t2, m);
  mod_init_noset0 (t3, m);
  mod_init_noset0 (u, m);
  mod_init_noset0 (uu, m);
  mod_init_noset0 (v, m);
  mod_init_noset0 (vv, m);
  mod_init_noset0 (vvv, m);
  mod_init_noset0 (r, m);
  mod_init_noset0 (A, m);

  mod_mul (t1, P->x, Q->z, m);
  mod_mul (t2, P->y, Q->z, m);
  mod_mul (t3, P->z, Q->z, m);

  mod_mul (u, Q->y, P->z, m);
  mod_sub (u, u, t2, m);
  mod_sqr (uu, u, m);

  mod_mul (v, Q->x, P->z, m);
  mod_sub (v, v, t1, m);
  mod_sqr (vv, v, m);
  mod_mul (vvv, vv, v, m);

  mod_mul (r, vv, t1, m);

  mod_mul (A, uu, t3, m);
  mod_sub (A, A, vvv, m);
  mod_sub (A, A, r, m);
  mod_sub (A, A, r, m);

  mod_mul (R->x, v, A, m);

  mod_sub (R->y, r, A, m);
  mod_mul (R->y, R->y, u, m);
  mod_mul (t2, t2, vvv, m);
  mod_sub (R->y, R->y, t2, m);

  mod_mul (R->z, vvv, t3, m);

  mod_clear (t1, m);
  mod_clear (t2, m);
  mod_clear (t3, m);
  mod_clear (u, m);
  mod_clear (uu, m);
  mod_clear (v, m);
  mod_clear (vv, m);
  mod_clear (vvv, m);
  mod_clear (r, m);
  mod_clear (A, m);
}

/* Computes P<-eP, with double-and-add algorithm.
 *    - m : modulus
 *    - a : curve coefficient
 */
static void
weierstrass_proj_smul_ui (ec_point_t P, const unsigned long e,
                          const residue_t a, const modulus_t m)
{
  if (e == 0)
    weierstrass_proj_point_set_zero (P, m);
  else if (e > 1)
  {
    unsigned long i;
    ec_point_t T;

    ec_point_init (T, m);

    i = ~(0UL);
    i -= i/2;   /* Now the most significant bit of i is set */
    while ((i & e) == 0)
      i >>= 1;

    ec_point_set (T, P, m, SHORT_WEIERSTRASS_proj);
    i >>= 1; /* skip most significant bit of e */

    for (; i > 0; i >>= 1)
    {
      weierstrass_proj_dbl (T, T, a, m);
      if (e & i)
        weierstrass_proj_add (T, T, P, m);
    }

    ec_point_set (P, T, m, SHORT_WEIERSTRASS_proj);
    ec_point_clear (T, m);
  }
  /* else do nothing for e == 1 */
}

#endif /* _EC_ARITH_WEIERSTRASS_H_ */
