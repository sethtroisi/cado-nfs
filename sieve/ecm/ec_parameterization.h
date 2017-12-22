#include "ec_arith_common.h"
#include "ec_arith_Edwards.h"
#include "ec_arith_Montgomery.h"
#include "ec_arith_Weierstrass.h"

/******************************************************************************/
/*********************** Brent--Suyama parameterization ***********************/
/******************************************************************************/

/* Produces curve in Montgomery form.
 *
 * Rational parameterization (parameter is called sigma).
 *
 * Valid parameter: sigma in Q \ { 0, -1, 1, -3, 3, -5, 5, -5/3, 5/3 }
 *
 * Note: sigma and -sigma give isomorphic curves.
 *
 * Parameterization:
      u = sigma^2-5
      v = 4*sigma
      x0 = u^3
      y0 = (sigma^2-1)*(sigma^2-25)*(sigma^4-25)
      z0 = v^3
      A = (v-u)^3*(3*u+v)/(4*u^3*v) - 2
      B = u/z0
      b = (v-u)^3*(3*u+v)/(16*u^3*v)                          # [ b = (A+2)/4 ]
 * The point of order 3 is given by:
      x3 = u
      y3 = 2*(sigma^2-25)*(sigma^2-1)*sigma/u
      z3 = v
 * To check with Sage:
      QQsigma.<sigma> = QQ[]
      # copy paste all the above formula
      (b-(A+2)/4).is_zero() # b is correct ?
      (B*y0^2*z0 - (x0^3+A*x0^2*z0+x0*z0^2)).is_zero() # P0 is on the curve ?
      (B*y3^2*z3 - (x3^3+A*x3^2*z3+x3*z3^2)).is_zero() # P3 is on the curve ?
      (4*A*(x3/z3)^3+3*(x3/z3)^4+6*(x3/z3)^2-1).is_zero() # P3 of order 3 ?
 */


static inline int
ec_parameterization_Brent_Suyama_is_valid (const unsigned long sigma)
{
  if (sigma == 0 || sigma == 1 || sigma == 3 || sigma == 5)
    return 0; /* invalid values */
  else
    return 1;
}

/* Return 1 if it worked, 0 if a modular inverse failed.
 * If modular inverse failed, return non-invertible value in P0->x.
 *
 * For the computation, we only need b and the coordinates (x0::z0) of the
 * starting point.
 */
static int
ec_parameterization_Brent_Suyama (residue_t b, ec_point_t P0,
                                  const unsigned long sigma, const modulus_t m)
{
  ASSERT_ALWAYS (ec_parameterization_Brent_Suyama_is_valid (sigma));

  residue_t s, u, v, t1, t2;
  int ret;

  mod_init_noset0 (s, m);
  mod_init_noset0 (u, m);
  mod_init_noset0 (v, m);
  mod_init_noset0 (t1, m);
  mod_init_noset0 (t2, m);

  mod_set_ul (s, sigma, m);

  mod_add (v, s, s, m);
  mod_add (v, v, v, m);             /* v = 4*sigma */
  mod_sqr (u, s, m);
  mod_set_ul (t1, 5, m);
  mod_sub (u, u, t1, m);            /* u = sigma^2 - 5 */
  mod_sqr (t1, u, m);
  mod_mul (P0->x, t1, u, m);        /* x0 = u^3 */
  mod_sqr (t1, v, m);
  mod_mul (P0->z, t1, v, m);        /* z0 = v^3 */
  mod_mul (t1, P0->x, v, m);
  mod_2pow_ul (t2, 4, m);
  mod_mul (t1, t1, t2, m);          /* t1 = 16*x0*v = 16*u^3*v */
  mod_add (b, u, u, m);
  mod_add (b, b, u, m);
  mod_add (b, b, v, m);             /* b = 3*u+v (so far) */
  mod_sub (t2, v, u, m);
  mod_mul (b, b, t2, m);            /* b = (v-u)*(3*u+v) (so far) */
  mod_sqr (t2, t2, m);
  mod_mul (b, b, t2, m);            /* b = (v-u)^3*(3*u+v) (so far) */

  ret = mod_inv (t2, t1, m);         /* t2 = 1/t1 */
  if (ret == 0) /* non-trivial gcd */
    mod_set (P0->x, t1, m);
  else
    mod_mul (b, b, t2, m);          /* b = (v-u)^3*(3*u+v)/(16*u^3*v) */

  mod_clear (s, m);
  mod_clear (u, m);
  mod_clear (v, m);
  mod_clear (t1, m);
  mod_clear (t2, m);

  return ret;
}

/******************************************************************************/
/********************* Montgomery Z/12Z parameterization **********************/
/******************************************************************************/

/* Produces curve in Montgomery form.
 *
 * Elliptic parameterization (parameter is called k).
 *    E: y^2 = x^3 - 12*x
 *      rank = 1
 *      P = (-2, 4) is a point of infinite order
 *      Torsion points: P2 = (0:0:1) of order 2
 *
 * Reference: Montgomery's thesis (Section 6.1)
 *
 * Valid parameter: k in Z \ {-1, 0, 1}
 *
 * Note: Using P2 does not produce a valid curve.
 *       A point and its opposite produce the same curve
 *       Two points that differ by P2 produce the same curve (but not the same
 *        starting point)
 *
 * Parameterization (where (u,v) := k*P):
      d1 = v^2+12*u^2
      d2 = v^2-4*u^2
      x0 = 2*(u*(u^2+12))^2
      y0 = d1*u*(u^2+12)
      z0 = 2*d1*d2
      A = (4*(32*v*u^3)^2-2*d1*d2^3)/(d1*d2^3)
      B = (16*u^2*(v^2+4*u^2))^2/(d1*d2^3)
      b = (32*v*u^3)^2/(d1*d2^3)
 * The point of order 12 is given by:
      x12 = 2*u*(2*u + v)^2
      y12 = v*(2*u + v)^2
      z12 = 2*u*(-2*u + v)^2
 * To check with Sage:
      QQuv.<u,v> = QQ[]
      I = QQuv.ideal (v^2 - (u^3-12*u))
      # copy paste all the above formula
      (b-(A+2)/4).is_zero() # b is correct ?
      # P0 is on the curve ?
      eq_P0 = (B*y0^2*z0 - (x0^3+A*x0^2*z0+x0*z0^2))
      eq_P0.numerator() in I and not eq_P0.denominator() in I
      # P12 is on the curve ?
      eq_P12 = (B*y12^2*z12 - (x12^3+A*x12^2*z12+x12*z12^2))
      eq_P12.numerator() in I and not eq_P12.denominator() in I
      # TODO check that 4*P12 is of order 3 and 3*P12 is of order 4
 */

static inline int
ec_parameterization_Montgomery12_is_valid (const unsigned long k)
{
  if (k <= 1)
    return 0; /* invalid values */
  else
    return 1;
}

/* Return 1 if it worked, 0 if a modular inverse failed.
 * If modular inverse failed, return non-invertible value in P0->x.
 *
 * For the computation, we only need b and the coordinates (x0::z0) of the
 * starting point.
 */
static int
ec_parameterization_Montgomery12 (residue_t b, ec_point_t P0,
                                  const unsigned long k, const modulus_t m)
{
  ASSERT_ALWAYS (ec_parameterization_Montgomery12_is_valid (k));

  residue_t u, v, a, uu, vv, d1, d2, t1, t2, one;
  ec_point_t T;
  int ret = 0;

  /* We want a multiple of the point (-2,4) on the curve y^2=x^3-12*x.
   * The curve has 2-torsion with torsion point (0,0), but adding it
   * does not seem to change the ECM curve we get out in the end.
   */
  mod_init_noset0 (u, m);
  mod_init_noset0 (v, m);
  mod_init_noset0 (a, m);
  mod_init_noset0 (uu, m);
  mod_init_noset0 (vv, m);
  mod_init_noset0 (d1, m);
  mod_init_noset0 (d2, m);
  mod_init_noset0 (t1, m);
  mod_init_noset0 (t2, m);
  mod_init_noset0 (one, m);

  ec_point_init (T, m);

  mod_set1 (one, m);

  /* Compute u and v */
  mod_add (T->y, one, one, m);
  mod_neg (T->x, T->y, m);              /* x = -2 */
  mod_add (T->y, T->y, T->y, m);        /* y = 4 */
  mod_add (a, T->y, T->y, m);
  mod_add (a, a, T->y, m);
  mod_neg (a, a, m);                    /* a = -12 */

  weierstrass_aff_smul_ui (T, k, a, m); // TODO check for failed inv here
  mod_set (u, T->x, m);
  mod_set (v, T->y, m);

  mod_sqr (uu, u, m);                   /* uu = u^2 */
  mod_sqr (vv, v, m);                   /* vv = v^2 */

  mod_add (t1, uu, uu, m);
  mod_add (t1, t1, t1, m);
  mod_sub (d2, vv, t1, m);              /* d2 = v^2-4*u^2 */
  mod_add (t2, t1, t1, m);
  mod_add (t2, t2, t1, m);
  mod_add (d1, vv, t2, m);              /* d1 = v^2+12*u^2 */

  mod_set_ul (t1, 12, m);               /* t1 = 12 */
  mod_add (P0->x, uu, t1, m);           /* x0 = (u^2+12) [so far] */
  mod_sqr (P0->x, P0->x, m);            /* x0 = (u^2+12)^2 [so far] */
  mod_mul (P0->x, P0->x, uu, m);        /* x0 = u^2*(u^2+12)^2 [so far] */
  mod_add (P0->x, P0->x, P0->x, m);     /* x0 = 2*u^2*(u^2+12)^2 */

  mod_mul (P0->z, d1, d2, m);           /* z0 = d1*d2 [so far] */
  mod_add (P0->z, P0->z, P0->z, m);     /* z0 = 2*d1*d2 */

  mod_sqr (t1, d2, m);
  mod_mul (t1, t1, d2, m);
  mod_mul (t1, t1, d1, m);              /* t1 = d1*d2^3 */

  ret = mod_inv (b, t1, m);             /* 1/t1 = 1/(d1*d2^3) (stored in b) */
  if (ret == 0) /* non-trivial gcd */
    mod_set (P0->x, t1, m);
  else
  {
    mod_sqr (t1, uu, m);
    mod_mul (t1, t1, uu, m);
    mod_mul (b, b, t1, m);              /* b = u^6/(d1*d2^3) [so far] */
    mod_mul (b, b, vv, m);              /* b = v^2*u^6/(d1*d2^3)  [so far] */
    mod_2pow_ul (t1, 10, m);            /* t1 = 2^10 = 1024 */
    mod_mul (b, b, t1, m);              /* b = 1024*v^2*u^6/(d1*d2^3) */
  }

  mod_clear (u, m);
  mod_clear (v, m);
  mod_clear (a, m);
  mod_clear (uu, m);
  mod_clear (vv, m);
  mod_clear (d1, m);
  mod_clear (d2, m);
  mod_clear (t1, m);
  mod_clear (t2, m);
  mod_clear (one, m);

  ec_point_clear (T, m);

  return ret;
}

/******************************************************************************/
/********************* Montgomery Z/16Z parameterization **********************/
/******************************************************************************/

/* Produces curve in Montgomery form.
 *
 * Finite number of curves.
 *
 * Reference: Montgomery's thesis (Section 6.2)
 *
 * Valid parameter: k in { 1 }
 */


static inline int
ec_parameterization_Montgomery16_is_valid (const unsigned long k)
{
  if (k != 1)
    return 0; /* invalid values */
  else
    return 1;
}


/* Return 1 if it worked, 0 if a modular inverse failed.
 * If modular inverse failed, return non-invertible value in P0->x.
 */
static int
ec_parameterization_Montgomery16 (residue_t b, ec_point_t P0,
                                  const unsigned long k, const modulus_t m)
{
  ASSERT_ALWAYS (ec_parameterization_Montgomery16_is_valid (k));

  /* Make curve corresponding to (a,b,c) = (8, 15, 17) in Montgomery's thesis,
   * equation (6.2.2).
   *    A = 54721/14400
   *    b = (A+2)/4 = (289/240)^2
   *    x0 = 8
   *    z0 = 15
   * [ see table 6.2.1 ]
   * We only need b and the coordinates of the starting point (x0::z0)
   * This curve is cheap to initialise: four div2, one div3, one div5,
   * three add, one mul.
   */
  residue_t t, one;

  mod_init (t, m);
  mod_init (one, m);

  mod_set1 (one, m);

  mod_set_ul (P0->x, 8UL, m);     /* x0 = 8 */
  mod_set_ul (P0->z, 15UL, m);    /* z0 = 15 */

  mod_div3 (t, one, m);
  mod_div2(t, t, m);
  mod_div2(t, t, m);
  mod_div2(t, t, m);
  mod_div2(t, t, m);              /* t = 1/48 */
  mod_add (t, t, one, m);         /* t = 49/48 */
  mod_div5 (t, t, m);             /* t = 49/240 */
  mod_add (t, t, one, m);         /* t = 289/240 */
  mod_sqr (b, t, m);              /* b = 83521/57600 */

#ifdef WANT_ASSERT_EXPENSIVE
  residue_t t2;
  mod_init (t2, m);
  mod_set_ul (t, 57600UL, m);
  mod_mul (t2, b, t, m);
  mod_set_ul (t, 83521UL, m);
  ASSERT (mod_equal (t, t2, m));
  mod_clear (t2, m);
#endif

  mod_clear (t, m);
  mod_clear (one, m);

  return 1; /* we assume gcd (m, 2*3*5) = 1 */
}

/******************************************************************************/
/******************* Z/6Z-rational-torsion parameterization *******************/
/******************************************************************************/


static inline int
ec_parameterization_Z6_is_valid (const unsigned long k)
{
  if (k == 0)
    return 0; /* invalid values */
  else
    return 1;
}

/* Produces curve
 *     - in "a=-1" Twisted Edwards form (in extended or projective coordinates)
 *     - or in Montgomery from
 * 
 * Elliptic parameterization (parameter is called k).
 *    E: y^2 = x^3 - 9747*x + 285714
 *      rank = 1
 *      P = (15, 378) is a point of infinite order
 *      Torsion points: P2 = (33:0:1) of order 2
 *                      Q2 = (78:0:1) of order 2
 *                      P2+Q2 = (-111:0:1) of order 2
 *
 * Reference: based on Theorem 5.4 of the article "Starfish on strike".
 *
 * Valid parameter: k in Z \ { 0 }
 *
 * Note: P2, Q2, P2+Q2 does not produce a valid curve.
 *       P+Q2, P+P2+Q2 does not produce a valid curve.
 *       Two points that differ by P2 produce the same curve.
 *       A point and its opposite produce the same curve (with opposite starting
 *        points)
 *       It generates curves isomorphic to Brent--Suyama curves with
 *        sigma = -(p-213)/(p+3)
 *       The value k=1 produces the same curve that the Brent-Suyama
 *        parameterization with sigma=11 (up to isomorphism, i.e. the values of
 *        B differ by a square factor).
 *
 * Parameterization (where (p,q) := k*P):
      alpha = (p-105)*(p+57)
      beta = q*(p-213)
      gamma = p^4 - 96*p^3 - 918*p^2 - 2808*p + 45346797
      delta = gamma + 2^5*3^7 *(p-33)^2
      epsilon =  p^2 - 66*p + 7569
      zeta = 2^2*3^6*q*(p-33)
      # Coeffs for "a=-1" Twisted Edwards curve
      a = -1
      d = -alpha^3*(p-51)*(p+327)/zeta^2
      xE0 = 2*3*(p+3)*alpha*beta*delta
      yE0 = (2*3^5)*(p-33)*alpha*gamma*epsilon
      zE0 = alpha^2*delta*epsilon
      tE0 = 2^2*3^6*(p-33)*(p+3)*beta*gamma
      # Coeffs for Montgomery curve
      A = 2*(a+d)/(a-d)
      B = 4/(a-d)
      b = zeta^2/(zeta^2-alpha^3*(p-51)*(p+327))               # [ b = (A+2)/4 ]
      # (xM0:yM0:zM0) = (xE0*(zE0+yE0):(zE0+yE0)*zE0:xE0*(zE0-yE0))
      xM0 = (p^2 + 114*p - 11331)^3
      yM0 = (alpha*epsilon*(epsilon+2^2*3^2*5*(p-105))^3)/(2*3*(p+3)*beta)
      zM0 = ((p-213)*(p+3))^3
 * The point of order 3 is given by:
      # Coeffs for "a=-1" Twisted Edwards curve
      xE3 = 2*3^2*q*alpha
      yE3 = 2*3^4*(p-33)*alpha
      zE3 = alpha^2
      tE3 = zeta
      # Coeffs for Montgomery curve
      xM3 = 2*3^2*q*(epsilon+2^2*3^2*5*(p-105))
      yM3 = alpha*(epsilon+2^2*3^2*5*(p-105))
      zM3 = 2*3^2*(p+3)*beta
 * To check with Sage:
      QQpq.<p,q> = QQ[]
      I = QQpq.ideal (q^2-(p^3-9747*p+285714))
      # copy paste all the above formula
      # P0 is in extended coordinates, we must have xE0*yE0 == zE0*tE0
      (xE0*yE0 - zE0*tE0).is_zero()
      # P0 is on the curve ?
      eq_P0 = -xE0^2+yE0^2 - (zE0^2+d*tE0^2)
      eq_P0.numerator() in I and not eq_P0.denominator() in I
      eq_P02 = B*yM0^2*zM0 - (xM0^3 + A*xM0^2*zM0 + xM0*zM0^2)
      eq_P02.numerator() in I and not eq_P02.denominator() in I
      # P3 is in extended coordinates, we must have xE3*yE3 == zE3*tE3
      (xE3*yE3 - zE3*tE3).is_zero()
      # P3 is on the curve ?
      eq_P3 = -xE3^2+yE3^2 - (zE3^2+d*tE3^2)
      eq_P3.numerator() in I and not eq_P3.denominator() in I
      eq_P32 = B*yM3^2*zM3 - (xM3^3 + A*xM3^2*zM3 + xM3*zM3^2)
      eq_P32.numerator() in I and not eq_P32.denominator() in I
      # P3 is of order 3 ?
      eq_order3 = 4*A*(xM3/zM3)^3+3*(xM3/zM3)^4+6*(xM3/zM3)^2-1
      eq_order3.numerator() in I and not eq_order3.denominator() in I
 */

/* Return 1 if it worked, 0 if a modular inverse failed.
 * If modular inverse failed, return non-invertible value in P0->x.
 *
 * We only need the starting point (xE0:yE0:zE0:tE0) (for Edwards extented) or
 * (xE0:yE0:zE0) (for Edwards projective) or (xM0::zM0) (for Montgomery) and
 * the curve coefficient b (if b is not NULL).
 */
static int
ec_parameterization_Z6 (residue_t b, ec_point_t P0, const unsigned long k,
                        ec_point_coord_type_t coord, const modulus_t m)
{
  ASSERT_ALWAYS (ec_parameterization_Z6_is_valid (k));
  ASSERT_ALWAYS (coord == MONTGOMERY_xz || coord == TWISTED_EDWARDS_ext
                                        || coord == TWISTED_EDWARDS_proj);

  residue_t a, p, q, alpha, beta, gamma, delta, epsilon, zeta, pm33, pp3, t;
  ec_point_t T;
  int ret = 0;

  ec_point_init (T, m);

  mod_init_noset0 (a, m);
  mod_init_noset0 (p, m);
  mod_init_noset0 (q, m);
  mod_init_noset0 (alpha, m);
  mod_init_noset0 (beta, m);
  mod_init_noset0 (gamma, m);
  mod_init_noset0 (delta, m);
  mod_init_noset0 (epsilon, m);
  mod_init_noset0 (zeta, m);
  mod_init_noset0 (pm33, m);
  mod_init_noset0 (pp3, m);
  mod_init_noset0 (t, m);

  /* Compute p and q */
  mod_set_ul (a, 9747UL, m);
  mod_neg (a, a, m);
  mod_set_ul (T->x, 15UL, m);
  mod_set_ul (T->y, 378UL, m);
  mod_set1 (T->z, m);

#if 0
  weierstrass_aff_smul_ui (T, k, a, m);
  mod_set (p, T->x, m);
  mod_set (q, T->y, m);
#else
  weierstrass_proj_smul_ui (T, k, a, m);
  mod_inv (a, T->z, m); // TODO check for failed inversion
  mod_mul (p, T->x, a, m);
  mod_mul (q, T->y, a, m);
#endif

  mod_set_ul (t, 3UL, m);
  mod_add (pp3, p, t, m);                 /* pm3 = p+3 */

  mod_set_ul (t, 33UL, m);
  mod_sub (pm33, p, t, m);                /* pm33 = p-33 */

  mod_set_ul (t, 105UL, m);
  mod_sub (alpha, p, t, m);
  mod_set_ul (t, 57UL, m);
  mod_add (t, p, t, m);
  mod_mul (alpha, alpha, t, m);           /* alpha = (p-105)*(p+57) */

  mod_set_ul (t, 213UL, m);
  mod_sub (t, p, t, m);
  mod_mul (beta, q, t, m);                /* beta = q*(p-213) */

  mod_set_ul (t, 96UL, m);
  mod_sub (gamma, p, t, m);
  mod_mul (gamma, gamma, p, m);
  mod_set_ul (t, 918UL, m);
  mod_sub (gamma, gamma, t, m);
  mod_mul (gamma, gamma, p, m);
  mod_set_ul (t, 2808UL, m);
  mod_sub (gamma, gamma, t, m);
  mod_mul (gamma, gamma, p, m);
  mod_set_ul (t, 45346797, m);
  mod_add (gamma, gamma, t, m); /* gamma = p^4-96*p^3-918*p^2-2808*p+45346797 */

  mod_set_ul (t, 69984UL, m); /* 2^5*3^7 */
  mod_sqr (delta, pm33, m);
  mod_mul (delta, delta, t, m);
  mod_add (delta, gamma, delta, m);       /* delta = gamma + 2^5*3^7*(p-33)^2 */

  mod_set_ul (t, 66UL, m);
  mod_sub (epsilon, p, t, m);
  mod_mul (epsilon, epsilon, p, m);
  mod_set_ul (t, 7569UL, m);
  mod_add (epsilon, epsilon, t, m);       /* epsilon =  p^2 - 66*p + 7569 */

  mod_set_ul (zeta, 2916UL, m);
  mod_mul (zeta, zeta, q, m);
  mod_mul (zeta, zeta, pm33, m);          /* zeta = 2^2*3^6*q*(p-33) */

  if (coord == MONTGOMERY_xz)
  {
    mod_set_ul (t, 114UL, m);
    mod_add (P0->x, p, t, m);
    mod_mul (P0->x, P0->x, p, m);
    mod_set_ul (t, 11331UL, m);
    mod_sub (P0->x, P0->x, t, m);
    mod_sqr (t, P0->x, m);
    mod_mul (P0->x, P0->x, t, m);         /* xM0 = (p^2 + 114*p - 11331)^3 */
    mod_set_ul (t, 213UL, m);
    mod_sub (P0->z, p, t, m);
    mod_mul (P0->z, P0->z, pp3, m);
    mod_sqr (t, P0->z, m);
    mod_mul (P0->z, P0->z, t, m);         /* zM0 = ((p-213)*(p+3))^3 */
  }
  else
  {
    mod_set_ul (t, 6UL, m);
    mod_mul (P0->x, t, pp3, m);
    mod_mul (P0->x, P0->x, alpha, m);
    mod_mul (P0->x, P0->x, beta, m);
    mod_mul (P0->x, P0->x, delta, m);     /* xE0 = 2*3*(p+3)*alpha*beta*delta */
    mod_set_ul (t, 486UL, m);
    mod_mul (P0->y, t, pm33, m);
    mod_mul (P0->y, P0->y, alpha, m);
    mod_mul (P0->y, P0->y, gamma, m);
    mod_mul (P0->y, P0->y, epsilon, m);
                                  /* yE0 = (2*3^5)*(p-33)*alpha*gamma*epsilon */
    mod_sqr (P0->z, alpha, m);
    mod_mul (P0->z, P0->z, delta, m);
    mod_mul (P0->z, P0->z, epsilon, m);   /* zE0 = alpha^2*delta*epsilon */
    if (coord == TWISTED_EDWARDS_ext)
    {
      mod_set_ul (t, 2916UL, m);
      mod_mul (P0->t, t, pm33, m);
      mod_mul (P0->t, P0->t, pp3, m);
      mod_mul (P0->t, P0->t, beta, m);
      mod_mul (P0->t, P0->t, gamma, m);
                                     /* tE0 = 2^2*3^6*(p-33)*(p+3)*beta*gamma */
    }
  }

  if (b != NULL)
  {
    mod_set_ul (t, 51UL, m);
    mod_sub (b, p, t, m);
    mod_set_ul (t, 327UL, m);
    mod_add (t, p, t, m);
    mod_mul (b, b, t, m);
    mod_sqr (t, alpha, m);
    mod_mul (b, b, t, m);
    mod_mul (b, b, alpha, m);
    mod_sqr (t, zeta, m);
    mod_sub (b, t, b, m);

    ret = mod_inv (a, b, m);
    if (ret == 0) /* non-trivial gcd */
      mod_set (P0->x, b, m);
    else
      mod_mul (b, a, t, m);     /* b = zeta^2/(zeta^2-alpha^3*(p-51)*(p+327)) */
  }

  mod_clear (p, m);
  mod_clear (q, m);
  mod_clear (a, m);
  mod_clear (alpha, m);
  mod_clear (beta, m);
  mod_clear (gamma, m);
  mod_clear (delta, m);
  mod_clear (epsilon, m);
  mod_clear (zeta, m);
  mod_clear (pm33, m);
  mod_clear (pp3, m);
  mod_clear (t, m);

  ec_point_clear (T, m);

  return ret;
}

int
ec_parameter_is_valid (ec_parameterization_t parameterization,
                       const unsigned long parameter)
{
  switch (parameterization)
  {
    case BRENT12:
      return ec_parameterization_Brent_Suyama_is_valid (parameter);
    case MONTY12:
      return ec_parameterization_Montgomery12_is_valid (parameter);
    case MONTY16:
      return ec_parameterization_Montgomery16_is_valid (parameter);
    case MONTYTWED12:
      return ec_parameterization_Z6_is_valid (parameter);
    default:
      printf ("Fatal error in %s at %s:%d -- unknown parameterization %u\n",
              __func__, __FILE__, __LINE__, parameterization);
      abort ();
  }
}
