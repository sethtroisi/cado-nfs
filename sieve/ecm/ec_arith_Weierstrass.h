#ifndef _EC_ARITH_WEIERSTRASS_H_
#define _EC_ARITH_WEIERSTRASS_H_

#include "ec_arith_common.h"

/*
  implemented functions:
    - weierstrass_add (R:weierstrass, P:weierstrass, Q:weierstrass)
        R <- P+Q

    - weierstrass_dbl (R:weierstrass, P:weierstrass)
        R <- 2*P
*/



/* R <- 2 * P for the curve y^2 = x^3 + a*x + b.

   For Weierstrass coordinates. Returns 1 if doubling worked normally,
   0 if the result is point at infinity.
*/
static int
weierstrass_dbl (ec_point_t R, const ec_point_t P,
		 const residue_t a, const modulus_t m)
{
  residue_t lambda, u, v;

  mod_init_noset0 (lambda, m);
  mod_init_noset0 (u, m);
  mod_init_noset0 (v, m);

  mod_sqr (u, P->x, m);
  mod_add (v, u, u, m);
  mod_add (v, v, u, m);
  mod_add (v, v, a, m); /* 3x^2 + a */
  mod_add (u, P->y, P->y, m);
  if (mod_inv (u, u, m) == 0)    /* 1/(2*y) */
  {
      mod_clear (v, m);
      mod_clear (u, m);
      mod_clear (lambda, m);
      return 0; /* y was 0  =>  result is point at infinity */
  }
  mod_mul (lambda, u, v, m);
  mod_sqr (u, lambda, m);
  mod_sub (u, u, P->x, m);
  mod_sub (u, u, P->x, m);    /* x3 = u = lambda^2 - 2*x */
  mod_sub (v, P->x, u, m);
  mod_mul (v, v, lambda, m);
  mod_sub (R->y, v, P->y, m);
  mod_set (R->x, u, m);

  mod_clear (v, m);
  mod_clear (u, m);
  mod_clear (lambda, m);
  return 1;
}



/* Adds two points P and Q on the curve y^2 = x^3 + a*x + b
   in Weierstrass coordinates and puts result in R.
   Returns 1 if the addition worked (i.e. the modular inverse existed)
   and 0 otherwise (resulting point is point at infinity) */
static int
weierstrass_add (ec_point_t R, const ec_point_t P, const ec_point_t Q,
		 const residue_t a, const modulus_t m)
{
  residue_t lambda, u, v;
  int r;

  mod_init_noset0 (u, m);
  mod_init_noset0 (v, m);
  mod_init_noset0 (lambda, m);

  mod_sub (u, Q->y, P->y, m);
  mod_sub (v, Q->x, P->x, m);
  if (mod_inv (v, v, m) == 0)
  {
      /* Maybe we were trying to add two identical points? If so,
         use the ellW_double() function instead */
      if (mod_equal (P->x, Q->x, m) && mod_equal (P->y, Q->y, m))
	  r = weierstrass_dbl (R, P, a, m);
      else
	{
	  /* Or maybe the points are negatives of each other? */
	  mod_neg (u, P->y, m);
	  if (mod_equal (P->x, Q->x, m) && mod_equal (u, Q->y, m))
	    r = 0; /* Signal point at infinity */
	  else
	    {
	      /* Neither identical, nor negatives (mod m). Looks like we
		 found a proper factor. FIXME: What do we do with it? */
	      r = 0;
	    }
	}
  }
  else
  {
      mod_mul (lambda, u, v, m);
      mod_sqr (u, lambda, m);
      mod_sub (u, u, P->x, m);
      mod_sub (u, u, Q->x, m);    /* x3 = u = lambda^2 - P->x - Q->x */
      mod_sub (v, P->x, u, m);
      mod_mul (v, v, lambda, m);
      mod_sub (R->y, v, P->y, m);
      mod_set (R->x, u, m);
      r = 1;
  }

  mod_clear (lambda, m);
  mod_clear (v, m);
  mod_clear (u, m);
  return r;
}


/* Convert the Montgomery curve B*Y^2*Z = X^3 + A*X^2*Z + X*Z^2 with a valid
 * point Pm into a affine Weierstrass curve y^2 = x^3 + a*x + b with a
 * valid point Pw.
 *
 * Return 1 if it worked, 0 if a modular inverse failed.
 *
 * The value of b will not be computed.
 * a and b can be the same variable.
 * Pm and Pw can be the same variable.
 */
static int
weierstrass_from_montgomery (residue_t a, ec_point_t Pw, const residue_t A,
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
    fprintf (stderr, "%s: could not invert Z\n", __func__);
  else
  {
    mod_mul (x, Pm->x, t, m); /* x = X/Z */
    mod_add (B, x, A, m);
    mod_mul (B, B, x, m);
    mod_add (B, B, one, m);
    mod_mul (B, B, x, m); /* B = x^3 + A*x^2 + x */

    /* Now (x,1) is on the curve B*y^2 = x^3 + A*x^2 + x. */
    ret = mod_inv (B, B, m);
    if (ret == 0)
      fprintf (stderr, "%s: could not invert B\n", __func__);
    else
    {
      mod_set (Pw->y, B, m);        /* y = 1/B */
      mod_div3 (t, A, m);
      mod_add (Pw->x, x, t, m);
      mod_mul (Pw->x, Pw->x, B, m); /* x = (X + A/3)/B */
      mod_mul (a, t, A, m);
      mod_sub (a, one, a, m);
      mod_mul (a, a, B, m);
      mod_mul (a, a, B, m);         /* a = (1 - (A^2)/3)/B^2 */
    }
  }

  mod_clear (one, m);
  mod_clear (B, m);
  mod_clear (t, m);
  mod_clear (x, m);

  return ret;
}


/* (x,y) <- e * (x,y) on the curve y^2 = x^3 + a*x + b (mod m) */
static int
weierstrass_smul_ui (ec_point_t P, const unsigned long e, const residue_t a,
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

  ec_point_set (T, P, m);
  tfinite = 1;
  i >>= 1;

  while (i > 0)
  {
      if (tfinite)
	tfinite = weierstrass_dbl (T, T, a, m);
      if (e & i)
      {
	  if (tfinite)
	      tfinite = weierstrass_add (T, T, P, a, m);
	  else
	  {
	      ec_point_set (T, P, m);
	      tfinite = 1;
	  }
      }
      i >>= 1;
  }

  if (tfinite)
    ec_point_set (P, T, m);

  ec_point_clear (T, m);

  return tfinite;
}

#endif /* _EC_ARITH_WEIERSTRASS_H_ */
