#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "portability.h"
#include "facul_ecm.h"
#include "ularith.h"
#include "getprime.h"

#include "ec_arith_common.h"
#include "ec_arith_Edwards.h"
#include "ec_arith_Montgomery.h"
#include "ec_arith_Weierstrass.h"


/* Do we want backtracking when processing factors of 2 in E? */
#ifndef ECM_BACKTRACKING
/* Default is "yes." Set to 0 for "no." */
#define ECM_BACKTRACKING 1
#endif

/* Define to 1 to make ellM_add() test if the two points are identical,
   and call ellM_double() if they are */
#ifndef ELLM_SAFE_ADD
#define ELLM_SAFE_ADD 0
#endif

#define COUNT_ELLM_OPS 0
#if COUNT_ELLM_OPS
static unsigned long ellM_dadd_count, ellM_double_count;
#endif

#define COUNT_ELLE_OPS 0
#if COUNT_ELLE_OPS
static unsigned long ellE_add_count, ellE_double_count, ellE_triple_count;
#endif

/******************************************************************************/
/*********************** bytecode interpreter functions ***********************/
/******************************************************************************/

/**** Interpret the double-base chain bytecode for Twisted Edwards curves ****/
/* Internal function.
 * Return a pointer to the last parsed byte.
 * Assume that the bytecode is correct (for example, assume that the size of the
 * array is correct).
 */
static bytecode_const
bytecode_dbchain_interpret_ec_edwards_internal (bytecode_const bc,
                                                ec_point_t *R,
                                                const modulus_t m)
{
  while (1)
  {
    uint8_t op, f, s, n, pow2 = 0, pow3 = 0;
    bytecode_elt_split_2_1_1_4 (&op, &f, &s, &n, *bc);

    if (op == DBCHAIN_OP_TPLADD || op == DBCHAIN_OP_TPLDBLADD)
    {
      bc++;
      pow3 = bytecode_elt_to_uint8 (*bc);
    }
    if (op == DBCHAIN_OP_DBLADD || op == DBCHAIN_OP_TPLDBLADD)
    {
      bc++;
      pow2 = bytecode_elt_to_uint8 (*bc);
    }

    for (uint8_t i = 0; i < pow3; i++)
      edwards_tpl (R[0], R[0], m, (i+1==pow3+pow2) ? EDW_ext : EDW_proj);
    for (uint8_t i = 0; i < pow2; i++)
      edwards_dbl (R[0], R[0], m, (i+1==pow2) ? EDW_ext : EDW_proj);

    ec_point_coord_type_t output;
    if (f)
    {
      uint8_t t;
      bytecode_elt_split_4_4 (&t, NULL, bc[1]);
      if (bc[1] == MISHMASH_FINAL || t == MISHMASH_PRAC_BLOCK)
        output = MONTG;
      else
        output = EDW_ext;
    }
    else
      output = EDW_proj;

    if (s == 0)
      edwards_add (R[f], R[0], R[n], m, output);
    else
      edwards_sub (R[f], R[0], R[n], m, output);

    if (f) /* is it finished ? */
      break;
    else
      bc++; /* go to next byte */
  }

  return bc;
}

/********* Interpret the precomp bytecode for Twisted Edwards curves **********/
/* Internal function.
 * Return a pointer to the last parsed byte.
 * Assume that the bytecode is correct (for example, assume that the size of the
 * array is correct).
 */
static bytecode_const
bytecode_precomp_interpret_ec_edwards_internal (bytecode_const bc,
                                                ec_point_t *R,
                                                const modulus_t m)
{
  while (*bc != PRECOMP_FINAL)
  {
    uint8_t op, a, s, i, j, k, pow2, pow3;
    bytecode_elt_split_2_1_1_4 (&op, &a, &s, &k, *bc);

    switch (op)
    {
      case PRECOMP_OP_ADD:
        bc++;
        bytecode_elt_split_4_4 (&i, &j, *bc);
        if (s == 0)
          edwards_add (R[k], R[i], R[j], m, a ? EDW_ext : EDW_proj);
        else
          edwards_sub (R[k], R[i], R[j], m, a ? EDW_ext : EDW_proj);
        break;
      case PRECOMP_OP_DBL:
        bc++;
        pow2 = bytecode_elt_to_uint8 (*bc);
        for (uint8_t i = 0; i < pow2; i++)
          edwards_dbl (R[0], R[0], m, (i+1==pow2 && a) ? EDW_ext : EDW_proj);
        if (k > 0)
          ec_point_set (R[k], R[0], m);
        break;
      case PRECOMP_OP_TPL:
        bc++;
        pow3 = bytecode_elt_to_uint8 (*bc);
        for (uint8_t i = 0; i < pow3; i++)
          edwards_tpl (R[0], R[0], m, (i+1==pow3 && a) ? EDW_ext : EDW_proj);
        if (k > 0)
          ec_point_set (R[k], R[0], m);
        break;
      default:
        printf ("Fatal error in %s at %s:%d -- unknown bytecode 0x%02x\n",
                __func__, __FILE__, __LINE__, *bc);
        abort ();
    }

    bc++; /* go to next byte */
  }

  return bc;
}

/************* Interpret the PRAC bytecode for Montgomery curves *************/
/* Internal function.
 * Return a pointer to the last parsed byte.
 * Assume that the bytecode is correct (for example, assume that the size of the
 * array R is >= 5).
 */
static bytecode_const
bytecode_prac_interpret_ec_montgomery_internal (bytecode_const bc,
                                                ec_point_t *R,
                                                const modulus_t m,
                                                const residue_t b)
{
  while (1)
  {
    int finished = 0;
    switch (*bc)
    {
      case PRAC_SWAP: /* [ = 's' ] Swap R[0], R[1] */
        ec_point_swap (R[0], R[1], m);
        break;
      case PRAC_SUBBLOCK_INIT: /* [ = 'i' ] Start of a sub-block */
        ec_point_set (R[1], R[0], m);
        ec_point_set (R[2], R[0], m);
        montgomery_dbl (R[0], R[0], m, b);
        break;
      case PRAC_SUBBLOCK_FINAL: /* [ = 'f' ] End of a sub-block */
        montgomery_dadd (R[0], R[0], R[1], R[2], b, m);
        break;
      case PRAC_BLOCK_FINAL: /* [ = 'F' ] End of the block */
        montgomery_dadd (R[1], R[0], R[1], R[2], b, m);
        finished = 1;
        break;
      case 1:
        montgomery_dadd (R[3], R[0], R[1], R[2], b, m);
        montgomery_dadd (R[4], R[3], R[0], R[1], b, m);
        montgomery_dadd (R[1], R[1], R[3], R[0], b, m);
        ec_point_set (R[0], R[4], m);
        break;
      case 2:
        montgomery_dadd (R[1], R[0], R[1], R[2], b, m);
        montgomery_dbl (R[0], R[0], m, b);
        break;
      case 3:
        montgomery_dadd (R[2], R[1], R[0], R[2], b, m);
        ec_point_swap (R[1], R[2], m);
        break;
      case 4:
        montgomery_dadd (R[1], R[1], R[0], R[2], b, m);
        montgomery_dbl (R[0], R[0], m, b);
        break;
      case 5:
        montgomery_dadd (R[2], R[2], R[0], R[1], b, m);
        montgomery_dbl (R[0], R[0], m, b);
        break;
      case 6:
        montgomery_dbl (R[3], R[0], m, b);
        montgomery_dadd (R[4], R[0], R[1], R[2], b, m);
        montgomery_dadd (R[0], R[3], R[0], R[0], b, m);
        montgomery_dadd (R[2], R[3], R[4], R[2], b, m);
        ec_point_swap (R[1], R[2], m);
        break;
      case 7:
        montgomery_dadd (R[3], R[0], R[1], R[2], b, m);
        montgomery_dadd (R[1], R[3], R[0], R[1], b, m);
        montgomery_dbl (R[3], R[0], m, b);
        montgomery_dadd (R[0], R[0], R[3], R[0], b, m);
        break;
      case 8:
        montgomery_dadd (R[3], R[0], R[1], R[2], b, m);
        montgomery_dadd (R[2], R[2], R[0], R[1], b, m);
        ec_point_swap (R[1], R[3], m);
        montgomery_dbl (R[3], R[0], m, b);
        montgomery_dadd (R[0], R[0], R[3], R[0], b, m);
        break;
      case 9:
        montgomery_dadd (R[2], R[2], R[1], R[0], b, m);
        montgomery_dbl (R[1], R[1], m, b);
        break;
      case 10:
        /* Combined final add of old subchain and init of new subchain [=fi] */
        montgomery_dadd (R[1], R[0], R[1], R[2], b, m);
        ec_point_set (R[2], R[1], m);
        montgomery_dbl (R[0], R[1], m, b);
        break;
      case 11:
        /* Combined rule 3 and rule 0 [=\x3s] */
        montgomery_dadd (R[2], R[1], R[0], R[2], b, m);
        /* (R[1],R[2],R[0]) := (R[0],R[1],R[2])  */
        ec_point_swap (R[1], R[2], m);
        ec_point_swap (R[0], R[1], m);
        break;
      case 12:
        /* Combined rule 3, then subchain end/start [=\x3fi] */
        montgomery_dadd (R[3], R[1], R[0], R[2], b, m);
        montgomery_dadd (R[2], R[0], R[3], R[1], b, m);
        ec_point_set (R[1], R[2], m);
        montgomery_dbl (R[0], R[2], m, b);
        break;
      case 13:
        /* Combined rule 3, swap, rule 3 and swap, merged a bit [=\x3s\x3s] */
        ec_point_set (R[3], R[1], m);
        montgomery_dadd (R[1], R[1], R[0], R[2], b, m);
        ec_point_set (R[2], R[0], m);
        montgomery_dadd (R[0], R[0], R[1], R[3], b, m);
        break;
      default:
        printf ("Fatal error in %s at %s:%d -- unknown bytecode 0x%02x\n",
                __func__, __FILE__, __LINE__, *bc);
        abort ();
    }

    if (finished) /* is it finished ? */
      break;
    else
      bc++; /* go to next byte */
  }

  return bc;
}

static void
bytecode_prac_interpret_ec_montgomery (ec_point_t P, bytecode_const bc,
                                       const modulus_t m, const residue_t b)
{
  ec_point_t *R = NULL;
  unsigned int R_nalloc;

  R_nalloc = 5; /* we need 5 points: 3 for PRAC + 2 temporary points */
  R = (ec_point_t *) malloc (R_nalloc * sizeof (ec_point_t));
  FATAL_ERROR_CHECK (R == NULL, "could not malloc R");
  for (unsigned int i = 0; i < R_nalloc; i++)
    ec_point_init (R[i], m);

  /* current point (here starting point) go into R[0] at init */
  ec_point_set (R[0], P, m);

  bytecode_prac_interpret_ec_montgomery_internal (bc, R, m, b);

  /* output is in R[1] */
  ec_point_set (P, R[1], m);

  for (unsigned int i = 0; i < R_nalloc; i++)
    ec_point_clear (R[i], m);
  free (R);
}

/** Interpret the MISHMASH bytecode on Twisted Edwards and Montgomery curves **/
MAYBE_UNUSED static void
bytecode_mishmash_interpret_ec_mixed_repr (ec_point_t P, bytecode_const bc,
                                           const modulus_t m, const residue_t b)
{
  ec_point_t *R = NULL;
  unsigned int R_nalloc;
  uint8_t n;

  /* parse first byte to alloc R and go to the next byte */
  bytecode_elt_split_4_4 (NULL, &n, *bc);
  bc++;

  R_nalloc = n + 2;
  R = (ec_point_t *) malloc (R_nalloc * sizeof (ec_point_t));
  FATAL_ERROR_CHECK (R == NULL, "could not malloc R");
  for (unsigned int i = 0; i < R_nalloc; i++)
    ec_point_init (R[i], m);

  /* starting point go into R[1] at init */
  ec_point_set (R[1], P, m);

  while (*bc != MISHMASH_FINAL)
  {
    uint8_t t, n;
    bytecode_elt_split_4_4 (&t, &n, *bc);

    if (n != 0)
      ec_point_set (R[0], R[n], m);

    if (t == MISHMASH_DBCHAIN_BLOCK)
      bc = bytecode_dbchain_interpret_ec_edwards_internal (++bc, R, m);
    else if (t == MISHMASH_PRECOMP_BLOCK)
      bc = bytecode_precomp_interpret_ec_edwards_internal (++bc, R, m);
    else if (t == MISHMASH_PRAC_BLOCK)
      bc = bytecode_prac_interpret_ec_montgomery_internal (++bc, R, m, b);
    else /* unexpected bytecode */
    {
      printf ("Fatal error in %s at %s:%d -- unexpected bytecode 0x%02x\n",
              __func__, __FILE__, __LINE__, *bc);
      abort ();
    }

    bc++; /* go to next byte */
  }

  /* output is in R[1] */
  ec_point_set (P, R[1], m);

  for (unsigned int i = 0; i < R_nalloc; i++)
    ec_point_clear (R[i], m);
  free (R);
}


/******************************************************************************/
/************************ parameterization functions **************************/
/******************************************************************************/

/* Produces curve in Montgomery form from 'sigma' value.
 * Also known as Suyama parameterization.
 * Return 1 if it worked, 0 if a modular inverse failed.
 * If modular inverse failed, return non-invertible value in P0->x.
 *
 * Parametrization:
      u = sigma^2-5
      v = 4*sigma
      x0 = u^3
      y0 = (sigma^2-1)*(sigma^2-25)*(sigma^4-25)
      z0 = v^3
      A = (v-u)^3*(3*u+v)/(4*u^3*v) - 2
      B = u/z0
      b = (v-u)^3*(3*u+v)/(16*u^3*v)                          # [ b = (A+2)/4 ]
 * We only need b and the coordinates of the starting point (x0::z0)
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
static int
Brent12_curve_from_sigma (residue_t b, ec_point_t P0, const unsigned long sigma,
                          const modulus_t m)
{
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

/* Produces curve in Montgomery form from k value, using parameterization for a
 * torsion 12 curve as in Montgomery's thesis (6.1).
 * Return 1 if it worked, 0 if a modular inverse failed.
 * If modular inverse failed, return non-invertible value in P0->x.
 *
 * Parametrization:
 *    (u,v) = k*P
 *      where P = (-2, 4) is a point on the Weierstrass curve y^2 = x^3 - 12*x
      d1 = v^2+12*u^2
      d2 = v^2-4*u^2
      x0 = 2*(u*(u^2+12))^2
      y0 = d1*u*(u^2+12)
      z0 = 2*d1*d2
      A = (4*(32*v*u^3)^2-2*d1*d2^3)/(d1*d2^3)
      B = (16*u^2*(v^2+4*u^2))^2/(d1*d2^3)
      b = (32*v*u^3)^2/(d1*d2^3)
 * We only need b and the coordinates of the starting point (x0::z0)
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
      # P12 of order 3 ?
      # TODO
      # P12 of order 4 ?
      # TODO
 */
static int
Montgomery12_curve_from_k (residue_t b, ec_point_t P0, const unsigned long k,
		                       const modulus_t m)
{
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

  weierstrass_smul_ui (T, k, a, m); // TODO check for failed inv here
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


/* Produces curve in Montgomery parameterization from k value, using
   parameters for a torsion 16 curve as in Montgomery's thesis.
   Return 1 if it worked, 0 if a modular inverse failed.
   Currently can produce only one, hard-coded curve that is cheap
   to initialise */
static int
Montgomery16_curve_from_k (residue_t b, ec_point_t P0, const unsigned long k,
		                       const modulus_t m)
{
  if (k == 1UL)
  {
    /* Make curve corresponding to (a,b,c) = (8, 15, 17) in Montgomery's thesis,
     * equation (6.2.2).
     *    A = 54721/14400
     *    b = (A+2)/4 = (289/240)^2
     *    x0 = 8
     *    z0 = 15
     * [ see table 6.2.1 ]
     * We only need b and the coordinates of the starting point (x0::z0)
     * This curve is cheap to initialise: four div2, one div3, two div5,
     * three add, one mul
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
  }
  else
  {
    abort();
  }

  return 1;
}

/* Produces curve
 *     - in "a=-1" Twisted Edwards form (in extended or projective coordinates)
 *     - or in Montgomery from
 * from a value k, using parameterization for a rational torsion Z/6Z (implies a
 * torsion of order 12 modulo every prime), based on Theorem 5.4 of the article
 * "Starfish on strike".
 * 'parameter' must be > 1
 * Return 1 if it worked, 0 if a modular inverse failed.
 * If modular inverse failed, return non-invertible value in P0->x.
 *
 * Parametrization:
 *    (p,q) = k*P
 *      where P = ?? is a point on the Weierstrass curve
 *        y^2 = x^3 - 9747*x + 285714
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
 * We only need the starting point (xE0:yE0:zE0:tE0) (for Edwards extented) or
 * (xE0:yE0:zE0) (for Edwards projective) or (xM0::zM0) (for Montgomery) and
 * maybe the curve coefficient b.
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
 * It corresponds to the Brent--Suyama parameterization with
      sigma = -(p-213)/(p+3)
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
 *
 * Remark: The value k=1 produces the same curve that Brent-Suyama
 * parameterization with sigma=11 (up to isomorphism, i.e. the values of B
 * differ by a square factor).
 */
static int
ec_parameterization_Z6 (residue_t b, ec_point_t P0, const unsigned long k,
                        ec_point_coord_type_t coord, const modulus_t m)
{
  ASSERT_ALWAYS (coord == MONTG || coord == EDW_ext || coord == EDW_proj);
  ASSERT_ALWAYS (k > 0);

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

  weierstrass_smul_ui (T, k, a, m); // TODO check for failed inv here
  mod_set (p, T->x, m);
  mod_set (q, T->y, m);

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

  if (coord == MONTG)
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
    if (coord == EDW_ext)
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

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/* Multiplies x[1] by z[2]*z[3]*z[4]...*z[n],
   x[2] by z[1]*z[3]*z[4]...*z[n] etc., generally
   x[i] by \prod_{1\leq j \leq n, j\neq i} z[j]
   Requires n > 1. Uses 4n-6 multiplications. */

MAYBE_UNUSED
static void ATTRIBUTE((__noinline__))
common_z (const int n1, residue_t *x1, residue_t *z1,
	  const int n2, residue_t *x2, residue_t *z2,
	  const modulus_t m)
{
  const int n = n1 + n2;
  int i, j;
  residue_t *t, p;
  const int verbose = 0;

  if (verbose)
    printf ("common_z: n1 = %d, n2 = %d, sum = %d, nr muls=%d\n",
            n1, n2, n1 + n2, 4*(n1 + n2) - 6);

  if (n < 2)
    return;

  t = (residue_t *) malloc (n * sizeof (residue_t));
  for (i = 0; i < n; i++)
    mod_init (t[i], m);

  /* Set t[i] = z_0 * z_1 * ... * z_i, where the z_i are taken
     from the two lists z1 and z2 */
  i = j = 0;
  if (n1 == 0)
    mod_set (t[0], z2[j++], m);
  else
    mod_set (t[0], z1[i++], m);

  for ( ; i < n1; i++)
    mod_mul (t[i], t[i - 1], z1[i], m);

  for ( ; j < n2; j++)
    mod_mul (t[j + n1], t[j + n1 - 1], z2[j], m);

  /* Now t[i] contains z_0 * ... * z_i */

#ifdef ECM15_UL
  if (verbose) {
    printf ("t[] = {");
    for (int verbose_i = 0; verbose_i < n; verbose_i++)
      printf ("{%lu,%lu}%s", t[verbose_i][0], t[verbose_i][1], (verbose_i < n-1) ? ", " : "");
    printf ("}\n");
  }
#endif
  mod_init (p, m);

  i = n - 1;
  if (i < n1) /* <==>  n2 == 0 */
    mod_mul (x1[i], x1[i], t[n - 2], m);
  else
    mod_mul (x2[i - n1], x2[i - n1], t[n - 2], m);

  /* Init the accumulator p, either to the last element of z2 if z2 is
     non-empty, or to the last element of z1 if z2 is empty */
  if (n2 > 0)
    mod_set (p, z2[n2 - 1], m);
  else
    mod_set (p, z1[n1 - 1], m);

  for (i = n2 - 2; i > -n1 && i >= 0; i--)
    {
      /* Here p = z_{i+1} * ... * z_{n-1} */
      mod_mul (x2[i], x2[i], p, m);
      mod_mul (x2[i], x2[i], t[i + n1 - 1], m);
      mod_mul (p, p, z2[i], m);
    }

  /* n1 = 0  =>  i = 0 */
  /* n1 > 0  =>  i = -1 or -2 */

  for (i = i + n1 ; i > 0; i--)
    {
      /* Here p = z_{i+1} * ... * z_{n-1} */
      mod_mul (x1[i], x1[i], p, m);
      mod_mul (x1[i], x1[i], t[i-1], m);
      mod_mul (p, p, z1[i], m);
    }

  if (n1 > 0)
    mod_mul (x1[0], x1[0], p, m);
  else
    mod_mul (x2[0], x2[0], p, m);

  mod_clear (p, m);

  for (i = 0; i < n; i++)
    mod_clear (t[i], m);
  free (t);
  t = NULL;
}


static int ATTRIBUTE((__noinline__))
ecm_stage2 (residue_t r, const ec_point_t P, const stage2_plan_t *plan,
	    const residue_t b, const modulus_t m)
{
  ec_point_t Pd, Pt; /* d*P, i*d*P, (i+1)*d*P and a temp */
  residue_t *Pid_x, *Pid_z, *Pj_x, *Pj_z; /* saved i*d*P, i0 <= i < i1,
              and jP, j in S_1, x and z coordinate stored separately */
  residue_t a, a_bk, t;
  unsigned int i, k, l;
  int bt = 0;
  const int verbose = 0;

  ec_point_init (Pt, m);
  mod_init_noset0 (t, m);
  mod_init_noset0 (a, m);
  Pj_x = (residue_t *) malloc (plan->s1 * sizeof(residue_t));
  ASSERT (Pj_x != NULL);
  Pj_z = (residue_t *) malloc (plan->s1 * sizeof(residue_t));
  ASSERT (Pj_z != NULL);
  ASSERT(plan->i0 < plan->i1);
  Pid_x = (residue_t *) malloc ((plan->i1 - plan->i0) * sizeof(residue_t));
  ASSERT (Pid_x != NULL);
  Pid_z = (residue_t *) malloc ((plan->i1 - plan->i0) * sizeof(residue_t));
  ASSERT (Pid_z != NULL);
  for (i = 0; i < plan->s1; i++)
    {
      mod_init_noset0 (Pj_x[i], m);
      mod_init_noset0 (Pj_z[i], m);
    }
  for (i = 0; i < plan->i1 - plan->i0; i++)
    {
      mod_init_noset0 (Pid_x[i], m);
      mod_init_noset0 (Pid_z[i], m);
    }

#if COUNT_ELLM_OPS
  ellM_add_count = ellM_double_count = 0;
#endif

  if (verbose)
    printf ("Stage 2: P = (%lu::%lu)\n",
            mod_get_ul (P[0].x, m), mod_get_ul (P[0].z, m));

  /* Compute jP for j in S_1. Compute all the j, 1 <= j < d/2, gcd(j,d)=1
     with two arithmetic progressions 1+6k and 5+6k (this assumes 6|d).
     We need two values of each progression (1, 7 and 5, 11) and the
     common difference 6. These can be computed with the Lucas chain
     1, 2, 3, 5, 6, 7, 11 at the cost of 6 additions and 1 doubling.
     For d=210, generating all 24 desired values 1 <= j < 210/2, gcd(j,d)=1,
     takes 6+16+15=37 point additions. If d=30, we could use
     1,2,3,4,6,7,11,13 which has 5 additions and 2 doublings */

  ASSERT (plan->d % 6 == 0);
  {
    ec_point_t ap1_0, ap1_1, ap5_0, ap5_1, P2, P6;
    int i1, i5;
    ec_point_init (ap1_0, m);
    ec_point_init (ap1_1, m);
    ec_point_init (ap5_0, m);
    ec_point_init (ap5_1, m);
    ec_point_init (P6, m);
    ec_point_init (P2, m);

    /* Init ap1_0 = 1P, ap1_1 = 7P, ap5_0 = 5P, ap5_1 = 11P
       and P6 = 6P */
    ec_point_set (ap1_0, P, m);            /* ap1_0 = 1*P */
    montgomery_dbl (P2, P, m, b);         /* P2 = 2*P */
    montgomery_dadd (P6, P2, P, P, b, m);     /* P6 = 3*P (for now) */
    montgomery_dadd (ap5_0, P6, P2, P, b, m); /* 5*P = 3*P + 2*P */
    montgomery_dbl (P6, P6, m, b);        /* P6 = 6*P = 2*(3*P) */
    montgomery_dadd (ap1_1, P6, P, ap5_0, b, m); /* 7*P = 6*P + P */
    montgomery_dadd (ap5_1, P6, ap5_0, P, b, m); /* 11*P = 6*P + 5*P */

    /* Now we generate all the j*P for j in S_1 */
    /* We treat the first two manually because those might correspond
       to ap1_0 = 1*P and ap5_0 = 5*P */
    k = 0;
    if (plan->s1 > k && plan->S1[k] == 1)
      {
        mod_set (Pj_x[k], ap1_0[0].x, m);
        mod_set (Pj_z[k], ap1_0[0].z, m);
        k++;
      }
    if (plan->s1 > k && plan->S1[k] == 5)
      {
        mod_set (Pj_x[k], ap5_0[0].x, m);
        mod_set (Pj_z[k], ap5_0[0].z, m);
        k++;
      }

    i1 = 7;
    i5 = 11;
    while (k < plan->s1)
      {
        if (plan->S1[k] == i1)
          {
            mod_set (Pj_x[k], ap1_1[0].x, m);
            mod_set (Pj_z[k], ap1_1[0].z, m);
	    k++;
	    continue;
          }
        if (plan->S1[k] == i5)
          {
            mod_set (Pj_x[k], ap5_1[0].x, m);
            mod_set (Pj_z[k], ap5_1[0].z, m);
	    k++;
	    continue;
          }
	
        montgomery_dadd (Pt, ap1_1, P6, ap1_0, b, m);
        ec_point_set (ap1_0, ap1_1, m);
        ec_point_set (ap1_1, Pt, m);
        i1 += 6;
	
        montgomery_dadd (Pt, ap5_1, P6, ap5_0, b, m);
        ec_point_set (ap5_0, ap5_1, m);
        ec_point_set (ap5_1, Pt, m);
        i5 += 6;
      }

    if ((unsigned) (i1 + i5) < plan->d)
      {
        if (i1 < i5)
          {
            montgomery_dadd (Pt, ap1_1, P6, ap1_0, b, m);
            ec_point_set (ap1_0, ap1_1, m);
            ec_point_set (ap1_1, Pt, m);
            i1 += 6;
          }
        else
          {
            montgomery_dadd (Pt, ap5_1, P6, ap5_0, b, m);
            ec_point_set (ap5_0, ap5_1, m);
            ec_point_set (ap5_1, Pt, m);
            i5 += 6;
          }
      }

    ec_point_init (Pd, m);
#if 0
    /* Also compute Pd = d*P while we've got 6*P */
    montgomery_mul_ul (Pd, P6, plan->d / 6, m, b); /* slow! */
#else
    ASSERT ((unsigned) (i1 + i5) == plan->d);
    if (i1 + 4 == i5)
      {
        montgomery_dbl (P2, P2, m, b); /* We need 4P for difference */
        montgomery_dadd (Pd, ap1_1, ap5_1, P2, b, m);
      }
    else if (i5 + 2 == i1)
      {
        montgomery_dadd (Pd, ap1_1, ap5_1, P2, b, m);
      }
    else
      abort ();
#endif

    ec_point_clear (ap1_0, m);
    ec_point_clear (ap1_1, m);
    ec_point_clear (ap5_0, m);
    ec_point_clear (ap5_1, m);
    ec_point_clear (P6, m);
    ec_point_clear (P2, m);

  }

  if (verbose)
    {
      printf ("Pj = [");
      for (i = 0; i < plan->s1; i++)
        printf ("%s(%lu::%lu)", (i>0) ? ", " : "",
                mod_get_ul (Pj_x[i], m), mod_get_ul (Pj_z[i], m));
      printf ("]\nPd = (%lu::%lu)\n",
                mod_get_ul (Pd[0].x, m), mod_get_ul (Pd[0].z, m));
    }

  /* Compute idP for i0 <= i < i1 */
  {
    ec_point_t Pid, Pid1;

    ec_point_init (Pid, m);
    ec_point_init (Pid1, m);
    k = 0; i = plan->i0;

    /* If i0 == 0, we simply leave the first point at (0::0) which is the
       point at infinity */
    if (plan->i0 == 0)
      {
	mod_set0 (Pid_x[k], m);
	mod_set0 (Pid_z[k], m);
	k++;
	i++;
      }
    /* XXX CB: I think this is buggy: with B1=17 and B2=42, i0=0 i1=1, Pid_x
     * has length 1, so k must be < 1.
     */

    /* Todo: do both Pid and Pid1 with one addition chain */
    montgomery_smul_ui (Pid, Pd, i, m, b); /* Pid = i_0 d * P */
    mod_set (Pid_x[k], Pid->x, m);
    mod_set (Pid_z[k], Pid->z, m);
    k++; i++;
    if (i < plan->i1)
      {
        montgomery_smul_ui (Pid1, Pd, i, m, b); /* Pid = (i_0 + 1) d * P */
        mod_set (Pid_x[k], Pid1[0].x, m);
        mod_set (Pid_z[k], Pid1[0].z, m);
        k++; i++;
      }
    while (i < plan->i1)
      {
        montgomery_dadd (Pt, Pid1, Pd, Pid, b, m);
        ec_point_set (Pid, Pid1, m);
        ec_point_set (Pid1, Pt, m);
        mod_set (Pid_x[k], Pt[0].x, m);
        mod_set (Pid_z[k], Pt[0].z, m);
        k++; i++;
      }

    ec_point_clear (Pd, m);
    ec_point_clear (Pid, m);
    ec_point_clear (Pid1, m);
  }

  if (verbose)
    {
      printf ("Pid = [");
      for (i = 0; i < plan->i1 - plan->i0; i++)
        printf ("%s(%lu:%lu)", (i>0) ? ", " : "",
	        mod_get_ul (Pid_x[i], m), mod_get_ul (Pid_z[i], m));
      printf ("]\n");
    }

  /* Now we've computed all the points we need, so multiply each by
     the Z-coordinates of all the others, using Zimmermann's
     two product-lists trick.
     If i0 == 0, then Pid[0] is the point at infinity (0::0),
     so we skip that one */
  {
    int skip = (plan->i0 == 0) ? 1 : 0;
#if defined(ECM15_UL)
    if (verbose)
      {
        residue_t *x1 = Pj_x, *z1 = Pj_z, *x2 = Pid_x + skip, *z2 = Pid_z + skip;
        int n1 = plan->s1, n2 = plan->i1 - plan->i0 - skip;

        printf ("Before common_z():\n");
        printf ("x1 = {");
        for (int verbose_i = 0; verbose_i < n1; verbose_i++)
          printf ("{%luUL,%luUL}%s", x1[verbose_i][0], x1[verbose_i][1], (verbose_i < n1-1) ? ", " : "}\n");
        printf ("z1 = {");
        for (int verbose_i = 0; verbose_i < n1; verbose_i++)
          printf ("{%luUL,%luUL}%s", z1[verbose_i][0], z1[verbose_i][1], (verbose_i < n1-1) ? ", " : "}\n");
        printf ("x2 = {");
        for (int verbose_i = 0; verbose_i < n2; verbose_i++)
          printf ("{%luUL,%luUL}%s", x2[verbose_i][0], x2[verbose_i][1], (verbose_i < n2-1) ? ", " : "}\n");
        printf ("z2 = {");
        for (int verbose_i = 0; verbose_i < n2; verbose_i++)
          printf ("{%luUL,%luUL}%s", z2[verbose_i][0], z2[verbose_i][1], (verbose_i < n2-1) ? ", " : "}\n");
      }
#endif
    common_z (plan->s1, Pj_x, Pj_z, plan->i1 - plan->i0 - skip,
	      Pid_x + skip, Pid_z + skip, m);
#if defined(ECM15_UL)
    if (verbose)
      {
        residue_t *x1 = Pj_x, *z1 = Pj_z, *x2 = Pid_x + skip, *z2 = Pid_z + skip;
        int n1 = plan->s1, n2 = plan->i1 - plan->i0 - skip;

        printf ("After common_z():\n");
        printf ("x1 = {");
        for (int verbose_i = 0; verbose_i < n1; verbose_i++)
          printf ("{%luUL,%luUL}%s", x1[verbose_i][0], x1[verbose_i][1], (verbose_i < n1-1) ? ", " : "}\n");
        printf ("z1 = {");
        for (int verbose_i = 0; verbose_i < n1; verbose_i++)
          printf ("{%luUL,%luUL}%s", z1[verbose_i][0], z1[verbose_i][1], (verbose_i < n1-1) ? ", " : "}\n");
        printf ("x2 = {");
        for (int verbose_i = 0; verbose_i < n2; verbose_i++)
          printf ("{%luUL,%luUL}%s", x2[verbose_i][0], x2[verbose_i][1], (verbose_i < n2-1) ? ", " : "}\n");
        printf ("z2 = {");
        for (int verbose_i = 0; verbose_i < n2; verbose_i++)
          printf ("{%luUL,%luUL}%s", z2[verbose_i][0], z2[verbose_i][1], (verbose_i < n2-1) ? ", " : "}\n");
      }
#endif
  }

  if (verbose)
    {
      printf ("After canonicalizing:\nPj = [");
      for (i = 0; i < plan->s1; i++)
        printf ("%s(%lu:%lu)", (i>0) ? ", " : "",
                mod_get_ul (Pj_x[i], m), mod_get_ul (Pj_z[i], m));
      printf ("]\n");

      printf ("Pid = [");
      for (i = 0; i < plan->i1 - plan->i0; i++)
        printf ("%s(%lu:%lu)", (i>0) ? ", " : "",
                mod_get_ul (Pid_x[i], m), mod_get_ul (Pid_z[i], m));
      printf ("]\n");
      fflush(stdout);
    }

  /* Now process all the primes p = id - j, B1 < p <= B2 and multiply
     (id*P)_x - (j*P)_x to the accumulator */

  /* Init the accumulator to Pj[0], which contains the product of
     the Z-coordinates of all the precomputed points, except Pj_z[0]
     which is equal to P, and we know that one is coprime to the modulus.
     Maybe one of the others was zero (mod p) for some prime factor p. */

  mod_set (a, Pj_x[0], m);
  mod_init_noset0 (a_bk, m); /* Backup value of a, in case we get a == 0 */
  mod_set (a_bk, a, m);

  i = 0;
  l = 0;
  unsigned char j = plan->pairs[0];
  while (j != NEXT_PASS)
    {
      __asm__ volatile ("# ECM stage 2 loop here");
      while (j < NEXT_D && j < NEXT_PASS)
	{
	  mod_sub (t, Pid_x[i], Pj_x[j], m);
	  j = plan->pairs[++l];
	  mod_mul (a, a, t, m);
	}

#if ECM_BACKTRACKING
      /* See if we got a == 0. If yes, restore previous a value and
	 end stage 2. Let's hope not all factors were found since
	 the last d increase. */
      if (mod_is0 (a, m))
	{
	  mod_set (a, a_bk, m);
	  bt = 1;
	  break;
	}
      mod_set (a_bk, a, m); /* Save new a value */
#endif

      if (j == NEXT_D)
	{
	  i++;
	  j = plan->pairs[++l];
	  ASSERT (i < plan->i1 - plan->i0);
	}
    }

  if (verbose)
    printf ("Accumulator = %lu\n", mod_get_ul (a, m));

#if COUNT_ELLM_OPS
  printf("Stage 2 used %lu point additions and %lu point doublings\n",
         ellM_add_count, ellM_double_count);
#endif

  mod_set (r, a, m);

  /* Clear everything */

  for (i = 0; i < plan->s1; i++)
    {
      mod_clear (Pj_x[i], m);
      mod_clear (Pj_z[i], m);
    }
  free (Pj_x);
  Pj_x = NULL;
  free (Pj_z);
  Pj_z = NULL;

  for (i = 0; i < plan->i1 - plan->i0; i++)
    {
      mod_clear (Pid_x[i], m);
      mod_clear (Pid_z[i], m);
    }
  free (Pid_x);
  Pid_x = NULL;
  free (Pid_z);
  Pid_z = NULL;

  ec_point_clear (Pt, m);
  mod_clear (t, m);
  mod_clear (a, m);
  mod_clear (a_bk, m);
  return bt;
}

/* Stores any factor found in f_out (1 if no factor found).
   If back-tracking was used, returns 1, otherwise returns 0. */

int
ecm (modint_t f, const modulus_t m, const ecm_plan_t *plan)
{
  residue_t u, b;
  ec_point_t P, Pt;

  unsigned int i;
  int r, bt = 0;

  ec_point_init (P, m);
  mod_init (b, m);

  mod_intset_ul (f, 1UL);

  if (plan->parameterization == BRENT12)
    r = Brent12_curve_from_sigma (b, P, plan->parameter, m);
  else if (plan->parameterization == MONTY12)
    r = Montgomery12_curve_from_k (b, P, plan->parameter, m);
  else if (plan->parameterization == MONTY16)
    r = Montgomery16_curve_from_k (b, P, plan->parameter, m);
  else if (plan->parameterization == MONTYTWED12)
    r = ec_parameterization_Z6 (b, P, plan->parameter, EDW_ext, m);
  else
  {
    fprintf (stderr, "%s: unknown parameterization\n", __func__);
    abort();
  }

  if (r == 0)
  {
    printf ("# factor found during parameterization\n"); // XXX
    mod_gcd (f, P->x, m);
    mod_clear (b, m);
    ec_point_clear (P, m);
    return 0;
  }

#ifdef TRACE
  if (plan->parameterization & FULLMONTY)
  {
    /* FIXME need multiple precision print */
    residue_t A;
    mod_init (A, m);
    montgomery_A_from_b (A, b, m);
    printf ("%s: curve: B*Y^2*Z = X^3 + (%lu*4-2)*X^2*Z + X*Z^2\n"
            "%s: starting point (%lu::%lu)\n", __func__, mod_get_ul (A, m),
            __func__, mod_get_ul (P->x, m), mod_get_ul (P->z, m));
    mod_clear (A, m);
  }
#endif

  /* now start ecm */

  /* Do stage 1 */
  if (plan->parameterization & FULLMONTY)
    bytecode_prac_interpret_ec_montgomery (P, plan->bc, m, b);
  else if (plan->parameterization & FULLMONTYTWED)
    bytecode_mishmash_interpret_ec_mixed_repr (P, plan->bc, m, b);
    /* input is in Edwards (extended), output is in Montgomery form */

  /* Add prime 2 in the desired power. If a zero residue for the Z-coordinate is
   * encountered, we backtrack to previous point and stop.
   * NOTE: This is not as effective as I hoped. It prevents trivial
   * factorizations only if after processing the odd part of the stage 1
   * multiplier, the resulting point has power-of-2 order on E_p for all p|N.
   * If that were to happen, the point probably had that (presumably small on
   * most E_p) power-of-2 order during the last couple of primes processed in
   * the precomputed Lucas chain, and then quite likely the Lucas chain
   * incorrectly used an addition of identical points, causing the Z-coordinate
   * to become zero, leading to 0 (mod N) before we even get here.  For example,
   * using 10^6 composites from an RSA155 sieving experiment, without
   * backtracking we get N as the factor 456 times, with backtracking still 360
   * times.
   * TODO: this could probably be fixed by treating 3 separately, too, instead
   * of putting it in the precomputed Lucas chain. Then the probability that a
   * point of very small order on all E_p is encountered during the Lucas chain
   * is reduced, and so the probability of using curve addition erroneously.
   *
   * The following code assume P is in Montgomery form.
   */
  ec_point_init (Pt, m);
  ec_point_set (Pt, P, m);
  for (i = 0; i < plan->exp2; i++)
  {
    montgomery_dbl (P, P, m, b);
#if ECM_BACKTRACKING
    if (mod_is0 (P->z, m))
    {
      ec_point_set (P, Pt, m);
      bt = 1;
      break;
    }
    ec_point_set (Pt, P, m);
#endif
  }
  mod_gcd (f, P->z, m);

#if 0
  printf ("After stage 1, P = (%lu: :%lu), bt = %d, i = %d, exp2 = %d\n",
          mod_get_ul (P->x, m), mod_get_ul (P->z, m), bt, i, plan->exp2);
#endif

  /* Do stage 2 (for P in Montgomery form) */
  mod_init (u, m);
  if (bt == 0 && mod_intcmp_ul (f, 1UL) == 0 && plan->B1 < plan->stage2.B2)
  {
    bt = ecm_stage2 (u, P, &(plan->stage2), b, m);
    mod_gcd (f, u, m);
  }

  mod_clear (u, m);
  mod_clear (b, m);
  ec_point_clear (P, m);
  ec_point_clear (Pt, m);

  return bt;
}


/* -------------------------------------------------------------------------- */

/* Determine order of a point P on a curve, both defined by the parameter value
   as in ECM.
   If the group order is known to be == r (mod m), this can be supplied in
   the variables "known_r" and" known_m".
   Looks for i in Hasse interval so that i*P = O, has complexity O(m^(1/4)). */

//TODO parameter is an unsigned long
unsigned long
ell_pointorder (const residue_t parameter,
                const ec_parameterization_t parameterization,
                const unsigned long known_m, const unsigned long known_r,
                const modulus_t m, const int verbose)
{
  ec_point_t P, Pi, Pg, Q, *baby;
  residue_t x, a, d;
  unsigned long min, max, i, j, order, cof, p;
  unsigned long giant_step, giant_min, baby_len;
  modint_t tm;

  ASSERT (known_r < known_m);

  mod_init (a, m);
  ec_point_init (P, m);

  mod_intinit (tm);
  mod_getmod_int (tm, m);

  if (parameterization & FULLMONTY)
  {
    residue_t A, b;
    int r = 1;

    mod_init (b, m);
    mod_init (A, m);

    if (parameterization == BRENT12)
      r = Brent12_curve_from_sigma (b, P, mod_get_ul (parameter, m), m);
    else if (parameterization == MONTY12)
      r = Montgomery12_curve_from_k (b, P, mod_get_ul (parameter, m), m);
    else if (parameterization == MONTY16)
      r = Montgomery16_curve_from_k (b, P, mod_get_ul (parameter, m), m);

    if (r)
    {
      montgomery_A_from_b (A, b, m);

      if (verbose >= 2)
      {
        /* FIXME need multiple precision print */
        printf ("%s: Montgomery curve: B * Y^2 = X^3 + %lu * X^2 * Z + X * Z^2 "
                "(mod %lu)\n%s:   with point: (X::Z) = (%lu::%lu)\n", __func__,
                mod_get_ul (A, m), mod_intget_ul (tm), __func__,
                mod_get_ul (P->x, m), mod_get_ul (P->x, m));
      }

      r = weierstrass_from_montgomery (a, P, A, P, m);
    }

    mod_clear (b, m);
    mod_clear (A, m);

    if (r == 0)
    {
      mod_clear (a, m);
      ec_point_clear (P, m);
      mod_intclear (tm);
      return 0UL;
    }
  }
  else
  {
    fprintf (stderr, "ecm: Unknown parameterization\n");
    abort();
  }

  if (verbose >= 2)
  {
    /* FIXME need multiple precision print */
    printf ("%s: Weierstrass curve: y^2 = x^3 + %lu * x + b (mod %lu)\n"
            "%s:   with point: (x,y) = (%lu,%lu)\n", __func__,
            mod_get_ul (a, m), mod_intget_ul (tm), __func__,
            mod_get_ul (P->x, m), mod_get_ul (P->y, m));
  }

  mod_init (x, m);
  mod_init (d, m);
  ec_point_init (Pi, m);
  ec_point_init (Pg, m);

  /* FIXME deal with multiple precision modulus */
  i = (unsigned long) (2. * sqrt((double) mod_intget_ul(tm)));
  min = mod_intget_ul(tm) - i;
  max = mod_intget_ul(tm) + i;

  /* Giant steps visit values == r (mod m), baby steps values == 0 (mod m) */
  giant_step = ceil(sqrt(2.*(double)i / (double) known_m));
  /* Round up to multiple of m */
  giant_step = ((giant_step - 1) / known_m + 1) * known_m;

  /* We test Pi +- Pj, where Pi = (giant_min + i*giant_step), i >= 0,
     and Pj = j*P, 0 <= j <= giant_step / 2.
     To ensure we can find all values >= min, ensure
     giant_min <= min + giant_step / 2.
     We also want giant_min == r (mod m) */
  giant_min = ((min + giant_step / 2) / known_m) * known_m + known_r;
  if (giant_min > min + giant_step / 2)
    giant_min -= known_m;
  if (verbose >= 2)
    printf ("known_m = %lu, known_r = %lu, giant_step = %lu, "
            "giant_min = %lu\n", known_m, known_r, giant_step, giant_min);

  baby_len = giant_step / known_m / 2 + 1;
  baby = (ec_point_t *) malloc (baby_len * sizeof (ec_point_t));
  for (i = 0; i < baby_len; i++)
    ec_point_init (baby[i], m);

  ec_point_set (Pg, P, m);
  i = known_m;
  if (weierstrass_smul_ui (Pg, i, a, m) == 0) /* Pg = m*P for now */
    goto found_inf;

  if (1 < baby_len)
    ec_point_set (baby[1], Pg, m);

  if (2 < baby_len)
    {
      if (weierstrass_dbl (Pi, Pg, a, m) == 0)
        {
          i = 2 * known_m;
          goto found_inf;
        }
      ec_point_set (baby[2], Pi, m);
    }

  for (i = 3; i < baby_len; i++)
    {
      if (weierstrass_add (Pi, Pi, Pg, a, m) == 0)
        {
          i *= known_m;
          goto found_inf;
        }
      ec_point_set (baby[i], Pi, m);
    }

  /* Now compute the giant steps in [giant_min, giant_max] */
  i = giant_step;
  ec_point_set (Pg, P, m);
  if (weierstrass_smul_ui (Pg, i, a, m) == 0)
    goto found_inf;

  i = giant_min;
  ec_point_set (Pi, P, m);
  if (weierstrass_smul_ui (Pi, i, a, m) == 0)
    goto found_inf;

  while (i <= max + giant_step - 1)
    {
      /* Compare x-coordinate with stored baby steps. This makes it
         O(sqrt(p)) complexity, strictly speaking. */
      for (j = 1; j < baby_len; j++)
        if (mod_equal (Pi[0].x, baby[j]->x, m))
          {
            if (mod_equal (Pi[0].y, baby[j]->y, m))
              i -= j * known_m; /* Equal, so iP = jP and (i-j)P = 0 */
            else
              {
                mod_neg (Pi[0].y, Pi[0].y, m);
                if (mod_equal (Pi[0].y, baby[j]->y, m))
                  i += j * known_m; /* Negatives, so iP = -jP and (i+j)P = 0 */
                else
                  {
                    fprintf (stderr, "Matching x-coordinates, but y neither "
                             "equal nor negatives\n");
                    abort();
                  }
              }
            goto found_inf;
          }

      i += giant_step;
      if (!weierstrass_add (Pi, Pi, Pg, a, m))
        goto found_inf;
    }

  if (i > max)
  {
      fprintf (stderr, "ec_order: Error, no match found for p = %lu, "
               "min = %lu, max = %lu, giant_step = %lu, giant_min = %lu\n",
               mod_intget_ul(tm), min, max, giant_step, giant_min);
      abort ();
  }

found_inf:
  /* Check that i is a multiple of the order */
  ec_point_set (Pi, P, m);
  if (weierstrass_smul_ui (Pi, i, a, m) != 0)
    {
      modint_t tx1, ty1;
      mod_intinit (tx1);
      mod_intinit (ty1);
      mod_get_int (tx1, P[0].x, m);
      mod_get_int (ty1, P[0].y, m);
#ifndef MODMPZ_MAXBITS
      fprintf (stderr, "ec_order: Error, %ld*(%ld, %ld) (mod %ld) is "
               "not the point at infinity\n",
               i, tx1[0], ty1[0], tm[0]);
#endif
      mod_intclear (tx1);
      mod_intclear (ty1);
      return 0UL;
    }

  /* Ok, now we have some i so that ord(P) | i. Find ord(P).
     We know that ord(P) > 1 since P is not at infinity */

  /* For each prime factor of the order, reduce the exponent of
     that prime factor as far as possible */

  cof = order = i;
  for (p = 2; p * p <= cof; p += 1 + p%2)
    if (cof % p == 0)
      {
        ASSERT (order % p == 0);
        /* Remove all factors of p */
        for (order /= p, cof /= p; order % p == 0; order /= p)
          {
            ASSERT(cof % p == 0);
            cof /= p;
          }
        ASSERT (cof % p != 0);

        /* Add factors of p again one by one, stopping when we hit
           point at infinity */
        ec_point_set (Pi, P, m);
        if (weierstrass_smul_ui (Pi, order, a, m) != 0)
          {
            order *= p;
            while (weierstrass_smul_ui (Pi, p, a, m) != 0)
              order *= p;
          }
      }
  /* Now cof is 1 or a prime */
  if (cof > 1)
    {
      ec_point_set (Pi, P, m);
      ASSERT (order % cof == 0);
      if (weierstrass_smul_ui (Pi, order / cof, a, m) == 0)
        order /= cof;
    }


  /* One last check that order divides real order */
  ec_point_set (Pi, P, m);
  if (weierstrass_smul_ui (Pi, order, a, m) != 0)
    {
      modint_t tx1, ty1;
      mod_intinit (tx1);
      mod_intinit (ty1);
      mod_get_int (tx1, P[0].x, m);
      mod_get_int (ty1, P[0].y, m);
      fprintf (stderr, "ec_order: Error, final order %ld is wrong\n",
               order);
      mod_intclear (tx1);
      mod_intclear (ty1);
      abort ();
    }

  for (i = 0; i < giant_step; i++)
    ec_point_clear (baby[i], m);
  free (baby);
  baby = NULL;
  mod_clear (x, m);
  mod_clear (a, m);
  mod_clear (d, m);
  mod_intclear (tm);
  ec_point_clear (P, m);
  ec_point_clear (Pi, m);
  ec_point_clear (Pg, m);
  ec_point_clear (Q, m);

  return order;
}


/* Count points on curve using the Jacobi symbol. This has complexity O(m). */

unsigned long
ellM_curveorder_jacobi (residue_t A, residue_t x, modulus_t m)
{
  residue_t t, one;
  unsigned long order, i;
  modint_t tm;
  int bchar;

  mod_init_noset0 (t, m);
  mod_init_noset0 (one, m);
  mod_set1 (one, m);

  /* Compute x^3 + A*x^2 + x and see if it is a square */
  mod_set (t, x, m);
  mod_add (t, t, A, m);
  mod_mul (t, t, x, m);
  mod_add (t, t, one, m);
  mod_mul (t, t, x, m);
  bchar = mod_jacobi (t, m);
  ASSERT (bchar != 0);

  order = 2; /* One for (0, 0, 1), one for the point at infinity */
  /* FIXME deal with multi-word modulus */
  mod_getmod_int (tm, m);
  for (i = 1; mod_intcmp_ul (tm, i) > 0; i++)
    {
      mod_set_ul (x, i, m);
      mod_set (t, x, m);
      mod_add (t, t, A, m);
      mod_mul (t, t, x, m);
      mod_add (t, t, one, m);
      mod_mul (t, t, x, m);
      if (bchar == 1)
	order = order + (unsigned long) (1L + (long) mod_jacobi (t, m));
      else
	order = order + (unsigned long) (1L - (long) mod_jacobi (t, m));
	/* Brackets put like this to avoid signedness warning */
    }

  mod_clear (one, m);
  mod_clear (t, m);

  return order;
}

unsigned long
ell_curveorder (const ec_parameterization_t parameterization,
                const unsigned long parameter, const unsigned long m_par)
{
  residue_t b, A;
  modint_t lm;
  modulus_t m;
  unsigned long order;
  ec_point_t P;

  mod_init_noset0 (b, m);
  mod_init_noset0 (A, m);
  ec_point_init (P, m);

  mod_intset_ul (lm, m_par);
  mod_initmod_int (m, lm);

  if (parameterization == BRENT12)
  {
    if (Brent12_curve_from_sigma (b, P, parameter, m) == 0)
      return 0UL;
    montgomery_A_from_b (A, b, m);
  }
  else if (parameterization == MONTY12)
  {
    if (Montgomery12_curve_from_k (b, P, parameter, m) == 0)
      return 0UL;
    montgomery_A_from_b (A, b, m);
  }
  else
  {
    fprintf (stderr, "ec_curveorder: Unknown parameterization\n");
    abort();
  }
  order = ellM_curveorder_jacobi (A, P->x, m);

  mod_clear (b, m);
  mod_clear (A, m);
  ec_point_clear (P, m);
  return order;
}

