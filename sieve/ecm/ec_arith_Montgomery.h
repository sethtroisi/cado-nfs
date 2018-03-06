#ifndef _EC_ARITH_MONTGOMERY_H_
#define _EC_ARITH_MONTGOMERY_H_

#include "ec_arith_common.h"

#ifdef ECM_COUNT_OPS
static unsigned int _count_montgomery_dadd, _count_montgomery_dbl;
#define MONTGOMERY_COUNT_OPS_M _count_montgomery_dadd*6+_count_montgomery_dbl*5
#define MONTGOMERY_COUNT_OPS_RESET() do {                   \
      _count_montgomery_dadd = _count_montgomery_dbl = 0;   \
    } while (0)
#endif

/************************ Montgomery elliptic curves **************************/

/* Equation:
 *    B*Y^2*Z = X^3 + A*X^2*Z + X*Z^2
 *
 * Constant needed in computation: b = (A+2)/4
 *
 * Implemented functions:
 *   - montgomery_dadd (R:MONTGOMERY_xz, P:MONTGOMERY_xz, Q:MONTGOMERY_xz,
 *                                                              D:MONTGOMERY_xz)
 *        R <- P+Q (where D=P-Q)
 *   - montgomery_dbl (R:MONTGOMERY_xz, P:MONTGOMERY_xz)
 *        R <- 2*P
 */

/* Compute A = 4*b-2. A and b can be the same variable. */
static inline void
montgomery_A_from_b (residue_t A, const residue_t b, const modulus_t m)
{
  mod_add (A, b, b, m);    /* A <- b+b = 2b */
  mod_add (A, A, A, m);    /* A <- 4b */
  mod_sub_ul (A, A, 2, m); /* A <- 4b-2 */
}

static inline void
montgomery_point_set_zero (ec_point_t P, const modulus_t m)
{
  mod_set0 (P->x, m);
  mod_set0 (P->z, m);
}

/* P and Q can be the same variables */
static inline int
montgomery_point_to_affine (ec_point_t Q, ec_point_t P, const modulus_t m)
{
  residue_t t;
  mod_init_noset0 (t, m);

  int ret = mod_inv (t, P->z, m);
  if (ret)
  {
    mod_mul (Q->x, P->x, t, m);
    mod_set1 (Q->z, m);
  }

  mod_clear (t, m);
  return ret;
}

static inline void
montgomery_point_fprintf_affine (FILE *out, ec_point_t P, const modulus_t m)
{
  ec_point_t Paff;
  ec_point_init (Paff, m);
  int ret = montgomery_point_to_affine (Paff, P, m);
  if (ret)
    ec_point_fprintf (out, Paff, MONTGOMERY_xz, m);
  else
    fprintf (out, "(not a affine point)");
  ec_point_clear (Paff, m);
}

static inline void
montgomery_point_from_edwards_point (ec_point_t PM, ec_point_t PE, int out_aff,
                                     const modulus_t m)
{
  mod_add (PM->x, PE->z, PE->y, m);
  mod_sub (PM->z, PE->z, PE->y, m);
  if (out_aff)
    montgomery_point_to_affine (PM, PM, m);
}

/* Computes Q=2P, with 5 muls (3 muls and 2 squares) and 4 add/sub.
 *    - m : modulus
 *    - b = (A+2)/4 mod m
 * It is permissible to let P and Q use the same memory.
 */
static inline void
montgomery_dbl (ec_point_t Q, const ec_point_t P, const modulus_t m,
             const residue_t b)
{
  residue_t u, v, w;

#if ELLM_SAFE_ADD
  if (mod_is0 (P->z, m))
    {
      ASSERT (mod_is0 (P->x, m));
    }
#endif

#ifdef ECM_COUNT_OPS
  _count_montgomery_dbl++;
#endif

  mod_init_noset0 (u, m);
  mod_init_noset0 (v, m);
  mod_init_noset0 (w, m);

  mod_add (u, P->x, P->z, m);
  mod_sqr (u, u, m);          /* u = (x + z)^2 */
  mod_sub (v, P->x, P->z, m);
  mod_sqr (v, v, m);          /* v = (x - z)^2 */
  mod_mul (Q->x, u, v, m);    /* x2 = (x^2 - z^2)^2 */
  mod_sub (w, u, v, m);       /* w = 4 * x * z */
  mod_mul (u, w, b, m);       /* u = x * z * (A + 2) */
  mod_add (u, u, v, m);       /* u = x^2 + x * z * A + z^2 */
  mod_mul (Q->z, w, u, m);    /* Q_z = (4xz) * (x^2 + xzA + z^2) */

#if ELLM_SAFE_ADD
  if (mod_is0 (Q->z, m))
    mod_set0 (Q->x, m);
#endif

  mod_clear (w, m);
  mod_clear (v, m);
  mod_clear (u, m);
}


/* adds P and Q and puts the result in R,
 * using 6 muls (4 muls and 2 squares), and 6 add/sub.
 * One assumes that Q-R=D or R-Q=D.
 * This function assumes that P !~= Q, i.e. that there is
 * no t!=0 so that P->x = t*Q->x and P->z = t*Q->z, for otherwise the result
 * is (0:0) although it shouldn't be (which actually is good for factoring!).
 *
 * R may be identical to P, Q and/or D.
 */
static inline void
montgomery_dadd (ec_point_t R, const ec_point_t P, const ec_point_t Q,
          const ec_point_t D, MAYBE_UNUSED const residue_t b,
          const modulus_t m)
{
  residue_t u, v, w;

#ifdef ECM_COUNT_OPS
  _count_montgomery_dadd++;
#endif

#if ELLM_SAFE_ADD
  /* Handle case where at least one input point is point at infinity */
  if (mod_is0 (P->z, m))
    {
      ASSERT (mod_is0 (P->x, m));
      ellM_set (R, Q, m);
      return;
    }
  if (mod_is0 (Q->z, m))
    {
      ASSERT (mod_is0 (Q->x, m));
      ellM_set (R, P, m);
      return;
    }
#endif

  mod_init_noset0 (u, m);
  mod_init_noset0 (v, m);
  mod_init_noset0 (w, m);

  mod_sub (u, P->x, P->z, m);
  mod_add (v, Q->x, Q->z, m);
  mod_mul (u, u, v, m);      /* u = (Px-Pz)*(Qx+Qz) */
  mod_add (w, P->x, P->z, m);
  mod_sub (v, Q->x, Q->z, m);
  mod_mul (v, w, v, m);      /* v = (Px+Pz)*(Qx-Qz) */
  mod_add (w, u, v, m);      /* w = 2*(Qx*Px - Qz*Pz)*/
  mod_sub (v, u, v, m);      /* v = 2*(Qz*Px - Qx*Pz) */
#if ELLM_SAFE_ADD
  /* Check if v == 0, which happens if P=Q or P=-Q.
     If P=-Q, set result to point at infinity.
     If P=Q, use ellM_double() instead.
     This test only works if P=Q on the pseudo-curve modulo N, i.e.,
     if N has several prime factors p, q, ... and P=Q or P=-Q on E_p but
     not on E_q, this test won't notice it. */
  if (mod_is0 (v, m))
    {
      mod_clear (w, m);
      mod_clear (v, m);
      mod_clear (u, m);
      /* Test if difference is point at infinity */
      if (mod_is0 (D->z, m))
        {
          ASSERT (mod_is0 (D->x, m));
          montgomery_dbl (R, P, m, b); /* Yes, points are identical, use doubling */
        }
      else
        montgomery_point_set_zero (R, m); /* Set result to point at infinity */
      return;
    }
#endif
  mod_sqr (w, w, m);          /* w = 4*(Qx*Px - Qz*Pz)^2 */
  mod_sqr (v, v, m);          /* v = 4*(Qz*Px - Qx*Pz)^2 */
  mod_set (u, D->x, m);       /* save D->x */
  mod_mul (R->x, w, D->z, m); /* may overwrite D->x */
  mod_mul (R->z, u, v, m);

  mod_clear (w, m);
  mod_clear (v, m);
  mod_clear (u, m);
}


/* (x:z) <- e*(x:z) (mod p) */
static inline void
montgomery_smul_ui (ec_point_t R, const ec_point_t P, unsigned long e,
		    const modulus_t m, const residue_t b)
{
  unsigned long l, n;
  ec_point_t t1, t2;

  if (e == 0UL)
    {
      mod_set0 (R->x, m);
      mod_set0 (R->z, m);
      return;
    }

  if (e == 1UL)
    {
      ec_point_set (R, P, m, MONTGOMERY_xz);
      return;
    }

  if (e == 2UL)
    {
      montgomery_dbl (R, P, m, b);
      return;
    }

  if (e == 4UL)
    {
      montgomery_dbl (R, P, m, b);
      montgomery_dbl (R, R, m, b);
      return;
    }

  ec_point_init (t1, m);

  if (e == 3UL)
    {
      montgomery_dbl (t1, P, m, b);
      montgomery_dadd (R, t1, P, P, b, m);
      ec_point_clear (t1, m);
      return;
    }

  ec_point_init (t2, m);
  e --;

  /* compute number of steps needed: we start from (1,2) and go from
     (i,i+1) to (2i,2i+1) or (2i+1,2i+2) */
  for (l = e, n = 0; l > 1; n ++, l /= 2) ;

  /* start from P1=P, P2=2P */
  ec_point_set (t1, P, m, MONTGOMERY_xz);
  montgomery_dbl (t2, t1, m, b);

  while (n--)
    {
      if ((e >> n) & 1) /* (i,i+1) -> (2i+1,2i+2) */
        {
          /* printf ("(i,i+1) -> (2i+1,2i+2)\n"); */
          montgomery_dadd (t1, t1, t2, P, b, m);
          montgomery_dbl (t2, t2, m, b);
        }
      else /* (i,i+1) -> (2i,2i+1) */
        {
          /* printf ("(i,i+1) -> (2i,2i+1)\n"); */
          montgomery_dadd (t2, t1, t2, P, b, m);
          montgomery_dbl (t1, t1, m, b);
        }
    }

  ec_point_set (R, t2, m, MONTGOMERY_xz);

  ec_point_clear (t1, m);
  ec_point_clear (t2, m);
}

#endif /* _EC_ARITH_MONTGOMERY_H_ */
