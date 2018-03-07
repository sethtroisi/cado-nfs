/* This file is not compilable as-is, it must be included by another file. */
#ifndef mod_init
  #error "One of the mod*_default.h headers must be included before this file"
#endif

#include "cado.h"
#include <math.h>
#include <inttypes.h>
#include "portability.h"
#include "facul_ecm.h"
#include "timing.h"

//#define ECM_COUNT_OPS /* define to print number of operations of ECM */
//#define ECM_TIMINGS /* define to print timings of ECM */

//#define ECM_TRACE

/* Define to 1 to make ellM_add() test if the two points are identical,
   and call ellM_double() if they are */
#ifndef ELLM_SAFE_ADD
#define ELLM_SAFE_ADD 0
#endif

/* Exported functions. */
#define ecm MOD_APPEND_TYPE(ecm)

/* Each static function are rename just before its declaration (there is no need
 * to rename them, but it's nice to have these functions distinguishable e.g. in
 * profiler output).
 */

#include "ec_arith_common.h"
#include "ec_arith_Edwards.h"
#include "ec_arith_Montgomery.h"
#include "ec_arith_Weierstrass.h"
#include "ec_parameterization.h"

/* Do we want backtracking when processing factors of 2 in E? */
#ifndef ECM_BACKTRACKING
/* Default is "yes." Set to 0 for "no." */
#define ECM_BACKTRACKING 1
#endif

#ifdef ECM_COUNT_OPS
static unsigned int _count_stage2_M;
#define ECM_COUNT_OPS_STAGE1_TOTAL_M EDWARDS_COUNT_OPS_M+MONTGOMERY_COUNT_OPS_M
#define ECM_COUNT_OPS_STAGE2_TOTAL_M MONTGOMERY_COUNT_OPS_M+_count_stage2_M
#define ECM_COUNT_OPS_RESET() do {                             \
      EDWARDS_COUNT_OPS_RESET(); MONTGOMERY_COUNT_OPS_RESET(); \
      _count_stage2_M = 0;                                     \
    } while (0)
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
#define bytecode_dbchain_interpret_edwards_internal MOD_APPEND_TYPE(bytecode_dbchain_interpret_edwards_internal)
static bytecode_const
bytecode_dbchain_interpret_edwards_internal (bytecode_const bc, ec_point_t *R,
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
      edwards_tpl (R[0], R[0], m, (i+1==pow3+pow2) ? TWISTED_EDWARDS_ext :
                                                          TWISTED_EDWARDS_proj);
    for (uint8_t i = 0; i < pow2; i++)
      edwards_dbl (R[0], R[0], m, (i+1==pow2) ? TWISTED_EDWARDS_ext :
                                                          TWISTED_EDWARDS_proj);

    ec_point_coord_type_t output;
    if (f)
    {
      uint8_t t;
      bytecode_elt_split_4_4 (&t, NULL, bc[1]);
      if (bc[1] == MISHMASH_FINAL || t == MISHMASH_PRAC_BLOCK)
        output = MONTGOMERY_xz;
      else
        output = TWISTED_EDWARDS_ext;
    }
    else
      output = TWISTED_EDWARDS_proj;

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
#define bytecode_precomp_interpret_edwards_internal MOD_APPEND_TYPE(bytecode_precomp_interpret_edwards_internal)
static bytecode_const
bytecode_precomp_interpret_edwards_internal (bytecode_const bc, ec_point_t *R,
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
          edwards_add (R[k], R[i], R[j], m, a ? TWISTED_EDWARDS_ext :
                                                          TWISTED_EDWARDS_proj);
        else
          edwards_sub (R[k], R[i], R[j], m, a ? TWISTED_EDWARDS_ext :
                                                          TWISTED_EDWARDS_proj);
        break;
      case PRECOMP_OP_DBL:
        bc++;
        pow2 = bytecode_elt_to_uint8 (*bc);
        for (uint8_t i = 0; i < pow2; i++)
          edwards_dbl (R[0], R[0], m, (i+1==pow2 && a) ? TWISTED_EDWARDS_ext :
                                                          TWISTED_EDWARDS_proj);
        if (k > 0)
          ec_point_set (R[k], R[0], m, a ? TWISTED_EDWARDS_ext :
                                                          TWISTED_EDWARDS_proj);
        break;
      case PRECOMP_OP_TPL:
        bc++;
        pow3 = bytecode_elt_to_uint8 (*bc);
        for (uint8_t i = 0; i < pow3; i++)
          edwards_tpl (R[0], R[0], m, (i+1==pow3 && a) ? TWISTED_EDWARDS_ext :
                                                          TWISTED_EDWARDS_proj);
        if (k > 0)
          ec_point_set (R[k], R[0], m, a ? TWISTED_EDWARDS_ext :
                                                          TWISTED_EDWARDS_proj);
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
#define bytecode_prac_interpret_montgomery_internal MOD_APPEND_TYPE(bytecode_prac_interpret_montgomery_internal)
static bytecode_const
bytecode_prac_interpret_montgomery_internal (bytecode_const bc, ec_point_t *R,
                                             const modulus_t m,
                                             const residue_t b)
{
  while (1)
  {
    int finished = 0;
    switch (*bc)
    {
      case PRAC_SWAP: /* [ = 's' ] Swap R[0], R[1] */
        ec_point_swap (R[0], R[1], m, MONTGOMERY_xz);
        break;
      case PRAC_SUBBLOCK_INIT: /* [ = 'i' ] Start of a sub-block */
        ec_point_set (R[1], R[0], m, MONTGOMERY_xz);
        ec_point_set (R[2], R[0], m, MONTGOMERY_xz);
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
        ec_point_set (R[0], R[4], m, MONTGOMERY_xz);
        break;
      case 2:
        montgomery_dadd (R[1], R[0], R[1], R[2], b, m);
        montgomery_dbl (R[0], R[0], m, b);
        break;
      case 3:
        montgomery_dadd (R[2], R[1], R[0], R[2], b, m);
        ec_point_swap (R[1], R[2], m, MONTGOMERY_xz);
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
        ec_point_swap (R[1], R[2], m, MONTGOMERY_xz);
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
        ec_point_swap (R[1], R[3], m, MONTGOMERY_xz);
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
        ec_point_set (R[2], R[1], m, MONTGOMERY_xz);
        montgomery_dbl (R[0], R[1], m, b);
        break;
      case 11:
        /* Combined rule 3 and rule 0 [=\x3s] */
        montgomery_dadd (R[2], R[1], R[0], R[2], b, m);
        /* (R[1],R[2],R[0]) := (R[0],R[1],R[2])  */
        ec_point_swap (R[1], R[2], m, MONTGOMERY_xz);
        ec_point_swap (R[0], R[1], m, MONTGOMERY_xz);
        break;
      case 12:
        /* Combined rule 3, then subchain end/start [=\x3fi] */
        montgomery_dadd (R[3], R[1], R[0], R[2], b, m);
        montgomery_dadd (R[2], R[0], R[3], R[1], b, m);
        ec_point_set (R[1], R[2], m, MONTGOMERY_xz);
        montgomery_dbl (R[0], R[2], m, b);
        break;
      case 13:
        /* Combined rule 3, swap, rule 3 and swap, merged a bit [=\x3s\x3s] */
        ec_point_set (R[3], R[1], m, MONTGOMERY_xz);
        montgomery_dadd (R[1], R[1], R[0], R[2], b, m);
        ec_point_set (R[2], R[0], m, MONTGOMERY_xz);
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

#define bytecode_prac_interpret_montgomery MOD_APPEND_TYPE(bytecode_prac_interpret_montgomery)
static void
bytecode_prac_interpret_montgomery (ec_point_t P, bytecode_const bc,
                                    const modulus_t m, const residue_t b)
{
  ec_point_t *R = NULL;
  unsigned int R_nalloc;

  if (bc == NULL)
    return;

  R_nalloc = 5; /* we need 5 points: 3 for PRAC + 2 temporary points */
  R = (ec_point_t *) malloc (R_nalloc * sizeof (ec_point_t));
  FATAL_ERROR_CHECK (R == NULL, "could not malloc R");
  for (unsigned int i = 0; i < R_nalloc; i++)
    ec_point_init (R[i], m);

  /* current point (here starting point) go into R[0] at init */
  ec_point_set (R[0], P, m, MONTGOMERY_xz);

  bytecode_prac_interpret_montgomery_internal (bc, R, m, b);

  /* output is in R[1] */
  ec_point_set (P, R[1], m, MONTGOMERY_xz);

  for (unsigned int i = 0; i < R_nalloc; i++)
    ec_point_clear (R[i], m);
  free (R);
}

/** Interpret the MISHMASH bytecode on Twisted Edwards and Montgomery curves **/
#define bytecode_mishmash_interpret_mixed_repr MOD_APPEND_TYPE(bytecode_mishmash_interpret_mixed_repr)
static void
bytecode_mishmash_interpret_mixed_repr (ec_point_t P, bytecode_const bc,
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
  ec_point_set (R[1], P, m, TWISTED_EDWARDS_ext);

  while (*bc != MISHMASH_FINAL)
  {
    uint8_t t, n;
    bytecode_elt_split_4_4 (&t, &n, *bc);

    if (n != 0)
      ec_point_set (R[0], R[n], m, TWISTED_EDWARDS_ext);

    if (t == MISHMASH_DBCHAIN_BLOCK)
      bc = bytecode_dbchain_interpret_edwards_internal (++bc, R, m);
    else if (t == MISHMASH_PRECOMP_BLOCK)
      bc = bytecode_precomp_interpret_edwards_internal (++bc, R, m);
    else if (t == MISHMASH_PRAC_BLOCK)
      bc = bytecode_prac_interpret_montgomery_internal (++bc, R, m, b);
    else /* unexpected bytecode */
    {
      printf ("Fatal error in %s at %s:%d -- unexpected bytecode 0x%02x\n",
              __func__, __FILE__, __LINE__, *bc);
      abort ();
    }

    bc++; /* go to next byte */
  }

  /* output (always in Montgomery form) is in R[1] */
  ec_point_set (P, R[1], m, MONTGOMERY_xz);

  for (unsigned int i = 0; i < R_nalloc; i++)
    ec_point_clear (R[i], m);
  free (R);
}


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/* Multiplies x[1] by z[2]*z[3]*z[4]...*z[n],
   x[2] by z[1]*z[3]*z[4]...*z[n] etc., generally
   x[i] by \prod_{1\leq j \leq n, j\neq i} z[j]
   Requires n > 1. Uses 4n-6 multiplications. */

#define common_z MOD_APPEND_TYPE(common_z)
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


#define ecm_stage2 MOD_APPEND_TYPE(ecm_stage_2)
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
    ec_point_set (ap1_0, P, m, MONTGOMERY_xz); /* ap1_0 = 1*P */
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
        ec_point_set (ap1_0, ap1_1, m, MONTGOMERY_xz);
        ec_point_set (ap1_1, Pt, m, MONTGOMERY_xz);
        i1 += 6;
	
        montgomery_dadd (Pt, ap5_1, P6, ap5_0, b, m);
        ec_point_set (ap5_0, ap5_1, m, MONTGOMERY_xz);
        ec_point_set (ap5_1, Pt, m, MONTGOMERY_xz);
        i5 += 6;
      }

    if ((unsigned) (i1 + i5) < plan->d)
      {
        if (i1 < i5)
          {
            montgomery_dadd (Pt, ap1_1, P6, ap1_0, b, m);
            ec_point_set (ap1_0, ap1_1, m, MONTGOMERY_xz);
            ec_point_set (ap1_1, Pt, m, MONTGOMERY_xz);
            i1 += 6;
          }
        else
          {
            montgomery_dadd (Pt, ap5_1, P6, ap5_0, b, m);
            ec_point_set (ap5_0, ap5_1, m, MONTGOMERY_xz);
            ec_point_set (ap5_1, Pt, m, MONTGOMERY_xz);
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
        ec_point_set (Pid, Pid1, m, MONTGOMERY_xz);
        ec_point_set (Pid1, Pt, m, MONTGOMERY_xz);
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
#ifdef ECM_TRACE
    mod_printf ("# TRACE: start %s with B1=%u and B2=%u for m = %" PRIMODu "\n",
            __func__, plan->B1, plan->stage2.B2, MOD_PRINT_MODULUS (m));
    printf ("# TRACE: using parameterization 0x%x with parameter %lu\n",
            plan->parameterization, plan->parameter);
#endif

  residue_t u, b;
  ec_point_t P, Pt;
  ec_point_coord_type_t param_output_type;

  unsigned int i;
  int r, bt = 0;

#ifdef ECM_TIMINGS
  uint64_t t0, param_dt, stage1_dt, stage2_dt;
  double ecm_starttime, ecm_dt;
  ecm_starttime = wct_seconds();
#endif

  ec_point_init (P, m);
  mod_init_noset0 (b, m);

  mod_intset_ul (f, 1UL);

  /**************************** parameterization ******************************/
#ifdef ECM_TIMINGS
  t0 = microseconds_thread ();
#endif

  if (plan->parameterization & FULLMONTY)
  {
    param_output_type = MONTGOMERY_xz;
    if (plan->parameterization == BRENT12)
      r = ec_parameterization_Brent_Suyama (b, P, plan->parameter, m);
    else if (plan->parameterization == MONTY12)
      r = ec_parameterization_Montgomery12 (b, P, plan->parameter, m);
    else /* if (plan->parameterization == MONTY16) */
      r = ec_parameterization_Montgomery16 (b, P, plan->parameter, m);
  }
  else if (plan->parameterization == MONTYTWED12)
  {
    param_output_type = plan->bc[1] & 0x80 ? MONTGOMERY_xz: TWISTED_EDWARDS_ext;
    r = ec_parameterization_Z6 (b, P, plan->parameter, param_output_type, m);
  }
  else
  {
    fprintf (stderr, "%s: unknown parameterization\n", __func__);
    abort();
  }

#ifdef ECM_TIMINGS
  param_dt = microseconds_thread() - t0;
#endif

  if (r == 0)
  {
    mod_gcd (f, P->x, m);
#ifdef ECM_TRACE
    mod_printf ("# TRACE: during parameterization, found factor %" PRIMODu "\n",
                MOD_PRINT_INT (f));
#endif
    mod_clear (b, m);
    ec_point_clear (P, m);
    return 0;
  }

#ifdef ECM_TRACE
  residue_t A;
  ec_point_t PM;
  ec_point_init (PM, m);
  mod_init_noset0 (A, m);
  montgomery_A_from_b (A, b, m);
  printf ("# TRACE: starting values:\n");

  if (param_output_type == MONTGOMERY_xz)
  {
    montgomery_curve_fprintf (stdout, "# TRACE:   ", A, P, m);
    printf ("# TRACE:                     = ");
    montgomery_point_fprintf_affine (stdout, P, m);
    fputc ('\n', stdout);
  }
  else
  {
    residue_t d;
    mod_init_noset0 (d, m);
    edwards_d_from_montgomery_A (d, A, m);
    edwards_ext_curve_fprintf (stdout, "# TRACE:   ", d, P, m);
    mod_clear (d, m);

    printf ("# TRACE:   Equivalent to Montgomery curve with:\n");
    montgomery_point_from_edwards_point (PM, P, 1, m);
    montgomery_curve_fprintf (stdout, "# TRACE:     ", A, PM, m);
  }
  mod_clear (A, m);
  ec_point_clear (PM, m);
#endif

  /******************************** stage 1 ***********************************/
#ifdef ECM_COUNT_OPS
  ECM_COUNT_OPS_RESET();
#endif
#ifdef ECM_TIMINGS
  t0 = microseconds_thread();
#endif

  /* output is always in Montgomery form */
  if (plan->parameterization & FULLMONTY)
    bytecode_prac_interpret_montgomery (P, plan->bc, m, b);
  else if (plan->parameterization & FULLMONTYTWED)
    bytecode_mishmash_interpret_mixed_repr (P, plan->bc, m, b);

#ifdef ECM_TRACE
  printf ("# TRACE: after stage 1 (without the powers of 2):\n");
  printf ("# TRACE:   output (X::Z) = ");
  ec_point_fprintf (stdout, P, MONTGOMERY_xz, m);
  printf ("\n# TRACE:                 = ");
  montgomery_point_fprintf_affine (stdout, P, m);
  fputc ('\n', stdout);
#endif
  
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
  ec_point_set (Pt, P, m, MONTGOMERY_xz);
  for (i = 0; i < plan->exp2; i++)
  {
    montgomery_dbl (P, P, m, b);
#if ECM_BACKTRACKING
    if (mod_is0 (P->z, m))
    {
      ec_point_set (P, Pt, m, MONTGOMERY_xz);
      bt = 1;
      break;
    }
    ec_point_set (Pt, P, m, MONTGOMERY_xz);
#endif
  }
  mod_gcd (f, P->z, m);

#ifdef ECM_TIMINGS
  stage1_dt = microseconds_thread() - t0;
#endif

#ifdef ECM_TRACE
  printf ("# TRACE: after stage 1 and %u power(s) of 2 (out of %u):\n", i,
          plan->exp2);
  printf ("# TRACE:   output (X::Z) = ");
  ec_point_fprintf (stdout, P, MONTGOMERY_xz, m);
  printf ("\n# TRACE:                 = ");
  montgomery_point_fprintf_affine (stdout, P, m);
#if ECM_BACKTRACKING
  mod_printf ("\n# TRACE:   backtracking was%s used", bt ? "" :" not");
#endif
  mod_printf ("\n# TRACE:   gcd = %" PRIMODu "\n", MOD_PRINT_INT (f));
#endif

#ifdef ECM_COUNT_OPS
  unsigned int tot_stage1_M = ECM_COUNT_OPS_STAGE1_TOTAL_M;
  fprintf (stderr, "COUNT OPS: for stage 1 with B1 = %u\n"
                   "COUNT OPS:   Montgomery: %3u dADD %3u DBL\n"
                   "COUNT OPS:      Edwards: %3u ADD  %3u DBL %3u TPL %3d M\n"
                   "COUNT OPS:        total: %5u M\n", plan->B1,
                   _count_montgomery_dadd, _count_montgomery_dbl,
                   _count_edwards_add, _count_edwards_dbl, _count_edwards_tpl,
                   _count_edwards_extraM, tot_stage1_M);
#endif

  /******************************** stage 2 ***********************************/
#ifdef ECM_COUNT_OPS
  ECM_COUNT_OPS_RESET();
#endif
#ifdef ECM_TIMINGS
  t0 = microseconds_thread();
#endif

  mod_init (u, m);
  if (bt == 0 && mod_intcmp_ul (f, 1UL) == 0 && plan->B1 < plan->stage2.B2)
  {
    bt = ecm_stage2 (u, P, &(plan->stage2), b, m);
    mod_gcd (f, u, m);
#ifdef ECM_TRACE
    mod_printf ("# TRACE: after stage 2, gcd=%" PRIMODu "\n", MOD_PRINT_INT(f));
#endif
  }
#ifdef ECM_TRACE
  else
    printf ("# TRACE: stage 2 not done\n");
#endif

#ifdef ECM_TIMINGS
  stage2_dt = microseconds_thread() - t0;
#endif
#ifdef ECM_COUNT_OPS
  unsigned int tot_stage2_M = ECM_COUNT_OPS_STAGE2_TOTAL_M;
  fprintf (stderr, "COUNT OPS: for stage 2 with B2 = %u\n"
                   "COUNT OPS:   Montgomery: %3u dADD %3u DBL\n"
                   "COUNT OPS:          mul: %3u M\n"
                   "COUNT OPS:        total: %5u M\n"
                   "COUNT OPS: ECM total: %5u M\n", plan->stage2.B2,
                   _count_montgomery_dadd, _count_montgomery_dbl,
                   _count_stage2_M, tot_stage2_M, tot_stage1_M + tot_stage2_M);
#endif

  mod_clear (u, m);
  mod_clear (b, m);
  ec_point_clear (P, m);
  ec_point_clear (Pt, m);

#ifdef ECM_TIMINGS
  ecm_dt = wct_seconds() - ecm_starttime;
  fprintf (stderr, "TIMINGS:   ECM took %.0f usec\n"
                   "TIMINGS:      param took %" PRIu64 "usec\n"
                   "TIMINGS:     stage1 took %" PRIu64 "usec\n"
                   "TIMINGS:     stage2 took %" PRIu64 "usec\n"
                   "TIMINGS:      extra time %.0fusec\n",
                   ecm_dt*1e6, param_dt, stage1_dt, stage2_dt,
                   ecm_dt*1e6-(double)(param_dt+stage1_dt+stage2_dt));
#endif

  return bt;
}
