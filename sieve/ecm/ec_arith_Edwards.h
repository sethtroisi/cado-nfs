#ifndef _EC_ARITH_EDWARDS_H_
#define _EC_ARITH_EDWARDS_H_

#include "ec_arith_common.h"

#ifdef ECM_COUNT_OPS
static unsigned int _count_edwards_add, _count_edwards_dbl, _count_edwards_tpl;
static int _count_edwards_extraM;
#define EDWARDS_COUNT_OPS_M _count_edwards_add*8 + _count_edwards_dbl*7 \
                            + _count_edwards_tpl*12 + _count_edwards_extraM
#define EDWARDS_COUNT_OPS_RESET() do { _count_edwards_extraM = 0;           \
      _count_edwards_add = _count_edwards_dbl = _count_edwards_tpl = 0;    \
    } while (0)
#endif

/* #define SAFE_TWISTED_EDWARDS_TO_MONTGOMERY */

/* Compute d = -(A-2)/(A+2). A and d can be the same variable. */
static inline void
edwards_d_from_montgomery_A (residue_t d, const residue_t A, const modulus_t m)
{
  residue_t (t);
  mod_init (t, m);

  mod_add_ul (t, A, 2UL, m);
  mod_inv (d, t, m);
  mod_sub_ul (t, t, 4UL, m); /* t = A-2 */
  mod_mul (d, d, t, m);
  mod_neg (d, d, m);

  mod_clear (t, m);
}

static inline void
edwards_neg (ec_point_t Q, const ec_point_t P, const modulus_t m)
{
  mod_neg (Q->x, P->x, m);
  mod_set (Q->y, P->y, m);
  mod_set (Q->z, P->z, m);
  mod_neg (Q->t, P->t, m);
}

/* - edwards_add (R:output_flag, P:edwards_ext, Q:edwards_ext, output_flag) */
/*     R <- P+Q */
/*     output_flag can be edwards_proj, edwards_ext or montgomery */
static inline void
edwards_add (ec_point_t R, const ec_point_t P, const ec_point_t Q,
             const modulus_t m, const ec_point_coord_type_t output_type)
{
  /* The "add-2008-hwcd-4" addition formulas */
  /* Cost: 8M + 8add + 2*2. */
  /* Cost: 8M + 6add dependent upon the first point. */
  /* Source: 2008 Hisil–Wong–Carter–Dawson */ 
  /* http://eprint.iacr.org/2008/522, Section 3.2. */
  
#ifdef ECM_COUNT_OPS
  _count_edwards_add++;
  if (output_type == TWISTED_EDWARDS_ext)
    _count_edwards_extraM++;
  else if (output_type == MONTGOMERY_xz)
    _count_edwards_extraM -= 4;
#endif

  residue_t t0, t1, A, B, C, D, E, F, G, H;
  
  /* ASSERT (output_flag != AFF); */

  mod_init_noset0 (t0, m);
  mod_init_noset0 (t1, m);
  mod_init_noset0 (A, m);
  mod_init_noset0 (B, m);
  mod_init_noset0 (C, m);
  mod_init_noset0 (D, m);
  mod_init_noset0 (E, m);
  mod_init_noset0 (F, m);
  mod_init_noset0 (G, m);
  mod_init_noset0 (H, m);

  mod_sub (t0, P->y, P->x, m);     // t0 := (Y1-X1)
  mod_add (t1, Q->y, Q->x, m);     // t1 := (Y2+X2)
  mod_mul (A, t0, t1, m);          // A := (Y1-X1)*(Y2+X2) 
  mod_add (t0, P->y, P->x, m);     // t0 := (X1+Y1)
  mod_sub (t1, Q->y, Q->x, m);     // t1 := (Y2-X2)
  mod_mul (B, t0, t1, m);          // B := (Y1+X1)*(Y2-X2)
  mod_mul (C, P->z, Q->t, m);      // C := Z1*T2
  mod_add (C, C, C, m);            // C := 2*Z1*T2
  mod_mul (D, P->t, Q->z, m);      // D := T1*Z2
  mod_add (D, D, D, m);            // D := 2*T1*Z2
  mod_add (E, D, C, m);            // E := D+C
  mod_sub (F, B, A, m);            // F := B-A
  mod_add (G, B, A, m);            // G := B+A
  mod_sub (H, D, C, m);            // H := D-C

  if (output_type != MONTGOMERY_xz)
  {
    mod_mul (R->x, E, F, m);     // X3 := E*F
    mod_mul (R->y, G, H, m);     // Y3 := G*H
    mod_mul (R->z, F, G, m);     // Z3 := F*G
    if (output_type == TWISTED_EDWARDS_ext)
      mod_mul (R->t, E, H, m);   // T3 := E*H
  }
  else
  {
#ifdef SAFE_TWISTED_EDWARDS_TO_MONTGOMERY
    mod_add (R->x, F, H, m);
    mod_mul (R->x, R->x, E, m);
    mod_sub (R->z, F, H, m);
    mod_mul (R->z, R->z, E, m);
#else
    /* CAUTION! */
    /* This may produce unstable results */
    /* But seems to "work" for our purpose */
    /* TODO: COMMENTS */
    mod_add (R->x, F, H, m);
    mod_sub (R->z, F, H, m);
#endif
  }
  
  mod_clear (t0, m);
  mod_clear (t1, m);
  mod_clear (A, m);
  mod_clear (B, m);
  mod_clear (C, m);
  mod_clear (D, m);
  mod_clear (E, m);
  mod_clear (F, m);
  mod_clear (G, m);
  mod_clear (H, m);  
}


/* - edwards_sub (R:output_flag, P:edwards_ext, Q:edwards_ext, output_flag) */
/*     R <- P-Q */
/*     output_flag can be edwards_proj, edwards_ext or montgomery */
static inline void
edwards_sub (ec_point_t R, const ec_point_t P, const ec_point_t Q,
	  const modulus_t m, const ec_point_coord_type_t output_type)
{
  ec_point_t QQ;
  ec_point_init (QQ, m);

  edwards_neg (QQ, Q, m);
  edwards_add (R, P, QQ, m, output_type);

  ec_point_clear (QQ, m);
}


/* - edwards_dbl (R:output_flag, P:edwards_proj, output_flag) */
/*     R <- 2*P */
/*     output_flag can be edwards_proj, edwards_ext */
static inline void
edwards_dbl (ec_point_t R, const ec_point_t P,
	  const modulus_t m, const ec_point_coord_type_t output_type)
{
  /* The "dbl-2008-hwcd" doubling formulas */
  /* Cost: 4M + 4S + 1*a + 6add + 1*2. */
  /* Source: 2008 Hisil–Wong–Carter–Dawson */
  /* http://eprint.iacr.org/2008/522, Section 3.3. */
    
  residue_t A, B, C, D, E, F, G, H;

#ifdef ECM_COUNT_OPS
  _count_edwards_dbl++;
  if (output_type == TWISTED_EDWARDS_ext)
    _count_edwards_extraM++;
#endif

  /* ASSERT (output_flag != AFF); */
  /* ASSERT (output_flag != MONTG); */
  
  mod_init_noset0 (A, m);
  mod_init_noset0 (B, m);
  mod_init_noset0 (C, m);
  mod_init_noset0 (D, m);
  mod_init_noset0 (E, m);
  mod_init_noset0 (F, m);
  mod_init_noset0 (G, m);
  mod_init_noset0 (H, m);

  mod_sqr (A, P->x, m);                // A := X1^2
  mod_sqr (B, P->y, m);                // B := Y1^2
  mod_sqr (C, P->z, m);                // C := Z1^2
  mod_add (C, C, C, m);                // C := 2*Z1^2
  mod_neg (D, A, m);                   // D := a*A  (a = -1)
  mod_add (E, P->x, P->y, m);          // E := (X1 + Y1)
  mod_sqr (E, E, m);                   // E := (X1 + Y1)^2
  mod_sub (E, E, A, m);                // E := (X1 + Y1)^2-A
  mod_sub (E, E, B, m);                // E := (X1 + Y1)^2-A-B
  mod_add (G, D, B, m);                // G := D+B
  mod_sub (F, G, C, m);                // F := G-C
  mod_sub (H, D, B, m);                // H := D-B
  mod_mul (R->x, E, F, m);             // X3 := E*F
  mod_mul (R->y, G, H, m);             // Y3 := G*H
  mod_mul (R->z, F, G, m);             // Z3 := F*G
  if (output_type == TWISTED_EDWARDS_ext)
    mod_mul (R->t, E, H, m);           // T3 := E*H
  
  mod_clear (A, m);
  mod_clear (B, m);
  mod_clear (C, m);
  mod_clear (D, m);
  mod_clear (E, m);
  mod_clear (F, m);
  mod_clear (G, m);
  mod_clear (H, m);  
}


/* - edwards_tpl (R:output_flag, P:edwards_proj, output_flag) */
/*     R <- 3*P */
/*     output_flag can be edwards_proj, edwards_ext */

/* The "tpl-2015-c" tripling formulas */
/* Cost: 11M + 3S + 1*a + 7add + 2*2. */
/* Source: 2015 Chuengsatiansup. */
/* https://hyperelliptic.org/EFD/g1p/auto-twisted-extended-1.html#tripling-tpl-2015-c */
static inline void
edwards_tpl (ec_point_t R, const ec_point_t P,
	  const modulus_t m, const ec_point_coord_type_t output_type)
{

  residue_t YY, aXX, Ap, B, xB, yB, AA, F, G, xE, yH, zF, zG;

#ifdef ECM_COUNT_OPS
  _count_edwards_tpl++;
  if (output_type == TWISTED_EDWARDS_ext)
    _count_edwards_extraM += 2;
#endif

  /* ASSERT (output_flag != AFF); */
  /* ASSERT (output_flag != MONTG); */
  
  mod_init_noset0 (YY, m);
  mod_init_noset0 (aXX, m);
  mod_init_noset0 (Ap, m);
  mod_init_noset0 (B, m);
  mod_init_noset0 (xB, m);
  mod_init_noset0 (yB, m);
  mod_init_noset0 (AA, m);
  mod_init_noset0 (F, m);
  mod_init_noset0 (G, m);
  mod_init_noset0 (xE, m);
  mod_init_noset0 (yH, m);
  mod_init_noset0 (zF, m);
  mod_init_noset0 (zG, m);

  mod_sqr (YY, P->y, m);                // YY := Y1^2
  mod_sqr (aXX, P->x, m);               // aXX := X1^2
  mod_neg (aXX, aXX, m);                // aXX := -X1^2
  mod_add (Ap, YY, aXX, m);             // Ap := YY+aXX 
  mod_sqr (B, P->z, m);                 // B := Z1^2
  mod_add (B, B, B, m);                 // B := 2*Z1^2
  mod_sub (B, B, Ap, m);                // B := 2*Z1^2-Ap
  mod_add (B, B, B, m);                 // B := 2*(2*Z1^2-Ap)
  mod_mul (xB, aXX, B, m);              // xB := aXX*B
  mod_mul (yB, YY, B, m);               // yB := YY*B
  mod_sub (AA, YY, aXX, m);             // AA := YY-aXX
  mod_mul (AA, Ap, AA, m);              // AA := Ap*(YY-aXX)
  mod_sub (F, AA, yB, m);               // F := AA-yB
  mod_add (G, AA, xB, m);               // G := AA+xB

  mod_add (xE, yB, AA, m);              // xE := yB+AA
  mod_mul (xE, P->x, xE, m);            // xE := X1*(yB+AA)
  mod_sub (yH, xB, AA, m);              // yH := xB-AA
  mod_mul (yH, P->y, yH, m);            // yH := Y1*(xB-AA)
  mod_mul (zF, P->z, F, m);             // zF := Z1*F
  
  switch (output_type)
    {
    case TWISTED_EDWARDS_proj:
      mod_mul (R->x, xE, F, m);         // X3 := xE*F
      mod_mul (R->y, yH, G, m);         // Y3 := yH*G
      mod_mul (R->z, zF, G, m);         // Z3 := zF*G
      break;
    case TWISTED_EDWARDS_ext:
      mod_mul (zG, P->z, G, m);         // zG := Z1*G
      mod_mul (R->x, xE, zF, m);        // X3 := xE*zF
      mod_mul (R->y, yH, zG, m);        // Y3 := yH*zG
      mod_mul (R->z, zF, zG, m);        // Z3 := zF*zG
      mod_mul (R->t, xE, yH, m);        // T3 := xE*yH
      break;
    default :
      break;
    }
  
  mod_clear (YY, m);
  mod_clear (aXX, m);
  mod_clear (Ap, m);
  mod_clear (B, m);
  mod_clear (xB, m);
  mod_clear (yB, m);
  mod_clear (AA, m);
  mod_clear (F, m);
  mod_clear (G, m);
  mod_clear (xE, m);
  mod_clear (yH, m);
  mod_clear (zF, m);
  mod_clear (zG, m);
}


/* ------------------------------------------------------------------------- */

/* Computes R = [e]P (mod m)  */
MAYBE_UNUSED 
static inline void
edwards_smul_ui (ec_point_t R, const ec_point_t P, const unsigned long e, 
		 const modulus_t m)
{
  unsigned long j;
  long k;
  ec_point_t T, Pe;
  
  if (e == 0UL)
    {
      mod_set0 (R->x, m);
      mod_set1 (R->y, m);
      mod_set1 (R->z, m);
      return;
    }
  
  if (e == 1UL)
    {
      ec_point_set (R, P, m);
      return;
    }
  
  if (e == 2UL)
    {
      edwards_dbl (R, P, m, TWISTED_EDWARDS_ext);
      return;
    }

  if (e == 4UL)
    {
      edwards_dbl (R, P, m, TWISTED_EDWARDS_ext);
      edwards_dbl (R, R, m, TWISTED_EDWARDS_ext);
      return;
    }

  ec_point_init (T, m);
  ec_point_init (Pe, m);
  ec_point_set (Pe, P, m);
  
  /* basic double-and-add */
  
  mod_set0 (T->x, m);
  mod_set0 (T->t, m);
  mod_set1 (T->y, m);
  mod_set1 (T->z, m);
  
  k = CHAR_BIT * sizeof(e) - 1;
  j = (1UL << k);
  
  while(k-- >= 0)
    {
      edwards_dbl (T, T, m, TWISTED_EDWARDS_ext);
      if (j & e)
	edwards_add (T, T, Pe, m, TWISTED_EDWARDS_ext);
      j >>= 1;
    }

  ec_point_set (R, T, m);
   
  ec_point_clear (T, m);
  ec_point_clear (Pe, m);
}

#endif /* _EC_ARITH_EDWARDS_H_ */
