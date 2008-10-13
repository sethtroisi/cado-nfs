#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ecm.h"
#if defined(MODREDCUL)
#include "modredc_ul.h"
#include "modredc_ul_default.h"
#define ecm ecm_ul
#define ell_pointorder ell_pointorder_ul
#define ellM_curveorder_jacobi ellM_curveorder_jacobi_ul
#define ell_curveorder ell_curveorder_ul
/* The next ones are static so there's no need to rename them, but it's 
   nice to have these functions distinguishable e.g. in profiler output */
#define ellM_double ellM_double_ul
#define ellM_add ellM_add_ul
#define ellM_interpret_bytecode ellM_interpret_bytecode_ul
#define ecm_stage2 ecm_stage2_ul
#elif defined(MODREDC15UL)
#include "modredc_15ul.h"
#include "modredc_15ul_default.h"
#define ecm ecm_15ul
#define ell_pointorder ell_pointorder_15ul
#define ellM_curveorder_jacobi ellM_curveorder_jacobi_15ul
#define ell_curveorder ell_curveorder_15ul
#define ellM_double ellM_double_15ul
#define ellM_add ellM_add_15ul
#define ellM_interpret_bytecode ellM_interpret_bytecode_15ul
#define ecm_stage2 ecm_stage2_15ul
#else
#error Please define MODREDCUL or MODREDC15UL
#endif

typedef struct {residue_t x, z;} __ellM_point_t;
typedef __ellM_point_t ellM_point_t[1];

static inline void
ellM_init (ellM_point_t P, const modulus_t m)
{
  mod_init (P->x, m);
  mod_init (P->z, m);
}

static inline void
ellM_clear (ellM_point_t P, const modulus_t m)
{
  mod_clear (P->x, m);
  mod_clear (P->z, m);
}

static inline void
ellM_set (ellM_point_t Q, const ellM_point_t P, const modulus_t m)
{
  mod_set (Q->x, P->x, m);
  mod_set (Q->z, P->z, m);
}

static inline void
ellM_swap (ellM_point_t Q, ellM_point_t P, const modulus_t m)
{
  mod_swap (Q->x, P->x, m);
  mod_swap (Q->z, P->z, m);
}

/* computes Q=2P, with 5 muls (3 muls and 2 squares) and 4 add/sub.
     - m : number to factor
     - b : (a+2)/4 mod n
  It is permissible to let x1, z1 and x2, z2 use the same memory. */

static void
ellM_double (ellM_point_t Q, const ellM_point_t P, const modulus_t m, 
             const residue_t b)
{
  residue_t u, v, w;

  mod_init_noset0 (u, m);
  mod_init_noset0 (v, m);
  mod_init_noset0 (w, m);

  mod_add (u, P->x, P->z, m);
  mod_mul (u, u, u, m);   /* u = (x1 + z1)^2 */
  mod_sub (v, P->x, P->z, m);
  mod_mul (v, v, v, m);   /* v = (x1 - z1)^2 */
  mod_mul (Q->x, u, v, m);  /* x2 = (x1^2 - z1^2)^2 */
  mod_sub (w, u, v, m);   /* w = 4 * x1 * z1 */
  mod_mul (u, w, b, m);   /* u = x1 * z1 * (A + 2) */
  mod_add (u, u, v, m);
  mod_mul (Q->z, w, u, m);

  mod_clear (w, m);
  mod_clear (v, m);
  mod_clear (u, m);
}


/* For Weierstrass coordinates. Returns 1 if doubling worked normally, 
   0 if the result is point at infinity */

static int
ellW_double (residue_t x3, residue_t y3, residue_t x1, residue_t y1,
	     residue_t a, const modulus_t m)
{
  residue_t lambda, u, v;

  mod_init_noset0 (lambda, m);
  mod_init_noset0 (u, m);
  mod_init_noset0 (v, m);

  mod_mul (u, x1, x1, m);
  mod_add (v, u, u, m);
  mod_add (v, v, u, m);
  mod_add (v, v, a, m); /* 3x^2 + a */
  mod_add (u, y1, y1, m);
  if (mod_inv (u, u, m) == 0)    /* 1/(2*y) */
  {
      mod_clear (v, m);
      mod_clear (u, m);
      mod_clear (lambda, m);
      return 0; /* y was 0  =>  result is point at infinity */
  }
  mod_mul (lambda, u, v, m);
  mod_mul (u, lambda, lambda, m);
  mod_sub (u, u, x1, m);
  mod_sub (u, u, x1, m);    /* x3 = u = lambda^2 - 2*x */
  mod_sub (v, x1, u, m);
  mod_mul (v, v, lambda, m);
  mod_sub (y3, v, y1, m);
  mod_set (x3, u, m);
  
  mod_clear (v, m);
  mod_clear (u, m);
  mod_clear (lambda, m);
  return 1;
}


/* adds P and Q and puts the result in R,
     using 6 muls (4 muls and 2 squares), and 6 add/sub.
   One assumes that Q-R=D or R-Q=D.
   This function assumes that P !~= Q, i.e. that there is 
   no t!=0 so that P->x = t*Q->x and P->z = t*Q->z, for otherwise the result 
   is (0:0) although it shouldn't be (which actually is good for factoring!).

   R may be identical to P, Q and/or D. */

static void
ellM_add (ellM_point_t R, const ellM_point_t P, const ellM_point_t Q, 
          const ellM_point_t D, const modulus_t m)
{
  residue_t u, v, w;

  mod_init_noset0 (u, m);
  mod_init_noset0 (v, m);
  mod_init_noset0 (w, m);

  mod_sub (u, P->x, P->z, m);
  mod_add (v, Q->x, Q->z, m);
  mod_mul (u, u, v, m);
  mod_add (w, P->x, P->z, m);
  mod_sub (v, Q->x, Q->z, m);
  mod_mul (v, w, v, m);
  mod_add (w, u, v, m);
  mod_sub (v, u, v, m);
  mod_mul (w, w, w, m);
  mod_mul (v, v, v, m);
  mod_set (u, D->x, m); /* save D->x */
  mod_mul (R->x, w, D->z, m);
  mod_mul (R->z, u, v, m);

  mod_clear (w, m);
  mod_clear (v, m);
  mod_clear (u, m);
}


/* Adds two points (x2, y2) and (x1, y1) on the curve y^2 = x^3 + a*x + b
   in Weierstrass coordinates and puts result in (x3, y3). 
   Returns 1 if the addition worked (i.e. the modular inverse existed) 
   and 0 otherwise (resulting point is point at infinity) */

static int
ellW_add3 (residue_t x3, residue_t y3, residue_t x2, residue_t y2, 
           residue_t x1, residue_t y1, residue_t a, const modulus_t m)
{
  residue_t lambda, u, v;
  int r;

  mod_init_noset0 (u, m);
  mod_init_noset0 (v, m);

  mod_sub (u, y2, y1, m);
  mod_sub (v, x2, x1, m);
  r = mod_inv (v, v, m);
  if (r == 0)
  {
      /* Maybe we were trying to add two identical points? If so,
         use the ellW_double() function instead */
      if (mod_equal (x1, x2, m) && mod_equal (y1, y2, m))
	  r = ellW_double (x3, y3, x1, y1, a, m);
      else
	  r = 0; /* No, the points were negatives of each other */
  }
  else
  {
      mod_mul (lambda, u, v, m);
      mod_mul (u, lambda, lambda, m);
      mod_sub (u, u, x1, m);
      mod_sub (u, u, x2, m);    /* x3 = u = lambda^2 - x1 - x2 */
      mod_sub (v, x1, u, m);
      mod_mul (v, v, lambda, m);
      mod_sub (y3, v, y1, m);
      mod_set (x3, u, m);
      r = 1;
  }

  mod_clear (v, m);
  mod_clear (u, m);
  return r;
}


/* (x:z) <- e*(x:z) (mod p)
   Assumes e >= 5.
*/
static void
ellM_mul_ul (ellM_point_t R, const ellM_point_t P, unsigned long e, 
             const modulus_t m, const residue_t b)
{
  unsigned long l, n;
  ellM_point_t t1, t2;

  if (e == 0UL)
    {
      mod_set0 (R[0].x, m);
      mod_set0 (R[0].z, m);
      return;
    }

  if (e == 1UL)
    {
      ellM_set (R, P, m);
      return;
    }
  
  if (e == 2UL)
    {
      ellM_double (R, P, m, b);
      return;
    }

  if (e == 4UL)
    {
      ellM_double (R, P, m, b);
      ellM_double (R, R, m, b);
      return;
    }

  ellM_init (t1, m);

  if (e == 3UL)
    {
      ellM_double (t1, P, m, b);
      ellM_add (R, t1, P, P, m);
      ellM_clear (t1, m);
      return;
    }

  ellM_init (t2, m);
  e --;

  /* compute number of steps needed: we start from (1,2) and go from
     (i,i+1) to (2i,2i+1) or (2i+1,2i+2) */
  for (l = e, n = 0; l > 1; n ++, l /= 2) ;

  /* start from P1=P, P2=2P */
  ellM_set (t1, P, m);
  ellM_double (t2, t1, m, b);

  while (n--)
    {
      if ((e >> n) & 1) /* (i,i+1) -> (2i+1,2i+2) */
        {
          /* printf ("(i,i+1) -> (2i+1,2i+2)\n"); */
          ellM_add (t1, t1, t2, P, m);
          ellM_double (t2, t2, m, b);
        }
      else /* (i,i+1) -> (2i,2i+1) */
        {
          /* printf ("(i,i+1) -> (2i,2i+1)\n"); */
          ellM_add (t2, t1, t2, P, m);
          ellM_double (t1, t1, m, b);
        }
    }
  
  ellM_set (R, t2, m);

  ellM_clear (t1, m);
  ellM_clear (t2, m);
}

static int
ellW_mul_ui (residue_t x, residue_t y, const unsigned long e, residue_t a, 
	     const modulus_t m)
{
  unsigned long i;
  residue_t xt, yt;
  int tfinite; /* Nonzero iff (xt, yt) is NOT point at infinity */

  if (e == 0)
    return 0; /* signal point at infinity */

  mod_init_noset0 (xt, m);
  mod_init_noset0 (yt, m);

  i = ~(0UL);
  i -= i/2;   /* Now the most significant bit of i is set */
  while ((i & e) == 0)
    i >>= 1;

  mod_set (xt, x, m);
  mod_set (yt, y, m);
  tfinite = 1;
  i >>= 1;

  while (i > 0)
  {
      if (tfinite)
        tfinite = ellW_double (xt, yt, xt, yt, a, m);
      if (e & i)
      {
	  if (tfinite)
	      tfinite = ellW_add3 (xt, yt, x, y, xt, yt, a, m);
	  else
	  {
	      mod_set (xt, x, m);
	      mod_set (yt, y, m);
	      tfinite = 1;
	  }
      }
      i >>= 1;
  }

  if (tfinite)
  {
      mod_set (x, xt, m);
      mod_set (y, yt, m);
  }
  mod_clear (yt, m);
  mod_clear (xt, m);

  return tfinite;
}


/* Interpret the "l" bytes of bytecode located at "code" and do the 
   corresponding elliptic curve operations on (x::z) */

static void
ellM_interpret_bytecode (ellM_point_t P, const char *code,
			 const unsigned long l, const modulus_t m, 
			 const residue_t b)
{
  unsigned long i;
  ellM_point_t A, B, C, t, t2;
  
  ellM_init (A, m);
  ellM_init (B, m);
  ellM_init (C, m);
  ellM_init (t, m);
  ellM_init (t2, m);

  ellM_set (A, P, m);

  for (i = 0; i < l; i++)
    {
      switch (code[i])
        {
          case 10: /* Init of subchain, B=A, C=A, A=2*A */
            ellM_set (B, A, m);
            ellM_set (C, A, m);
            ellM_double (A, A, m, b);
            break;
          case 0: /* Swap A, B */
            ellM_swap (A, B, m);
            break;
          case 1:
            ellM_add (t, A, B, C, m);
            ellM_add (t2, t, A, B, m);
            ellM_add (B, B, t, A, m);
            ellM_set (A, t2, m);
            break;
          case 2:
            ellM_add (B, A, B, C, m);
            ellM_double (A, A, m, b);
            break;
          case 3:
            ellM_add (C, B, A, C, m);
            ellM_swap (B, C, m);
            break;
          case 4:
            ellM_add (B, B, A, C, m);
            ellM_double (A, A, m, b);
            break;
          case 5:
            ellM_add (C, C, A, B, m);
            ellM_double (A, A, m, b);
            break;
          case 6:
            ellM_double (t, A, m, b);
            ellM_add (t2, A, B, C, m);
            ellM_add (A, t, A, A, m);
            ellM_add (C, t, t2, C, m);
            ellM_swap (B, C, m);
            break;
          case 7:
            ellM_add (t, A, B, C, m);
            ellM_add (B, t, A, B, m);
            ellM_double (t, A, m, b);
            ellM_add (A, A, t, A, m);
            break;
          case 8:
            ellM_add (t, A, B, C, m);
            ellM_add (C, C, A, B, m);
            ellM_swap (B, t, m);
            ellM_double (t, A, m, b);
            ellM_add (A, A, t, A, m);
            break;
          case 9:
            ellM_add (C, C, B, A, m);
            ellM_double (B, B, m, b);
            break;
          case 11:
            ellM_add (A, A, B, C, m); /* Final add */
            break;
#if 0
	    /* p=2 is handled outside of the byte code interpreter now */
          case 12:
            ellM_double (A, A, m, b); /* For p=2 */
            break;
#endif
          default:
            abort ();
        }
    }

  ellM_set (P, A, m);

  ellM_clear (A, m);
  ellM_clear (B, m);
  ellM_clear (C, m);
  ellM_clear (t, m);
  ellM_clear (t2, m);
}


/* Produces curve in Montgomery form from sigma value.
   Return 1 if it worked, 0 if a modular inverse failed */

static int
Brent12_curve_from_sigma (residue_t A, residue_t x, const residue_t sigma, 
			  const modulus_t m)
{
  residue_t u, v, t, b, z;
  int r;

  mod_init_noset0 (u, m);
  mod_init_noset0 (v, m);
  mod_init_noset0 (t, m);
  mod_init_noset0 (b, m);
  mod_init_noset0 (z, m);

  /* compute b, x */
  mod_add (v, sigma, sigma, m);
  mod_add (v, v, v, m); /* v = 4*sigma */
  mod_mul (u, sigma, sigma, m);
  mod_set1 (b, m);
  mod_add (t, b, b, m);
  mod_add (t, t, t, m);
  mod_add (t, t, b, m); /* t = 5 */
  mod_sub (u, u, t, m); /* u = sigma^2 - 5 */
  mod_mul (t, u, u, m);
  mod_mul (x, t, u, m);
  mod_mul (t, v, v, m);
  mod_mul (z, t, v, m);
  mod_mul (t, x, v, m);
  mod_add (b, t, t, m);
  mod_add (b, b, b, m); /* b = 4 * t */
  mod_add (t, u, u, m);
  mod_add (t, t, u, m); /* t = 3 * u */
  mod_sub (u, v, u, m);
  mod_add (v, t, v, m);
  mod_mul (t, u, u, m);
  mod_mul (u, t, u, m);
  mod_mul (A, u, v, m);
  mod_mul (v, b, z, m);

  r = mod_inv (u, v, m);
  if (r) /* non trivial gcd */
  {
      mod_mul (v, u, b, m);
      mod_mul (x, x, v, m);
      mod_mul (v, u, z, m);
      mod_mul (t, A, v, m);
      mod_set1 (u, m);
      mod_add (u, u, u, m);
      mod_sub (A, t, u, m);
  }

  mod_clear (z, m);
  mod_clear (b, m);
  mod_clear (t, m);
  mod_clear (v, m);
  mod_clear (u, m);

  return r;
}

/* Produces curve in Montgomery parameterization from n value, using
   parameters for a torsion 12 curve as in Montgomery's thesis.
   Return 1 if it worked, 0 if a modular inverse failed */

static int
Monty12_curve_from_k (residue_t A, residue_t x, unsigned long n, 
		      const modulus_t m)
{
  residue_t u, v, u0, v0, a, t2, one;
  
  /* We want a multiple of the point (-2,4) on the curve Y^2=X^3-12*X */
  mod_init (a, m);
  mod_init_noset0 (u, m);
  mod_init_noset0 (v, m);
  mod_init (u0, m);
  mod_init (v0, m);
  mod_init_noset0 (one, m);

  mod_set1 (one, m);
  mod_add (v, one, one, m);
  mod_neg (u, v, m); /* u = -2 */
  mod_add (v, v, v, m); /* v = 4 */
  mod_sub (a, a, v, m);
  mod_sub (a, a, v, m);
  mod_sub (a, a, v, m); /* a = -12 */
  ellW_mul_ui (u, v, n/2, a, m);
  if (n % 2 == 1)
    ellW_add3 (u, v, u, v, u0, v0, a, m);
  /* Now we have a $u$ so that $u^3-12u$ is a square */
  mod_clear (u0, m);
  mod_clear (v0, m);
  /* printf ("Monty12_curve_from_k: u = %lu\n", mod_get_ul (u)); */
  
  mod_init_noset0 (t2, m);
  mod_div2 (v, u, m);
  mod_mul (t2, v, v, m); /* u^2/4 */
  mod_sub (t2, t2, one, m);
  mod_sub (t2, t2, one, m);
  mod_sub (t2, t2, one, m); /* u^2/4 - 3 */
  if (mod_inv (u, u, m) == 0)
  {
    fprintf (stderr, "Monty12_curve_from_k: u = 0\n");
    mod_clear (t2, m);
    mod_clear (v, m);
    mod_clear (u, m);
    mod_clear (a, m);
    return 0;
  }
  mod_mul (t2, t2, u, m); /* t^2 = (u^2/4 - 3)/u = (u^2 - 12)/4u */

  mod_sub (u, t2, one, m);
  mod_add (v, t2, one, m);
  mod_add (v, v, one, m);
  mod_add (v, v, one, m);
  mod_mul (a, u, v, m);
  if (mod_inv (a, a, m) == 0) /* a  = 1/(uv), I want u/v and v/u */
  {
    fprintf (stderr, "Monty12_curve_from_k: (t^2 - 1)(t^2 + 3) = 0\n");
    mod_clear (t2, m);
    mod_clear (v, m);
    mod_clear (u, m);
    mod_clear (a, m);
    return 0;
  }
  mod_mul (u, u, u, m); /* u^2 */
  mod_mul (v, v, v, m); /* v^2 */
  mod_mul (v, v, a, m); /* v^2 * (1/(uv)) = v/u = 1/a*/
  mod_mul (a, a, u, m); /* u^2 * (1/(uv)) = u/v = a*/

  mod_mul (u, a, a, m); /* a^2 */
  mod_add (A, u, one, m);
  mod_add (A, A, one, m); /* a^2 + 2 */
  mod_add (t2, A, A, m);
  mod_add (A, A, t2, m); /* 3*(a^2 + 2) */
  mod_mul (t2, A, a, m);
  mod_set (A, v, m);
  mod_sub (A, A, t2, m); /* 1/a - 3 a (a^2 + 2) */
  mod_div2 (v, v, m); /* v = 1/(2a) */
  mod_mul (t2, v, v, m); /* t2 = 1/(2a)^2 */
  mod_mul (A, A, t2, m);

  mod_add (x, u, u, m);
  mod_add (x, x, u, m); /* 3*a^2 */
  mod_add (x, x, one, m); /* 3*a^2 + 1 */
  mod_div2 (v, v, m); /* v = 1/(4a) */
  mod_mul (x, x, v, m);
  
  mod_clear (one, m);
  mod_clear (t2, m);
  mod_clear (v, m);
  mod_clear (u, m);
  mod_clear (a, m);
  return 1;
}


/* Make a curve of the form y^2 = x^3 + a*x^2 + b with a valid point
   (x, y) from a curve Y^2 = X^3 + A*X^2 + X. The value of b will not
   be computed. 
   x and X may be the same variable. */

static int
curveW_from_Montgomery (residue_t a, residue_t x, residue_t y,
			residue_t X, residue_t A, const modulus_t m)
{
  residue_t g, one;
  int r;

  mod_init_noset0 (g, m);
  mod_init_noset0 (one, m);
  mod_set1 (one, m);

  mod_add (g, X, A, m);
  mod_mul (g, g, X, m);
  mod_add (g, g, one, m);
  mod_mul (g, g, X, m); /* G = X^3 + A*X^2 + X */
  /* printf ("curveW_from_Montgomery: Y^2 = %lu\n", g[0]); */

  /* Now (x,1) is on the curve G*Y^2 = X^3 + A*X^2 + X. */
  r = mod_inv (g, g, m);
  if (r != 0)
  {
      mod_set (y, g, m);       /* y = 1/G */
      mod_div3 (a, A, m);
      mod_add (x, X, a, m);
      mod_mul (x, x, g, m); /* x = (X + A/3)/G */
      mod_mul (a, a, A, m);
      mod_sub (a, one, a, m);
      mod_mul (a, a, g, m);
      mod_mul (a, a, g, m); /* a = (1 - (A^2)/3)/G^2 */
  }
  else
    fprintf (stderr, "curveW_from_Montgomery: r = 0\n");

  mod_clear (one, m);
  mod_clear (g, m);

  return r;
}

static int 
ecm_stage2 (residue_t r, const ellM_point_t P, const stage2_plan_t *plan, 
	    const residue_t b, const modulus_t m)
{
  ellM_point_t Pd, Pt; /* d*P, i*d*P, (i+1)*d*P and a temp */
  residue_t *Pid_x, *Pid_z, *Pj_x, *Pj_z; /* saved i*d*P, i0 <= i < i1, 
              and jP, j in S_1, x and z coordinate stored separately */
  residue_t a, a_bk, t;
  unsigned int i, k, l;
  int bt = 0;
  
  ellM_init (Pt, m);
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
    ellM_point_t ap1_0, ap1_1, ap5_0, ap5_1, P2, P6;
    int i1, i5;
    ellM_init (ap1_0, m);
    ellM_init (ap1_1, m);
    ellM_init (ap5_0, m);
    ellM_init (ap5_1, m);
    ellM_init (P6, m);
    ellM_init (P2, m);

    /* Init ap1_0 = 1P, ap1_1 = 7P, ap5_0 = 5P, ap5_1 = 11P
       and P6 = 6P */
    ellM_set (ap1_0, P, m);         /* ap1_0 = 1*P */
    ellM_double (P2, P, m, b);      /* P2 = 2*P */
    ellM_add (P6, P2, P, P, m);     /* P6 = 3*P (for now) */
    ellM_add (ap5_0, P6, P2, P, m); /* 5*P = 3*P + 2*P */
    ellM_double (P6, P6, m, b);     /* P6 = 6*P = 2*(3*P) */
    ellM_add (ap1_1, P6, P, ap5_0, m); /* 7*P = 6*P + P */
    ellM_add (ap5_1, P6, ap5_0, P, m); /* 11*P = 6*P + 5*P */
    
    ellM_clear (P2, m);

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
	
        ellM_add (Pt, ap1_1, P6, ap1_0, m);
        ellM_set (ap1_0, ap1_1, m);
        ellM_set (ap1_1, Pt, m);
        i1 += 6;
	
        ellM_add (Pt, ap5_1, P6, ap5_0, m);
        ellM_set (ap5_0, ap5_1, m);
        ellM_set (ap5_1, Pt, m);
        i5 += 6;
      }
    
    /* Also compute Pd = d*P while we've got 6*P */
    ellM_mul_ul (Pd, P6, plan->d / 6, m, b);

    ellM_clear (ap1_0, m);
    ellM_clear (ap1_1, m);
    ellM_clear (ap5_0, m);
    ellM_clear (ap5_1, m);
    ellM_clear (P6, m);
    ellM_clear (P2, m);

#if 0
    printf ("Pj = [");
    for (i = 0; i < plan->s1; i++)
      printf ("%s(%lu:%lu)", (i>0) ? ", " : "", 
	      mod_get_ul (Pj_x[i], m), mod_get_ul (Pj_z[i], m));
    printf ("]\n");
#endif
  }
  
  /* Compute idP for i0 <= i < i1 */
  {
    ellM_point_t Pid, Pid1;
    
    ellM_init (Pid, m);
    ellM_init (Pid1, m);
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

    /* Todo: do both Pid and Pid1 with one addition chain */
    ellM_mul_ul (Pid, Pd, i, m, b); /* Pid = i_0 d * P */
    mod_set (Pid_x[k], Pid[0].x, m);
    mod_set (Pid_z[k], Pid[0].z, m);
    k++; i++;
    if (i < plan->i1)
      {
        ellM_mul_ul (Pid1, Pd, i, m, b); /* Pid = (i_0 + 1) d * P */
        mod_set (Pid_x[k], Pid1[0].x, m);
        mod_set (Pid_z[k], Pid1[0].z, m);
        k++; i++;
      }
    while (i < plan->i1)
      {
        ellM_add (Pt, Pid1, Pd, Pid, m);
        ellM_set (Pid, Pid1, m);
        ellM_set (Pid1, Pt, m);
        mod_set (Pid_x[k], Pt[0].x, m);
        mod_set (Pid_z[k], Pt[0].z, m);
        k++; i++;
      }

    ellM_clear (Pid, m);
    ellM_clear (Pid1, m);
#if 0
    printf ("Pid = [");
    for (i = 0; i < plan->i1 - plan->i0; i++)
      printf ("%s(%lu:%lu)", (i>0) ? ", " : "", 
	      mod_get_ul (Pid_x[i], m), mod_get_ul (Pid_z[i], m));
    printf ("]\n");
#endif
  }

  /* Now we've computed all the points we need, so normalize them,
     using Montgomery's batch inversion trick */
  {
    /* If i0 == 0, we skip that point in the inversion */
    /* Let a_0, ..., a_n = Pj_z[0], ..., Pj_z[s1 - 1], 
                           Pid_z[0], ..., Pid_z[i1 - i0 - 1]
       be the list of residues we want to invert. */
    const unsigned int skip_i0 = (plan->i0 == 0) ? 1 : 0;
    const unsigned long n = plan->s1 + plan->i1 - plan->i0 - skip_i0;
    residue_t *invarray;

    invarray = (residue_t *) malloc (n * sizeof (residue_t));
    ASSERT (invarray != NULL);
    for (i = 0; i < n; i++)
      mod_init (invarray[i], m);
    /* Set invarray[k] = prod_{j=0}^{k} a_j */
    mod_set (invarray[0], Pj_z[0], m);
    for (i = 1, k = 1; i < plan->s1; i++, k++)
      mod_mul (invarray[k], invarray[k - 1], Pj_z[i], m);
    for (i = skip_i0; i < plan->i1 - plan->i0; i++, k++)
      mod_mul (invarray[k], invarray[k - 1], Pid_z[i], m);
    
    /* Set a_n := 1 / (a_0, ..., a_n) */
    if (! mod_inv (a, invarray[k - 1], m))
      {
        mod_set (r, invarray[k - 1], m);
	for (i = 0; i < n; i++)
	  mod_clear (invarray[i], m);
	free (invarray);
	invarray = NULL;
	goto clear_and_exit;
      }

    /* Compute 1/a_i and normalize the Pid points */
    mod_set (invarray[k - 1], a, m);
    for (i = plan->i1 - plan->i0; i > skip_i0; i--, k--)
      {
	/* Here, invarray[k - 1] = 1 / (a_0, ..., a_{k-1}) and
	   invarray[k - 2] = (a_0, ..., a_{k-2}) */
	mod_mul (a, invarray[k - 1], invarray[k - 2], m);
	/* a = 1/a_{k-1} = 1/Pid_z[i - 1] */
#ifdef WANT_ASSERT_EXPENSIVE
	mod_mul (a, a, Pid_z[i - 1], m);
	ASSERT (mod_get_ul (a, m) == 1UL);
	mod_mul (a, invarray[k - 1], invarray[k - 2], m);
#endif
	/* Normalize point */
	mod_mul (Pid_x[i - 1], Pid_x[i - 1], a, m);
	/* Set invarray[k - 2] = 1 / (a_0, ..., a_{k-2}) = 
	   1 / (a_0, ..., a_{k-1}) * a_{k-1} */
	mod_mul (invarray[k - 2], invarray[k - 1],  Pid_z[i - 1], m);
      }
    
    /* Compute 1/a_i and normalize the Pj points */
    for (i = plan->s1; i > 1; i--, k--)
      {
	/* Here, invarray[k - 1] = 1 / (a_0, ..., a_{k-1}) and
	   invarray[k - 2] = (a_0, ..., a_{k-2}) */
	mod_mul (a, invarray[k - 1], invarray[k - 2], m);
	/* a = 1/a_{k-1} = 1/Pj_z[i - 1] */
#ifdef WANT_ASSERT_EXPENSIVE
	mod_mul (a, a, Pj_z[i - 1], m);
	ASSERT (mod_get_ul (a, m) == 1UL);
	mod_mul (a, invarray[k - 1], invarray[k - 2], m);
#endif
	/* Normalize point */
	mod_mul (Pj_x[i - 1], Pj_x[i - 1], a, m);
	/* Set invarray[k - 2] = 1 / (a_0, ..., a_{k-2}) = 
	   1 / (a_0, ..., a_{k-1}) * a_{k-1} */
	mod_mul (invarray[k - 2], invarray[k - 1],  Pj_z[i - 1], m);
      }
    /* Do the first Pj point. Here k = 1 and invarray[k] = 1/a_1 */
    ASSERT (k == 1);
#ifdef WANT_ASSERT_EXPENSIVE
    mod_mul (a, invarray[k - 1], Pj_z[i - 1], m);
    ASSERT (mod_get_ul (a, m) == 1UL);
#endif
    /* Normalize point */
    mod_mul (Pj_x[i - 1], Pj_x[i - 1], invarray[k - 1], m);

    for (i = 0; i < n; i++)
      mod_clear (invarray[i], m);
    free (invarray);
    invarray = NULL;
  }
  
  /* Now process all the primes p = id - j, B1 < p <= B2 and multiply
     (id*P)_x - (j*P)_x to the accumulator */
  mod_set1 (a, m);
  mod_init (a_bk, m); /* Backup value of a, in case we get a == 0 */
  mod_set (a_bk, a, m);
  i = 0;
  l = 0;
  while (plan->pairs[l] != NEXT_PASS) 
    {
      while (plan->pairs[l] < NEXT_D && plan->pairs[l] < NEXT_PASS)
	{
	  mod_sub (t, Pid_x[i], Pj_x[plan->pairs[l]], m);
	  mod_mul (a, a, t, m);
	  l++;
	}
      
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
      
      if (plan->pairs[l] == NEXT_D)
	{
	  i++; l++;
	  ASSERT (i < plan->i1 - plan->i0);
	}
    }
  
  mod_set (r, a, m);
  
  /* At end of stage 2 */

 clear_and_exit:
  /* Clear everything */
  for (i = 0; i < plan->s1; i++)
    mod_clear (Pj_z[i], m);
  free (Pj_z);
  Pj_z = NULL;
  for (i = 0; i < plan->i1 - plan->i0; i++)
    mod_clear (Pid_z[i], m);
  free (Pid_z);
  Pid_z = NULL;
  
   for (i = 0; i < plan->s1; i++)
    {
      mod_clear (Pj_x[k], m);
      mod_clear (Pj_z[k], m);
    }
 for (i = 0; i < plan->i1 - plan->i0; i++)
    {
      mod_clear (Pid_x[k], m);
      mod_clear (Pid_z[k], m);
    }
  free (Pj_x);
  free (Pid_x);
  
  ellM_clear (Pt, m);
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
  residue_t u, A, b, r;
  ellM_point_t P, Pt;
  unsigned int i;
  int bt = 0;

  mod_init (u, m);
  mod_init (A, m);
  mod_init (b, m);
  mod_init (r, m);
  ellM_init (P, m);
  ellM_init (Pt, m);

  mod_intset_ul (f, 1UL);

  if (plan->parameterization == BRENT12)
  {
    residue_t s;
    mod_init_noset0 (s, m);
    mod_set_ul (s, plan->sigma, m);
    if (Brent12_curve_from_sigma (A, P->x, s, m) == 0)
      {
	mod_clear (u, m);
	mod_clear (A, m);
	mod_clear (b, m);
	ellM_clear (P, m);
	ellM_clear (Pt, m);
	mod_clear (s, m);
	return 0;
      }
    mod_clear (s, m);
    mod_set1 (P->z, m);
  }
  else if (plan->parameterization == MONTY12)
  {
    if (Monty12_curve_from_k (A, P->x, plan->sigma, m) == 0)
      {
	mod_clear (u, m);
	mod_clear (A, m);
	mod_clear (b, m);
	ellM_clear (P, m);
	ellM_clear (Pt, m);
	return 0;
      }
    mod_set1 (P->z, m);
  }
  else
  {
    fprintf (stderr, "ecm: Unknown parameterization\n");
    abort();
  }

#ifdef TRACE
  printf ("starting point: (%lu::%lu)\n", 
	  mod_get_ul (P->x, m), mod_get_ul (P->z, m));
#endif

  mod_set1 (b, m);
  mod_add (b, b, b, m);
  mod_add (b, b, A, m); /* b = A + 2 */
  mod_div2 (b, b, m);
  mod_div2 (b, b, m);

  /* now start ecm */

  /* Do stage 1 */
  ellM_interpret_bytecode (P, plan->bc, plan->bc_len, m, b);

  /* Add prime 2 in the desired power. If a zero residue for the Z-coordinate 
     is encountered, we go back to previous point and stop */
  ellM_set (Pt, P, m);
  for (i = 0; i < plan->exp2; i++)
    {
      ellM_double (P, P, m, b);
      if (mod_is0 (P[0].z, m))
	{
	  ellM_set (P, Pt, m);
	  bt = 1;
	  break;
	}
      ellM_set (Pt, P, m);
    }

  mod_gcd (f, P[0].z, m); /* FIXME: skip this gcd and let the extgcd
			     in stage 2 init find factors? */
  if (bt == 0 && f[0] == 1UL && plan->B1 < plan->stage2.B2)
    {
      bt = ecm_stage2 (r, P, &(plan->stage2), b, m);
      mod_gcd (f, r, m);
    }
  
  mod_clear (u, m);
  mod_clear (A, m);
  mod_clear (b, m);
  mod_clear (r, m);
  ellM_clear (P, m);
  ellM_clear (Pt, m);

  return bt;
}


/* Determine order of a point P on a curve, both defined by the sigma value
   as in ECM. Looks for i in Hasse interval so that i*P = O, has complexity
   O(sqrt(m)). */

unsigned long
ell_pointorder (const residue_t sigma, const int parameterization, 
		const modulus_t m, const int verbose)
{
  residue_t A, x, a, xi, yi, x1, y1;
  unsigned long min, max, i, order, p;

  mod_init (A, m);

  if (parameterization == BRENT12)
    {
      if (Brent12_curve_from_sigma (A, x, sigma, m) == 0)
	return 0;
    }
  else if (parameterization == MONTY12)
  {
    if (Monty12_curve_from_k (A, x, mod_get_ul (sigma, m), m) == 0)
      return 1;
  }
  else
  {
    fprintf (stderr, "ecm: Unknown parameterization\n");
    abort();
  }
  
  if (verbose >= 2)
    printf ("Curve parameters: A = %ld, x = %ld (mod %ld)\n", 
            mod_get_ul (A, m), mod_get_ul (x, m), mod_getmod_ul (m));

  if (curveW_from_Montgomery (a, x1, y1, x, A, m) == 0)
    return 0UL;

  if (verbose >= 2)
    printf ("Finding order of point (%ld, %ld) on curve "
	    "y^2 = x^3 + %ld * x + b (mod %ld)\n", 
	    mod_get_ul (x1, m), mod_get_ul (y1, m), mod_get_ul (a, m), 
	    mod_getmod_ul (m));

  i = 2 * (unsigned long) sqrt((double) mod_getmod_ul (m));
  min = mod_getmod_ul (m) - i + 1;
  max = mod_getmod_ul (m) + i + 1;
  mod_set (xi, x1, m);
  mod_set (yi, y1, m);
  if (ellW_mul_ui (xi, yi, min, a, m) == 0)
  {
      i = min;
  }
  else
  {
      for (i = min + 1; i <= max; i++)
      {
	  if (!ellW_add3 (xi, yi, xi, yi, x1, y1, a, m))
	      break;
      }
      
      if (i > max)
      {
	  fprintf (stderr, "ell_order: Error, point at infinity not "
		   "reached with i*(x0, z0), i in [%ld, %ld]\n", min, max);
	  return 0UL;
      }

      /* Check that this is the correct order */
      mod_set (xi, x1, m);
      mod_set (yi, y1, m);
      if (ellW_mul_ui (xi, yi, i, a, m) != 0)
      {
	  fprintf (stderr, "ell_order: Error, %ld*(%ld, %ld) (mod %ld) is "
		   "not the point at infinity\n", 
		   i, mod_get_ul (x1, m), mod_get_ul (y1, m), 
		   mod_getmod_ul (m));
	  return 0UL;
      }
  }
  
  /* Ok, now we have some i so that ord(P) | i. Find ord(P).
     We know that ord(P) > 1 since P is not at infinity */

  order = i;
  for (p = 2; p * p <= order; p++)
  {
      mod_set (xi, x1, m);
      mod_set (yi, y1, m);
      while (order % p == 0 && ellW_mul_ui (xi, yi, order / p, a, m) == 0)
	  order /= p;
  }

  return order;
}


/* Count points on curve using the Jacobi symbol. This has complexity O(m). */

unsigned long 
ellM_curveorder_jacobi (residue_t A, residue_t x, modulus_t m)
{
  residue_t t, one;
  unsigned long order, i;
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
  for (i = 1; i < mod_getmod_ul(m); i++)
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
ell_curveorder (const unsigned long sigma_par, int parameterization, 
		const unsigned long m_par)
{
  residue_t sigma, A, X;
  modulus_t m;
  unsigned long order;

  mod_initmod_ul (m, m_par);
  mod_set_ul (sigma, sigma_par, m);

  if (parameterization == BRENT12)
  {
    if (Brent12_curve_from_sigma (A, X, sigma, m) == 0)
      return 0UL;
  }
  else if (parameterization == MONTY12)
  {
    if (Monty12_curve_from_k (A, X, sigma_par, m) == 0)
      return 0UL;
  }
  else
  {
    fprintf (stderr, "ell_curveorder: Unknown parameterization\n");
    abort();
  }
  order = ellM_curveorder_jacobi (A, X, m);

#ifndef NDEBUG
  ASSERT (parameterization != BRENT12 || order == ell_pointorder (sigma, parameterization, m, 0));
#endif

  return order;
}
