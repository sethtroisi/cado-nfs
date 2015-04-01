/*
 * Prototypes for fixmp.h
 *
 * The functions are implemented in a separate .h file that depends on
 * the architecture.
 *
 * NB: For all these functions, overlaps are not allowed, but pointer
 * aliasing is ok.  (i.e. x == z, or even x==y==z)
 *
 * Exceptions: mul, shortmul, and sqr do not support any input pointer to
 * alias the output pointer.
 */

#ifndef MPFQ_FIXMP_H_
#define MPFQ_FIXMP_H_

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <assert.h>

/* We need umul_ppmm. If available from gmp, take it, otherwise redefine it. */
/* #undef HAVE_GMP_LONGLONG_AND_IMPL */
#ifdef HAVE_GMP_LONGLONG_AND_IMPL
#include "gmp-impl.h"
#include "longlong.h"
#define mpfq_umul_ppmm(w1, w0, u, v) umul_ppmm(w1, w0, u, v)
#else
/* This function is adapted from gmp/longlong.h, copyright FSF. */
static inline void 
mpfq_umul_ppmm_func(mp_limb_t *pw1,
    mp_limb_t *pw0,
    mp_limb_t u,
    mp_limb_t v)
{                   
    mp_limb_t x0, x1, x2, x3;
    mp_limb_t ul, vl, uh, vh;                                   
    mp_limb_t mask_high = (1UL<<(GMP_NUMB_BITS/2))-1;

    ul = u & mask_high; 
    uh = u >> (GMP_NUMB_BITS/2);                                       
    vl = v & mask_high;                                        
    vh = v >> (GMP_NUMB_BITS/2);                                       
                                                                      
    x0 = ul * vl;                                      
    x1 = ul * vh;                                      
    x2 = uh * vl;                                      
    x3 = uh * vh;                                      
                                                                      
    x1 += x0 >> (GMP_NUMB_BITS/2);/* this can't give carry */          
    x1 += x2;               /* but this indeed can */             
    if (x1 < x2)            /* did we get it? */                  
      x3 += 1UL << (GMP_NUMB_BITS/2);  /* yes, add it in the proper pos. */  
                                                                      
    (*pw1) = x3 + (x1 >> (GMP_NUMB_BITS/2));                               
    (*pw0) = (x1 << GMP_NUMB_BITS/2) + (x0 & mask_high);
}
#define mpfq_umul_ppmm(w1, w0, u, v) mpfq_umul_ppmm_func(&(w1), &(w0), u, v)
#endif

#include "mpfq.h"

/* Load architecture-specific definitions */
#if defined(__x86_64__) && defined(__GNUC__)
# include "mpfq_fixmp_x86_64.h"
#endif

#if defined (__i386__) && defined(__GNUC__)
# include "mpfq_fixmp_i386.h"
#endif

/* Complete missing definitions with generic functions based on longlong.h */

#include "mpfq_fixmp_longlong.h"


/* The rest of this file is only documentation, for convenience.
 * Authoritative documentation and prototype information is to be found
 * in fixmp/gen_fixmp.pl
 */

/*----------------------- Addition - Subtraction ---------------------*/

/* Each function has a _nc variant, in which the carry/borrow is lost */

/*
FUNCTION: mp_limb_t add_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t *y)
  Adds y to x and put the result in z.
  x,y and z have k limbs.
  The potential carry is returned.
*/

/*
FUNCTION: mp_limb_t sub_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t *y)
  Subtracts y to x and put the result in z.
  x,y and z have k limbs.
  The potential borrow is returned
*/

/*
FUNCTION: mp_limb_t add_ui_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t y)
  Adds y to x and put the result in z.
  x and z have k limbs.
  The potential carry is returned.
*/

/*
FUNCTION: mp_limb_t sub_ui_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t y)
  Subtracts y to x and put the result in z.
  x and z have k limbs.
  The potential borrow is returned
*/


/*
FUNCTION: void add_nc_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t *y)
  Adds y to x and put the result in z.
  x,y and z have k limbs.
  The potential carry is lost.
*/

/*
FUNCTION: void sub_nc_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t *y)
  Subtracts y to x and put the result in z.
  x,y and z have k limbs.
  The potential borrow is lost
*/

/*
FUNCTION: void add_ui_nc__k(mp_limb_t *z, mp_limb_t *x, mp_limb_t y)
  Adds y to x and put the result in z.
  x and z have k limbs.
  The potential carry is lost.
*/

/*
FUNCTION: void sub_ui_nc_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t y)
  Subtracts y to x and put the result in z.
  x and z have k limbs.
  The potential borrow is lost
*/

/*--------------------------- Multiplication ---- ---------------------*/

/*
FUNCTION: void addmul1_nc_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t c)
  Multiplies the limb c to x and adds the result in z.
  x has k limbs,
  z has k+1 limbs.
  The potential carry is lost (better have z[k]=0 !!!)
*/


/*
FUNCTION: mp_limb_t addmul1_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t c)
  Multiplies the limb c to x and adds the result in z.
  x has k limbs,
  z has k+1 limbs.
  The potential carry is returned
*/

/*
FUNCTION: void mul1_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t c)
  Multiplies the limb c to x and puts the result in z.
  x has k limbs,
  z has k+1 limbs.
*/

/*
FUNCTION: void mul_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t *y)
  Multiply x and y and put the result in z.
  x and y  have k limbs.
  z has 2*k limbs.
*/

/*
FUNCTION: void sqr_k(mp_limb_t *z, mp_limb_t *x)
  Square x and put the result in z.
  x has k limbs.
  z has 2*k limbs.
*/

/*
FUNCTION: void shortmul_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t *y)
  Multiply x and y and put the truncated result in z.
  x and y  have k limbs.
  z has k limbs.
*/

/*---------------------- Other ---------------------------------------*/

/*
FUNCTION: void mod_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t *p)
  Reduce x modulo p and put the result in z.
  p has k limbs, and x has 2*k limbs. Furthermore, p[k-1]!=0.
  z has k limbs.
*/

/*
FUNCTION: int invmod_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t *p)
  Put in z the inverse of x modulo p if it exists (and then return 1)
  If x is 0 modulo p, abort.
  If x is non invertible, put a factor of p in z and return 0.
*/

/*
FUNCTION: int cmp_k(mp_limb_t *x, mp_limb_t *y)
  Compares x and y and returns the result:
    0   if x==y
    1   if x>y
    -1  if x<y
  x and y have k limbs. 
*/

/*
FUNCTION: int cmp_ui_k(mp_limb_t *x, mp_limb_t y)
  Compares x and y and returns the result:
    0   if x==y
    1   if x>y
    -1  if x<y
  x has k limbs. 
*/


/*
FUNCTIONS:
  void lshift_k(mp_limb_t *z, int cnt)
  void rshift_k(mp_limb_t *z, int cnt)
  void long_lshift_k(mp_limb_t *z, int offset, int cnt)
  void long_rshift_k(mp_limb_t *z, int offset, int cnt)
  z has k limbs. cnt is in [0, GMP_NUMB_BITS[.
  shifts z by (offset*GMP_NUMB_BITS)+cnt bits, in place.
  The result is put in the k limbs of z, and the bits that are shifted out are lost.
*/




/*--------------------- Montgomery rep ------------------------------*/

/*
FUNCTION: void redc_k(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p)
  REDC computation: divide x by R modulo p, where R is 2^(k*GMP_LIMB_BITS).
  p has k limbs, and x has 2*k limbs. Furthermore, p[k-1]!=0.
  z has k limbs.
  mip has k limbs and is -1/p mod R.
  Main Redc condition: 0 <= x <= R*p.
*/

/* 
FUNCTION: void mgy_encode_k(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p)
  Put in z the encoded representation of x mod p:
    z := x*R mod p, where R is 2^(k*GMP_LIMB_BITS)
  z, x and p have k limbs.
*/

/* 
FUNCTION: void mgy_decode_k(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *invR, const mp_limb_t *p)
  Put in z the decoded representation of x mod p:
    z := x/R mod p, where R is 2^(k*GMP_LIMB_BITS)
  z, x, invR and p have k limbs.
  invR should contain 1/R mod p.
*/

#endif /* MPFQ_FIXMP_H_ */
