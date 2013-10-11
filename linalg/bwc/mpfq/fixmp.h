/*
  Prototypes for fixmp.h
  The functions are implemented in a separate .h file that depends on the architecture.

  NB: For all these functions, overlaps are not allowed, but pointer aliasing is ok.
  (i.e. x == z, or even x==y==z)
*/

#ifndef MPFQ_FIXMP_H_
#define MPFQ_FIXMP_H_

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <assert.h>

/* We need umul_ppmm. If available from gmp, take it, otherwise redefine
 * it.  */
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

/*----------------------- Addition - Subtraction ---------------------*/

/* Each function has a _nc variant, in which the carry/borrow is lost */

/*
FUNCTION: mp_limb_t add_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t *y)
  Adds y to x and put the result in z.
  x,y and z have k limbs.
  The potential carry is returned.
*/
static mp_limb_t add_1(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static mp_limb_t add_2(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static mp_limb_t add_3(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static mp_limb_t add_4(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static mp_limb_t add_5(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static mp_limb_t add_6(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static mp_limb_t add_7(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static mp_limb_t add_8(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static mp_limb_t add_9(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;

/*
FUNCTION: mp_limb_t sub_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t *y)
  Subtracts y to x and put the result in z.
  x,y and z have k limbs.
  The potential borrow is returned
*/
static mp_limb_t sub_1(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static mp_limb_t sub_2(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static mp_limb_t sub_3(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static mp_limb_t sub_4(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static mp_limb_t sub_5(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static mp_limb_t sub_6(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static mp_limb_t sub_7(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static mp_limb_t sub_8(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static mp_limb_t sub_9(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;

/*
FUNCTION: mp_limb_t add_ui_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t y)
  Adds y to x and put the result in z.
  x and z have k limbs.
  The potential carry is returned.
*/
static mp_limb_t add_ui_1(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static mp_limb_t add_ui_2(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static mp_limb_t add_ui_3(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static mp_limb_t add_ui_4(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static mp_limb_t add_ui_5(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static mp_limb_t add_ui_6(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static mp_limb_t add_ui_7(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static mp_limb_t add_ui_8(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static mp_limb_t add_ui_9(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;

/*
FUNCTION: mp_limb_t sub_ui_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t y)
  Subtracts y to x and put the result in z.
  x and z have k limbs.
  The potential borrow is returned
*/
static mp_limb_t sub_ui_1(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static mp_limb_t sub_ui_2(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static mp_limb_t sub_ui_3(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static mp_limb_t sub_ui_4(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static mp_limb_t sub_ui_5(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static mp_limb_t sub_ui_6(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static mp_limb_t sub_ui_7(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static mp_limb_t sub_ui_8(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static mp_limb_t sub_ui_9(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;


/*
FUNCTION: void add_nc_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t *y)
  Adds y to x and put the result in z.
  x,y and z have k limbs.
  The potential carry is lost.
*/
static void add_nc_1(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void add_nc_2(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void add_nc_3(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void add_nc_4(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void add_nc_5(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void add_nc_6(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void add_nc_7(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void add_nc_8(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void add_nc_9(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;

/*
FUNCTION: void sub_nc_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t *y)
  Subtracts y to x and put the result in z.
  x,y and z have k limbs.
  The potential borrow is lost
*/
static void sub_nc_1(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void sub_nc_2(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void sub_nc_3(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void sub_nc_4(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void sub_nc_5(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void sub_nc_6(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void sub_nc_7(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void sub_nc_8(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void sub_nc_9(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;

/*
FUNCTION: void add_ui_nc__k(mp_limb_t *z, mp_limb_t *x, mp_limb_t y)
  Adds y to x and put the result in z.
  x and z have k limbs.
  The potential carry is lost.
*/
static void add_ui_nc_1(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static void add_ui_nc_2(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static void add_ui_nc_3(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static void add_ui_nc_4(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static void add_ui_nc_5(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static void add_ui_nc_6(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static void add_ui_nc_7(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static void add_ui_nc_8(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static void add_ui_nc_9(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;

/*
FUNCTION: void sub_ui_nc_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t y)
  Subtracts y to x and put the result in z.
  x and z have k limbs.
  The potential borrow is lost
*/
static void sub_ui_nc_1(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static void sub_ui_nc_2(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static void sub_ui_nc_3(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static void sub_ui_nc_4(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static void sub_ui_nc_5(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static void sub_ui_nc_6(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static void sub_ui_nc_7(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static void sub_ui_nc_8(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static void sub_ui_nc_9(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;

/*--------------------------- Multiplication ---- ---------------------*/

/*
FUNCTION: void addmul1_nc_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t c)
  Multiplies the limb c to x and adds the result in z.
  x has k limbs,
  z has k+1 limbs.
  The potential carry is lost (better have z[k]=0 !!!)
*/
static void addmul1_nc_1(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void addmul1_nc_2(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void addmul1_nc_3(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void addmul1_nc_4(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void addmul1_nc_5(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void addmul1_nc_6(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void addmul1_nc_7(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void addmul1_nc_8(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void addmul1_nc_9(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;


/*
FUNCTION: mp_limb_t addmul1_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t c)
  Multiplies the limb c to x and adds the result in z.
  x has k limbs,
  z has k+1 limbs.
  The potential carry is returned
*/
static mp_limb_t addmul1_1(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static mp_limb_t addmul1_2(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static mp_limb_t addmul1_3(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static mp_limb_t addmul1_4(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static mp_limb_t addmul1_5(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static mp_limb_t addmul1_6(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static mp_limb_t addmul1_7(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static mp_limb_t addmul1_8(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static mp_limb_t addmul1_9(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;

/*
FUNCTION: void mul1_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t c)
  Multiplies the limb c to x and puts the result in z.
  x has k limbs,
  z has k+1 limbs.
*/
static void mul1_1(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void mul1_2(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void mul1_3(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void mul1_4(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void mul1_5(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void mul1_6(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void mul1_7(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void mul1_8(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void mul1_9(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;

/*
FUNCTION: void mul_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t *y)
  Multiply x and y and put the result in z.
  x and y  have k limbs.
  z has 2*k limbs.
*/
static void mul_1(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void mul_2(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void mul_3(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void mul_4(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void mul_5(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void mul_6(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void mul_7(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void mul_8(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void mul_9(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;

/*
FUNCTION: void sqr_k(mp_limb_t *z, mp_limb_t *x)
  Square x and put the result in z.
  x has k limbs.
  z has 2*k limbs.
*/
static void sqr_1(mp_limb_t *z, const mp_limb_t *x) MAYBE_UNUSED;
static void sqr_2(mp_limb_t *z, const mp_limb_t *x) MAYBE_UNUSED;
static void sqr_3(mp_limb_t *z, const mp_limb_t *x) MAYBE_UNUSED;
static void sqr_4(mp_limb_t *z, const mp_limb_t *x) MAYBE_UNUSED;
static void sqr_5(mp_limb_t *z, const mp_limb_t *x) MAYBE_UNUSED;
static void sqr_6(mp_limb_t *z, const mp_limb_t *x) MAYBE_UNUSED;
static void sqr_7(mp_limb_t *z, const mp_limb_t *x) MAYBE_UNUSED;
static void sqr_8(mp_limb_t *z, const mp_limb_t *x) MAYBE_UNUSED;
static void sqr_9(mp_limb_t *z, const mp_limb_t *x) MAYBE_UNUSED;

/*
FUNCTION: void shortmul_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t *y)
  Multiply x and y and put the truncated result in z.
  x and y  have k limbs.
  z has k limbs.
*/
static void shortmul_1(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void shortmul_2(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void shortmul_3(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void shortmul_4(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void shortmul_5(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void shortmul_6(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void shortmul_7(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void shortmul_8(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void shortmul_9(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;

/*---------------------- Other ---------------------------------------*/

/*
FUNCTION: void mod_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t *p)
  Reduce x modulo p and put the result in z.
  p has k limbs, and x has 2*k limbs. Furthermore, p[k-1]!=0.
  z has k limbs.
*/
static void mod_1(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mod_2(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mod_3(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mod_4(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mod_5(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mod_6(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mod_7(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mod_8(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mod_9(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;

/*
FUNCTION: int invmod_k(mp_limb_t *z, mp_limb_t *x, mp_limb_t *p)
  Put in z the inverse of x modulo p if it exists (and then return 1)
  If x is 0 modulo p, abort.
  If x is non invertible, put a factor of p in z and return 0.
*/
static int invmod_1(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static int invmod_2(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static int invmod_3(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static int invmod_4(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static int invmod_5(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static int invmod_6(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static int invmod_7(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static int invmod_8(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static int invmod_9(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;

/*
FUNCTION: int cmp_k(mp_limb_t *x, mp_limb_t *y)
  Compares x and y and returns the result:
    0   if x==y
    1   if x>y
    -1  if x<y
  x and y have k limbs. 
*/
static int cmp_1(const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static int cmp_2(const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static int cmp_3(const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static int cmp_4(const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static int cmp_5(const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static int cmp_6(const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static int cmp_7(const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static int cmp_8(const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static int cmp_9(const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;

/*
FUNCTION: int cmp_ui_k(mp_limb_t *x, mp_limb_t y)
  Compares x and y and returns the result:
    0   if x==y
    1   if x>y
    -1  if x<y
  x has k limbs. 
*/
static int cmp_ui_1(const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static int cmp_ui_2(const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static int cmp_ui_3(const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static int cmp_ui_4(const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static int cmp_ui_5(const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static int cmp_ui_6(const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static int cmp_ui_7(const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static int cmp_ui_8(const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;
static int cmp_ui_9(const mp_limb_t *x, const mp_limb_t y) MAYBE_UNUSED;


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
static void lshift_1(mp_limb_t *z, int cnt) MAYBE_UNUSED;
static void lshift_2(mp_limb_t *z, int cnt) MAYBE_UNUSED;
static void lshift_3(mp_limb_t *z, int cnt) MAYBE_UNUSED;
static void lshift_4(mp_limb_t *z, int cnt) MAYBE_UNUSED;
static void lshift_5(mp_limb_t *z, int cnt) MAYBE_UNUSED;
static void lshift_6(mp_limb_t *z, int cnt) MAYBE_UNUSED;
static void lshift_7(mp_limb_t *z, int cnt) MAYBE_UNUSED;
static void lshift_8(mp_limb_t *z, int cnt) MAYBE_UNUSED;
static void lshift_9(mp_limb_t *z, int cnt) MAYBE_UNUSED;

static void rshift_1(mp_limb_t *z, int cnt) MAYBE_UNUSED;
static void rshift_2(mp_limb_t *z, int cnt) MAYBE_UNUSED;
static void rshift_3(mp_limb_t *z, int cnt) MAYBE_UNUSED;
static void rshift_4(mp_limb_t *z, int cnt) MAYBE_UNUSED;
static void rshift_5(mp_limb_t *z, int cnt) MAYBE_UNUSED;
static void rshift_6(mp_limb_t *z, int cnt) MAYBE_UNUSED;
static void rshift_7(mp_limb_t *z, int cnt) MAYBE_UNUSED;
static void rshift_8(mp_limb_t *z, int cnt) MAYBE_UNUSED;
static void rshift_9(mp_limb_t *z, int cnt) MAYBE_UNUSED;

static void long_lshift_1(mp_limb_t *z, int offset, int cnt) MAYBE_UNUSED;
static void long_lshift_2(mp_limb_t *z, int offset, int cnt) MAYBE_UNUSED;
static void long_lshift_3(mp_limb_t *z, int offset, int cnt) MAYBE_UNUSED;
static void long_lshift_4(mp_limb_t *z, int offset, int cnt) MAYBE_UNUSED;
static void long_lshift_5(mp_limb_t *z, int offset, int cnt) MAYBE_UNUSED;
static void long_lshift_6(mp_limb_t *z, int offset, int cnt) MAYBE_UNUSED;
static void long_lshift_7(mp_limb_t *z, int offset, int cnt) MAYBE_UNUSED;
static void long_lshift_8(mp_limb_t *z, int offset, int cnt) MAYBE_UNUSED;
static void long_lshift_9(mp_limb_t *z, int offset, int cnt) MAYBE_UNUSED;

static void long_rshift_1(mp_limb_t *z, int offset, int cnt) MAYBE_UNUSED;
static void long_rshift_2(mp_limb_t *z, int offset, int cnt) MAYBE_UNUSED;
static void long_rshift_3(mp_limb_t *z, int offset, int cnt) MAYBE_UNUSED;
static void long_rshift_4(mp_limb_t *z, int offset, int cnt) MAYBE_UNUSED;
static void long_rshift_5(mp_limb_t *z, int offset, int cnt) MAYBE_UNUSED;
static void long_rshift_6(mp_limb_t *z, int offset, int cnt) MAYBE_UNUSED;
static void long_rshift_7(mp_limb_t *z, int offset, int cnt) MAYBE_UNUSED;
static void long_rshift_8(mp_limb_t *z, int offset, int cnt) MAYBE_UNUSED;
static void long_rshift_9(mp_limb_t *z, int offset, int cnt) MAYBE_UNUSED;

/*--------------------- Montgomery rep ------------------------------*/

/*
FUNCTION: void redc_k(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p)
  REDC computation: divide x by R modulo p, where R is 2^(k*GMP_LIMB_BITS).
  p has k limbs, and x has 2*k limbs. Furthermore, p[k-1]!=0.
  z has k limbs.
  mip has k limbs and is -1/p mod R.
  Main Redc condition: 0 <= x <= R*p.
*/
static void redc_1(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;
static void redc_2(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;
static void redc_3(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;
static void redc_4(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;
static void redc_5(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;
static void redc_6(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;
static void redc_7(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;
static void redc_8(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;
static void redc_9(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;

/* 
FUNCTION: void mgy_encode_k(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p)
  Put in z the encoded representation of x mod p:
    z := x*R mod p, where R is 2^(k*GMP_LIMB_BITS)
  z, x and p have k limbs.
*/
static void mgy_encode_1(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_encode_2(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_encode_3(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_encode_4(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_encode_5(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_encode_6(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_encode_7(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_encode_8(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_encode_9(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;

/* 
FUNCTION: void mgy_decode_k(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *invR, const mp_limb_t *p)
  Put in z the decoded representation of x mod p:
    z := x/R mod p, where R is 2^(k*GMP_LIMB_BITS)
  z, x, invR and p have k limbs.
  invR should contain 1/R mod p.
*/
static void mgy_decode_1(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *invR, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_decode_2(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *invR, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_decode_3(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *invR, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_decode_4(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *invR, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_decode_5(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *invR, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_decode_6(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *invR, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_decode_7(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *invR, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_decode_8(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *invR, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_decode_9(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *invR, const mp_limb_t *p) MAYBE_UNUSED;



static mp_limb_t addmul1_1hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static mp_limb_t addmul1_2hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static mp_limb_t addmul1_3hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static mp_limb_t addmul1_4hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static mp_limb_t addmul1_5hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static mp_limb_t addmul1_6hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static mp_limb_t addmul1_7hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static mp_limb_t addmul1_8hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static mp_limb_t addmul1_9hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;


static mp_limb_t addmul1_smallz_1hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static mp_limb_t addmul1_smallz_2hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static mp_limb_t addmul1_smallz_3hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static mp_limb_t addmul1_smallz_4hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static mp_limb_t addmul1_smallz_5hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static mp_limb_t addmul1_smallz_6hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static mp_limb_t addmul1_smallz_7hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static mp_limb_t addmul1_smallz_8hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static mp_limb_t addmul1_smallz_9hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;

static void mul1_1hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void mul1_2hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void mul1_3hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void mul1_4hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void mul1_5hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void mul1_6hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void mul1_7hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void mul1_8hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void mul1_9hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;

static void mul_1hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void mul_2hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void mul_3hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void mul_4hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void mul_5hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void mul_6hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void mul_7hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void mul_8hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;
static void mul_9hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y) MAYBE_UNUSED;

static void sqr_1hw(mp_limb_t *z, const mp_limb_t *x) MAYBE_UNUSED;
static void sqr_2hw(mp_limb_t *z, const mp_limb_t *x) MAYBE_UNUSED;
static void sqr_3hw(mp_limb_t *z, const mp_limb_t *x) MAYBE_UNUSED;
static void sqr_4hw(mp_limb_t *z, const mp_limb_t *x) MAYBE_UNUSED;
static void sqr_5hw(mp_limb_t *z, const mp_limb_t *x) MAYBE_UNUSED;
static void sqr_6hw(mp_limb_t *z, const mp_limb_t *x) MAYBE_UNUSED;
static void sqr_7hw(mp_limb_t *z, const mp_limb_t *x) MAYBE_UNUSED;
static void sqr_8hw(mp_limb_t *z, const mp_limb_t *x) MAYBE_UNUSED;
static void sqr_9hw(mp_limb_t *z, const mp_limb_t *x) MAYBE_UNUSED;

static void mod_1hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mod_2hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mod_3hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mod_4hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mod_5hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mod_6hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mod_7hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mod_8hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mod_9hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;

static void redc_1hw(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;
static void redc_2hw(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;
static void redc_3hw(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;
static void redc_4hw(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;
static void redc_5hw(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;
static void redc_6hw(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;
static void redc_7hw(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;
static void redc_8hw(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;
static void redc_9hw(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;

static void mgy_encode_1hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_encode_2hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_encode_3hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_encode_4hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_encode_5hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_encode_6hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_encode_7hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_encode_8hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_encode_9hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *p) MAYBE_UNUSED;

static void mgy_decode_1hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *invR, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_decode_2hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *invR, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_decode_3hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *invR, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_decode_4hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *invR, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_decode_5hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *invR, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_decode_6hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *invR, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_decode_7hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *invR, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_decode_8hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *invR, const mp_limb_t *p) MAYBE_UNUSED;
static void mgy_decode_9hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *invR, const mp_limb_t *p) MAYBE_UNUSED;

static void addmul1_nc_1hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void addmul1_nc_2hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void addmul1_nc_3hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void addmul1_nc_4hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void addmul1_nc_5hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void addmul1_nc_6hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void addmul1_nc_7hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void addmul1_nc_8hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;
static void addmul1_nc_9hw(mp_limb_t *z, const mp_limb_t *x, const mp_limb_t c) MAYBE_UNUSED;

static void redc_ur_1(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;
static void redc_ur_2(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;
static void redc_ur_3(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;
static void redc_ur_4(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;
static void redc_ur_5(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;
static void redc_ur_6(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;
static void redc_ur_7(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;
static void redc_ur_8(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;
static void redc_ur_9(mp_limb_t *z, mp_limb_t *x, const mp_limb_t *mip, const mp_limb_t *p) MAYBE_UNUSED;

/* Load architecture-specific definitions */
#ifdef __x86_64__
# include "fixmp_x86_64.h"
#endif

#if (defined (__i386__) || defined (__i486__))
# include "fixmp_x86_32.h"
#endif

/* Complete missing definitions with generic functions based on longlong.h */

#include "fixmp_longlong.h"

#endif /* MPFQ_FIXMP_H_ */
