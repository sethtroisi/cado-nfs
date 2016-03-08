// Always include "fppol.h" first so that everything gets loaded in the correct
// order.
#include "fppol.h"

#ifndef __FP_H__
#define __FP_H__



/* Base field elements.
 *****************************************************************************/

// Set to zero.
static inline
void fp_set_zero(fp_ptr r)
{ for (unsigned k = 0; k < __FP_BITS; ++k) r[k] = 0; }


// Set to one.
static inline
void fp_set_one(fp_ptr r)
{ __FP_ONE(8, r); }


// Set to another element.
static inline
void fp_set(fp_ptr r, fp_srcptr p)
{ for (unsigned k = 0; k < __FP_BITS; ++k) r[k] = p[k]; }


// Swap two base field elements.
static inline
void fp_swap(fp_ptr x, fp_ptr y)
{ fp_t t; fp_set(t, x); fp_set(x, y); fp_set(y, t); }


// Set from an integer.
static inline
void fp_set_z(fp_ptr r, int x)
{ x %= __FP_CHAR; if (x < 0) x += __FP_CHAR; __FP_SET_Z(8, r, x); }


// Test if zero.
static inline
int fp_is_zero(fp_srcptr p)
{ uint8_t t = p[0];
  for (unsigned k = 1; k < __FP_BITS; ++k) t |= p[k];
  return !t; }


// Test if one.
static inline
int fp_is_one(fp_srcptr p)
{ uint8_t t = 0;
  for (unsigned k = 1; k < __FP_BITS; ++k) t |= p[k];
  return p[0] == 1 && !t; }


// Opposite.
static inline
void fp_opp(fp_ptr r, fp_srcptr p)
{ __FP_OPP(8, r, p); }


// Addition.
static inline
void fp_add(fp_ptr r, fp_srcptr p, fp_srcptr q)
{ __FP_ADD(8, r, p, q); }


// Subtraction.
static inline
void fp_sub(fp_ptr r, fp_srcptr p, fp_srcptr q)
{ __FP_SUB(8, r, p, q); }


// Inverse.
static inline
void fp_inv(fp_ptr r, fp_srcptr p)
{ __FP_SINV(8, r, p); }


// Multiplication.
static inline
void fp_mul(fp_ptr r, fp_srcptr p, fp_srcptr q)
{ __FP_SMUL(8, r, p, q); }


// Division.
static inline
void fp_div(fp_ptr r, fp_srcptr p, fp_srcptr q)
{ __FP_SDIV(8, r, p, q); }

#endif
