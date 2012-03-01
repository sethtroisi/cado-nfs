#ifndef __FPPOL_H__
#define __FPPOL_H__

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif



/* Load field-specific data.
 *****************************************************************************/

// Include the f<q>pol.h file corresponding to the specified base field
// GF(<q>).
#if   defined(USE_F2)
# include "f2pol.h"
#elif defined(USE_F3)
# include "f3pol.h"
#endif

// Number of elements in the field.
#define FP_SIZE __FP_SIZE

// Field characteristic.
#define FP_CHAR __FP_CHAR



/* Base field elements.
 *****************************************************************************/
  
// Type and pointer shorthands.
// - fp_t:      type of an element of the base field.
// - fp_ptr:    read/write pointer (internal type)
// - fp_srcptr: read-only  pointer (internal type)
typedef       uint8_t  fp_t[__FP_BITS];
typedef       uint8_t *fp_ptr;
typedef const uint8_t *fp_srcptr;



/* Fixed-size polynomials.
 *****************************************************************************/

// Fixed-size polynomial of at most <sz> terms.
// - fppol<sz>_t:      type of <sz>-term polynomials.
// - fppol<sz>_ptr:    read/write pointer (internal type)
// - fppol<sz>_srcptr: read-only  pointer (internal type)
#define __DECL_FPPOLxx_T(sz)                            \
  typedef       uint##sz##_t  fppol##sz##_t[__FP_BITS]; \
  typedef       uint##sz##_t *fppol##sz##_ptr;          \
  typedef const uint##sz##_t *fppol##sz##_srcptr;

__DECL_FPPOLxx_T(16)
__DECL_FPPOLxx_T(32)
__DECL_FPPOLxx_T(64)

#undef __DECL_FPPOLxx_T



/* Multiprecision polynomials.
 *****************************************************************************/

// By convention, deg(0) = -1.
// All the exported functions must ensure that they return fppol such
// that "deg" contains the true degree of the polynomial represented
// by "limbs".
typedef struct {
  int        deg;
  unsigned   alloc;
  fppol64_t *limbs;
} __fppol_struct;

// Type and pointer shorthands.
// - fppol_t:      type of multiprecision polynomials.
// - fppol_ptr:    read/write pointer (internal type)
// - fppol_srcptr: read-only  pointer (internal type)
typedef       __fppol_struct  fppol_t[1];
typedef       __fppol_struct *fppol_ptr;
typedef const __fppol_struct *fppol_srcptr;

// Initialize polynomial.
void fppol_init(fppol_ptr r);

// Initialize a NULL-terminated list of polynomials.
void fppol_inits(fppol_ptr r, ...);

// Initialize polynomial with space for n terms.
void fppol_init2(fppol_ptr r, unsigned n);

// Free polynomial.
void fppol_clear(fppol_ptr r);

// Free a NULL-terminated list of polynomials.
void fppol_clears(fppol_ptr r, ...);

// Reallocate polynomial so as to free unused limbs.
void fppol_trim(fppol_ptr r);



/* Include the other fppol_* headers.
 *****************************************************************************/

#include "fppol_arith.h"
#include "fppol_mul.h"
#include "fppol_div.h"
#include "fppol_mod.h"
#include "fppol_gcd.h"
#include "fppol_io.h"

#ifdef __cplusplus
}
#endif
#endif  /* __FPPOL_H__ */
