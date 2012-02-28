// Always include "fppol.h" first so that everything gets loaded in the correct
// order.
#include "fppol.h"

#ifndef __FPPOL_MUL_H__
#define __FPPOL_MUL_H__



/* Macros for internal use.
 *****************************************************************************/

// Operations for addmul and submul variants.
#define __add_OP(r, p) , add(r, r, p)
#define __sub_OP(r, p) , sub(r, r, p)

// Size of multiplication.
#define __MUL_16x16_SIZE ,  32
#define __MUL_32x16_SIZE ,  64
#define __MUL_32x32_SIZE ,  64
#define __MUL_SIZE(sp, sq) SWITCH(MUL_##sp##x##sq, SIZE, 128)



/* Fixed-size polynomials.
 *****************************************************************************/

// <sz> <- <sz> x <sz> <op>-multiplication, where <op> = [empty], add, or sub.
// Generic prototype:
//   void fppol<sz>_<op>mul(fppol<sz>_ptr    r,
//                          fppol<sz>_srcptr p,
//                          fppol<sz>_srcptr q);
#define __DEF_FPPOLxx_opMUL(sz, op)                \
  static inline                                    \
  void fppol##sz##_##op##mul(fppol##sz##_ptr    r, \
                             fppol##sz##_srcptr p, \
                             fppol##sz##_srcptr q) \
  { __FP_MUL_##sz##_##sz##x##sz(r, p, q, op); }


// <sz> <- <sz> x <sq> <op>-multiplication, where <op> = [empty], add, or sub.
// Generic prototype:
//   void fppol<sz>_<op>mul_<sq>(fppol<sz>_ptr    r,
//                               fppol<sz>_srcptr p,
//                               fppol<sq>_srcptr q);
#define __DEF_FPPOLxx_opMUL_yy(sz, sq, op)              \
  static inline                                         \
  void fppol##sz##_##op##mul_##sq(fppol##sz##_ptr    r, \
                                  fppol##sz##_srcptr p, \
                                  fppol##sq##_srcptr q) \
  { __FP_MUL_##sz##_##sz##x##sq(r, p, q, op); }


// <sz> <- <sp> x <sq> <op>-multiplication, where <op> = [empty], add, or sub,
// and <sp> >= <sq>.
// Generic prototype:
//   void fppol<sz>_<op>mul_<sp>x<sq>(fppol<sz>_ptr    r,
//                                    fppol<sp>_srcptr p,
//                                    fppol<sq>_srcptr q);
#define __DEF_FPPOLxx_opMUL_yyxzz(sz, sp, sq, op)              \
  static inline                                                \
  void fppol##sz##_##op##mul_##sp##x##sq(fppol##sz##_ptr    r, \
                                         fppol##sp##_srcptr p, \
                                         fppol##sq##_srcptr q) \
  { __FP_MUL_##sz##_##sp##x##sq(r, p, q, op); }


// <sz> <- <sp> x <sq> <op>-multiplication, where <op> = [empty], add, or sub,
// and <sp> < <sq>.
// Generic prototype:
//   void fppol<sz>_<op>mul_<sp>x<sq>(fppol<sz>_ptr    r,
//                                    fppol<sp>_srcptr p,
//                                    fppol<sq>_srcptr q);
#define __DEF_FPPOLxx_opMUL_zzxyy(sz, sp, sq, op)              \
  static inline                                                \
  void fppol##sz##_##op##mul_##sp##x##sq(fppol##sz##_ptr    r, \
                                         fppol##sp##_srcptr p, \
                                         fppol##sq##_srcptr q) \
  { fppol##sz##_##op##mul_##sq##x##sp(r, q, p); }


// All definitions bundled up into a single macro.
#define __DEF_FPPOLxx_opMUL_ALL(op)               \
        __DEF_FPPOLxx_opMUL      (16,         op) \
        __DEF_FPPOLxx_opMUL      (32,         op) \
        __DEF_FPPOLxx_opMUL      (64,         op) \
        __DEF_FPPOLxx_opMUL_yy   (16, 16,     op) \
        __DEF_FPPOLxx_opMUL_yy   (32, 16,     op) \
        __DEF_FPPOLxx_opMUL_yy   (32, 32,     op) \
        __DEF_FPPOLxx_opMUL_yy   (64, 16,     op) \
        __DEF_FPPOLxx_opMUL_yy   (64, 32,     op) \
        __DEF_FPPOLxx_opMUL_yy   (64, 64,     op) \
        __DEF_FPPOLxx_opMUL_yyxzz(16, 16, 16, op) \
        __DEF_FPPOLxx_opMUL_yyxzz(32, 16, 16, op) \
        __DEF_FPPOLxx_opMUL_yyxzz(32, 32, 16, op) \
        __DEF_FPPOLxx_opMUL_yyxzz(32, 32, 32, op) \
        __DEF_FPPOLxx_opMUL_yyxzz(64, 32, 16, op) \
        __DEF_FPPOLxx_opMUL_yyxzz(64, 32, 32, op) \
        __DEF_FPPOLxx_opMUL_yyxzz(64, 64, 16, op) \
        __DEF_FPPOLxx_opMUL_yyxzz(64, 64, 32, op) \
        __DEF_FPPOLxx_opMUL_yyxzz(64, 64, 64, op) \
        __DEF_FPPOLxx_opMUL_zzxyy(32, 16, 32, op) \
        __DEF_FPPOLxx_opMUL_zzxyy(64, 16, 32, op) \
        __DEF_FPPOLxx_opMUL_zzxyy(64, 16, 64, op) \
        __DEF_FPPOLxx_opMUL_zzxyy(64, 32, 64, op)

__DEF_FPPOLxx_opMUL_ALL()
__DEF_FPPOLxx_opMUL_ALL(add)
__DEF_FPPOLxx_opMUL_ALL(sub)

#undef __DEF_FPPOLxx_opMUL
#undef __DEF_FPPOLxx_opMUL_yy
#undef __DEF_FPPOLxx_opMUL_yyxzz
#undef __DEF_FPPOLxx_opMUL_zzxyy
#undef __DEF_FPPOLxx_opMUL_ALL



/* Multiprecision polynomials.
 *****************************************************************************/

// <op>-multiplication, where <op> = [empty], add, or sub.
// Generic prototype:
//   void fppol_<op>mul(fppol_ptr r, fppol_srcptr p, fppol_srcptr q);
#define __DECL_FPPOL_opMUL(op)                                       \
  void fppol_##op##mul(fppol_ptr r, fppol_srcptr p, fppol_srcptr q);


// <op>-multiplication by an <sq>-term polynomial,
// where <op> = [empty], add, or sub.
// Generic prototype:
//   void fppol_<op>mul_<sq>(fppol_ptr        r,
//                           fppol_srcptr     p,
//                           fppol<sq>_srcptr q);
#define __DECL_FPPOL_opMUL_xx(sq, op)              \
  void fppol_##op##mul_##sq(fppol_ptr          r,  \
                            fppol_srcptr       p,  \
                            fppol##sq##_srcptr q);


// <op>-multiplication of <sp> x <sq>-term polynomials,
// where <op> = [empty], add, or sub, and <sp> >= <sq>.
// Generic prototype:
//   void fppol_<op>mul_<sp>x<sq>(fppol_ptr        r,
//                                fppol<sp>_srcptr p,
//                                fppol<sq>_srcptr q);
#define __DECL_FPPOL_opMUL_xxxyy(sp, sq, op)              \
  void fppol_##op##mul_##sp##x##sq(fppol_ptr          r,  \
                                   fppol##sp##_srcptr p,  \
                                   fppol##sq##_srcptr q);


// <op>-multiplication of <sp> x <sq>-term polynomials,
// where <op> = [empty], add, or sub, and <sp> < <sq>.
// Generic prototype:
//   void fppol_<op>mul_<sp>x<sq>(fppol_ptr        r,
//                                fppol<sp>_srcptr p,
//                                fppol<sq>_srcptr q);
#define __DECL_FPPOL_opMUL_zzxyy(sp, sq, op)             \
  static inline                                          \
  void fppol_##op##mul_##sp##x##sq(fppol_ptr          r, \
                                   fppol##sp##_srcptr p, \
                                   fppol##sq##_srcptr q) \
  { fppol_##op##mul_##sq##x##sp(r, q, p); }


// All declarations bundled up into a single macro.
#define __DECL_FPPOL_opMUL_ALL(op)           \
        __DECL_FPPOL_opMUL      (        op) \
        __DECL_FPPOL_opMUL_xx   (16,     op) \
        __DECL_FPPOL_opMUL_xx   (32,     op) \
        __DECL_FPPOL_opMUL_xx   (64,     op) \
        __DECL_FPPOL_opMUL_xx   (16,     op) \
        __DECL_FPPOL_opMUL_xx   (32,     op) \
        __DECL_FPPOL_opMUL_xx   (64,     op) \
        __DECL_FPPOL_opMUL_xxxyy(16, 16, op) \
        __DECL_FPPOL_opMUL_xxxyy(32, 16, op) \
        __DECL_FPPOL_opMUL_xxxyy(32, 32, op) \
        __DECL_FPPOL_opMUL_xxxyy(64, 16, op) \
        __DECL_FPPOL_opMUL_xxxyy(64, 32, op) \
        __DECL_FPPOL_opMUL_xxxyy(64, 64, op) \
        __DECL_FPPOL_opMUL_zzxyy(16, 32, op) \
        __DECL_FPPOL_opMUL_zzxyy(16, 64, op) \
        __DECL_FPPOL_opMUL_zzxyy(32, 64, op)

__DECL_FPPOL_opMUL_ALL()
__DECL_FPPOL_opMUL_ALL(add)
__DECL_FPPOL_opMUL_ALL(sub)

#undef __DECL_FPPOL_opMUL
#undef __DECL_FPPOL_opMUL_xx
#undef __DECL_FPPOL_opMUL_xxxyy
#undef __DECL_FPPOL_opMUL_ALL

#endif  /* __FPPOL_MUL_H__ */
