//
// Copyright (C) 2006, 2007 INRIA (French National Institute for Research in
// Computer Science and Control)
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
//

/**
 * \file    macros.h
 * \author  Jerome Milan
 * \date    Thu Nov 23 2006
 * \version 1.0
 *
 * \brief Various CPP macros.
 *
 * Defines some C preprocessor macros that should be kept internal to
 * the TIFA library to avoid poluting client code.
 */

 /*
  *  License: GNU Lesser General Public License (LGPL)
  *  History:
  */

#if !defined(_TIFA_MACROS_H_)
   /**
    * \def _TIFA_MACROS_H_
    * Standard include guard.
    */
#define _TIFA_MACROS_H_

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__GMP_IMPL_H__)
    //
    // The following macros are defined in the gmp-impl.h header. Undefine
    // them to make sure we're always using TIFA's versions.
    //
    #undef MPN_NORMALIZE
    #undef SIZ
    #undef ABSIZ
    #undef PTR
    #undef ALLOC
    #undef MAX
    #undef MIN
    #undef ABS
#endif

    //
    // The following few macros are from the GMP library, in the gmp-impl.h
    // header. Since this header is not installed by "make install", chances
    // are that GMP users don't have access to it. Consequently, some macros
    // from gmp-impl.h are replicated here. These macros are only useful if
    // one wants to manipulate multi-precision integers with internal
    // functions from the mpn layer.
    //

   /**
    * \def MPN_NORMALIZE(dest, nlimbs)
    * Macro from the GMP library: Computes the effective size of an MPN number.
    *
    * Given \c dest, a pointer to an array of \c nlimbs \c mp_limbs_t
    * integers giving the representation of a multi-precision integer n,
    * computes the effective size of n, i.e the number of significant
    * \c mp_size_t integers needed to represent n and modifies the value
    * of \c nlimbs accordingly.
    *
    * \note This macro is originally the MPN_NORMALIZE macro from the GMP
    * library. It has been slightly modified.
    */
#define MPN_NORMALIZE(dest, nlimbs)                          \
do {                                                         \
    while (((nlimbs) > 0) && ((dest)[(nlimbs) - 1] == 0)) {  \
        (nlimbs)--;                                          \
    }                                                        \
} while (0)

   /**
    * \def SIZ(x)
    * Macro from the GMP library: Returns the \c _mp_size field of an \c mpz_t
    * integer.
    *
    * Returns the \c _mp_size field of the variable \c x of type <tt>mpz_t</tt>,
    * that is to say the number of \c mp_limbs_t integers needed to represent
    * the value of <tt>x</tt>. The sign of the returned value is given by the
    * sign of <tt>x</tt>'s value.
    *
    * \note This macro is the SIZ macro from the GMP library. It is
    * redistributed under the GNU LGPL license.
    */
#define SIZ(x) ((x)->_mp_size)

   /**
    * \def ABSIZ(x)
    * Macro from the GMP library: Returns the absolute value of
    * <tt>SIZ(x)</tt>.
    *
    * Returns the absolute value of <tt>SIZ(x)</tt> that is to say the number
    * of \c mp_limbs_t integers needed to represent the value of <tt>x</tt>.
    *
    * \note This macro is the ABSIZ macro from the GMP library. It is
    * redistributed under the GNU LGPL license.
    */
#define ABSIZ(x) (ABS(SIZ(x)))

   /**
    * \def PTR(x)
    * Macro from the GMP library: Returns the \c _mp_d field of an \c mpz_t
    * integer.
    *
    * Returns the \c _mp_d field of an \c mpz_t integer, that is to say a
    * pointer to an array of \c mp_limbs_t integers giving the representation
    * of the value of <tt>x</tt>.
    *
    * \note This macro is the PTR macro from the GMP library. It is
    * redistributed under the GNU LGPL license.
    */
#define PTR(x) ((x)->_mp_d)

   /**
    * \def ALLOC(x)
    * Macro from the GMP library: Returns the \c _mp_alloc field of an \c mpz_t
    * integer.
    *
    * Returns the \c _mp_alloc field of an \c mpz_t integer, that is to say
    * the size (in units of \c mp_limb_t) of the \c x->_mp_d array.
    *
    * \note This macro is the ALLOC macro from the GMP library. It is
    * redistributed under the GNU LGPL license.
    */
#define ALLOC(x) ((x)->_mp_alloc)

   /**
    * \def MAX(a, b)
    * Standard macro returning the maximum of a and b.
    *
    * \note As usual, be careful of possible side effects when using this kind
    * of macro. The standard disclaimers apply.
    *
    */
#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )

   /**
    * \def MIN(a, b)
    * Standard macro returning the minimum of a and b.
    *
    * \note As usual, be careful of possible side effects when using this kind
    * of macro. The standard disclaimers apply.
    *
    */
#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )

   /**
    * \def ABS(a)
    * Standard macro returning the absolute value of a.
    *
    * \note As usual, be careful of possible side effects when using this
    * kind of macro. The standard disclaimers apply.
    *
    */
#define ABS(a) ( ((a) < 0) ? (-(a)) : (a) )

   /**
    * \def IS_POWER_OF_2_UI(ui)
    * Macro returning a non-zero value if the unsigned integer \c ui is a
    * power of 2.
    *
    * \note As usual, be careful of possible side effects when using this
    * kind of macro. The standard disclaimers apply.
    *
    */
#define IS_POWER_OF_2_UI(ui) ( ((ui) & ((ui) - 1)) == 0 )

   /**
    * \def IS_EVEN(ui)
    * Macro returning True if the unsigned integer \c ui is even,
    * False otherwise.
    */
#define IS_EVEN(ui) (((ui) & 1) == 0)

   /**
    * \def IS_ODD(ui)
    * Macro returning True if the unsigned integer \c ui is odd,
    * False otherwise.
    */
#define IS_ODD(ui) (((ui) & 1) != 0)

   /**
    * \def ARE_EVEN(uia, uib)
    * Macro returning True if both of the unsigned integers \c uia and \c uib
    * are even, False otherwise.
    */
#define ARE_EVEN(uia, uib) ((((uia) | (uib)) & 1) == 0)

   /**
    * \def ARE_ODD(uia, uib)
    * Macro returning True if both of the unsigned integers \c uia and \c uib
    * are odd, False otherwise.
    */
#define ARE_ODD(uia, uib) ((((uia) | (uib)) & 1) != 0)

   /**
    * \def MPN_ADD(A, B, C)
    * 
    * Syntaxic sugar macro wrapping a call to <tt>mpn_add</tt>. Performs
    * size normalization on the result and takes care of the possible carry
    * out. Does not perform any reallocation: the user should make sure
    * the result has enough space to accomodate the possible carry out.
    *
    * Takes as parameters three \c mpz_t (and not arrays of <tt>mp_limb_t</tt>).
    * \c SIZ(B) should be greater than or equal to <tt>SIZ(C)</tt>.
    *
    * \see The GMP documentation for more information on the \c mpn_add
    * function.
    */
#define MPN_ADD(A, B, C)                                            \
    do {                                                            \
       if (mpn_add(PTR(A), PTR(B), SIZ(B), PTR(C), SIZ(C))) {       \
           SIZ(A) = SIZ(B);                                         \
           MPN_NORMALIZE(PTR(A), SIZ(A));                           \
           PTR(A)[SIZ(A)] = 1;                                      \
           SIZ(A)++;                                                \
       } else {                                                     \
           SIZ(A) = SIZ(B);                                         \
           MPN_NORMALIZE(PTR(A), SIZ(A));                           \
       }                                                            \
    } while (0)

   /**
    * \def MPN_ADD_CS(A, B, C)
    * 
    * Syntaxic sugar macro wrapping a call to <tt>mpn_add</tt>. Performs
    * size normalization on the result, takes care of the possible carry
    * out, and Checks the Sizes of the operands to call \c mpn_add with the
    * proper parameters' order. However, does not perform any reallocation: the 
    * user should make sure the result has enough space to accomodate the 
    * possible carry out.
    *
    * Takes as parameters three \c mpz_t (and not arrays of <tt>mp_limb_t</tt>).
    * \c B and \c C can be used interchangeably.
    *
    * \see The GMP documentation for more information on the \c mpn_add
    * function.
    */
#define MPN_ADD_CS(A, B, C)                                         \
    do {                                                            \
       if (SIZ(B) > SIZ(C)) {                                       \
           MPN_ADD(A, B, C);                                        \
       } else {                                                     \
           MPN_ADD(A, C, B);                                        \
       }                                                            \
    } while (0)

   /**
    * \def MPN_SUB(A, B, C)
    *
    * Syntaxic sugar macro wrapping a call to <tt>mpn_sub</tt>. Performs
    * size normalization on the result but does not take care of the possible 
    * borrow out.
    *
    * Takes as parameters three \c mpz_t (and not arrays of <tt>mp_limb_t</tt>).
    * \c B should be greater than or equal to <tt>C</tt>.
    *
    * \see The GMP documentation for more information on the \c mpn_sub
    * function.
    */
#define MPN_SUB(A, B, C)                                            \
    do {                                                            \
       mpn_sub(PTR(A), PTR(B), SIZ(B), PTR(C), SIZ(C));             \
       SIZ(A) = SIZ(B);                                             \
       MPN_NORMALIZE(PTR(A), SIZ(A));                               \
    } while (0)

   /**
    * \def MPN_SUB(A, B, C)
    * 
    * Syntaxic sugar macro wrapping a call to <tt>mpn_sub_n</tt>. Performs
    * size normalization on the result but does not take care of the possible 
    * borrow out.
    *
    * Takes as parameters three \c mpz_t (and not arrays of <tt>mp_limb_t</tt>).
    * \c B should be greater than or equal to <tt>C</tt>.
    *
    * \see The GMP documentation for more information on the \c mpn_sub_n
    * function.
    */
#define MPN_SUB_N(A, B, C)                                          \
    do {                                                            \
       mpn_sub_n(PTR(A), PTR(B), PTR(C), SIZ(B));                   \
       SIZ(A) = SIZ(B);                                             \
       MPN_NORMALIZE(PTR(A), SIZ(A));                               \
    } while (0)

   /**
    * \def MPN_TDIV_QR(Q, R, N, D)
    *
    * Syntaxic sugar macro wrapping a call to <tt>mpn_tdiv_qr</tt>. Performs
    * size normalization on both the quotient and remainder.
    *
    * Takes as parameters four \c mpz_t (and not arrays of <tt>mp_limb_t</tt>).
    *
    * \see The GMP documentation for more information on the \c mpn_tdiv_qr
    * function.
    */
#define MPN_TDIV_QR(Q, R, N, D)                                             \
    do {                                                                    \
       mpn_tdiv_qr(PTR(Q), PTR(R), 0,  PTR(N), SIZ(N), PTR(D), SIZ(D));     \
       SIZ(Q) = SIZ(N) - SIZ(D) + 1;                                        \
       MPN_NORMALIZE(PTR(Q), SIZ(Q));                                       \
       SIZ(R) = SIZ(D);                                                     \
       MPN_NORMALIZE(PTR(R), SIZ(R));                                       \
    } while (0)

   /**
    * \def MPN_MUL(A, B, C)
    *
    * Syntaxic sugar macro wrapping a call to <tt>mpn_mul</tt>. Performs
    * size normalization on the result.
    *
    * Takes as parameters three \c mpz_t (and not arrays of <tt>mp_limb_t</tt>).
    * \c SIZ(B) should be greater than or equal to <tt>SIZ(C)</tt>.
    *
    * \see The GMP documentation for more information on the \c mpn_mul
    * function.
     
    */
#define MPN_MUL(A, B, C)                                            \
    do {                                                            \
       mpn_mul(PTR(A), PTR(B), SIZ(B), PTR(C), SIZ(C));             \
       SIZ(A) = SIZ(B) + SIZ(C);                                    \
       MPN_NORMALIZE(PTR(A), SIZ(A));                               \
    } while (0)

   /**
    * \def MPN_MUL_N(A, B, C)
    * 
    * Syntaxic sugar macro wrapping a call to <tt>mpn_mul_n</tt>. Performs
    * size normalization on the result.
    * 
    * Takes as parameters three \c mpz_t (and not arrays of <tt>mp_limb_t</tt>).
    * \c SIZ(B) and \c SIZ(C) should be the same.
    *
    * \see The GMP documentation for more information on the \c mpn_mul_n
    * function.
    */
#define MPN_MUL_N(A, B, C)                                          \
    do {                                                            \
       mpn_mul_n(PTR(A), PTR(B), PTR(C), SIZ(B));                   \
       SIZ(A) = SIZ(B) << 1;                                        \
       MPN_NORMALIZE(PTR(A), SIZ(A));                               \
    } while (0)

   /**
    * \def MPN_MUL_CS(A, B, C)
    * 
    * Syntaxic sugar macro wrapping a call to <tt>mpn_mul</tt>. Performs
    * size normalization on the result, and Checks the Sizes of the operands to 
    * call \c mpn_mul with the proper parameters' order.
    * 
    * Takes as parameters three \c mpz_t (and not arrays of <tt>mp_limb_t</tt>).
    * \c B and \c C can be used interchangeably.
    *
    * \see The GMP documentation for more information on the \c mpn_mul
    * function.
    */
#define MPN_MUL_CS(A, B, C)                                         \
    do {                                                            \
       if (SIZ(B) > SIZ(C)) {                                       \
           mpn_mul(PTR(A), PTR(B), SIZ(B), PTR(C), SIZ(C));         \
       } else {                                                     \
           mpn_mul(PTR(A), PTR(C), SIZ(C), PTR(B), SIZ(B));         \
       }                                                            \
       SIZ(A) = SIZ(B) + SIZ(C);                                    \
       MPN_NORMALIZE(PTR(A), SIZ(A));                               \
    } while (0)

   /**
    * \def DECLARE_MPZ_SWAP_VARS
    * 
    * Macro declaring local variables needed by the \c MPZ_SWAP macro.
    * Should be called \emph once prior to any use of the \c MPZ_SWAP macro.
    *
    * \warning
    * Declares the variables \c __TMPPTR__MACROS_H__a9b3c01__ and
    * <tt>__TMPSIZ__MACROS_H__a9b3c01__<\tt>. Hopefully their names are fancy 
    * enough to avoid any local conflict.
    */
#define DECLARE_MPZ_SWAP_VARS                                       \
    mp_ptr __TMPPTR__MACROS_H__a9b3c01__;                           \
    mp_size_t __TMPSIZ__MACROS_H__a9b3c01__;

   /**
    * \def MPZ_SWAP(A, B)
    * 
    * Macro swapping the values of the two \c mpz_t \c A and <tt>B</tt>.
    *
    * \warning
    * The macro \c DECLARE_MPZ_SWAP_VARS should be called \emph once before 
    * using <tt>MPZ_SWAP</tt>.
    */
#define MPZ_SWAP(A, B)                                              \
    do {                                                            \
        __TMPPTR__MACROS_H__a9b3c01__ = PTR(A);                     \
        __TMPSIZ__MACROS_H__a9b3c01__ = SIZ(A);                     \
        PTR(A) = PTR(B);                                            \
        SIZ(A) = SIZ(B);                                            \
        PTR(B) = __TMPPTR__MACROS_H__a9b3c01__;                     \
        SIZ(B) = __TMPSIZ__MACROS_H__a9b3c01__;                     \
    } while (0)

#ifdef __cplusplus
}
#endif

#endif
