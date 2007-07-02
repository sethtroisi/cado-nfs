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
    // them to make sure we're always using the TIFA's versions.
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
    * library and it has been slightly modified.
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

#ifdef __cplusplus
}
#endif

#endif
