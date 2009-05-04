//
// Copyright (C) 2006, 2007, 2008 INRIA (French National Institute for Research
// in Computer Science and Control)
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
 * \file    sqrt_cont_frac.h
 * \author  Jerome Milan
 * \date    Mon Mar 10 2008
 * \version 1.0.2
 *
 * \brief Continued fraction expansion for square root of integers.
 *
 * Defines the continued fraction expansion for the square root of non-perfect
 * square integers.
 *
 * The expansion is computed via an iterative process, each step giving the
 * value of a new numerator. All the variables needed to perform this
 * computation is stored in an ad-hoc structure called
 * <tt>struct_cont_frac_state_t</tt>.
 *
 * \note Since the denominator of the continued fraction is not
 * used in the CFRAC algorithm, it is not computed here. Also, the numerator
 * of the fraction is only given modulo n. These restrictions are completely
 * trivial to fix should one need the complete approximation a/b of a square
 * root.
 */

/*
 * History:
 * --------
 *  1.0.2: Mon Mar 10 2008 by JM:
 *         - Inlined step_cont_frac_state(...) function.
 *  1.0.1: Mon Dec  4 2006 by JM:
 *         - Got rid of useless variables in the struct_cont_state_tstructure.
 *    1.0: Thu Mar 2 2006 by JM:
 *         - Initial version.
 */

#if !defined(_TIFA_SQRT_CONT_FRAC_H_)
   /**
    * \def _TIFA_SQRT_CONT_FRAC_H_
    * Standard include guard.
    */
#define _TIFA_SQRT_CONT_FRAC_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <inttypes.h>
#include <gmp.h>

   /**
    * \struct struct_cont_frac_state_t sqrt_cont_frac.h
    *                                  lib/utils/include/sqrt_cont_frac.h
    * \brief  An ad-hoc structure for the computation of the
    * continued fraction of a square root.
    *
    * This ad-hoc structure defines the variables needed for the computation
    * of the expansion of the continued fraction of the square root of a non
    * perfect square.
    *
    * \note Since the denominators of the fraction are not needed in the CFRAC
    * algorithm, they are not computed in this particular implementation.
    */
struct struct_cont_frac_state_t {
        //
        // _NOTE_: In the CFRAC algorithm, we don't need the value of the
        //         denominator, b.
        //
       /**
        * Current numerator of the continued fraction approximating \c sqrtn,
        * given modulo \c n.
        */
    mpz_t a;
        //
        // The following variable is not needed in the CFRAC algorithm.
        // Uncomment this declaration to also compute the denominator b.
        //
    //mpz_t b;
       /**
        * (Used in the computation of <tt>a</tt>).
        */
    mpz_t p;
       /**
        * One has: \c (-1)^(nsteps_performed) q = <tt>a*a - b*b*n</tt>
        */
    mpz_t q;
       /**
        * The current term of the expansion of the continued fraction.
        */
    mpz_t t;
       /**
        * The (truncated) square root of which to compute the continued
        * fraction.
        */
    mpz_t sqrtn;
       /**
        * The non perfect square integer whose square root will be approximated
        * by the computation of a continued fraction.
        */
    mpz_t n;
       /**
        * The number of non trivial terms of the continued fraction already
        * computed.
        */
    uint32_t nsteps_performed;
        //
        // "Private" temporary variables.
        // These are declared in the structure just to avoid redundant
        // initializations and deallocations.
        //
       /**
        * \internal
        * (Temporary variable declared in the structure just to avoid redundant
        * initializations and deallocations).
        */
    mpz_t _ztmp_q_;
       /**
        * \internal
        * (Temporary variable declared in the structure just to avoid redundant
        * initializations and deallocations).
        */
    mpz_t _ztmp_;
       /**
        * \internal
        * (Used in the computation of <tt>a</tt>).
        */
    mpz_t _a_old_;
        //
        // The following variables are not needed in the CFRAC algorithm.
        // Uncomment these declarations to also compute the denominator b.
        //
    //mpz_t _b_old_;
       /**
        * \internal
        * Pointer to the function used to compute convergents. This pointer
        * is initialized in the \c init_cont_frac_state function according
        * to the size of the number <tt>n</tt> whose square root has to be
        * approximated by a continued fraction.
        */
    void (*step_function)(struct struct_cont_frac_state_t* const, uint32_t);
};

   /**
    * \typedef cont_frac_state_t
    * \brief Equivalent to <tt>struct struct_cont_frac_state_t</tt>.
    */
typedef struct struct_cont_frac_state_t cont_frac_state_t;

   /**
    * \brief Initializes a <tt>cont_frac_state_t</tt>.
    *
    * Initializes a <tt>cont_frac_state_t</tt> to begin the computation of
    * a continued fraction. After invocation of this function, the fields
    * of \c state corresponds to the calculation of the second term of the
    * computed fraction, the first term beeing of course <tt>ceil(sqrt(n))</tt>.
    *
    *
    * \param[in] state A pointer to the <tt>cont_frac_state_t</tt> to
    *                  initialize.
    * \param[in] n     The non perfect square integer whose square root will
    *                  be approximated by the computation of a continued
    *                  fraction.
    */
void init_cont_frac_state(cont_frac_state_t* const state, const mpz_t n);

   /**
    * \brief Clears a <tt>cont_frac_state_t</tt>.
    *
    * Clears a <tt>cont_frac_state_t</tt>.
    *
    * \param[in] state A pointer to the <tt>cont_frac_state_t</tt> to clear.
    */
void clear_cont_frac_state(cont_frac_state_t* const state);

   /**
    * \brief Computes another term of a continued fraction.
    *
    * Computes another coefficient in the expansion of a continued fraction and
    * updates the structure <tt>state</tt>. The parameter \c nsteps gives the
    * number of iteration to perform. A new term is computed at each iteration.
    *
    * \note This function is actually given by <tt>state->step_function</tt>.
    *
    * \param[in] state A pointer to the <tt>cont_frac_state_t</tt>.
    * \param[in] nsteps Number of steps to perfom.
    */
inline static void step_cont_frac_state(cont_frac_state_t* const state,
                                        uint32_t nsteps) {
    state->step_function(state, nsteps);
}

#ifdef __cplusplus
}
#endif

#endif
