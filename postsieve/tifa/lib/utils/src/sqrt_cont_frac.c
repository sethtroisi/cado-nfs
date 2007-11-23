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
 * \file    sqrt_cont_frac.c
 * \author  Jerome Milan
 * \date    Wed July 4 2007
 * \version 1.1
 */

 /*
  *  Copyright (C) 2006, 2007 INRIA
  *  License: GNU Lesser General Public License (LGPL)
  *  History:
  *
  *  1.1: Wed July 4 2007 by JM:
  *         - Partly rewrote the step functions to add minor optimizations.
  *
  *  1.0.2: Tue Dec 19 2006 by JM:
  *         - Choose the step function according to the size of n, the number
  *           whose square root has to be expanded. Also added an mpn version
  *           of the general case step function.
  *
  *  1.0.1: Mon Dec  4 2006 by JM:
  *         - Got rid of the floating-point operations (how did I make such
  *           a blunder?) This is not only performance-related though. Indeed,
  *           floating point operations gave wrong results if the precision
  *           wasn't set high enough. /me still wonders how he could have been
  *           so stupid to use completely useless floating-point operations...
  *
  *  1.0.0: Thu Mar 2 2006 by JM:
  *         - Initial version.
  */

#include <stdlib.h>

#include "tifa_config.h"
#include <gmp.h>
#if TIFA_USE_GMP_INTERNAL_FUNCS
    #include "gmp-impl.h"
#endif
#include "gmp_utils.h"
#include "sqrt_cont_frac.h"
#include "macros.h"

#include <stdio.h>

//------------------------------------------------------------------------------
#define VAL(X) PTR(state->X)[0]
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void step_function_mpz(cont_frac_state_t* const state, uint32_t nb_steps);
//------------------------------------------------------------------------------
void step_function_mpn(cont_frac_state_t* const state, uint32_t nb_steps);
//------------------------------------------------------------------------------
void step_function_mpn_ui(cont_frac_state_t* const state, uint32_t nb_steps);
//------------------------------------------------------------------------------
void step_function_ui(cont_frac_state_t* const state, uint32_t nb_steps);
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void init_cont_frac_state(cont_frac_state_t* const state, const mpz_t n) {

    uint32_t size_n      = mpz_sizeinbase(n, 2);
    uint32_t size_sqrtn  = 1 + size_n / 2;
    uint32_t size_sqrt2n = 1 + size_sqrtn;

    mpz_init_set(state->n, n);

    mpz_init2(state->sqrtn, size_sqrtn);
    mpz_sqrt(state->sqrtn, n);

    mpz_init2(state->p, size_sqrt2n);
    mpz_init2(state->q, size_sqrt2n);
    //
    // The way it is used in the step functions, state->_ztmp_ is bounded by
    // n, but it is used in the MPN step function to hold the result of
    // the function mpn_mul to get state->p * state->p. Since mpn_mul(X, A, A)
    // require that ALLOC(X) >= 2 * SIZ(A) even if the result fits in
    // 2 * SIZ(A) - 1 limbs, we should make sure to allocate enough limbs for
    // state->_ztmp_...
    //
    mpz_init2(state->_ztmp_,  2 * ALLOC(state->p) * GMP_NUMB_BITS);

    mpz_set_ui(state->p, 0);
    mpz_set_ui(state->q, 1);

    mpz_init2(state->t, size_sqrt2n);
    mpz_set(state->t, state->sqrtn);

    state->nsteps_performed = 1;

    mpz_init2(state->_ztmp_q_, 2 * size_sqrt2n);
    //
    // state->_a_old_ is bounded by n but the step functions use mpz_swap's 
    // instead of mpz_set's, so _a_old_ should have as many limbs as
    // state->a...
    //
    mpz_init2(state->_a_old_, 2 * size_n);
    mpz_init2(state->a, 2 * size_n);
    //
    // Initial values for the a's and b's. These are the values after the
    // computation of one coefficient in the expansion of the continued
    // fraction. So we have to perform this first step for the other
    // variables in order to keep a coherent state...
    //
    mpz_set_ui(state->_a_old_, 1);
    mpz_set(state->a, state->t);

    //
    // Uncomment the folowing block to compute the b coefficient. This is of
    // no use for the CFRAC algorithm so it is not performed here...
    //
    //mpz_init2(state->_b_old_, 2 * size_n);
    //mpz_init2(state->b, 2 * size_n);
    //mpz_set_ui(state->_b_old_, 1);
    //mpz_set(state->b, state->t);

    //
    // Perform a first step for p, q and t to complete initialization...
    //
    mpz_add_ui(state->t, state->sqrtn, 0);
    mpz_tdiv_q_ui(state->t, state->t, 1);
    //
    // Now, t = (unsigned int)((sqrtn + p)/q);
    //
    mpz_addmul_ui(state->p, state->t, 1);
    //
    // Now, p = t*q - p"old";
    //
    mpz_set(state->q, state->n);
    mpz_submul(state->q, state->p, state->p);
    mpz_divexact_ui(state->q, state->q, 1);
    //
    // Now, q = (n - p*p) / q"old";
    //

    //
    // Choose step function according to maximum size of the variables (i.e.
    // 2*sqrt(N))
    //
    uint32_t bits_per_word = GMP_NUMB_BITS;
    uint32_t maxsize       = mpz_sizeinbase(state->sqrtn, 2);
    maxsize++;

    if (maxsize > bits_per_word) {
        state->step_function = &step_function_mpn;
        //
        // The faint of heart can use the following MPZ step function to
        // be on the safe side...
        //
        //state->step_function = &step_function_mpz;
    } else {
        if (SIZ(n) > 1) {
            state->step_function = &step_function_mpn_ui;
        } else {
            state->step_function = &step_function_ui;
        }
    }
    return;
}
//------------------------------------------------------------------------------
void clear_cont_frac_state(cont_frac_state_t* const state) {
    mpz_clear(state->a);
    mpz_clear(state->_a_old_);
    //
    // Uncomment the following block if the b coefficient is computed.
    //
    //mpz_clear(state->b);
    //mpz_clear(state->_b_old_);

    mpz_clear(state->p);
    mpz_clear(state->q);
    mpz_clear(state->t);
    mpz_clear(state->n);
    mpz_clear(state->sqrtn);

    mpz_clear(state->_ztmp_q_);
    mpz_clear(state->_ztmp_);

    return;
}
//------------------------------------------------------------------------------
inline void
step_cont_frac_state(cont_frac_state_t* const state, uint32_t nb_steps) {
    state->step_function(state, nb_steps);
}
//------------------------------------------------------------------------------

    //
    //                     "Private" functions
    //

//------------------------------------------------------------------------------
void step_function_mpz(cont_frac_state_t* const state, uint32_t nb_steps) {

    for (uint32_t i = 0; i < nb_steps; i++) {
        mpz_add(state->t, state->sqrtn, state->p);
        mpz_tdiv_qr(state->t, state->_ztmp_q_, state->t, state->q);
        mpz_sub(state->p, state->sqrtn, state->_ztmp_q_);
        //
        // Now, t = floor((sqrtn + p) / q)
        //      p = t*q - p = sqrtn - ((sqrtn + p) % q)
        //
        mpz_set(state->_ztmp_q_, state->n);
        mpz_submul(state->_ztmp_q_, state->p, state->p);
        mpz_divexact(state->q, state->_ztmp_q_, state->q);
        //
        // Now, q = (n - p*p) / q
        //
        mpz_swap(state->a, state->_a_old_);
        mpz_addmul(state->a, state->t, state->_a_old_);
        mpz_mod(state->a, state->a, state->n);
        //
        // Now, a"new" = t*a"old" + _a_old_ and _a_old_ = a"old"
        //
        // _WARNING_: The numerator 'a' of the continued fraction expansion
        //            is computed modulo n since this is what is of interest
        //            in the CFRAC algorithm.
        //
        // Uncomment the folowing instructions if the b coefficient should be
        // computed. Note that b will be computed modulo n which may not be
        // what you want in the general case (i.e. not CFRAC related).
        //
        //mpz_swap(state->b, state->_b_old_);
        //mpz_addmul(state->b, state->t, state->_b_old_);
        //mpz_mod(state->b, state->b, state->n);

        state->nsteps_performed++;
    }
    return;
}
//------------------------------------------------------------------------------
void step_function_mpn(cont_frac_state_t* const state, uint32_t nb_steps) {

    DECLARE_MPZ_SWAP_VARS;

    //
    // Mostly MPN version of the step function. Some MPZ functions remain.
    //
    for (uint32_t i = 0; i < nb_steps; i++) {

        MPN_ADD(state->_ztmp_q_, state->sqrtn, state->p);
        MPN_TDIV_QR(state->t, state->_ztmp_, state->_ztmp_q_, state->q);
        MPN_SUB(state->p, state->sqrtn, state->_ztmp_);
        //
        // Now, t = floor((sqrtn + p) / q)
        //      p = t*q - p = sqrtn - ((sqrtn + p) % q)
        //

        //
        // Computes q = (n - p*p) / q
        //
        MPN_MUL_N(state->_ztmp_, state->p, state->p);
        MPN_SUB(state->_ztmp_q_, state->n, state->_ztmp_);
        //
        // _NOTE_: What? An mpz function in the middle of this mpn paradise?
        //         Well, actually, mpn_tdiv_qr is slower according to my
        //         tedious tests... It is certainly possible to speed this
        //         up by delving into GMP's internal to expose some non public
        //         function, but it's probably not worth it.
        //
        // _NOTE_: There is a way to express these recurrence relations which
        //         can bypass the following division. However, benchmarks
        //         showed that avoiding this division lead to an overhead
        //         actually resulting in a slower code on Opteron!
        //
        // _TO_DO_: Dive into GMP's wonderful innerworld and speed this up...
        //
#if 1
        mpz_divexact(state->q, state->_ztmp_q_, state->q);
#else
        MPN_TDIV_QR(state->q, state->_ztmp_, state->_ztmp_q_, state->q);
#endif
        //
        // Now, q = (n - p*p) / q"old"
        //
        MPZ_SWAP(state->a, state->_a_old_);
        //
        // _NOTE_: What? mpz functions _again_? Yup! And again, faster than
        //         one call to mpn_add and mpn_mul. Maybe we could have a look
        //         at mpz_addmul's code to see how this is performed but again,
        //         this part is not really a bottleneck anyway...
        //
        // _TO_DO_: Find out why the mpn_add and mpn_mul is slower and come up
        //          with the good replacements from the mpn layer...
        //
#if 1
        mpz_addmul(state->a, state->t, state->_a_old_);
#else
        MPN_MUL_CS(state->_ztmp_, state->t, state->_a_old_);
        MPN_ADD_CS(state->a, state->_ztmp_, state->a);
#endif

        if (SIZ(state->a) >= SIZ(state->n)) {
            MPN_TDIV_QR(state->_ztmp_q_, state->a, state->a, state->n);
        }
        //
        // Now, a"new" = t*a"old" + _a_old_ and _a_old_ = a"old"
        //
        // _WARNING_: The numerator a of the continued fraction expansion
        //            is computed modulo n since this is what is of interest
        //            in the CFRAC algorithm.
        //
        // Uncomment the folowing block if the b coefficient is computed.
        // Note that b will be computed modulo n which may not be what you
        // want in the general case (i.e. not CFRAC related).
        //
        // _NOTE_: The following commented block has not be double-check by
        //         running it. Maintainers are advised to double check it
        //         before going full speed ahead and uncommenting it...
        //
        //MPZ_SWAP(state->b, state->_b_old_);
        //mpz_addmul(state->b, state->t, state->_b_old_);
        //if (SIZ(state->b) >= SIZ(state->n)) {
        //    MPN_TDIV_QR(state->_ztmp_q_, state->b, state->b, state->n);
        //}

        state->nsteps_performed++;
   }
   return;
}
//------------------------------------------------------------------------------
void step_function_mpn_ui(cont_frac_state_t* const state, uint32_t nb_steps) {

    //
    // Mostly MPN version of the step function, with a few single precision
    // computations when possible.
    //
    DECLARE_MPZ_SWAP_VARS;

    for (uint32_t i = 0; i < nb_steps; i++) {

        VAL(t) = (VAL(sqrtn) + VAL(p)) / VAL(q);
        VAL(p) = VAL(t) * VAL(q) - VAL(p);
        //
        // Now, t = floor((sqrtn + p) / q)
        //      p = t*q - p = sqrtn - ((sqrtn + p) % q)
        //

        //
        // Computes q = (n - p*p) / q
        //
        MPN_MUL_N(state->_ztmp_, state->p, state->p);
        MPN_SUB(state->_ztmp_q_, state->n, state->_ztmp_);

        //
        // _NOTE_: What? An mpz function in the middle of this mpn paradise?
        //         Well, actually, mpn_tdiv_qr is slower according to my
        //         tedious tests... It is certainly possible to speed this
        //         up by delving into GMP's internal to expose some non public
        //         function, but it's probably not worth it.
        //
        // _NOTE_: There is a way to express these recurrence relations which
        //         can bypass the following division. However, benchmarks
        //         showed that avoiding this division lead to an overhead
        //         actually resulting in a slower code on Opteron!
        //
        // _TO_DO_: Dive into GMP's wonderful innerworld and speed this up...
        //
#if 1
        mpz_divexact(state->q, state->_ztmp_q_, state->q);
#else
        MPN_TDIV_QR(state->q, state->_ztmp_, state->_ztmp_q_, state->q);
#endif
        //
        // Now, q = (n - p*p) / q
        //
        MPZ_SWAP(state->a, state->_a_old_);
        //
        // _NOTE_: What? mpz functions _again_? Yup! And again, faster than
        //         one call to mpn_add and mpn_mul. Maybe we could have a look
        //         at mpz_addmul's code to see how this is performed but again,
        //         this part is not really a bottleneck anyway...
        //
        // _TO_DO_: Find out why the mpn_add and mpn_mul is slower and come up
        //          with the good replacements from the mpn layer...
        //
#if 1
        mpz_addmul(state->a, state->t, state->_a_old_);
#else
        MPN_MUL_CS(state->_ztmp_, state->t, state->_a_old_);
        MPN_ADD_CS(state->a, state->_ztmp_, state->a);
#endif
        if (SIZ(state->a) >= SIZ(state->n)) {
            MPN_TDIV_QR(state->_ztmp_q_, state->a, state->a, state->n);
        }
        //
        // Now, a"new" = t*a"old" + _a_old_ and _a_old_ = a"old"
        //
        // _WARNING_: The numerator a of the continued fraction expansion
        //            is computed modulo n since this is what is of interest
        //            in the CFRAC algorithm.
        //
        // Uncomment the folowing block if the b coefficient is computed.
        // Note that b will be computed modulo n which may not be what you
        // want in the general case (i.e. not CFRAC related).
        //
        // _NOTE_: The following commented block has not be double-check by
        //         running it. Maintainers are advised to double check it
        //         before going full speed ahead and uncommenting it...
        //
        //MPZ_SWAP(state->b, state->_b_old_);
        //mpz_addmul(state->b, state->t, state->_b_old_);
        //if (SIZ(state->b) >= SIZ(state->n)) {
        //    MPN_TDIV_QR(state->_ztmp_q_, state->b, state->b, state->n);
        //}

        state->nsteps_performed++;
    }
    return;
}
//------------------------------------------------------------------------------
void step_function_ui(cont_frac_state_t* const state, uint32_t nb_steps) {

    mp_limb_t _tmp_;
    //
    // Version of the step function with mostly single precision computations.
    //
    for (uint32_t i = 0; i < nb_steps; i++) {

        VAL(t) = (VAL(sqrtn) + VAL(p)) / VAL(q);
        //
        // Now, t = floor((sqrtn + p) / q);
        //
        VAL(p) = VAL(t) * VAL(q) - VAL(p);
        //
        // Now, p = t*q - p"old";
        //
        VAL(q) = (VAL(n) - VAL(p) * VAL(p)) / VAL(q);
        //
        // Now, q = (n - p*p) / q"old";
        //
        _tmp_        = VAL(_a_old_);
        VAL(_a_old_) = VAL(a);
        VAL(a)       = _tmp_;
        //
        // On 32 bit architectures, the following call to a GMP function can
        // be avoided by casting the operands to uint64_t integers, taking
        // the modulo and casting back to a mp_limb_t (32 bits in that case).
        // This is not quite doable on 64 bit architectures where most of
        // the time sizeof(long long unsigned int) is 8 bytes, just as on 32
        // bit systems.
        //
        mpz_addmul(state->a, state->t, state->_a_old_);
        mpz_mod_ui(state->a, state->a, VAL(n));
        //
        // Now, a"new" = t*a"old" + _a_old_ and _a_old_ = a"old"
        //
        // _WARNING_: The numerator a of the continued fraction expansion
        //            is computed modulo n since this is what is of interest
        //            in the CFRAC algorithm.
        //
        // Uncomment the folowing block if the b coefficient is computed.
        // Note that b will be computed modulo n which may not be what you
        // want in the general case (i.e. not CFRAC related).
        //
        //_tmp_        = VAL(_b_old_);
        //VAL(_b_old_) = VAL(b);
        //VAL(b)       = _tmp_;
        //mpz_addmul(state->b, state->t, state->_b_old_);
        //mpz_mod_ui(state->b, state->b, VAL(n));

        state->nsteps_performed++;
    }
    return;
}
//------------------------------------------------------------------------------
