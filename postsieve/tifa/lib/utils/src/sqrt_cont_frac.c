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
 * \date    Tue Dec 19 2006
 * \version 1.0.2
 */

 /*
  *  Copyright (C) 2006, 2007 INRIA
  *  License: GNU Lesser General Public License (LGPL)
  *  History:
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

#include "gmp.h"
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

    uint32_t bits_in_n = mpz_sizeinbase(n, 2);

    mpz_init_set(state->n, n);

    mpz_init2(state->sqrtn, 1 + bits_in_n / 2);
    mpz_sqrt(state->sqrtn, n);

    mpz_init2(state->p, bits_in_n);
    mpz_init2(state->q, bits_in_n);

    mpz_set_ui(state->p, 0);
    mpz_set_ui(state->q, 1);

    mpz_init2(state->t, bits_in_n);
    mpz_set(state->t, state->sqrtn);

    state->nb_step_performed = 1;

    mpz_init2(state->_ztmp_q_, bits_in_n);
    //
    // Initial values for the a's and b's. These are the values after the
    // computation of one coefficient in the expansion of the continued
    // fraction. So we have to perform this first step for the other
    // variables in order to keep a coherent state...
    //
    mpz_init2(state->_a0_, bits_in_n);
    mpz_init2(state->_a1_, bits_in_n);
    //
    // _NOTE_: I ran into a hard to pin down bug under my Opteron machine.
    //         If state->a is initialized with bits_in_n bits I get a
    //         "realloc() : invalid next size" error for some numbers. Of
    //         course, it works without any problem on my PowerPC laptop (that
    //         would be too easy!).
    //         I need to spend more time on this... In the meanwhile,
    //         state->a is initialized with 2*bits_in_n bits...
    //         Ugly workaround, I know...
    //
    // _TO_DO_: Important: correct this realloc bug!!
    //
    mpz_init2(state->a, 2*bits_in_n);

    mpz_set_ui(state->_a0_, 1);
    mpz_set(state->_a1_, state->t);
    mpz_set(state->a, state->_a1_);

    //
    // Uncomment the folowing block to compute the b coefficient. This is of
    // no use for the CFRAC algorithm so it is not performed here...
    //
    //mpz_init_set_ui(state->_b0_, 0);
    //mpz_init_set_ui(state->_b1_, 1);
    //mpz_init_set(state->b, state->_b1_);

    //
    // Perform a first step for p,q and t to complete initialization...
    //
    mpz_add_ui(state->t, state->sqrtn, 0);
    mpz_tdiv_q_ui(state->t, state->t, 1);
    //
    // Now, t = (unsigned int)((sqrtn + p)/q);
    //
    mpz_addmul_ui(state->p, state->t, 1);
    //
    // Now, p = t*q - p_old;
    //
    mpz_set(state->q, state->n);
    mpz_submul(state->q, state->p, state->p);
    mpz_divexact_ui(state->q, state->q, 1);
    //
    // Now, q = (n - p*p)/q_old;
    //

    //
    // Choose step function according to maximum size of the variables
    // (i.e. 2*sqrt(N))
    //
    uint32_t bits_per_word = sizeof(mp_limb_t) << 3;
    uint32_t maxsize       = mpz_sizeinbase(state->sqrtn, 2);
    maxsize++;

    if (maxsize > bits_per_word) {
        //
        // _TO_DO_: I'm still not 100% confident in the MPN version of the
        //          step function. A thourough review should be made... 
        //          sometime...
        //
        state->step_function = &step_function_mpz;
        //
        // The faint of heart can still use the following MPZ step function...
        //
        //state->step_function = &step_function_mpn;
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
    mpz_clear(state->_a0_);
    mpz_clear(state->_a1_);

    //
    // Uncomment the folowing block if the b coefficient is computed.
    //
    //mpz_clear(state->b);
    //mpz_clear(state->_b0_);
    //mpz_clear(state->_b1_);

    mpz_clear(state->p);
    mpz_clear(state->q);
    mpz_clear(state->t);
    mpz_clear(state->n);
    mpz_clear(state->sqrtn);

    mpz_clear(state->_ztmp_q_);

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
        mpz_tdiv_q(state->t, state->t, state->q);
        //
        // Now, t = (unsigned int)((sqrtn + p)/q);
        //
        mpz_neg(state->p, state->p);
        mpz_addmul(state->p, state->t, state->q);
        //
        // Now, p = t*q - p_old;
        //
        mpz_set(state->_ztmp_q_, state->n);
        mpz_submul(state->_ztmp_q_, state->p, state->p);
        mpz_divexact(state->q, state->_ztmp_q_, state->q);
        //
        // Now, q = (n - p*p)/q_old;
        //
        mpz_set(state->a, state->_a0_);
        mpz_addmul(state->a, state->t, state->_a1_);
        //
        // Now, a = t*_a1_ + _a0_;
        //
        // _WARNING_: The numerator a of the continued fraction expansion
        //            is computed modulo n since this is what is of interest
        //            in the CFRAC algorithm.
        //
        mpz_mod(state->_a0_, state->_a1_, state->n);
        mpz_mod(state->_a1_, state->a,    state->n);
        mpz_set(state->a, state->_a1_);
        //
        // Uncomment the folowing block if the b coefficient is computed.
        // Note that b will be computed modulo n which may not be what you
        // want in the general case (i.e. not CFRAC related).
        //
        //mpz_set(state->b, state->_b0_);
        //mpz_addmul(state->b, state->t, state->_b1_);
        //mpz_mod(state->_b0_, state->_b1_, state->n);
        //mpz_mod(state->_b1_, state->b,    state->n);
        //mpz_set(state->b, state->_b1_);

        state->nb_step_performed++;
    }
    return;
}
//------------------------------------------------------------------------------
void step_function_mpn(cont_frac_state_t* const state, uint32_t nb_steps) {
    //
    // Mostly MPN version of the step function. Some MPZ functions remain.
    //
    // _NOTE_: Pure mpn code can quickly become un-understandable. Maintainers
    //         can still have a look at the mpz version of this function
    //         (i.e step_function_all_gmp) which does exactly the same thing...
    //
    mp_limb_t carry = 0;
    for (uint32_t i = 0; i < nb_steps; i++) {

       //
       // Computes t = (unsigned int)((sqrtn + p)/q);
       //
       carry = mpn_add(
                  PTR(state->t),
                  PTR(state->sqrtn), SIZ(state->sqrtn),
                  PTR(state->p),     SIZ(state->p)
              );

       SIZ(state->t) = SIZ(state->sqrtn);
       MPN_NORMALIZE(PTR(state->t), SIZ(state->t));

       if (carry != 0) {
           PTR(state->t)[SIZ(state->t)] = 1;
           SIZ(state->t)++;
       }

       mpn_tdiv_qr(PTR(state->t), PTR(state->a), 0,
                   PTR(state->t), SIZ(state->t),
                   PTR(state->q), SIZ(state->q));

       SIZ(state->t) = SIZ(state->t) - SIZ(state->q) + 1;
       MPN_NORMALIZE(PTR(state->t), SIZ(state->t));

        //
        // Computes a = t*q. Here a is just used as a temporary variable
        //
        if (SIZ(state->t) < SIZ(state->q)) {
            mpn_mul(PTR(state->a),
                    PTR(state->q), ABSIZ(state->q),
                    PTR(state->t), ABSIZ(state->t) );
        } else {
            mpn_mul(PTR(state->a),
                    PTR(state->t), ABSIZ(state->t),
                    PTR(state->q), ABSIZ(state->q) );
        }
        SIZ(state->a) = SIZ(state->t) + SIZ(state->q);
        MPN_NORMALIZE(PTR(state->a), SIZ(state->a));
        //
        // Computes p = a - p_old
        //
        mpn_sub_n(PTR(state->p),
                  PTR(state->a), PTR(state->p), SIZ(state->a));

        SIZ(state->p) = SIZ(state->a);
        MPN_NORMALIZE(PTR(state->p), SIZ(state->p));
        //
        // Now, p = t*q - p_old;
        //
        // Computes a = p*p. Again, this is not the "real" value of a. a is
        // just recycled as a temporary...
        //
        mpn_mul_n(PTR(state->a),
                  PTR(state->p), PTR(state->p), ABSIZ(state->p) );

        SIZ(state->a) = SIZ(state->p) << 1;
        MPN_NORMALIZE(PTR(state->a), SIZ(state->a));

        //
        // Computes _ztmp_q_ = n - a
        //
        if (SIZ(state->n) > SIZ(state->a)) {
            //
            // Since state->n and state->a can have different sizes, add
            // extra 0-bits padding in PTR(state->a) to make sure no leftover
            // values are taken into account in the substraction...
            //
            // Alternatively, mpn_sub could be used instead of mpn_sub_n but
            // we'll then have to propagate the borrow by hand...
            //
            int tmp_size = SIZ(state->n);
            do {
                PTR(state->a)[tmp_size-1] = 0;
                tmp_size--;
            } while (tmp_size != SIZ(state->a));
        }
        mpn_sub_n(PTR(state->_ztmp_q_),
                  PTR(state->n), PTR(state->a), SIZ(state->n));

        SIZ(state->_ztmp_q_) = SIZ(state->n);
        MPN_NORMALIZE(PTR(state->_ztmp_q_), SIZ(state->_ztmp_q_));

        //
        // Computes q = _ztmp_q_ / q_old
        //
        // _NOTE_: What? An mpz function in the middle of this mpn paradise?
        //         Well, actually, mpn_tdiv_qr is slower according to my
        //         tedious tests... It is certainly possible to speed this
        //         up by delving into GMP's internal to expose some non public
        //         function, but it's probably not worth it.
        //
        // _TO_DO_: Dive into GMP's wonderful innerworld and speed this up...
        //
        mpz_divexact(state->q, state->_ztmp_q_, state->q);

        //
        // Now, q = (n - p*p)/q_old;
        //
        // _NOTE_: What? mpz functions _again_? Yup! And again, faster than
        //         one call to mpn_add and mpn_mul. Maybe we could have a look
        //         at mpz_addmul's code to see how this is performed but again,
        //         this part is not really a bottleneck anyway...
        //
        // _TO_DO_: Find out why the mpn_add and mpn_mul is slower and come up
        //          with the good replacements from the mpn layer...
        //
        mpz_set(state->a, state->_a0_);
        //
        // _TO_DO_: Here the nasty realloc bug that has to be found and
        //          corrected... It is triggered by this call to mpz_addmul...
        //
        mpz_addmul(state->a, state->t, state->_a1_);
        //
        // Now, a = t*_a1_ + _a0_;
        //
        // _WARNING_: The numerator a of the continued fraction expansion
        //            is computed modulo n since this is what is of interest
        //            in the CFRAC algorithm.
        //
        if (SIZ(state->_a1_) >= SIZ(state->n)) {
            mpn_tdiv_qr(PTR(state->_ztmp_q_), PTR(state->_a0_), 0,
                        PTR(state->_a1_), SIZ(state->_a1_),
                        PTR(state->n), SIZ(state->n));

            SIZ(state->_a0_) = SIZ(state->n);
            MPN_NORMALIZE(PTR(state->_a0_), SIZ(state->_a0_));
        } else {
            mpz_swap(state->_a0_, state->_a1_);
        }
        if (SIZ(state->a) >= SIZ(state->n)) {
            mpn_tdiv_qr(PTR(state->_ztmp_q_), PTR(state->_a1_), 0,
                        PTR(state->a), SIZ(state->a),
                        PTR(state->n), SIZ(state->n));

            SIZ(state->_a1_) = SIZ(state->n);
            MPN_NORMALIZE(PTR(state->_a1_), SIZ(state->_a1_));
            mpz_set(state->a, state->_a1_);
        } else {
            mpz_set(state->_a1_, state->a);
        }
        //
        // Uncomment the folowing block if the b coefficient is computed.
        // Note that b will be computed modulo n which may not be what you
        // want in the general case (i.e. not CFRAC related).
        //
        // _NOTE_: The following commented block has not be double-check by
        //         running it. Maintainers are advised to double check it
        //         before going full speed ahead and uncommenting it...
        //
        //mpz_set(state->b, state->_b0_);
        //mpz_addmul(state->b, state->t, state->_b1_);
        //if (SIZ(state->_b1_) >= SIZ(state->n)) {
        //    mpn_tdiv_qr(PTR(state->_ztmp_q_), PTR(state->_b0_), 0,
        //                PTR(state->_b1_), SIZ(state->_b1_),
        //                PTR(state->n), SIZ(state->n));
        //    SIZ(state->_b0_) = SIZ(state->n);
        //    MPN_NORMALIZE(PTR(state->_b0_), SIZ(state->_b0_));
        //} else {
        //    mpz_swap(state->_b0_, state->_b1_);
        //}
        //if (SIZ(state->b) >= SIZ(state->n)) {
        //    mpn_tdiv_qr(PTR(state->_ztmp_q_), PTR(state->_b1_), 0,
        //                PTR(state->b), SIZ(state->b),
        //                PTR(state->n), SIZ(state->n));
        //    SIZ(state->_b1_) = SIZ(state->n);
        //    MPN_NORMALIZE(PTR(state->_b1_), SIZ(state->_b1_));
        //    mpz_set(state->b, state->_b1_);
        //} else {
        //    mpz_set(state->_b1_, state->b);
        //}

        state->nb_step_performed++;    
    }
    
    return;
}
//------------------------------------------------------------------------------
void step_function_mpn_ui(cont_frac_state_t* const state, uint32_t nb_steps) {
    //
    // Mostly MPN version of the step function, with a few single precision
    // computations when possible.
    //
    for (uint32_t i = 0; i < nb_steps; i++) {
        //
        // Computes t = (unsigned int)((sqrtn + p)/q);
        //
        VAL(t) = (VAL(sqrtn) + VAL(p)) / VAL(q);
        //
        // Now, t = (unsigned int)((sqrtn + p)/q);
        //
        VAL(p) = VAL(t) * VAL(q) - VAL(p);
        //
        // Now, p = t*q - p_old;
        //
        // Computes a = p*p. Again, this is not the "real" value of a. a is
        // just recycled as a temporary...
        //
        mpn_mul_n(PTR(state->a),
                  PTR(state->p), PTR(state->p), ABSIZ(state->p) );

        SIZ(state->a) = SIZ(state->p) << 1;
        MPN_NORMALIZE(PTR(state->a), SIZ(state->a));
        //
        // Computes _ztmp_q_ = n - a
        //
        mpn_sub_n(PTR(state->_ztmp_q_),
                  PTR(state->n), PTR(state->a), SIZ(state->n));

        SIZ(state->_ztmp_q_) = SIZ(state->n);
        MPN_NORMALIZE(PTR(state->_ztmp_q_), SIZ(state->_ztmp_q_));
        //
        // Computes q = _ztmp_q_ / q_old
        //
        // _NOTE_: What? An mpz function in the middle of this mpn paradise?
        //         Well, actually, mpn_tdiv_qr is slower according to my
        //         tedious tests... It is certainly possible to speed this
        //         up by delving into GMP's internal to expose some non public
        //         function, but it's probably not worth it.
        //
        // _TO_DO_: Dive into GMP's wonderful innerworld and speed this up...
        //
        mpz_divexact(state->q, state->_ztmp_q_, state->q);
        //
        // Now, q = (n - p*p)/q_old;
        //
        // _NOTE_: What? mpz functions _again_? Yup! And again, faster than
        //         one call to mpn_add and mpn_mul. Maybe we could have a look
        //         at mpz_addmul's code to see how this is performed but again,
        //         this part is not really a bottleneck anyway...
        //
        // _TO_DO_: Find out why the mpn_add and mpn_mul is slower and come up
        //          with the good replacements from the mpn layer...
        //
        mpz_set(state->a, state->_a0_);
        mpz_addmul(state->a, state->t, state->_a1_);
        //
        // Now, a = t*_a1_ + _a0_;
        //
        // _WARNING_: The numerator a of the continued fraction expansion
        //            is computed modulo n since this is what is of interest
        //            in the CFRAC algorithm.
        //
        if (SIZ(state->_a1_) >= SIZ(state->n)) {
            mpn_tdiv_qr(PTR(state->_ztmp_q_), PTR(state->_a0_), 0,
                        PTR(state->_a1_), SIZ(state->_a1_),
                        PTR(state->n), SIZ(state->n));

            SIZ(state->_a0_) = SIZ(state->n);
            MPN_NORMALIZE(PTR(state->_a0_), SIZ(state->_a0_));
        } else {
            mpz_swap(state->_a0_, state->_a1_);
        }
        if (SIZ(state->a) >= SIZ(state->n)) {
            mpn_tdiv_qr(PTR(state->_ztmp_q_), PTR(state->_a1_), 0,
                        PTR(state->a), SIZ(state->a),
                        PTR(state->n), SIZ(state->n));

            SIZ(state->_a1_) = SIZ(state->n);
            MPN_NORMALIZE(PTR(state->_a1_), SIZ(state->_a1_));
            mpz_set(state->a, state->_a1_);
        } else {
            mpz_set(state->_a1_, state->a);
        }
        //
        // Uncomment the folowing block if the b coefficient is computed.
        // Note that b will be computed modulo n which may not be what you
        // want in the general case (i.e. not CFRAC related).
        //
        // _NOTE_: The following commented block has not be double-check by
        //         running it. Maintainers are advised to double check it
        //         before going full speed ahead and uncommenting it...
        //
        //mpz_set(state->b, state->_b0_);
        //mpz_addmul(state->b, state->t, state->_b1_);
        //if (SIZ(state->_b1_) >= SIZ(state->n)) {
        //    mpn_tdiv_qr(PTR(state->_ztmp_q_), PTR(state->_b0_), 0,
        //                PTR(state->_b1_), SIZ(state->_b1_),
        //                PTR(state->n), SIZ(state->n));
        //    SIZ(state->_b0_) = SIZ(state->n);
        //    MPN_NORMALIZE(PTR(state->_b0_), SIZ(state->_b0_));
        //} else {
        //    mpz_swap(state->_b0_, state->_b1_);
        //}
        //if (SIZ(state->b) >= SIZ(state->n)) {
        //    mpn_tdiv_qr(PTR(state->_ztmp_q_), PTR(state->_b1_), 0,
        //                PTR(state->b), SIZ(state->b),
        //                PTR(state->n), SIZ(state->n));
        //    SIZ(state->_b1_) = SIZ(state->n);
        //    MPN_NORMALIZE(PTR(state->_b1_), SIZ(state->_b1_));
        //    mpz_set(state->b, state->_b1_);
        //} else {
        //    mpz_set(state->_b1_, state->b);
        //}

        state->nb_step_performed++;
    }
    return;
}
//------------------------------------------------------------------------------
void step_function_ui(cont_frac_state_t* const state, uint32_t nb_steps) {
    //
    // Version of the step function with mostly single precision computations.
    //
    for (uint32_t i = 0; i < nb_steps; i++) {

        VAL(t) = (VAL(sqrtn) + VAL(p)) / VAL(q);
        //
        // Now, t = (unsigned int)((sqrtn + p)/q);
        //
        VAL(p) = VAL(t) * VAL(q) - VAL(p);
        //
        // Now, p = t*q - p_old;
        //
        VAL(q) = (VAL(n) - VAL(p) * VAL(p)) / VAL(q);
        //
        // Now, q = (n - p*p)/q_old;
        //
        VAL(a) = VAL(_a0_);
        //
        // On 32 bit architectures, the following call to a GMP function can
        // be avoided by casting the operands to uint64_t integers, taking
        // the modulo and casting back to a mp_limb_t (32 bits in that case).
        // This is not quite doable on 64 bit architectures where most of
        // the time sizeof(long long int) is 8 bytes, just like on 32 bit
        // systems.
        //
        mpz_addmul(state->a, state->t, state->_a1_);
        //
        // Now, a = t*_a1_ + _a0_;
        //
        // _WARNING_: The numerator a of the continued fraction expansion
        //            is computed modulo n since this is what is of interest
        //            in the CFRAC algorithm.
        //
        VAL(_a0_) = VAL(_a1_);
        mpz_mod_ui(state->_a1_, state->a, VAL(n));
        VAL(a)         = VAL(_a1_);
        SIZ(state->a) = 1;
        //
        // Uncomment the folowing block if the b coefficient is computed.
        // Note that b will be computed modulo n which may not be what you
        // want in the general case (i.e. not CFRAC related).
        //
        //VAL(b)  = VAL(_b0_);
        //mpz_addmul(state->a, state->t, state->_b1_);
        //VAL(_b0_) = VAL(_b1_);
        //mpz_mod_ui(state->_b1_, state->b, VAL(n));
        //VAL(a)         = VAL(_a1_);
        //SIZ(state->a) = 1;

        state->nb_step_performed++;
    }
    return;
}
//------------------------------------------------------------------------------
