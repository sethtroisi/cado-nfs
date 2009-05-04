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
 * \file    squfof.c
 * \author  Jerome Milan
 * \date    Mon Mar 3 2008
 * \version 1.4
 */

 /*
  * History:
  *     1.4: Mon Mar 3 2008 by JM
  *        - Replaced ugly macros by inline code. Unsure if that affect 
  *          performances though...
  *   1.3.2: Wed Feb 27 2008 by JM
  *        - Ugly, ugly hack in LARGE_STEP macro to avoid infinite loop.
  *   1.3.1: Wed Dec 19 2007 by JM
  *        - Fixed two mighty bugs in LARGE_STEP macro.
  *   1.3:   Thu Aug 2 2007 by JM
  *        - Implemented Williams & Wunderlich's "large step" algorithm for
  *          a "fast return" in step 4) of the SQUFOF algorithm.
  *   1.2.1: Tue Jun 27 2007 by JM
  *        - Suppressed user parameters. Now performs a multiplier race only
  *          if needed.
  *   1.2: Fri Jun  1 2007 by JM
  *        - Added multiplier race.
  *   1.1: Wed May 30 2007 by JM
  *        - Fixed (stupid) bug to factor (roughly) double precision integers.
  *   1.0: Wed Jan 24 2007 by JM
  *        - Initial version.
  */

  //------------------------------------------------------------------------
  // References:
  //------------------------------------------------------------------------
  // - For a thorough explanation of the SQUFOF algorithm:
  //         "Square Form Factorization",
  //         Jason E. Gowen, Samuel S. Wagstaff Jr.,
  //         Mathematics of Computation,
  //         S 0025-5718(07)02010-8,
  //         Article electronically published on May 14, 2007.
  //
  // - For an alternate explanation, including Shanks' very notes on SQUFOF,
  //   see David Joyner and Stephen McMath's work at:
  //       http://cadigweb.ew.usna.edu/~wdj/mcmath
  //   URL last retrieved on Wed Aug 1 2007.
  //
  // - For a description of the "large step algorithm" used to quickly
  //   jump over several forms:
  //         "On the Parallel Generation of the Residues for the Continued
  //         Fraction Factoring Algorithm",
  //         Hugh C. Williams, Marvin C. Wunderlich,
  //         Mathematics of Computation, Volume 48, Number 177,
  //         January 1987, pages 405-423
  //------------------------------------------------------------------------

#include "tifa_config.h"

#include <stdlib.h>
#include <stdbool.h>
#include <inttypes.h>
#include <gmp.h>

#include "funcs.h"
#include "macros.h"
#include "linked_list.h"
#include "factoring_machine.h"
#include "squfof.h"
#include "tifa_factor.h"

#define SQUFOF_DEBUG 0
#if SQUFOF_DEBUG
  #include <stdio.h>
  #define SQUFOF_GMP_PRINT(...) gmp_printf(__VA_ARGS__)
  #define SQUFOF_PRINT(...)     printf(__VA_ARGS__)
#else
  #define SQUFOF_GMP_PRINT(...) /* intentionally left empty */
  #define SQUFOF_PRINT(...)     /* intentionally left empty */
#endif

#define __PREFIX__  "squfof: "
#define __VERBOSE__ TIFA_VERBOSE_SQUFOF
#define __TIMING__  TIFA_TIMING_SQUFOF

#include "messages.h"

//-----------------------------------------------------------------------------
//                       GENERAL NON PUBLIC DEFINE(S)
//-----------------------------------------------------------------------------
//
// Maximum size of the queues.
//
#define MAX_QUEUE_ALLOC  64
//
// Number of multipliers used if a multiplier race is performed.
//
#define NMULTIPLIERS     15
//
// Number of forward cycle iterations to perform in a row for each racer.
//
#define NSTEP_CYCLE_FWD  2000


//-----------------------------------------------------------------------------
//          SYMBOLS AND FUNCTIONS RELATED TO OUR CHEAP "DIVISION TRICK"
//-----------------------------------------------------------------------------
//
// A switch toggling our cheap division trick on and off
//
#define USE_DIVISION_TRICK 1
//
// Switch from the standard division to our cheap trick here, to avoid
// disseminating #if's in the code...
//
#if USE_DIVISION_TRICK
    #define divide(N, D)  divide_trick((N), (D))
#else
    #define divide(N, D)  ((N) / (D))
#endif
//-----------------------------------------------------------------------------
inline static unsigned long int
divide_trick(const unsigned long int N, const unsigned long int D);  
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//           SYMBOLS AND FUNCTIONS RELATED TO THE LARGE STEP ALGORITHM
//-----------------------------------------------------------------------------

//
// _NOTE_: We only implement here a simplified version of Williams &
//         Wunderlich's large step algorithm. This is used in step 4 of the
//         SQUFOF algorithm (using Gowen & Wagstaff's description) to achieve
//         what Shanks calls the "fast return". Note that Gowen & Wagstaff's
//         description does not include this "fast return".
//

//
// Should we perform the "fast return" variant?
//
// _NOTE_: The "fast return" variant is now enabled by default since the
//         mighty infinite loop bug has been solved (hopefully for good!).
//
#define PERFORM_FAST_RETURN 1

//
// Only use the large step algorithm if the number of forms to compute (from
// the inverse square root of a found square form) is larger than a certain
// threshold.
//
#if PERFORM_FAST_RETURN
    #define LARGE_STEP_THRESHOLD 4096
#else
    #define LARGE_STEP_THRESHOLD TIFA_LONG_MAX
#endif

//
// Mininum number of forms to step through sequentially before using the
// large step algorithm. If MIN_NSTEP_LS is too small, the large step's
// overhead will degrade performances. If set too high, we'll be less
// precise in the length's estimation of the large step.
//
#define MIN_NSTEP_LS 128
//-----------------------------------------------------------------------------
inline static void
single_step(
    unsigned long int* const P,
    unsigned long int* const Q,
    unsigned long int* const QQ,
    const unsigned long int S
);
//-----------------------------------------------------------------------------
inline static unsigned long int
large_step(
    unsigned long int* const P,
    unsigned long int* const Q,
    unsigned long int* const QQ,
    unsigned long int PA,
    unsigned long int QA,
    unsigned long int QQA,
    unsigned long int PB,
    unsigned long int QB,
    unsigned long int QQB,
    const unsigned long int S,
    const mpz_t D,
    const mpz_t N
);
//-----------------------------------------------------------------------------
inline static unsigned long int
doubling_algo(
    unsigned long int* const P2,
    unsigned long int* const Q2,
    unsigned long int* const QQ2,
    unsigned long int P,
    unsigned long int Q,
    unsigned long int QQ,
    const unsigned long int S,
    const mpz_t D,
    const mpz_t N
);
//-----------------------------------------------------------------------------
inline static unsigned long int
jump_over_k_steps(
    unsigned long int* const P,
    unsigned long int* const Q,
    unsigned long int* const QQ,
    unsigned long int P0,
    unsigned long int Q0,
    unsigned long int QQ0,
    unsigned long int S,
    const mpz_t D,
    unsigned long int K,
    const mpz_t N
);
//-----------------------------------------------------------------------------
inline static unsigned long int
jump_near_sym_point(
    unsigned long int* const P,
    unsigned long int* const Q,
    unsigned long int* const QQ,
    unsigned long int P0,
    unsigned long int Q0,
    unsigned long int QQ0,
    unsigned long int S,
    const mpz_t D,
    unsigned long int K,
    const mpz_t N
);
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                           NON PUBLIC DATA
//-----------------------------------------------------------------------------

//
// List of multipliers used in case of a multiplier race. These are the 15
// squarefree divisors of 1155 = 3 * 5 * 7 * 11.
//
static unsigned int multipliers[NMULTIPLIERS] = {
       3,   5,   7,  11,   15,
      21,  33,  35,  55,   77,
     105, 165, 231, 385, 1155
};

//-----------------------------------------------------------------------------
//                   NON PUBLIC STRUCTURE(S) AND TYPEDEF(S)
//-----------------------------------------------------------------------------

//
// An enumeration of the possible status of a squfof racer.
//
enum racer_status_enum {
    RUNNING,
    STOPPED
};
//-----------------------------------------------------------------------------
typedef enum racer_status_enum racer_status_t;
//-----------------------------------------------------------------------------
//
// Ad-hoc 'pair' structure used for elements placed in a queue.
//
struct struct_ulint_pair_t {
    unsigned long int x;
    unsigned long int y;
};
typedef struct struct_ulint_pair_t ulint_pair_t;
//-----------------------------------------------------------------------------
//
// Ad-hoc structure holding variables needed for each SQUFOF racer
//
struct struct_squfof_racer_t {
    //
    // Multiplier used.
    //
    unsigned long int multiplier;
    //
    // A basic queue implemented as a singly linked list.
    //
    linked_list_t queue;
    //
    // Nodes of the aforementionned linked list, "statically" allocated (so
    // to speak) since we know in advance the maximum number of nodes...
    //
    linked_list_node_t queue_node[MAX_QUEUE_ALLOC];
    //
    // Nodes' data of the linked list...
    //
    ulint_pair_t       queue_data[MAX_QUEUE_ALLOC];
    //
    // SQUFOF's working variables: notations follow the one given in Gowen and
    // Wagstaff's paper.
    //
    mpz_t D;
    mpz_t D_minus_P2;
    mpz_t P2;
    mpf_t sqrtD;

    unsigned long int S;

    unsigned long int P;
    unsigned long int Q;
    unsigned long int QQ;
    unsigned long int PP;

    unsigned long int P0;
    unsigned long int Q0;
    unsigned long int QQ0;

    unsigned long int L;
    unsigned long int B;
    unsigned long int i;

    unsigned long int q;
    unsigned long int t;
    unsigned long int r;

    //
    // Misc. working variables...
    //
    unsigned long int inode;

    linked_list_node_t* current;
    ulint_pair_t* current_pair;

    unsigned long int g;
    uint64_t istep;

    racer_status_t status;
    bool number_is_too_large;
};
//-----------------------------------------------------------------------------
typedef struct struct_squfof_racer_t squfof_racer_t;
//-----------------------------------------------------------------------------
//
// An ad-hoc structure holding the data variables used by the SQUFOF
// implementation.
//
struct struct_squfof_context_t {
    //
    // SQUFOF's parameters. Actually none for the time being...
    //
    squfof_params_t* params;
    //
    // The number to factor.
    //
    mpz_ptr n;
    //
    // Number of racers used (relevant only if a race is performed).
    //
    unsigned long int nracers;
    //
    // Array of SQUFOF racers used (relevant only if a race is performed).
    //
    squfof_racer_t* racers;
    //
    // Whether the context has been updated or not
    //
    bool updated;
    //
    // Whether the number to factor is too large
    //
    bool integer_too_large;
};
//-----------------------------------------------------------------------------
typedef struct struct_squfof_context_t squfof_context_t;
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//               PROTOTYPES OF SQUFOF'S NON PUBLIC FUNCTION(S)
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
static ecode_t init_squfof_context(factoring_machine_t* const machine);
//-----------------------------------------------------------------------------
static ecode_t clear_squfof_context(factoring_machine_t* const machine);
//-----------------------------------------------------------------------------
static ecode_t update_squfof_context(factoring_machine_t* const machine);
//-----------------------------------------------------------------------------
static ecode_t perform_squfof(factoring_machine_t* const machine);
//-----------------------------------------------------------------------------
static ecode_t perform_squfof_no_race(factoring_machine_t* const machine);
//-----------------------------------------------------------------------------
static ecode_t perform_squfof_race(factoring_machine_t* const machine);
//-----------------------------------------------------------------------------
static ecode_t init_squfof_racer(
    squfof_racer_t* const racer,
    const mpz_t n,
    unsigned int multiplier
);
//-----------------------------------------------------------------------------
static ecode_t clear_squfof_racer(squfof_racer_t* const racer);
//-----------------------------------------------------------------------------
static ecode_t cycle_forward(squfof_racer_t* const racer);
//-----------------------------------------------------------------------------
void add_pair_in_queue(
    unsigned long int x,
    unsigned long int y,
    linked_list_t* const queue
);
//-----------------------------------------------------------------------------
#if SQUFOF_DEBUG
void print_queue(linked_list_t* queue);
#endif
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                 "PUBLIC" FUNCTION(S) IMPLEMENTATION
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void set_squfof_params_to_default(squfof_params_t* const params
                                  MAYBE_UNUSED) {
    //
    // While there is no user parameters for this SQUFOF implementation, this
    // function is kept as a placeholder: it might be needed in future code
    // revisions.
    //
    return;
}
//-----------------------------------------------------------------------------
static ecode_t recurse(mpz_array_t* const factors, uint32_array_t* const multis,
                       const mpz_t n, factoring_mode_t mode) {

    return tifa_factor(factors, multis, n, mode);

    //
    // The following code should be used to factor a number using _only_ SQUFOF.
    //
    //squfof_params_t params;
    //set_squfof_params_to_default(&params);
    //return squfof(factors, multis, n, &params, mode);
}
//-----------------------------------------------------------------------------
ecode_t squfof(mpz_array_t* const factors, uint32_array_t* const multis,
               const mpz_t n, const squfof_params_t* const params,
               const factoring_mode_t mode) {

    PRINT_INTRO_MSG("square form");
    INIT_TIMER;
    START_TIMER;

    ecode_t ecode = NO_FACTOR_FOUND;

    factoring_machine_t machine;

    mpz_init_set(machine.n, n);

    machine.mode                = mode;
    machine.params              = (void*) params;
    machine.init_context_func   = init_squfof_context;
    machine.perform_algo_func   = perform_squfof;
    machine.update_context_func = update_squfof_context;
    machine.clear_context_func  = clear_squfof_context;
    machine.recurse_func        = recurse;
    machine.factors             = factors;
    machine.multis              = multis;

    ecode = run_machine(&machine);

    mpz_clear(machine.n);

    STOP_TIMER;
    PRINT_STATUS(machine.success, ecode);

    return ecode;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                  "PRIVATE" FUNCTION(S) IMPLEMENTATION
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                  Functions used by factoring machine
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
static ecode_t init_squfof_context(factoring_machine_t* const machine) {
    //
    // Init the SQUFOF implementation specific variables
    //
    PRINT_INIT_MSG;
    INIT_TIMER;
    START_TIMER;

    mpf_set_default_prec(128);

    machine->context          = (void*) malloc(sizeof(squfof_context_t));
    squfof_context_t* context = (squfof_context_t*) machine->context;

    context->params  = machine->params;
    context->n       = machine->n;
    context->updated = false;
    context->integer_too_large = false;

    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//-----------------------------------------------------------------------------
static ecode_t clear_squfof_context(factoring_machine_t* const machine) {

    PRINT_CLEAN_MSG;
    INIT_TIMER;
    START_TIMER;

    squfof_context_t* context = (squfof_context_t*) machine->context;

    if (context->updated) {
        for (uint32_t i = 0; i < context->nracers; i++) {
            clear_squfof_racer(&context->racers[i]);
        }
        free(context->racers);
    }
    free(context);

    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//-----------------------------------------------------------------------------
static ecode_t update_squfof_context(factoring_machine_t* const machine) {

    INIT_TIMER;
    START_TIMER;

    squfof_context_t* const context = (squfof_context_t*) machine->context;

    if (context->integer_too_large) {
        PRINT_UPDATE_GIVEUP_MSG;
        STOP_TIMER;
        PRINT_TIMING;
        return INTEGER_TOO_LARGE;
    }

    ecode_t ecode = SUCCESS;

    if (!context->updated) {
        //
        // Initialize all the SQUFOF racers needed. For each racer, this will
        // perform the step 1) of the SQUFOF algorithm as described by Gowen and
        // Wagstaff in their paper.
        //
        PRINT_UPDATE_RACE_MSG;

        context->nracers = NMULTIPLIERS;
        context->racers  = malloc(
                               context->nracers * sizeof(squfof_racer_t)
                           );
        for (uint32_t i = 0; i < context->nracers; i++) {
            init_squfof_racer(
                &context->racers[i], context->n, multipliers[i]
            );
        }
        context->updated = true;
        ecode = SUCCESS;
    } else {
        //
        // The context has already been updated and a race has been performed.
        // For the time being, just give up the factorization. It is hard to
        // see what can be done in such case anyway...
        //
        PRINT_UPDATE_GIVEUP_MSG;
        ecode = GIVING_UP;
    }
    STOP_TIMER;
    PRINT_TIMING;

    return ecode;
}
//-----------------------------------------------------------------------------
static ecode_t perform_squfof(factoring_machine_t* const machine) {
    squfof_context_t* const context = (squfof_context_t*) machine->context;
    //
    // Begin by performing SQUFOF without multiplier. This will fail to find
    // a factor approximately 5% of the time. Should it fail and if the
    // factoring mode on the machine is different from SINGLE_RUN, continue
    // the factorization attempt by using a multiplier race.
    //
    if (!context->updated) {
        return perform_squfof_no_race(machine);
    } else {
        return perform_squfof_race(machine);
    }
}
//-----------------------------------------------------------------------------
static ecode_t perform_squfof_no_race(factoring_machine_t* const machine) {

    INIT_TIMER;
    squfof_context_t* context = (squfof_context_t*) machine->context;

    ecode_t ecode = SOME_FACTORS_FOUND;

    mpf_set_default_prec(128);

    //
    // Step 1: Initialization
    //
    mpz_t D;
    mpz_t D_minus_P2;

    mpz_init_set(D, context->n);

    uint8_t n_mod_4 = mpz_get_ui(context->n) & 3;
    if (n_mod_4 == 1) {
        mpz_mul_2exp(D, D, 1);
    }

    //
    // We need a floating point computation here. See next comment.
    //
    mpf_t sqrtD;
    mpf_init(sqrtD);
    mpf_set_z(sqrtD, D);
    mpf_sqrt(sqrtD, sqrtD);

    if (0 == mpf_fits_ulong_p(sqrtD)) {
        mpz_clear(D);
        mpf_clear(sqrtD);
        context->integer_too_large = true;
        return INTEGER_TOO_LARGE;
    }
    unsigned long int S = mpf_get_ui(sqrtD);

    if (S > TIFA_ULONG_MAX_DIVIDED_BY_6) {
        mpz_clear(D);
        mpf_clear(sqrtD);
        context->integer_too_large = true;
        return INTEGER_TOO_LARGE;
    }

    //
    // Set up the linked list used as a queue...
    //
    linked_list_node_t* current = NULL;
    linked_list_t queue;
    //
    // Nodes of the aforementionned linked list, "statically" allocated (so
    // to speak) since we know in advance the maximum number of nodes...
    //
    linked_list_node_t queue_node[MAX_QUEUE_ALLOC];
    //
    // Nodes' data of the linked list...
    //
    ulint_pair_t       queue_data[MAX_QUEUE_ALLOC];

    for (uint32_t i = 0; i < MAX_QUEUE_ALLOC-1; i++) {
        queue_node[i].next = &(queue_node[i+1]);
        queue_node[i].data = (void*)&(queue_data[i]);
    }
    queue_node[MAX_QUEUE_ALLOC-1].next = &queue_node[0];
    queue_node[MAX_QUEUE_ALLOC-1].data = (void*)&queue_data[MAX_QUEUE_ALLOC-1];
    queue.head   = &queue_node[0];
    queue.tail   = queue.head;
    queue.length = 0;

    unsigned long int QQ = 1;
    unsigned long int P  = S;

    mpz_t P2;
    mpz_init_set_ui(P2, P);
    mpz_mul(P2, P2, P2);

    mpz_init_set(D_minus_P2, D);
    mpz_sub(D_minus_P2, D_minus_P2, P2);

    unsigned long int Q  = mpz_get_ui(D_minus_P2);

    mpf_mul_2exp(sqrtD, sqrtD, 1);
    mpf_sqrt(sqrtD, sqrtD);
    mpf_mul_2exp(sqrtD, sqrtD, 1);
    //
    // Now sqrtD = 2*sqrt(2*sqrt(D)).
    //
    // The use of floating point operation is made necessary to properly compute
    // L = floor(2*sqrt(2*sqrt(D))) without introducing approximations that
    // would lead to L = 2*floor(sqrt(2*floor(sqrt(D)))). However, such an
    // approximation would probably have no impact whatsoever...
    //
    unsigned long int L  = mpf_get_ui(sqrtD);
    unsigned long int B  = L << 1;
    unsigned long int i  = 0;

    unsigned long int q  = 0;
    //
    // _NOTE_: An amusing remark: the gmpl-impl.h header already defines
    //         a symbol named PP. So make sure that gmp-impl.h is not included!
    //
    unsigned long int PP = 0;

    unsigned long int t  = 0;
    unsigned long int r  = 0;

    unsigned long int inode = 0;

    ulint_pair_t* current_pair = NULL;

    PRINT_FWD_CYCL_MSG;
    START_TIMER;

    SQUFOF_GMP_PRINT("\nD = %Zd\n", D);
    SQUFOF_PRINT("S = %lu\n", S);
    SQUFOF_PRINT("B = %lu\n", B);

    unsigned long int P0  = P;
    unsigned long int Q0  = Q;
    unsigned long int QQ0 = 1;

    //
    // Step 2: Cycle forward to find a proper square form
    //
    while (true) {
        SQUFOF_PRINT("\n----------------------------\n");
        SQUFOF_PRINT("Step 2: Iteration %lu\n", i);
        SQUFOF_PRINT("----------------------------\n");
        //
        // Step 2 a)
        //
        // Compute q  = (S + P)/Q;
        //
        q = divide(S + P, Q);

        //
        // Compute PP = q*Q - P;
        //
        PP  = q;
        PP *= Q;
        PP -= P;

        SQUFOF_PRINT("step 2a):\n");
        SQUFOF_PRINT("-----------------------------------\n");
        SQUFOF_PRINT("2a)     S  = %lu\n", S);
        SQUFOF_PRINT("2a)     Q  = %lu\n", Q);
        SQUFOF_PRINT("2a)     P  = %lu\n", P);
        SQUFOF_PRINT("2a)     QQ = %lu\n", QQ);
        SQUFOF_PRINT("2a)     PP = %lu\n", PP);
        SQUFOF_PRINT("2a)     t  = %lu\n", t);
        SQUFOF_PRINT("2a)     q  = %lu\n", q);
        //
        // Step 2 b)
        //
        if (Q <= L) {

            if (IS_EVEN(Q)) {
                if (queue.length == MAX_QUEUE_ALLOC) {
                    ecode = QUEUE_OVERFLOW;
                    PRINT_FAILURE_NL;
                    goto clean_and_return;
                }
                add_pair_in_queue(Q >> 1, P % (Q>>1), &queue);
            } else {
                if (Q <= (L >> 1)) {
                    if (queue.length == MAX_QUEUE_ALLOC) {
                        ecode = QUEUE_OVERFLOW;
                        PRINT_FAILURE_NL;
                        goto clean_and_return;
                    }
                    add_pair_in_queue(Q, P % Q, &queue);
                }
            }
        }
        //
        // Step 2 c)
        //
        // Compute t  = QQ + q * (P - PP);
        //
        t  = P;
        t -= PP;
        t *= q;
        t += QQ;

        QQ = Q;
        Q  = t;
        P  = PP;

        SQUFOF_PRINT("step 2c):\n");
        SQUFOF_PRINT("-----------------------------------\n");
        SQUFOF_PRINT("2c)     Q  = %lu\n", Q);
        SQUFOF_PRINT("2c)     P  = %lu\n", P);
        SQUFOF_PRINT("2c)     QQ = %lu\n", QQ);
        SQUFOF_PRINT("2c)     PP = %lu\n", PP);
        SQUFOF_PRINT("2c)     t  = %lu\n", t);
        SQUFOF_PRINT("2c)     q  = %lu\n", q);

        //
        // Step 2 d)
        //
        if (IS_EVEN(i) && ((r = is_square(Q)) != 0)) {

            inode   = 0;
            current = queue.head;

            bool pair_found = false;

            while (inode < queue.length) {

                current_pair = (ulint_pair_t*)(current->data);

                if (current_pair->x == r) {

                    if (r == 1) {
                        ecode = NO_PROPER_FORM_FOUND;
                        PRINT_FAILURE_NL;
                        goto clean_and_return;
                    }
                    if (r > 1) {
                        if ((P - current_pair->y) % current_pair->x == 0) {
                            //
                            // Remove all pairs from the beginning of the queue
                            // up to (and including) this pair
                            //
                            queue.head    = current->next;
                            queue.length -= (inode + 1);
                            pair_found    = true;

                            break;
                        }
                    }
                }
                current = current->next;
                inode++;
            }

            if (pair_found == false) {
                //
                // No pairs (r, t) in the queue for which r divides (P - t):
                // exit the loop and go to Step 3.
                //
                break;
            }
        }
        //
        // Step 2 e)
        //
        i++;
        if (i > B) {
            ecode = GIVING_UP;
            PRINT_FAILURE_NL;
            goto clean_and_return;
        }
    }
    STOP_TIMER;
    PRINT_TIMING;
    PRINT_INV_SQRT_MSG;
    RESET_TIMER;
    START_TIMER;

    //
    // Step 3: Compute an inverse square root of the square form
    //
    QQ = r;
    P += r * (unsigned long int)((S - P)/r);

    mpz_set_ui(P2, P);
    mpz_mul(P2, P2, P2);

    mpz_sub(D_minus_P2, D, P2);
    mpz_divexact_ui(D_minus_P2, D_minus_P2, QQ);

    Q  = mpz_get_ui(D_minus_P2);

    SQUFOF_PRINT("after step 3):\n");
    SQUFOF_PRINT("-----------------------------------\n");
    SQUFOF_PRINT("3)      Q  = %lu\n", Q);
    SQUFOF_PRINT("3)      P  = %lu\n", P);
    SQUFOF_PRINT("3)      QQ = %lu\n", QQ);
    SQUFOF_PRINT("3)      PP = %lu\n", PP);
    SQUFOF_PRINT("3)      S  = %lu\n", S);
    SQUFOF_PRINT("3)      t  = %lu\n", t);
    SQUFOF_PRINT("3)      q  = %lu\n", q);

    STOP_TIMER;
    PRINT_TIMING;
    RESET_TIMER;
    START_TIMER;

    unsigned long int factor = 0;

    if ((i / 2) < LARGE_STEP_THRESHOLD) {
        PRINT_REV_CYCL_MSG;
        //
        // Step 4: Cycle in the reverse direction to find a factor of N
        //
        while (true) {

            q = divide(S + P, Q);

            PP  = q;
            PP *= Q;
            PP -= P;
            if (P == PP) {
                factor = Q;
                break;
            }
            t  = P;
            t -= PP;
            t *= q;
            t += QQ;
            QQ = Q;
            Q  = t;
            P  = PP;

            SQUFOF_PRINT("after iteration of step 4):\n");
            SQUFOF_PRINT("-----------------------------------\n");
            SQUFOF_PRINT("4)      Q  = %lu\n", Q);
            SQUFOF_PRINT("4)      P  = %lu\n", P);
            SQUFOF_PRINT("4)      QQ = %lu\n", QQ);
            SQUFOF_PRINT("4)      PP = %lu\n", PP);
            SQUFOF_PRINT("4)      t  = %lu\n", t);
        }
    } else {
        PRINT_FR_REV_CYCL_MSG;
        //
        // Use the large step algorithm to get closer to the point of symmetry.
        //
    
        factor = jump_near_sym_point(&P, &Q, &QQ, P0, Q0,
                                     QQ0, S, D, i / 2, context->n);

        if (factor == 0) {
            //
            // Step 4: Cycle in _both_ directions to find a factor of N
            //
            unsigned long int PR  = P;
            unsigned long int QR  = QQ;
            unsigned long int QQR = Q;
            unsigned long int PPR = 0;
            unsigned long int qr  = 0;
            
            while (true) {
                //
                // Cycle in the "forward" direction
                //
                q = divide(S + P, Q);
                PP  = q;
                PP *= Q;
                PP -= P;
                if (P == PP) {
                    factor = Q;
                    break;
                }
                //
                // Cycle in the "reverse" direction
                //
                qr = divide(S + PR, QR);
                PPR  = qr;
                PPR *= QR;
                PPR -= PR;
                if (PR == PPR) {
                    factor = QR;
                    break;
                }
                //
                // Complete the "forward" direction cycle
                //
                t  = P;
                t -= PP;
                t *= q;
                t += QQ;
                QQ = Q;
                Q  = t;
                P  = PP;
                //
                // Complete the "reverse" direction cycle
                //
                t  = PR;
                t -= PPR;
                t *= qr;
                t += QQR;
                QQR = QR;
                QR  = t;
                PR  = PPR;
            
                SQUFOF_PRINT("after forward iteration of step 4):\n");
                SQUFOF_PRINT("-----------------------------------\n");
                SQUFOF_PRINT("4)      Q  = %lu\n", Q);
                SQUFOF_PRINT("4)      P  = %lu\n", P);
                SQUFOF_PRINT("4)      QQ = %lu\n", QQ);
                SQUFOF_PRINT("4)      PP = %lu\n", PP);
            
                SQUFOF_PRINT("\tafter reverse iteration of step 4):\n");
                SQUFOF_PRINT("\t-----------------------------------\n");
                SQUFOF_PRINT("\t4)      QR  = %lu\n", QR);
                SQUFOF_PRINT("\t4)      PR  = %lu\n", PR);
                SQUFOF_PRINT("\t4)      QQR = %lu\n", QQR);
                SQUFOF_PRINT("\t4)      PPR = %lu\n", PPR);
            }
        }
    }
    STOP_TIMER;
    PRINT_TIMING;

    //
    // Step 5: We found a factor of n
    //
    if (IS_EVEN(factor)) {
        factor >>= 1;
    }
    if (factor != 1) {
        mpz_t factor_z;
        mpz_init_set_ui(factor_z, factor);
        append_mpz_to_array(machine->factors, factor_z);
        mpz_divexact(factor_z, context->n, factor_z);
        append_mpz_to_array(machine->factors, factor_z);
        mpz_clear(factor_z);
    } else {
        ecode = NO_FACTOR_FOUND;
    }

  clean_and_return:

    mpf_clear(sqrtD);
    mpz_clear(D);
    mpz_clear(D_minus_P2);
    mpz_clear(P2);

    return ecode;
}
//-----------------------------------------------------------------------------
static ecode_t perform_squfof_race(factoring_machine_t* const machine) {

    INIT_TIMER;

    squfof_context_t* const context = (squfof_context_t*) machine->context;

    bool winner_found      = false;
    squfof_racer_t *winner = NULL;

    ecode_t rcode = 0;

    PRINT_FWD_CYCL_MSG;
    START_TIMER;
    //
    // Step 2: cycle forward to find a proper form.
    //
    // This is step 2) as described by Gowen and Wagstaff in their paper. Note
    // that the step 1) has already been performed via the init_squfof_context
    // function...
    //

    //
    // First get the number of racers still racing: indeed some may have stopped
    // during their initialization if the number to factor multiplied by the
    // multiplier is too large for this single-precision implementation.
    //
    unsigned int nrunning_racers = 0;

    for (uint32_t iracer = 0; iracer < context->nracers; iracer++) {
        if (RUNNING == context->racers[iracer].status) {
            nrunning_racers++;
        }
    }

    if (nrunning_racers == 0) {
        //
        // The only way this can happen is if the number given to factorize
        // was too large...
        //
        PRINT_FAILURE_NL;
        return INTEGER_TOO_LARGE;
    }

    while (! winner_found) {
        //
        // Make each racer cycle forward for NSTEP_CYCLE_FWD iterations...
        //
        for (uint32_t iracer = 0; iracer < context->nracers; iracer++) {

            if (RUNNING == context->racers[iracer].status) {

                rcode = cycle_forward(&context->racers[iracer]);

                if (rcode == SUCCESS) {
                    winner       = &context->racers[iracer];
                    winner_found = true;
                    break;
                }
                if (STOPPED == context->racers[iracer].status) {
                    nrunning_racers--;
                }
            }
        }
        if (nrunning_racers == 0) {
            //
            // All the racers have stopped and none of them gave a proper
            // square form: abort the factorization process and admit defeat.
            // Note that we return as exit code the code returned by the last
            // failing racer in the race...
            //
            PRINT_FAILURE_NL;
            return rcode;
        }
    }

    STOP_TIMER;
    PRINT_TIMING;
    PRINT_INV_SQRT_MSG;
    RESET_TIMER;
    START_TIMER;

    //
    // Step 3: Compute an inverse square root of the square form
    //
    winner->QQ = winner->r;
    winner->P += winner->r * (unsigned long int)(
                                 (winner->S - winner->P) /winner->r
                             );

    mpz_set_ui(winner->P2, winner->P);
    mpz_mul(winner->P2, winner->P2, winner->P2);

    mpz_sub(winner->D_minus_P2, winner->D, winner->P2);
    mpz_divexact_ui(winner->D_minus_P2, winner->D_minus_P2, winner->QQ);

    winner->Q = mpz_get_ui(winner->D_minus_P2);

    STOP_TIMER;
    PRINT_TIMING;
    RESET_TIMER;
    START_TIMER;

    unsigned long int factor = 0;
    //
    // Use local variables to bypass readings of a non constant pointer. It is
    // not clear to me if this makes for a measurable difference though...
    //
    unsigned long int S  = winner->S;
    unsigned long int P  = winner->P;
    unsigned long int Q  = winner->Q;
    unsigned long int QQ = winner->QQ;
    unsigned long int PP = 0;
    unsigned long int q  = 0;
    unsigned long int t  = 0;

    if ((winner->i / 2) < LARGE_STEP_THRESHOLD) {
        PRINT_REV_CYCL_MSG;
        //
        // Step 4: Cycle in the reverse direction to find a factor of N
        //
        while (true) {
            q = divide(S + P, Q);
            PP  = q;
            PP *= Q;
            PP -= P;
            if (P == PP) {
                factor = Q;
                break;
            }
            t   = P;
            t  -= PP;
            t  *= q;
            t  += QQ;
            QQ  = Q;
            Q   = t;
            P   = PP;
        }
    } else {
        PRINT_FR_REV_CYCL_MSG;

        unsigned long int P0  = winner->P0;
        unsigned long int Q0  = winner->Q0;
        unsigned long int QQ0 = winner->QQ0;
        //
        // Use the large step algorithm to get closer to the point of symmetry.
        //
        factor = jump_near_sym_point(&P, &Q, &QQ, P0, Q0, QQ0,
                                     S, winner->D, winner->i / 2, context->n);

        if (factor == 0) {
            //
            // Step 4: Cycle in _both_ directions to find a factor of N
            //
            unsigned long int PR  = P;
            unsigned long int QR  = QQ;
            unsigned long int QQR = Q;
            unsigned long int PPR = 0;
            unsigned long int qr  = 0;
            
            while (true) {
                //
                // Cycle in the "forward" direction
                //
                q = divide(S + P, Q);
                PP  = q;
                PP *= Q;
                PP -= P;
                if (P == PP) {
                    factor = Q;
                    break;
                }
                //
                // Cycle in the "reverse" direction
                //
                qr = divide(S + PR, QR);
                PPR  = qr;
                PPR *= QR;
                PPR -= PR;
                if (PR == PPR) {
                    factor = QR;
                    break;
                }
                //
                // Complete the "forward" direction cycle
                //
                t   = P;
                t  -= PP;
                t  *= q;
                t  += QQ;
                QQ  = Q;
                Q   = t;
                P   = PP;
                //
                // Complete the "reverse" direction cycle
                //
                t  = PR;
                t -= PPR;
                t *= qr;
                t += QQR;
                QQR = QR;
                QR  = t;
                PR  = PPR;
            }
        }
    }
    STOP_TIMER;
    PRINT_TIMING;

    //
    // Step 5: We found a factor of n
    //
    factor /= gcd_ulint(factor, 2 * winner->multiplier);
    
    if (factor ==  1) {
        return NO_FACTOR_FOUND;
    }
    
    mpz_t factor_z;
    
    mpz_init_set_ui(factor_z, factor);
    append_mpz_to_array(machine->factors, factor_z);
    mpz_divexact(factor_z, context->n, factor_z);
    append_mpz_to_array(machine->factors, factor_z);
    
    mpz_clear(factor_z);
    
    return SOME_FACTORS_FOUND;
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *              Step 2) of the SQUFOF algorithm: forward cycle
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
static ecode_t cycle_forward(squfof_racer_t* const racer) {
    //
    // Cycle forward the given racer to find a proper form, but stop after
    // NSTEP_CYCLE_FWD iterations.
    //
    if (racer->status == STOPPED) {
        return FAILURE;
    }
    racer->istep = 0;
    
    while (true) {
        racer->istep++;

        //
        // Step 2 a)
        //
        // Compute q  = (S + P)/Q;
        //
        racer->q = divide(racer->S + racer->P, racer->Q);
        //
        // Compute PP = q*Q - P;
        //
        racer->PP  = racer->q;
        racer->PP *= racer->Q;
        racer->PP -= racer->P;

        //
        // Step 2 b)
        //
        racer->g = racer->Q / gcd_ulint(racer->Q, 2 * racer->multiplier);

        if (racer->g <= racer->L) {

            if (racer->queue.length == MAX_QUEUE_ALLOC) {
                racer->status = STOPPED;
                return QUEUE_OVERFLOW;
            }
            add_pair_in_queue(racer->g, racer->P % racer->g, &(racer->queue));
        }
        //
        // Step 2 c)
        //
        // Compute t  = QQ + q * (P - PP);
        //
        racer->t   = racer->P;
        racer->t  -= racer->PP;
        racer->t  *= racer->q;
        racer->t  += racer->QQ;
        racer->QQ  = racer->Q;
        racer->Q   = racer->t;
        racer->P   = racer->PP;

        //
        // Step 2 d)
        //
        if (IS_EVEN(racer->i) && ((racer->r = is_square(racer->Q)) != 0)) {

            racer->inode   = 0;
            racer->current = racer->queue.head;

            bool pair_found = false;

            while (racer->inode < racer->queue.length) {

                racer->current_pair = (ulint_pair_t*)(racer->current->data);

                if (racer->current_pair->x == racer->r) {
                    if (racer->r == 1) {
                        racer->status = STOPPED;
                        return NO_PROPER_FORM_FOUND;
                    }
                    if (racer->r > 1) {
                        if ( (racer->P - racer->current_pair->y)
                            % racer->current_pair->x == 0) {
                            //
                            // Remove all pairs from the beginning of the queue
                            // up to (and including) this pair
                            //
                            racer->queue.head    = racer->current->next;
                            racer->queue.length -= (racer->inode + 1);
                            pair_found           = true;
                            break;
                        }
                    }
                }
                racer->current = racer->current->next;
                racer->inode++;
            }
            if (pair_found == false) {
                //
                // No pairs (r, t) in the queue for which r divides (P - t):
                // exit the loop and go to Step 3 (i.e. we have a winner!).
                //
                break;
            }
        }
        //
        // Step 2 e)
        //
        racer->i++;
        if (racer->i > racer->B) {
            //
            // Giving up: stop the race as the abort limit has been reached
            //
            racer->status = STOPPED;
            return GIVING_UP;
        }
        if (racer->istep == NSTEP_CYCLE_FWD) {
            //
            // The racer failed to find a proper form for now, but keep this
            // racer in the race...
            //
            return FAILURE;
        }
    }
    //
    // We found a proper form. We have our winning racer!
    //
    return SUCCESS;
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *              SQUFOF racer's initialization and clean-up
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
static ecode_t init_squfof_racer(squfof_racer_t* const racer, const mpz_t n,
                                 unsigned int multiplier) {
    //
    // Initialize the SQUFOF 'racer' and execute the Step 1) of the SQUFOF
    // algorithm as described by Gowen and Wagstaff in their paper "Square Form
    // Factorization".
    //

    //
    // Step 1: Initialization
    //
    mpz_init_set(racer->D, n);

    racer->multiplier = multiplier;

    mpz_mul_ui(racer->D, racer->D, racer->multiplier);

    uint8_t n_mod_4 = mpz_get_ui(racer->D) & 3;
    if (n_mod_4 == 1) {
        mpz_mul_2exp(racer->D, racer->D, 1);
    }
    //
    // We need a floating point computation here. See next comment.
    //
    mpf_init(racer->sqrtD);
    mpf_set_z(racer->sqrtD, racer->D);
    mpf_sqrt(racer->sqrtD, racer->sqrtD);

    if (0 == mpf_fits_ulong_p(racer->sqrtD)) {
        racer->status = STOPPED;
        racer->number_is_too_large = true;
        return SUCCESS;
    }
    racer->S = mpf_get_ui(racer->sqrtD);

    if (racer->S > TIFA_ULONG_MAX_DIVIDED_BY_6) {
        racer->status = STOPPED;
        racer->number_is_too_large = true;
        return SUCCESS;
    }

    racer->status = RUNNING;
    racer->number_is_too_large = false;

    racer->QQ = 1;
    racer->P  = racer->S;

    mpz_init_set_ui(racer->P2, racer->P);
    mpz_mul(racer->P2, racer->P2, racer->P2);

    mpz_init_set(racer->D_minus_P2, racer->D);
    mpz_sub(racer->D_minus_P2, racer->D_minus_P2, racer->P2);

    racer->Q = mpz_get_ui(racer->D_minus_P2);

    mpf_mul_2exp(racer->sqrtD, racer->sqrtD, 1);
    mpf_sqrt(racer->sqrtD, racer->sqrtD);
    mpf_mul_2exp(racer->sqrtD, racer->sqrtD, 1);
    //
    // Now sqrtD = 2*sqrt(2*sqrt(D)).
    //
    // The use of floating point operation is made necessary to properly compute
    // L = floor(2*sqrt(2*sqrt(D))) without introducing approximations that
    // would lead to L = 2*floor(sqrt(2*floor(sqrt(D)))).
    //
    racer->L  = mpf_get_ui(racer->sqrtD);
    racer->B  = racer->L << 1;
    racer->i  = 0;
    racer->q  = 0;
    //
    // _NOTE_: An amusing remark: the gmpl-impl.h header already defines a
    //         symbol named PP. So make sure that gmp-impl.h is not included!
    //
    racer->PP = 0;
    racer->t  = 0;
    racer->r  = 0;

    racer->P0  = racer->P;
    racer->Q0  = racer->Q;
    racer->QQ0 = 1;

    //
    // Set up the linked list used as the racer's queue and make it circular...
    //
    for (uint32_t i = 0; i < MAX_QUEUE_ALLOC-1; i++) {
        racer->queue_node[i].next = &(racer->queue_node[i+1]);
        racer->queue_node[i].data = (void*)&(racer->queue_data[i]);
    }
    racer->queue_node[MAX_QUEUE_ALLOC-1].next = &(racer->queue_node[0]);
    racer->queue_node[MAX_QUEUE_ALLOC-1].data =
                                (void*)&(racer->queue_data[MAX_QUEUE_ALLOC-1]);
    racer->queue.head   = &(racer->queue_node[0]);
    racer->queue.tail   = racer->queue.head;
    racer->queue.length = 0;

    racer->inode        = 0;
    racer->current_pair = NULL;
    racer->current      = NULL;

    racer->g     = 0;
    racer->istep = 0;

    return SUCCESS;
}
//-----------------------------------------------------------------------------
static ecode_t clear_squfof_racer(squfof_racer_t* const racer) {
    mpz_clear(racer->D);
    mpz_clear(racer->D_minus_P2);
    if (false == racer->number_is_too_large) {
        mpz_clear(racer->P2);
        mpf_clear(racer->sqrtD);
    }
    return SUCCESS;
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                      Queue management functions
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
void add_pair_in_queue(unsigned long int x, unsigned long int y,
                       linked_list_t* const queue) {

    ulint_pair_t* current_pair = NULL;
    if (queue->length == 0) {
        queue->tail = queue->head;
    } else {
        queue->tail = queue->tail->next;
    }
    current_pair = (ulint_pair_t*)(queue->tail->data);
    current_pair->x = x;
    current_pair->y = y;

    queue->length++;
}
//-----------------------------------------------------------------------------
#if SQUFOF_DEBUG
void print_queue(linked_list_t* queue) {
    //
    // Only used for debugging purpose
    //
    linked_list_node_t* current_node = queue->head;
    ulint_pair_t*       current_pair = NULL;

    uint32_t inode = 0;
    printf("\t-------- Queue --------\n");
    while (inode < queue->length) {
        current_pair = (ulint_pair_t*)(current_node->data);
        printf("\t(%ld, %ld)\n", current_pair->x, current_pair->y);
        current_node = current_node->next;
        inode++;
    }
    printf("\t-----------------------\n");
}
#endif
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                 Functions used for our cheap division trick
 *-----------------------------------------------------------------------------
 */

inline static unsigned long int divide_trick(const unsigned long int N,
                                             const unsigned long int D) {
    //
    // A really quick and dirty trick to (hopefully) speed up divisions if
    // we know with a relatively high probability that the result will be one
    // or two _and_ that the result *won't be zero*. In our case, the quotient
    // distribution has been experimentally determined to be roughly:
    //
    //      Quotient Percentage
    //          1      41.50 %
    //          2      17.00 %
    //          3       9.30 %
    //          4       5.88 %
    //          5       4.06 %
    //          6       2.97 %
    //          7       2.27 %
    //          8       1.79 %
    //          9       1.44 %
    //         10       1.19 %
    //
    unsigned long int __D__ = D << 1;
    if (N < __D__) {
        return 1;
    } else {
        __D__ += D;
        if (N < __D__) {
            return 2;
        } else {
            __D__ += D;
            if (N < __D__) {
                return 3;
            } else {
                __D__ += D;
                if (N < __D__) {
                    return 4;
                } else {
                    __D__ += D;
                    if (N < __D__) {
                        return 5;
                    } else {
                        __D__ += D;
                        if (N < __D__) {
                            return 6;
                        } else {
                            __D__ += D;
                            if (N < __D__) {
                                return 7;
                            } else {
                                __D__ += D;
                                if (N < __D__) {
                                    return 8;
                                } else {
                                    __D__ += D;
                                    if (N < __D__) {
                                        return 9;
                                    } else {
                                        __D__ += D;
                                        if (N < __D__) {
                                            return 10;
                                        } else {
                                            return (N / D);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                 Functions used for the large step algorithm
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
inline static void single_step(unsigned long int* const P,
                               unsigned long int* const Q,
                               unsigned long int* const QQ,
                               const unsigned long int S) {
    //
    // Performs the single step algorithm. This is actually just a step in
    // a continued fraction expansion and is used to step from a binary
    // quadratic form to the next adjacent one (notations follow Gower's & 
    // Wagstaff's).
    //
    unsigned long int  q  = 0;
    unsigned long int  PP = 0;
    unsigned long int  t  = 0;

    q = divide(S + *P, *Q);
    PP   =  q;
    PP  *= *Q;
    PP  -= *P;
    t    = *P;
    t   -=  PP;
    t   *=  q;
    t   += *QQ;
    *QQ  = *Q;
    *Q   =  t;
    *P   =  PP;
}
//-----------------------------------------------------------------------------
inline static unsigned long int large_step(unsigned long int* const P,
                                           unsigned long int* const Q,
                                           unsigned long int* const QQ,
                                           unsigned long int PA,
                                           unsigned long int QA,
                                           unsigned long int QQA,
                                           unsigned long int PB,
                                           unsigned long int QB,
                                           unsigned long int QQB,
                                           const unsigned long int S,
                                           const mpz_t D,
                                           const mpz_t N) {
    //
    // Performs the large step algorithm, "composing" the form given by
    // (-QQA, PA, QA) with the one given by (-QQB, PB, QB) (or a nearby form)
    // and stores the resulting form in (-QQ, P, Q).
    // Notations mostly follow Williams & Wunderlich's with S = sqrt(D).
    //
    unsigned long int gcd = 0;
    unsigned long int qk  = 0;
    unsigned long int PP  = 0;
    
    mpz_t P0;
    mpz_t Q0;
    mpz_t qa_z;
    mpz_t qb_z;
    mpz_t M;
    mpz_t u;
    mpz_t rem;
    mpz_t quot;
    
    mpz_init(Q0);
    mpz_init(P0);
    mpz_init(qa_z);
    mpz_init(qb_z);
    mpz_init(M);
    mpz_init(u);
    mpz_init(rem);
    mpz_init(quot);

    unsigned long int factor = 0;

    gcd = gcd_ulint(QA, QB);
    while (gcd != 1) {
        single_step(&PB, &QB, &QQB, S);
        gcd = gcd_ulint(QA, QB);
    }
    //
    // Step 1: Compute Q0 = QA * QB and its inverse mod D
    //
    mpz_set_ui(Q0, QA);
    mpz_mul_ui(Q0, Q0, QB);
    //
    // We require Q0 to be inversible in Z/DZ and gcd(QA, QB) = 1,
    // so we cycle though the nearby forms until we find a suitable
    // Q pairs...
    //
    while (mpz_invert(M, Q0, D) == 0) {
        mpz_gcd(M, Q0, N);
        //
        // _NOTE_: If the number to factor is not squarefree the current loop
        //         will (most of the time, but not always!) never terminates!
        //         (the gcd M will give a factor of n or k times that factor
        //         where k is the multiplier used: this has been determined
        //         empirically, theorists --or anyone using their brain, unlike
        //         yours truly--, please correct me!). The good news is that we
        //         can escape from that loop and extract that factor...
        //
        if (mpz_cmp_ui(M, 1) != 0) {
            factor = mpz_get_ui(M);
            goto clean_and_return;
        }
        single_step(&PA, &QA, &QQA, S);
        single_step(&PB, &QB, &QQB, S);
        gcd = gcd_ulint(QA, QB);
        while (gcd != 1) {
            single_step(&PB, &QB, &QQB, S);
            gcd = gcd_ulint(QA, QB);
        }
        mpz_set_ui(Q0, QA);
        mpz_mul_ui(Q0, Q0, QB);
    }
    //
    // Solve X*QA - Y*QB = PA - PB for X
    //
    mpz_set_ui(qa_z, QA);
    mpz_set_ui(qb_z, QB);
    mpz_gcdext(M, u, NULL, qa_z, qb_z);
    mpz_set_ui(P0, PB);
    mpz_sub_ui(P0, P0, PA);
    mpz_mul(P0, u, P0);
    mpz_mul_ui(P0, P0, QA);
    mpz_add_ui(P0, P0, PA);
    mpz_mod(P0, P0, Q0);
    //
    // Step 2: Find Q and P' congruent to P mod Q
    //
    mpz_add_ui(quot, P0, S);
    mpz_fdiv_qr(quot, rem, quot, Q0);
    while (true) {
        mpz_ui_sub(P0, S, rem);
        mpz_set(M, D);
        mpz_submul(M, P0, P0);
        mpz_divexact(Q0, M, Q0);
        mpz_add_ui(quot, P0, S);
        mpz_fdiv_qr(quot, rem, quot, Q0);

        if ((mpz_sgn(Q0) > 0) && (mpz_cmp_ui(Q0, S) <= 0)){
            break;
        }
    }
    //
    // Step 3: Keep Q and P mod Q
    //
    *Q = mpz_get_ui(Q0);
    //
    // Remark: P0 may be negative and/or not fit in a unsigned long
    // int, so the mpz_mod is necessary...
    //
    mpz_mod_ui(P0, P0, *Q);
    *P = mpz_get_ui(P0);
    //
    // While not stricly needed, we compute the actual value of P
    // (and not just P mod Q) together with QQ (which is the previous
    // or the next value of Q depending on whether or not we've gone
    // beyond the period of the form cycle)
    //
    qk = divide(S + *P, *Q);
    PP = qk * (*Q) - *P;
    if (S >= *P) {
        //
        // Remember that the 'division_trick' function potentially used by
        // the 'divide' macro doesn't work if the result is zero, so we
        // handle this case separately...
        //
        if ((S - *P) >= *Q) {
            unsigned long int tmp = 0;
            tmp = divide(S - *P, *Q);
            qk += tmp;
        }
    } else {
        if ((*P - S) >= *Q) {
            unsigned long int tmp = 0;
            tmp = divide(*P - S, *Q);
            qk -= tmp;
            if ((*P - S) != tmp * (*Q) ) {
                qk--;
            }
        } else {
            if (*P != S) {
                qk--;
            }
        }
    }
    *P = qk * (*Q) - PP;

    mpz_set_ui(M, *P);
    mpz_mul_ui(M, M, *P);
    mpz_sub(M, D, M);
    mpz_fdiv_q_ui(M, M, *Q);

    *QQ = mpz_get_ui(M);

  clean_and_return:

    mpz_clear(Q0);
    mpz_clear(P0);
    mpz_clear(qa_z);
    mpz_clear(qb_z);
    mpz_clear(u);
    mpz_clear(M);
    mpz_clear(rem);
    mpz_clear(quot);

    return factor;
}
//-----------------------------------------------------------------------------
inline static unsigned long int doubling_algo(unsigned long int* const P2,
                                              unsigned long int* const Q2,
                                              unsigned long int* const QQ2,
                                              unsigned long int P,
                                              unsigned long int Q,
                                              unsigned long int QQ,
                                              const unsigned long int S,
                                              const mpz_t D,
                                              const mpz_t N) {
    //
    // Performs the doubling algorithm, composing the form given by (P, Q, QQ)
    // with the nearest suitable form and stores the result in (P2, Q2, QQ2).
    //
    unsigned long int PB  = P;
    unsigned long int QB  = Q;
    unsigned long int QQB = QQ;
    return large_step(P2, Q2, QQ2, P, Q, QQ, PB, QB, QQB, S, D, N);
}
//-----------------------------------------------------------------------------
inline static unsigned long int jump_over_k_steps(unsigned long int* const P,
                                                  unsigned long int* const Q,
                                                  unsigned long int* const QQ,
                                                  unsigned long int P0,
                                                  unsigned long int Q0,
                                                  unsigned long int QQ0,
                                                  unsigned long int S,
                                                  const mpz_t D,
                                                  unsigned long int K,
                                                  const mpz_t N) {
    //
    // Starting with the principal form (1, P0, Q0), jump over (approximately!)
    // K forms to get the form (QQ, P, Q) using Williams & Wunderlich's large 
    // step algorithm.
    //
    unsigned long int factor = 0;
    if (K <= MIN_NSTEP_LS) {
        //
        // Just perform K single steps sequentially since we're under the
        // threshold for large step use.
        //
        *P  = P0;
        *Q  = Q0;
        *QQ = QQ0;
        for (uint32_t i = 0; i < K; i++) {
            single_step(P, Q, QQ, S);
        }
    } else {
        unsigned long int PMIN  = P0;
        unsigned long int QMIN  = Q0;
        unsigned long int QQMIN = QQ0;

        for (uint32_t i = 0; i < MIN_NSTEP_LS; i++) {
            single_step(&PMIN, &QMIN, &QQMIN, S);
        }
        uint32_t n = K / MIN_NSTEP_LS;

        *P  = PMIN;
        *Q  = QMIN;
        *QQ = QQMIN;

        uint32_t f = most_significant_bit(n);

        while (f != 0) {
            f--;
            factor = doubling_algo(P, Q, QQ, *P, *Q, *QQ, S, D, N);
            if (factor != 0) {
                return factor;
            }
            if (BIT(n, f) == 1) {
                factor = large_step(P, Q, QQ, *P, *Q, *QQ,
                                    PMIN, QMIN, QQMIN, S, D, N);
                if (factor != 0) {
                    return factor;
                }
            }
        }
    }
    return 0;
}
//-----------------------------------------------------------------------------
inline static unsigned long int jump_near_sym_point(unsigned long int* const P,
                                                    unsigned long int* const Q,
                                                    unsigned long int* const QQ,
                                                    unsigned long int P0,
                                                    unsigned long int Q0,
                                                    unsigned long int QQ0,
                                                    unsigned long int S,
                                                    const mpz_t D,
                                                    unsigned long int K,
                                                    const mpz_t N) {
    //
    // Given the principal form (1, P0, Q0) and the inverse square root of a
    // found square form given by (QQ, P, Q), jump through K forms starting
    // from (QQ, P, Q) to get in the vicinity of a point of symmetry of the
    // form cycle of (QQ, P, Q). The resulting form is stored in (QQ, P, Q).
    //
    unsigned long int PH  = 0;
    unsigned long int QH  = 0;
    unsigned long int QQH = 0;
    
    unsigned long int factor = 0;
    
    factor = jump_over_k_steps(&PH, &QH, &QQH, P0, Q0, QQ0, S, D, K, N);
    if (factor != 0) {
        return factor;
    }
    factor = large_step(P, Q, QQ, *P, *Q, *QQ, PH, QH, QQH, S, D, N);
    if (factor != 0) {
        return factor;
    }
    return 0;
}
//-----------------------------------------------------------------------------
