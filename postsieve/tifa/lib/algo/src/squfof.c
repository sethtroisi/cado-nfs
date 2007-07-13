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
 * \file    squfof.c
 * \author  Jerome Milan
 * \date    Tue Jun 27 2007
 * \version 1.2.1
 */

 /*
  * History:
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

#define SQUFOF_DEBUG 0
#if SQUFOF_DEBUG
  #define SQUFOF_PRINT(...) printf(__VA_ARGS__)
#else
  #define SQUFOF_PRINT(...) /* intentionally left empty */
#endif

  //------------------------------------------------------------------------
  // Ref. "Square Form Factorization",
  //       Jason E. Gowen, Samuel S. Wagstaff Jr.,
  //       accepted for publication in Mathematics of Computation
  //------------------------------------------------------------------------

#include "tifa_config.h"

#include <stdlib.h>
#include <stdbool.h>
#include <inttypes.h>
#include <gmp.h>

#include "squfof.h"
#include "funcs.h"
#include "macros.h"
#include "linked_list.h"

#define __PREFIX__  "squfof: "
#define __VERBOSE__ TIFA_VERBOSE_SQUFOF
#define __TIMING__  TIFA_TIMING_SQUFOF

#include "messages.h"

//-----------------------------------------------------------------------------
//                         NON PUBLIC DEFINE(S)
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
//
// A really quick and dirty trick to (hopefully) speed up divisions if
// we know with a relatively high probability that the result will be one or 
// two _and_ that the result won't be zero. In our case, the quotient 
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
// First, a macro to declare needed local variables...
//
#define INIT_DIVIDE_TRICK_VARS      \
    unsigned long int __D__ = 0;    \
    unsigned long int __N__ = 0
//
// Then the division macro strictly speaking...
//
#define DIVIDE_TRICK(Q, N, D)                                   \
    __D__ = (D) << 1;                                           \
    __N__ = (N);                                                \
    if (__N__ < __D__) {                                        \
        Q = 1;                                                  \
    } else {                                                    \
        __D__ += (D);                                           \
        if (__N__ < __D__) {                                    \
            Q = 2;                                              \
        } else {                                                \
            __D__ += (D);                                       \
            if (__N__ < __D__) {                                \
                Q = 3;                                          \
            } else {                                            \
                __D__ += (D);                                   \
                if (__N__ < __D__) {                            \
                    Q = 4;                                      \
                } else {                                        \
                    __D__ += (D);                               \
                    if (__N__ < __D__) {                        \
                        Q = 5;                                  \
                    } else {                                    \
                        __D__ += (D);                           \
                        if (__N__ < __D__) {                    \
                            Q = 6;                              \
                        } else {                                \
                            __D__ += (D);                       \
                            if (__N__ < __D__) {                \
                                Q = 7;                          \
                            } else {                            \
                                __D__ += (D);                   \
                                if (__N__ < __D__) {            \
                                    Q = 8;                      \
                                } else {                        \
                                    __D__ += (D);               \
                                    if (__N__ < __D__) {        \
                                        Q = 9;                  \
                                    } else {                    \
                                        __D__ += (D);           \
                                        if (__N__ < __D__) {    \
                                            Q = 10;             \
                                        } else {                \
                                            Q = __N__ / (D);    \
                                        }                       \
                                    }                           \
                                }                               \
                            }                                   \
                        }                                       \
                    }                                           \
                }                                               \
            }                                                   \
        }                                                       \
    }
    
//
// Finally a switch toggling this cheap trick on and off
//
#define USE_DIVISION_TRICK 1
//
// Switch from the standard division to our cheap trick here, to avoid
// disseminating #if's in the code...
//
#if USE_DIVISION_TRICK
    #define DIVIDE(Q, N, D) DIVIDE_TRICK(Q, N, D)
    #define INIT_DIVIDE_VARS INIT_DIVIDE_TRICK_VARS
#else
    #define DIVIDE(Q, N, D) Q = ((N) / (D))
    #define INIT_DIVIDE_VARS /* intentionally left empty */
#endif
   
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
    unsigned long int QQ;
    unsigned long int P;

    unsigned long int Q;
    unsigned long int L;
    unsigned long int B;
    unsigned long int i;

    unsigned long int q;
    unsigned long int PP;

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
//                 PROTOTYPES OF NON PUBLIC FUNCTION(S)
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
                                  __attribute__ ((unused))) {
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

    squfof_params_t params;
    set_squfof_params_to_default(&params);

    return squfof(factors, multis, n, &params, mode);
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
    
    if (S > TIFA_ULONG_MAX_DIVIDED_BY_3) {
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

    INIT_DIVIDE_VARS;

    PRINT_FWD_CYCL_MSG;
    START_TIMER;

    SQUFOF_PRINT("S = %lu\n", S);
    SQUFOF_PRINT("B = %lu\n", B);

    //
    // Step 2: Cycle forward to find a proper square form
    //
    unsigned long int istep = 0;
    while (true) {
        istep++;

        SQUFOF_PRINT("\n----------------------------\n");
        SQUFOF_PRINT("Step 2: Iteration %lu\n", istep);
        SQUFOF_PRINT("----------------------------\n");

        //
        // Step 2 a)
        //
        // Compute q  = (S + P)/Q;
        //
        DIVIDE(q, S + P, Q);
        
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

    printf("after step 3):\n");
    printf("-----------------------------------\n");
    printf("3)      Q  = %lu\n", Q);
    printf("3)      P  = %lu\n", P);
    printf("3)      QQ = %lu\n", QQ);
    printf("3)      PP = %lu\n", PP);
    printf("3)      S  = %lu\n", S);
    printf("3)      t  = %lu\n", t);
    printf("3)      q  = %lu\n", q);

    STOP_TIMER;
    PRINT_TIMING;
    PRINT_REV_CYCL_MSG;
    START_TIMER;
    START_TIMER;

    //
    // Step 4: Cycle in the reverse direction to find a factor of N
    //
    while (true) {
        DIVIDE(q, S + P, Q);
        
        PP  = q;
        PP *= Q;
        PP -= P;
        if (P == PP) {
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

    STOP_TIMER;
    PRINT_TIMING;

    //
    // Step 5: We found a factor of n
    //
    if (IS_EVEN(Q)) {
        Q >>= 1;
    }
    unsigned long int factor = Q;

    mpz_init_set_ui(machine->factors->data[machine->factors->length], factor);
    machine->factors->length++;

    mpz_init_set_ui(machine->factors->data[machine->factors->length], factor);
    mpz_fdiv_ui(machine->factors->data[machine->factors->length], factor);
    machine->factors->length++;

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

    INIT_DIVIDE_VARS;

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
    PRINT_REV_CYCL_MSG;
    START_TIMER;
    START_TIMER;

    //
    // Step 4: Cycle in the reverse direction to find a factor of N
    //
    while (true) {    
        DIVIDE(winner->q, winner->S + winner->P, winner->Q);
        
        winner->PP  = winner->q;
        winner->PP *= winner->Q;
        winner->PP -= winner->P;

        if (winner->P == winner->PP) {
            break;
        }
        winner->t   = winner->P;
        winner->t  -= winner->PP;
        winner->t  *= winner->q;
        winner->t  += winner->QQ;
        winner->QQ  = winner->Q;
        winner->Q   = winner->t;
        winner->P   = winner->PP;
    }

    STOP_TIMER;
    PRINT_TIMING;

    //
    // Step 5: We found a factor of n
    //
    if (IS_EVEN(winner->Q)) {
        winner->Q >>= 1;
    }

    mpz_t factor;

    mpz_init_set_ui(
        factor,
        winner->Q / gcd_ulint(winner->Q, 2 * winner->multiplier)
    );
    append_mpz_to_array(machine->factors, factor);

    mpz_divexact(factor, context->n, factor);
    append_mpz_to_array(machine->factors, factor);

    mpz_clear(factor);

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
    INIT_DIVIDE_VARS;
    while (true) {
        racer->istep++;
        //
        // Step 2 a)
        //
        // Compute q  = (S + P)/Q;
        //
        DIVIDE(racer->q, racer->S + racer->P, racer->Q);
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
    
    if (racer->S > TIFA_ULONG_MAX_DIVIDED_BY_3) {
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

