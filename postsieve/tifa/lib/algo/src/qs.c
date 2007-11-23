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
 * \file    qs.c
 * \author  Jerome Milan
 * \date    Late October / early November 2007
 * \version 1.2
 */

/*
 * History:
 *
 * 1.2: Late October / early November 2007 by JM
 *      - Some performance enhancements (finally!) such as: better size of
 *        sieve chunks, rewritten factor extraction, new multipliers, etc.
 *
 * 1.1: Late March / early April 2007 by JM
 *      - Completely refactored to use a factoring_machine_t.
 *
 * 1.0: Circa July 2006 by JM
 *      - Initial version.
 */

#include "tifa_config.h"

#include <stdlib.h>
#include <inttypes.h>

#if TIFA_VERBOSE_QS || TIFA_TIMING_QS
    #include <stdio.h>
#endif

#if TIFA_USE_CALLOC_MEMSET
    #include <string.h>
#endif

#include <math.h>

#include "first_primes.h"
#include "gmp_utils.h"
#include "funcs.h"
#include "macros.h"
#include "hashtable.h"
#include "x_tree.h"
#include "bernsteinisms.h"
#include "qs.h"

#define __PREFIX__  "qs: "
#define __VERBOSE__ TIFA_VERBOSE_QS
#define __TIMING__  TIFA_TIMING_QS

#include "messages.h"

//-----------------------------------------------------------------------------
//                         NON PUBLIC DEFINE(S)
//-----------------------------------------------------------------------------

//
// Multiplicative correction factor applied to the computed threshold to
// account for the fact that we do not sieve with powers of primes.
//
#define SIEVE_THRESHOLD_MULTIPLIER 0.7
//
// Largest size of a sieve chunk. The optimal value is archecture dependant, so
// this value should be tweaked according to their target processor.
//
#define SIEVE_CHUNK_SIZE           65535
//
// Begin to sieve with the (NPRIMES_TO_SKIP + 1)-th prime in the factor base
// (the so-called small-primes variation).
//
#define NPRIMES_TO_SKIP            10
#if (NPRIMES_TO_SKIP == 0)
    //
    // Should we sieve with the prime 2?
    //
    #define SIEVE_WITH_2           1
#else
    #define SIEVE_WITH_2           0
#endif
//
// In order to not degrade too much the performance of the hashtable used
// for the large prime variation (due to too many collisions), we set its size
// to the size of the factor base multiplied by HTABLE_SIZE_MULTIPLIER.
//
#define HTABLE_SIZE_MULTIPLIER     16
//
// For the positions in the sieve lesser than (or equal to)
// MAX_X_FOR_FULL_LOG_APPROX, compute the sieving criterion (about the logarithm
// of g(x)) for each x.
// For positions greater than MAX_X_FOR_FULL_LOG_APPROX, just use the logarithm
// of g(avx) where avx is the average value of x in the sieve chunk to scan.
//
#define MAX_X_FOR_FULL_LOG_APPROX 2048
//
// Minimum number of candidates for smoothness test
//
#define MIN_SIZE_SMOOTHNESS_BATCH 64
//
// Number of ranges for the size of the number to factor. Each range is
// define by the optimal value of the factor base's size to use.
//
#define NRANGES 18
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                   NON PUBLIC STRUCTURE(S) AND TYPEDEF(S)
//-----------------------------------------------------------------------------

//
// An ad-hoc structure holding the data variables used by the QS implementation.
//
struct struct_qs_context_t {
    //
    // Number to factor.
    //
    mpz_t n;
    //
    // Multiplier to use.
    //
    uint32_t multiplier;
    //
    // Number to factor x multiplier.
    //
    mpz_t kn;
    //
    // b = ceil(sqrt(n)).
    //
    mpz_t b;
    //
    // Product tree of primes in factor base.
    //
    mpz_tree_t *ptree;
    //
    // Number of columns of the constructed matrix.
    //
    uint32_t matrix_ncols;
    //
    // Maximum number of rows of the constructed matrix.
    //
    uint32_t matrix_nrows;
    //
    // Matrix generated after collecting relations
    //
    binary_matrix_t *matrix;
    //
    // Relations to find: (xi+b)^2 = g(xi) (mod n) with g(xi) smooth.
    // These are the g(xi) to check for smoothness together with their
    // corresponding ui = xi+b.
    //
    mpz_array_t *cand_gx;
    mpz_array_t *cand_u;
    //
    // All accepted relations, i.e. the g(xi) and their corresponding ui = xi+b,
    // such that ui^2 = g(xi) (mod n) with g(xi) smooth.
    //
    mpz_array_t *smooth_gx;
    mpz_array_t *u;
    //
    // Factor base used. Only primes pi such that t^2 = N mod (pi) has a
    // solution for t are kept.
    //
    uint32_array_t* factor_base;
    //
    // Logarithms of the elements in the factor base.
    //
    byte_array_t *log_factor_base;
    //
    // Stores sol1i = (ti-b) mod (pi) where ti^2 = N (mod pi) with pi in the
    // factor base.
    //
    uint32_array_t *sol1;
    //
    // Stores sol2i = (-ti-b) mod (pi) where ti^2 = N (mod pi) with pi in the
    // factor base.
    //
    uint32_array_t *sol2;
    //
    // Pool keeping the xi values that has survived the sieving stage, i.e.
    // the xi such that g(xi) _may_ be smooth.
    //
    int32_array_t *xpool;
    //
    // The sieve is implemented as a byte_array_t. Once the sieve is filled
    // and if we need further x values the same array is re-used. The array
    // sieve is consequently only a chunk of the conceptual global sieve. This
    // chunk is given by sieve_begin so that sieve->data[x] really means
    // the (x + sieve_begin)-th position in the complete sieve.
    //
    byte_array_t *sieve;
    //
    // True iff the current sieve chunk has been fully scanned for smooth
    // candidate residues.
    //
    bool chunk_fully_scanned;
    //
    // Hashtable used in large primes variation. It will hold mpz_pair_t's
    // such as (xi+b, g(xi)) with an index in first_primes_array as key.
    //
    hashtable_t *htable;
    //
    // Floor of the logarithm in base 2 of sqrt(n).
    //
    uint32_t log_sqrtn;
    //
    // Number of relations (i.e. of pairs (xi+b, g(xi)) with g(xi) smooth)
    // we want to collect.
    //
    uint32_t npairs_wanted;
    //
    // Number of pairs (xi+b, g(xi)) to collect before proceeding to a
    // smoothness test.
    //
    uint32_t to_collect;
    //
    // First x of the current sieve slice.
    //
    long int sieve_begin;
    //
    // Size (i.e. length) of a sieve slice.
    //
    uint32_t size_sieve;
    //
    // True if we are sieving / scanning in the negative x side.
    //
    bool sieve_negative_side;
    //
    // The position in the sieve slice where we should begin the next scanning
    // phase (that will lead to keep only the xi that satisfy the sieving
    // criterion). The equality (scan_begin == sieve->alloced) means that all
    // positions in the sieve array were checked and that another chunk
    // of the sieve should be computed if more values of x are needed.
    //
    uint32_t scan_begin;
    //
    // Ceil of the logarithm in base 2 of the least prime greater than the
    // largest prime in the factor base.
    //
    uint32_t loglast;
};
//-----------------------------------------------------------------------------
typedef struct struct_qs_context_t qs_context_t;
//-----------------------------------------------------------------------------
struct struct_u32pair_t {
    //
    // If size of number to factor greater than or equal to min_size_n...
    //
    uint32_t min_size_n;
    //
    // use size_base as size of factor base.
    //
    uint32_t size_base;
};
//------------------------------------------------------------------------------
typedef struct struct_u32pair_t u32pair_t;
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                         NON PUBLIC DATA
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
static const u32pair_t optimal_base_sizes[NRANGES] = {
   //
   // These "optimal" values were computed by experimental determination.
   // Actually, they should be regarded more as "good enough" values than
   // "optimal" values as they are still quite crude and not as fine-tuned as
   // CFRAC's.
   //
   {.min_size_n = 0	 , .size_base = 92   },
   {.min_size_n = 65 , .size_base = 128  },
   {.min_size_n = 80 , .size_base = 192  },
   {.min_size_n = 85 , .size_base = 320  },
   {.min_size_n = 90 , .size_base = 448  },
   {.min_size_n = 95 , .size_base = 512  },
   {.min_size_n = 100, .size_base = 512  },
   {.min_size_n = 110, .size_base = 768  },
   {.min_size_n = 120, .size_base = 1024 },
   {.min_size_n = 126, .size_base = 1280 },
   {.min_size_n = 132, .size_base = 1536 },
   {.min_size_n = 150, .size_base = 2048 },
   {.min_size_n = 155, .size_base = 3072 },
   {.min_size_n = 160, .size_base = 4096 }, //      Experimental determination
   {.min_size_n = 165, .size_base = 5120 }, // <--- up to here.
   {.min_size_n = 180, .size_base = 6656 },
   {.min_size_n = 200, .size_base = 8192 },
   {.min_size_n = 220, .size_base = 12288}
};
//------------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                 PROTOTYPES OF NON PUBLIC FUNCTION(S)
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
static ecode_t init_qs_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t clear_qs_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t update_qs_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t perform_qs(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t recurse(
    mpz_array_t* const,
    uint32_array_t* const,
    const mpz_t,
    factoring_mode_t
);
//-----------------------------------------------------------------------------
static void init_startup_data(qs_context_t* const context);
//-----------------------------------------------------------------------------
static void compute_polynomial_values(
    mpz_array_t* const cand_u,
    mpz_array_t* const cand_gx,
    const mpz_t n,
    const mpz_t b,
    const int32_array_t* const x
);
//-----------------------------------------------------------------------------
static void fill_sieve(qs_context_t* const sieve);
//-----------------------------------------------------------------------------
static void scan_sieve(qs_context_t* const context);
//-----------------------------------------------------------------------------
inline uint32_t get_sieve_threshold(const uint32_t log_sqrtn, const uint32_t x);
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                 "PUBLIC" FUNCTION(S) IMPLEMENTATION
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
ecode_t qs(mpz_array_t* const factors, uint32_array_t* const multis,
           const mpz_t n, const qs_params_t* const params,
           const factoring_mode_t mode) {

    PRINT_INTRO_MSG("quadratic sieve");
    INIT_TIMER;
    START_TIMER;

    ecode_t ecode = NO_FACTOR_FOUND;

    factoring_machine_t machine;

    mpz_init_set(machine.n, n);

    machine.mode                = mode;
    machine.params              = (void*) params;
    machine.init_context_func   = init_qs_context;
    machine.perform_algo_func   = perform_qs;
    machine.update_context_func = update_qs_context;
    machine.clear_context_func  = clear_qs_context;
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
void set_qs_params_to_default(const mpz_t n, qs_params_t* const params) {
    //
    // _NOTE_: For the time being these default values are a bit crude.
    //
    // _TO_DO_: Gather more data to provide better parameter values depending
    //          on the size of the number to factor.
    //

    //
    // Size of number to factor = size(n),
    //
    // If:    optimal_base_sizes[i].min_size_n <= size(n)
    // and:   size(n) < optimal_base_sizes[i+1].min_size_n
    // Then:  use optimal_base_sizes[i].size_base as size of factor base.
    //
    uint32_t size_n = mpz_sizeinbase(n, 2);
    uint32_t i      = NRANGES - 1;

    while (size_n < optimal_base_sizes[i].min_size_n) {
        i--;
    }
    params->nprimes_in_base    = optimal_base_sizes[i].size_base;
    params->nprimes_tdiv       = params->nprimes_in_base;
    params->nrelations         = QS_DFLT_NRELATIONS;
    params->linalg_method      = QS_DFLT_LINALG_METHOD;
    params->use_large_primes   = QS_DFLT_USE_LARGE_PRIMES;
}
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                  "PRIVATE" FUNCTION(S) IMPLEMENTATION
//------------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                  Functions used by factoring machine
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
static ecode_t init_qs_context(factoring_machine_t* const machine) {
    //
    // Initialize the QS implementation specific variables...
    //
    PRINT_INIT_MSG;
    INIT_TIMER;
    START_TIMER;

    machine->context      = (void*) malloc(sizeof(qs_context_t));
    qs_context_t* context = (qs_context_t*) machine->context;
    qs_params_t*  params  = (qs_params_t*)  machine->params;

    mpz_init_set(context->n, machine->n);

    //
    // Allocates the needed arrays
    //
    context->factor_base     = alloc_uint32_array(params->nprimes_in_base);
    context->log_factor_base = alloc_byte_array(context->factor_base->alloced);
    context->sol1 = alloc_uint32_array(context->factor_base->alloced);
    context->sol2 = alloc_uint32_array(context->factor_base->alloced);

    //
    // Use Silverman's modified Knuth-Schroeppel function to determine
    // which multiplier to use.
    //
    context->multiplier = ks_multiplier(
                              context->n,
                              context->factor_base->alloced
                          );

    mpz_init(context->kn);
    mpz_mul_ui(context->kn, context->n, context->multiplier);
    //
    // Set b to floor(sqrt(b)) + 1
    //
    mpz_init(context->b);

    mpz_sqrt(context->b, context->kn);
    mpz_add_ui(context->b, context->b, 1);

    //
    // Compute the factor base, its logarithms and the solutions arrays.
    //
    init_startup_data(context);
    //
    // Compute the product tree of all primes in the base once and for all...
    //
    context->ptree = prod_tree_ui(context->factor_base);
    //
    // The number of columns is indeed factor_base->length + 1 since the first
    // row is used to keep track of the sign of the residue.
    //
    context->matrix_ncols  = params->nprimes_in_base + 1;
    context->matrix_nrows  = context->matrix_ncols;
    context->matrix_nrows += params->nrelations;

    context->matrix = alloc_binary_matrix(
                          context->matrix_nrows,
                          context->matrix_ncols
                      );
    context->matrix->nrows = 0;

    context->npairs_wanted = context->matrix_nrows;
    context->to_collect    = context->npairs_wanted;

    context->sieve_negative_side = false;

    context->xpool = alloc_int32_array(context->to_collect);

    context->sieve_begin = 0;
    context->size_sieve  = SIEVE_CHUNK_SIZE;

    context->sieve         = alloc_byte_array(context->size_sieve);
    context->sieve->length = context->sieve->alloced;

    context->cand_gx = alloc_mpz_array(context->to_collect);
    context->cand_u  = alloc_mpz_array(context->to_collect);

    for (uint32_t i = 0; i < context->to_collect; i++) {
        //
        // We fully initialize these arrays since they will be used over
        // and over again. Of course, we should remember to reset their
        // lengths to their alloced values before clearing them!
        //
        mpz_init(context->cand_gx->data[i]);
        mpz_init(context->cand_u->data[i]);
    }
    context->cand_u->length  = context->cand_u->alloced;
    context->cand_gx->length = context->cand_gx->alloced;

    context->smooth_gx = alloc_mpz_array(context->matrix_nrows);
    context->u         = alloc_mpz_array(context->matrix_nrows);

    context->scan_begin = 0;

    context->htable = NULL;
    if (params->use_large_primes) {
        uint32_t size   = HTABLE_SIZE_MULTIPLIER * context->smooth_gx->alloced;
        context->htable = alloc_init_hashtable(
                              size, uint32_cmp_func, hash_rj_32
                          );
    }
    uint32_t lastp   = context->factor_base->data[params->nprimes_in_base - 1];
    context->loglast = ceil_log2(lastp);

    context->log_sqrtn = mpz_sizeinbase(context->n, 2) / 2;
    //
    // Init to true to force the filling of the first sieve chunk.
    //
    context->chunk_fully_scanned = true;

    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//-----------------------------------------------------------------------------
static ecode_t clear_qs_context(factoring_machine_t* const machine) {

    PRINT_CLEAN_MSG;
    INIT_TIMER;
    START_TIMER;

    qs_context_t* context = (qs_context_t*) machine->context;

    if (context != NULL) {
        //
        // Since the folowing arrays were fully initialized, don't forget to
        // set their current length to their alloced length before clearing
        // them...
        //
        context->cand_u->length  = context->cand_u->alloced;
        context->cand_gx->length = context->cand_gx->alloced;
        context->xpool->length   = context->xpool->alloced;

        clear_mpz_array(context->cand_u);
        free(context->cand_u);

        clear_mpz_array(context->u);
        free(context->u);

        clear_mpz_array(context->cand_gx);
        free(context->cand_gx);

        clear_mpz_array(context->smooth_gx);
        free(context->smooth_gx);

        clear_mpz_tree(context->ptree);
        free(context->ptree);

        clear_int32_array(context->xpool);
        free(context->xpool);

        clear_uint32_array(context->factor_base);
        free(context->factor_base);

        clear_byte_array(context->log_factor_base);
        free(context->log_factor_base);

        clear_uint32_array(context->sol1);
        free(context->sol1);

        clear_uint32_array(context->sol2);
        free(context->sol2);

        clear_binary_matrix(context->matrix);
        free(context->matrix);

        clear_byte_array(context->sieve);
        free(context->sieve);

        mpz_clear(context->b);
        mpz_clear(context->n);
        mpz_clear(context->kn);

        if (context->htable != NULL) {
            clear_mpzpair_htable(context->htable);
            free(context->htable);
        }
        free(context);
    }
    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//-----------------------------------------------------------------------------
static ecode_t update_qs_context(factoring_machine_t* const machine) {
    //
    // _TO_DO_: No fallback strategy is implemented yet. We could for example
    //          try to gather more relations...
    //
    qs_context_t* context = (qs_context_t*) machine->context;
    context = NULL;
    return GIVING_UP;
}
//------------------------------------------------------------------------------
static ecode_t perform_qs(factoring_machine_t* const machine) {

    INIT_TIMER;

    qs_context_t* context = (qs_context_t*) machine->context;
    qs_params_t*  params  = (qs_params_t*)  machine->params;

    PRINT_COLLECT_RELS_MSG;
    START_TIMER;

    //
    // Collect relations (xi+b)^2 = g(xi) (mod n) with g(xi) smooth until we
    // have enough to fill the matrix.
    //
    while (context->smooth_gx->length != context->smooth_gx->alloced) {
        //
        // Note that we may collect more than the remaining number of relations
        // needed before attempting a batch smoothness test a la Bernstein.
        //
        // _NOTE_: The value of to_collect can indeed impact performance.
        //         Taking it too small can induce an overhead due to numerous
        //         batches while taking it too big is an overkill that could
        //         lead to useless computation. A careful assessment of the
        //         "optimal" (or rather, of a "not too bad") value must take
        //         into account the probability that a number will indeed be
        //         smooth (based on the other parameters' values... not an easy
        //         task) and cannot rely only on the theoretical complexity of
        //         the batch algorithm used.
        //
        uint32_t remtc = context->npairs_wanted - context->smooth_gx->length;

        context->to_collect = MAX(MIN_SIZE_SMOOTHNESS_BATCH, remtc);

        //
        // We wait until we have filled xpool before attempting to find
        // relations via a batch smoothness test.
        //
        while (context->xpool->length != context->to_collect) {

            if (context->chunk_fully_scanned) {
                //
                // The sieve slice was either completely scanned for candidate
                // xi values or not filled at all. In either case, we need to
                // fill a (new) sieve slice.
                //
                fill_sieve(context);
            }
            //
            // Scan the yet unscanned part of the slice to keep xi values for
            // which g(xi) _may_ be smooth...
            //
            scan_sieve(context);

            if (context->chunk_fully_scanned) {
                if (context->sieve_negative_side) {
                    //
                    // We were sieving for x < 0. Prepare to fill the next
                    // chunk for x > 0.
                    //
                    context->sieve_begin = -(context->sieve_begin);
                    context->scan_begin  = 0;
                    context->sieve_negative_side = false;
                } else {
                    //
                    // We were sieving for x >= 0. Prepare to fill the next
                    // chunk for x < 0.
                    //
                    context->sieve_begin += (int32_t)context->sieve->length;
                    context->sieve_begin  = -context->sieve_begin;
                    context->scan_begin   = context->sieve->length - 1;
                    context->sieve_negative_side = true;
                }
            }
        }
        //
        // We're now ready for a batch smoothness test.
        // First, compute the g(xi)...
        //
        compute_polynomial_values(
            context->cand_u,
            context->cand_gx,
            context->kn,
            context->b,
            context->xpool
        );
        //
        // ... and keep the ones that are smooth...
        //
        if (params->use_large_primes) {
            bern_21_rt_pairs_lp(
                context->n,
                context->htable,
                context->u,
                context->smooth_gx,
                context->cand_u,
                context->cand_gx,
                context->ptree->data[0]
            );
        } else {
            bern_21_rt_pairs(
                context->u,
                context->smooth_gx,
                context->cand_u,
                context->cand_gx,
                context->ptree->data[0]
            );
        }

        PRINT_NRELS_FOUND(
            context->smooth_gx->length,
            context->smooth_gx->alloced
        );

        context->xpool->length   = 0;
    }
    //
    // We have now all the relations we need to fill the matrix and solve the
    // resulting linear system. The rest of the algorithm is now completely
    // standard to all congruence of squares methods.
    //
    STOP_TIMER;
    PRINT_COLLECT_RELS_DONE_MSG;
    PRINT_TIMING;
    PRINT_FACTOR_RES_MSG;
    RESET_TIMER;
    START_TIMER;

    //
    // Cofactors of generated smooth g(xi) after trial divisions.
    //
    mpz_array_t* partial_gx_array = alloc_mpz_array(context->smooth_gx->length);

    uint32_array_t primes_array;
    primes_array.alloced = 0;
    primes_array.length  = params->nprimes_tdiv;
    primes_array.data    = (uint32_t*)(context->factor_base->data);

    byte_matrix_t* decomp_matrix = alloc_byte_matrix(
                                      context->matrix_nrows,
                                      context->matrix_ncols
                                   );

    //
    // Fill the matrix: first by trial divisions...
    //
    fill_trial_div_decomp(
        context->matrix,
        decomp_matrix,
        partial_gx_array,
        context->smooth_gx,
        &primes_array
    );

    if (params->nprimes_tdiv < params->nprimes_in_base) {
        //
        // Finish to fill the matrix with one of Bernstein's batch algo...
        //
        uint32_array_list_t* decomp_list = alloc_uint32_array_list(
                                               context->smooth_gx->length
                                           );
        primes_array.length = params->nprimes_in_base - params->nprimes_tdiv;
        primes_array.data   = (uint32_t*) &(
                                context->factor_base->data[params->nprimes_tdiv]
                              );

        bern_71(decomp_list, partial_gx_array, &primes_array);

        fill_matrix_from_list_decomp(
            context->matrix,
            decomp_matrix,
            partial_gx_array,
            decomp_list,
            context->factor_base
        );

        clear_uint32_array_list(decomp_list);
        free(decomp_list);
    }
    STOP_TIMER;
    PRINT_TIMING;
    PRINT_LIN_ALG_MSG;
    RESET_TIMER;
    START_TIMER;

    //
    // The matrix is now completely filled. Solve the linear system to find
    // relations.
    //
    uint32_array_list_t *relations;
    relations = find_dependencies(context->matrix, params->linalg_method);

    STOP_TIMER;
    PRINT_TIMING;
    PRINT_DED_FACTORS_MSG;
    RESET_TIMER;
    START_TIMER;

    //
    // Deduce factors of n using the found relations.
    //
    ecode_t ecode = find_factors_decomp(
                        machine->factors,
                        context->n,
                        context->u,
                        decomp_matrix,
                        relations,
                        context->factor_base
                    );

    STOP_TIMER;
    PRINT_TIMING;

    clear_byte_matrix(decomp_matrix);
    free(decomp_matrix);

    clear_uint32_array_list(relations);
    free(relations);

    clear_mpz_array(partial_gx_array);
    free(partial_gx_array);

    return ecode;
}
//-----------------------------------------------------------------------------
static ecode_t recurse(mpz_array_t* const factors, uint32_array_t* const multis,
                       const mpz_t n, factoring_mode_t mode) {

    qs_params_t params;
    set_qs_params_to_default(n, &params);

    return qs(factors, multis, n, &params, mode);
}
//------------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                     Sieving utilitarian functions
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
static void fill_sieve(qs_context_t* const context) {
    //
    // Fill the sieve array sieve.
    //
    // The sieve is computed by chunks of size sieve->alloced and the same
    // array sieve is constantly re-used to store the next chunk.
    //
    // The current chunk is given by sieve_begin such that sieve->data[0]
    // really correspond to the (sieve_begin)-th position of the real sieve.
    //
    // _NOTE_: For the time being, we only sieve with primes, not primes
    //         powers. It could be worthwhile to sieve also with small powers
    //         of the few smallest primes in the factor base.
    //
    // _TO_DO_: Also sieve with small powers of the few smallest primes in the
    //          factor base.
    //
    byte_array_t* const sieve                 = context->sieve;
    const long int sieve_begin                = context->sieve_begin;
    const uint32_array_t* const factor_base   = context->factor_base;
    const byte_array_t* const log_factor_base = context->log_factor_base;
    const uint32_array_t* const sol1          = context->sol1;
    const uint32_array_t* const sol2          = context->sol2;

    uint32_t logp     = 0;
    long int imax     = 0;
    long int imin     = 0;
    uint32_t curprime = 0;
    uint32_t cursol   = 0;
    long int sindex   = 0;
    long int smin     = 0;
    long int smax     = 0;

    const long int sieve_end = sieve_begin + sieve->length - 1;

    //
    // (Re-)initialize the sieve array with zeroes...
    //
#if TIFA_USE_CALLOC_MEMSET
    memset(sieve->data, 0x00, sieve->length * sizeof(unsigned char));
#else
    for (uint32_t i = 0; i < sieve->length; i++) {
        sieve->data[i] = 0;
    }
#endif

    //
    // Sieve with the prime 2 only with sol1 and only if SIEVE_WITH_2 != 0.
    //
#if SIEVE_WITH_2
    //
    // imax (respectively imin) is the maximum (resp. minimum) value of i
    // for which (sol1->data[0] + 2*i) >= sieve_begin and
    //           (sol1->data[0] + 2*i) < sieve_begin + sieve->alloced
    //
    cursol = sol1->data[0];

    imax  = sieve_end - cursol;
    imin  = sieve_begin - cursol;
    imax /= 2;
    imin /= 2;

    smin = cursol + (2 * imin) - sieve_begin;
    smax = cursol + (2 * imax) - sieve_begin;

    if (smin < 0) {
        smin += 2;
        imin++;
    }
    if ((uint32_t)smax >= sieve->length) {
        smax -= 2;
        imax--;
    }

    if (smin <= smax) {
        sindex = smin;
        for (int i = imin; i <= imax; i++) {
            sieve->data[sindex] += 1;
            sindex += 2;
        }
    }
#endif

#if (NPRIMES_TO_SKIP == 0)
    uint32_t p = 1;
#else
    uint32_t p = NPRIMES_TO_SKIP;
#endif

    //
    // Sieve with primes > 2 with sol1 and sol2
    //
    for ( ; p < factor_base->length; p++) {

        curprime = factor_base->data[p];
        logp     = log_factor_base->data[p];
        //
        // Sieve with sol1
        //
        // imax (respectively imin) is the maximum (resp. minimum) value of i
        // for which:
        //   (sol1->data[p] + factor_base->data[p]*i) >= sieve_begin
        // and
        //   (sol1->data[p] + factor_base->data[p]*i) < sieve_begin +
        //                                              sieve->alloced
        //
        cursol = sol1->data[p];

        imax  = sieve_end - cursol;
        imin  = sieve_begin - cursol;
        imax /= (long int)curprime;
        imin /= (long int)curprime;

        smin = cursol + (curprime * imin) - sieve_begin;
        smax = cursol + (curprime * imax) - sieve_begin;

        if (smin < 0) {
            smin += curprime;
            imin++;
        }
        if ((uint32_t)smax >= sieve->length) {
            smax -= curprime;
            imax--;
        }
        if (smin <= smax) {
            sindex = smin;
            for (int i = imin; i <= imax; i++) {
                sieve->data[sindex] += logp;
                sindex += curprime;
            }
        }
        //
        // Sieve with sol2
        //
        // imax (respectively imin) is the maximum (resp. minimum) value of i
        // for which:
        //   (sol2->data[p] + factor_base->data[p]*i) >= sieve_begin
        // and
        //   (sol2->data[p] + factor_base->data[p]*i) < sieve_begin +
        //                                              sieve->alloced
        //
        cursol = sol2->data[p];

        imax  = sieve_end - cursol;
        imin  = sieve_begin - cursol;
        imax /= (long int)curprime;
        imin /= (long int)curprime;

        smin = cursol + (curprime * imin) - sieve_begin;
        smax = cursol + (curprime * imax) - sieve_begin;

        if (smin < 0) {
            smin += curprime;
            imin++;
        }
        if ((uint32_t)smax >= sieve->length) {
            smax -= curprime;
            imax--;
        }
        if (smin <= smax) {
            sindex = smin;
            for (int i = imin; i <= imax; i++) {
                sieve->data[sindex] += logp;
                sindex += curprime;
            }
        }
    }
}
//-----------------------------------------------------------------------------
static void scan_sieve(qs_context_t* const context) {
    //
    // _NOTE_: Use local variables to avoid numerous readings of non constant
    //         pointed data...
    //
    int32_array_t* const xpool      = context->xpool;
    int32_t* const xpool_data       = xpool->data;
    uint32_t xpool_length           = xpool->length;
    const byte_array_t* const sieve = context->sieve;
    const long int sieve_begin      = context->sieve_begin;
    uint32_t scan_begin             = context->scan_begin;
    const uint32_t log_sqrtn        = context->log_sqrtn;
    const uint32_t loglast          = context->loglast;
    const uint32_t to_collect       = context->to_collect;

    //
    // Scans the sieve array 'sieve' beginning at the position 'scan_begin'
    // and stores in 'xpool' only the locations x for which:
    //
    //    sieve->data[x] >= 'threshold'
    //
    // 'threshold' is computed in the following way:
    //
    //    - If ABS('scan_begin' + 'sieve_begin') <= MAX_X_FOR_FULL_LOG_APPROX
    //      (that is for the beginning of the sieve) we compute 'threshold' for
    //      every position x in the sieve. 'threshold' is taken as an
    //      approximation (minus a correction factor since we don't sieve with
    //      prime powers) of the value of the polynomial g(x).
    //
    //    - Otherwise, 'threshold' is computed as the aforementionned log
    //      approximate taken for x being the average value of the x index
    //      in the portion of the sieve we want to scan. Indeed, as soon as
    //      x is in the order of the sieve's length, log(g(x)) does not change
    //      much in the sieving interval.
    //
    // _NOTE_: scan_begin is given as a _relative_ value in the given chunk
    //         i.e. 0 <= scan_begin < sieve->alloced. If we are scanning on
    //         the positive side (i.e. context->sieve_negative_side==false)
    //         we start at the beginning of the array (or at the next unscanned
    //         position) and scan towards +infinity. If we are scanning on the
    //         negative side, we start at the end of the array and scan towards
    //         -infinity.
    //
    context->chunk_fully_scanned = false;

    if (! context->sieve_negative_side) {
        //
        // Scanning on the positive side (with x >= 0)
        //
        if (scan_begin + sieve_begin <= MAX_X_FOR_FULL_LOG_APPROX) {

            uint32_t global_x  = scan_begin + sieve_begin;
            uint32_t threshold = 0;
            uint32_t scan_end  = MIN(
                                     sieve->length - 1,
                                     MAX_X_FOR_FULL_LOG_APPROX
                                 );
            
            while (scan_begin <= scan_end) {
                //
                // We compute a different 'threshold' for each x as it can
                // significantly vary for small x values.
                //
                threshold = get_sieve_threshold(log_sqrtn, global_x);
                
                
                if (loglast < threshold) {
                    threshold -= loglast;
                    threshold  = (uint32_t)(
                                   (float)threshold * SIEVE_THRESHOLD_MULTIPLIER
                                 );
                } else {
                    threshold = 1;
                }
                
                if (sieve->data[scan_begin] >= threshold) {
                    //
                    // Check passed: store this x in the xpool array...
                    //
                    xpool_data[xpool_length] = global_x;
                    xpool_length++;
                    if (xpool_length == to_collect) {
                        if (scan_begin == scan_end) {
                            context->chunk_fully_scanned = true;
                        }
                        scan_begin++;
                        goto update_context_and_return;
                    }
                }
                scan_begin++;
                global_x++;
            }
        }

        if (scan_begin < sieve->length) {
            //
            // Just compute 'threshold' for the average value of x in the sieve
            // chunk to scan. In practice, this won't make any difference and
            // can save a few cycles.
            //
            uint32_t average_x = scan_begin + sieve_begin;
            average_x         += (sieve->length - scan_begin) / 2;

            uint32_t threshold = get_sieve_threshold(log_sqrtn, average_x);

            if (loglast < threshold) {
                threshold -= loglast;
                threshold  = (uint32_t)(
                                 (float)threshold
                               * SIEVE_THRESHOLD_MULTIPLIER
                             );
            } else {
                threshold = 1;
            }

            while (scan_begin < sieve->length) {
                if (sieve->data[scan_begin] >= threshold) {
                    //
                    // Check passed: store this x in the xpool array...
                    //
                    xpool_data[xpool_length] = scan_begin + sieve_begin;
                    xpool_length++;
                    if (xpool_length == to_collect) {
                        if (scan_begin == sieve->length - 1) {
                            context->chunk_fully_scanned = true;
                        }
                        scan_begin++;
                        goto update_context_and_return;
                    }
                }
                scan_begin++;
            }
        }
        context->chunk_fully_scanned = true;

    } else {
        //
        // Scanning on the negative side (with x < 0)
        //
        if (sieve_begin + (long int)scan_begin >= -MAX_X_FOR_FULL_LOG_APPROX) {
            uint32_t global_x  = -(sieve_begin + (long int)scan_begin);
            uint32_t threshold = 0;
            uint32_t scan_end  = MAX(
                                    0,
                                    -sieve_begin - MAX_X_FOR_FULL_LOG_APPROX
                                 );
            scan_begin++;

            while (scan_begin != scan_end) {
                scan_begin--;
                //
                // We compute a different 'threshold' for each x as it can
                // significantly vary for small x values.
                //
                threshold = get_sieve_threshold(log_sqrtn, global_x);

                if (loglast < threshold) {
                    threshold -= loglast;
                    threshold  = (uint32_t)(
                                   (float)threshold * SIEVE_THRESHOLD_MULTIPLIER
                                 );
                } else {
                    threshold = 1;
                }

                if (sieve->data[scan_begin] >= threshold) {
                    //
                    // Check passed: store this x in the xpool array...
                    //
                    xpool_data[xpool_length] = -((long int)global_x);
                    xpool_length++;
                    if (xpool_length == to_collect) {
                        if (scan_begin == 0) {
                            context->chunk_fully_scanned = true;
                        }
                        goto update_context_and_return;
                    }
                }
                global_x++;
            }
        }

        if (scan_begin != 0) {
            //
            // Just compute 'threshold' for the average value of x in the sieve
            // chunk to scan. In practice, this won't make any difference and
            // can save a few (not really significant) cycles.
            //
            uint32_t average_x = ABS(sieve_begin) - scan_begin / 2;
            uint32_t threshold = get_sieve_threshold(log_sqrtn, average_x);

            if (loglast < threshold) {
                threshold -= loglast;
                threshold  = (uint32_t)(
                                 (float)threshold
                               * SIEVE_THRESHOLD_MULTIPLIER
                             );
            } else {
                threshold = 1;
            }
            scan_begin++;

            while (scan_begin != 0) {
                scan_begin--;

                if (sieve->data[scan_begin] >= threshold) {
                    //
                    // Check passed: store this x in the xpool array...
                    //
                    xpool_data[xpool_length] = scan_begin + sieve_begin;
                    xpool_length++;
                    if (xpool_length == to_collect) {
                        if (scan_begin == 0) {
                            context->chunk_fully_scanned = true;
                        }
                        goto update_context_and_return;
                    }
                }
            }
        }
        context->chunk_fully_scanned = true;
    }

  update_context_and_return:

    context->scan_begin    = scan_begin;
    context->xpool->length = xpool_length;
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                Determination of factor base and misc. data
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
static void init_startup_data(qs_context_t* const context) {
    //
    // Initializes the factor base factor_base by keeping only those primes
    // pi such that the number to factor n has a square root mod pi. Also,
    // keeps (t - b) mod pi and (-t - b) mod pi in the arrays sol1 and sol2
    // respectively, where t is the modular square root of n mod pi and
    // b = ceil(sqrt(n)).
    //
    // Also precomputes the logarithms of the primes in the factor base and
    // stores them in log_factor_base.
    //
    uint32_array_t* const factor_base   = context->factor_base;
    byte_array_t* const log_factor_base = context->log_factor_base;
    uint32_array_t* const sol1          = context->sol1;
    uint32_array_t* const sol2          = context->sol2;
    const mpz_srcptr kn                 = context->kn;

    mpz_t b;
    mpz_init(b);
    mpz_sqrt(b, kn);
    mpz_add_ui(b, b, 1);
    //
    // Now b = floor(sqrt(b)) + 1
    //

    mpz_t negb;
    mpz_init(negb);
    mpz_neg(negb, b);
    //
    // Now negb = -b
    //

    mpz_t tmp;
    mpz_init(tmp);

    //
    // Always put 2 in the factor base...
    //
    factor_base->data[0] = 2;
    factor_base->length++;

    log_factor_base->data[0] = 1;
    log_factor_base->length++;

    sol1->data[sol1->length] = 0;
    sol1->length++;

    sol2->data[sol2->length] = 0;
    sol2->length++;

    uint32_t i = 0;
    uint32_t t = 0;
    uint32_t nmodp = 0;
    //
    // Determine the factor base and the solutions to the equations
    // t^2 = N (mod p) where the p's are the primes in the factor base
    //
    while (factor_base->length != factor_base->alloced) {

        i++;
        nmodp = mpz_fdiv_ui(kn, first_primes[i]);
        t     = sqrtm(nmodp, first_primes[i]);

        if (t != NO_SQRT_MOD_P) {
            //
            // i.e. if n has a square root mod first_primes[i]...
            //
            factor_base->data[factor_base->length] = first_primes[i];
            factor_base->length++;

            mpz_add_ui(tmp, negb, t);
            sol1->data[sol1->length] = mpz_fdiv_ui(tmp, first_primes[i]);
            sol1->length++;

            mpz_sub_ui(tmp, negb, t);
            sol2->data[sol2->length] = mpz_fdiv_ui(tmp, first_primes[i]);
            sol2->length++;
        }
    }
    context->loglast = ceil_log2(first_primes[i+1]);
    //
    // Precompute the logarithms of the primes in the factor base...
    //
    for (uint32_t j = 0; j < factor_base->length; j++) {
        log_factor_base->data[j] = most_significant_bit(factor_base->data[j]);
    }
    log_factor_base->length = factor_base->length;

    mpz_clear(b);
    mpz_clear(negb);
    mpz_clear(tmp);
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                  Miscellaneous utilitarian functions.
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
static void compute_polynomial_values(mpz_array_t* const cand_u,
                                      mpz_array_t* const cand_gx,
                                      const mpz_t n,
                                      const mpz_t b,
                                      const int32_array_t* const x) {
    //
    // Computes the values of the polynomial g(xi) = (xi + b)^2 - n for
    // all xi values in the array x and stores them in the array cand_gx.
    // Also keeps the (xi + b) in the array cand_u. Previous contents of the
    // cand_u and cand_gx array is overridden.
    //
    mpz_t polgx;
    mpz_init(polgx);

    for (uint32_t i = 0; i < x->length; i++) {

        if (x->data[i] > 0) {
            mpz_add_ui(polgx, b, x->data[i]);
        } else {
            mpz_sub_ui(polgx, b, -(x->data[i]));
        }
        mpz_set(cand_u->data[i], polgx);

        mpz_mul(polgx, polgx, polgx);
        mpz_sub(polgx, polgx, n);

        mpz_set(cand_gx->data[i], polgx);
    }
    mpz_clear(polgx);

    cand_u->length  = x->length;
    cand_gx->length = x->length;
}
//-----------------------------------------------------------------------------
inline uint32_t get_sieve_threshold(const uint32_t log_sqrtn,
                                    const uint32_t x) {
    //
    // Returns the threshold to be used as a criterion (before applying
    // a correction factor) during the sieve scanning phase. It is more or less
    // the value of the logarithm of g(x) which is about the logarithm of
    // 2 * x * sqrt(n).
    //
    return (floor_log2(x) + 1 + log_sqrtn);
};
//-----------------------------------------------------------------------------
