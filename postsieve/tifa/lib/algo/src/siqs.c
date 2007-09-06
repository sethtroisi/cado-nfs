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
 * \file    siqs.c
 * \author  Jerome Milan
 * \date    Late March / early April 2007
 * \version 1.1
 */

/*
 * History:
 *
 * 1.1: Late March / early April 2007 by JM
 *      - Completely refactored to use a factoring_machine_t.
 *
 * 1.0: Circa late summer/early fall 2006 by JM
 *      - Initial version.
 */

#include "tifa_config.h"

#include <stdlib.h>
#include <stdio.h>

#if TIFA_TIMING_SIQS
    #include <time.h>
#endif

#if TIFA_USE_CALLOC_MEMSET
    #include <string.h>
#endif

#include <math.h>

#include "first_primes.h"
#include "gmp_utils.h"
#include "funcs.h"
#include "macros.h"
#include "x_tree.h"
#include "hashtable.h"
#include "bernsteinisms.h"
#include "siqs.h"

#define __PREFIX__  "siqs: "
#define __VERBOSE__ TIFA_VERBOSE_SIQS
#define __TIMING__  TIFA_TIMING_SIQS

#include "messages.h"

//-----------------------------------------------------------------------------
//                         NON PUBLIC DEFINE(S)
//-----------------------------------------------------------------------------

//
// Multiplicative correction factor applied to the computed threshold to
// account for the fact that we do not sieve with powers of primes.
//
#define SIEVE_THRESHOLD_MULTIPLIER 0.80
//
// Largest size of a sieve chunk. The optimal value is architecture dependant,
// so this value should be tweaked according to the target processor.
//
#define SIEVE_CHUNK_MAX_SIZE       4096
//
// In order to not degrade too much the performance of the hashtable used
// for the large prime variation (due to too many collisions), we set its size
// to the size of the factor base multiplied by HTABLE_SIZE_MULTIPLIER.
//
#define HTABLE_SIZE_MULTIPLIER     8
//
// The first parameter of the polynomial used (i.e. 'a') is taken to be a
// product of primes in the factor base. The divisors of 'a' are choosen to be
// either, greater than MIN_PRIME_DECOMP...
//
#define MIN_PRIME_DECOMP           2000
//
// or than the NTH_ROOT_MIN_PRIME_DECOMP-th root of sqrt(2*kn)/sieve_half_width
// whichever is the lesser...
//
#define NTH_ROOT_MIN_PRIME_DECOMP  3
//
// Initial size of the array holding the different values of the paramaters 'a'
// used.
//
#define MAX_NA                     16
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                   NON PUBLIC STRUCTURE(S) AND TYPEDEF(S)
//-----------------------------------------------------------------------------

//
// An ad-hoc structure holding the data variables used by the
// SIQS implementation.
//
struct struct_siqs_context_t {
    //
    // Number to factor.
    //
    mpz_t n;
    //
    // A pointer to the SIQS parameters used.
    //
    siqs_params_t* params;
    //
    // Multiplier to use.
    //
    uint32_t multiplier;
    //
    // Number to factor x multiplier.
    //
    mpz_t kn;
    //
    // Floor of the logarithm in base 2 of sqrt(kn).
    //
    uint32_t log_sqrtkn;
    //
    // Factor base used. Only primes pi such that t^2 = N mod (pi) has a
    // solution for t are kept.
    //
    uint32_array_t *factor_base;
    //
    // Logarithms of the elements in the factor base.
    //
    uint32_array_t *log_factor_base;
    //
    // For each prime pi in the base, n has a square root mod pi. Keep
    // these modular square roots in the array sqrtm_pi such that
    // sqrtm_pi->data[i]^2 = n mod factor_base->data[i].
    //
    uint32_array_t *sqrtm_pi;
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
    // Matrix generated after collecting relations.
    //
    binary_matrix_t *matrix;
    //
    // Coefficients a and b of polynomial g_{a,b} used for sieving.
    // g_{a,b}(x) = (a.x + b)^2 - n
    //
    mpz_t a;
    mpz_t b;
    //
    // b^2 - kn = a.c, thus g_{a,b}(xi) = a.(a.xi^2 + 2.b.xi + c)
    //
    mpz_t c;
    //
    // Array keeping the indexes of the relations for which the polynomial
    // parameter 'a' changes. E.g, if min_index_of_a_array->data[i] = k1
    // and min_index_of_a_array->data[i+1] = k2, then the collected relations
    // numbered k1 to (k2-1) are associated to the 'a' value given by
    // a_array->data[i].
    //
    uint32_array_t* min_index_of_a_array;
    //
    // Number of relations (i.e. of pairs (a.xi+b, g_{a,b}(xi)) with
    // g_{a,b}(xi) smooth) we want to collect.
    //
    uint32_t npairs_wanted;
    //
    // Number of pairs (a.xi+b, g_{a,b}(xi)) to collect before proceeding to
    // a smoothness test.
    //
    uint32_t to_collect;
    //
    // Pool keeping the xi values that has survived the sieving stage, i.e.
    // the xi such that g_{a,b}(xi) _may_ be smooth.
    //
    int32_array_t *xpool;
    //
    // Size (i.e. length) of a sieve slice.
    //
    uint32_t size_sieve;
    //
    // The sieve is implemented as a uint32_array_t. Once the sieve is filled
    // and if we need further x values the same array is re-used. The array
    // sieve is consequently only a chunk of the conceptual global sieve. This
    // chunk is given by sieve_begin so that sieve->data[x] really means
    // the (x + sieve_begin)-th position in the complete sieve.
    //
    uint32_array_t *sieve;
    //
    // Relations to find: (a.xi+b)^2 = g_{a,b}(xi) (mod kn) with
    // g_{a,b}(xi) smooth. These are the g_{a, b}(xi)/a to check for
    // smoothness (cand_redgx) together with their corresponding ui = a.xi+b
    // and the values of the 'a' parameter...
    //
    mpz_array_t *cand_redgx;
    mpz_array_t *cand_u;
    mpz_array_t* cand_a_array;
    //
    // All accepted relations, i.e. the g_{a,b}(xi)/a and their
    // corresponding ui = a.xi+b, such that ui^2 = g_{a, b}(xi) (mod n) with
    // g_{a,b}(xi) smooth.
    //
    mpz_array_t *smooth_redgx;
    mpz_array_t *u;
    //
    // Values of the 'a' parameter associated to each smooth g_{a,b}(xi)/a in
    // the smooth_redgx array.
    //
    mpz_array_t *a_for_smooth_redgx;
    //
    // Hashtable used in large primes variation. It will hold mpz_pair_t's
    // such as ((a.xi + b), g_{a,b}(xi)) with an index in first_primes_array
    // as key.
    //
    hashtable_t *htable;
    //
    // First x of the current sieve slice. Since the sieve in an interval
    // [-M, M], sieve_begin can be negative in contrast to the QS case.
    //
    int32_t sieve_begin;
    //
    // Last x of the current sieve slice. Since the sieve in an interval
    // [-M, M], sieve_begin can be negative in contrast to the QS case.
    // In general sieve_end will be equal to sieve_begin + sieve_>alloced - 1
    // unless we have reach the end of the sieving interval.
    //
    int32_t sieve_end;
    //
    // The position in the _global_ sieve interval where we should begin the
    // next scanning phase (that will lead to keep only the xi that satisfy the
    // sieving criterion). The equality (scan_begin == sieve_end + 1) means
    // that all positions in the sieve array were checked and that another chunk
    // of the sieve should be computed if more values of x are needed.
    //
    int32_t scan_begin;
    //
    // Stores sol1i = inv(a).(ti-b) mod (pi) where ti^2 = N (mod pi) with
    // pi in the factor base.
    //
    uint32_array_t *sol1;
    //
    // Stores sol2i = inv(a).(-ti-b) mod (pi) where ti^2 = N (mod pi) with
    // pi in the factor base.
    //
    uint32_array_t *sol2;
    //
    // The number 'a' should ideally be as close as possible as 'to_approx',
    // that is to say sqrt(2*kn)/params->sieve_half_width.
    //
    mpz_t to_approx;
    //
    // The number of the current polynomial g_{a,b}(x) used.
    //
    uint32_t polynomial_number;
    //
    // The number of primes in the decomposition of 'a', 'a' being a product of
    // primes.
    //
    uint32_t nprimes_in_a;
    //
    // Index of the smallest prime in the prime decomposition of 'a'
    //
    uint32_t imin;
    //
    // Index of the largest prime in the prime decomposition of 'a'.
    //
    uint32_t imax;
    //
    // The number of polynomials g_{a,b}(x) that can be used to sieve, with 'a'
    // fixed. In other words, this is the number of 'b' values (with 'a' being
    // fixed) yielding a suitable g_{a,b}(x) polynomial).
    //
    uint32_t nb_poly_same_a;
    //
    // The bi such that bi^2 = N mod pi and bi = 0 mod pj where the pi's are
    // the primes involved in the decomposition of 'a'. These bi's values can
    // be used to deduce several values of the polynomial's 'b' parameter from
    // a single 'a' value.
    //
    mpz_array_t* Bl;
    //
    // Bainv2[i][p] = 2 x (Bl->data[i]) x (inv(a) mod p)
    //
    uint32_t** Bainv2;
    //
    // The threshold defined the sieving criteria. Position xi in the sieve
    // will survive (and g(xi + sieve_begin) tested for smoothness) if
    // sieve->data[xi] >= threshold.
    //
    uint32_t threshold;
    //
    // Floor of the logarithm in base 2 of pimax where pimax is the greatest
    // prime dividing 'a'.
    //
    uint32_t loglast;
    //
    // Floor of the logarithm in base 2 of 'a', 'b' and 'c'.
    //
    uint32_t loga;
    uint32_t logb;
    uint32_t logc;
    //
    // Average value of x in the sieve chunk to scan.
    //
    int32_t  av_x;
};
//-----------------------------------------------------------------------------
typedef struct struct_siqs_context_t siqs_context_t;
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                 PROTOTYPES OF NON PUBLIC FUNCTION(S)
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
static ecode_t init_siqs_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t clear_siqs_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t update_siqs_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t perform_siqs(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t recurse(
    mpz_array_t* const,
    uint32_array_t* const,
    const mpz_t,
    factoring_mode_t
);
//-----------------------------------------------------------------------------
static void init_startup_data(siqs_context_t* const context);
//-----------------------------------------------------------------------------
static void compute_reduced_polynomial_values(siqs_context_t* const context);
//-----------------------------------------------------------------------------
static void fill_sieve_array(
    uint32_array_t* const sieve,
    const int32_t sieve_begin,
    const int32_t sieve_end,
    const uint32_array_t* const factor_base,
    const uint32_array_t* const log_factor_base,
    const uint32_array_t* const sol1,
    const uint32_array_t* const sol2,
    const uint32_t a_pmin,
    const uint32_t a_pmax
);
//-----------------------------------------------------------------------------
static void scan_sieve_array(
    int32_array_t* const xpool,
    const uint32_array_t* const sieve,
    const int32_t sieve_begin,
    const uint32_t sieve_half_width,
    int32_t* scan_begin,
    const uint32_t threshold
);
//-----------------------------------------------------------------------------
static void determine_first_a(siqs_context_t* const context);
//-----------------------------------------------------------------------------
static void determine_next_a(siqs_context_t* const context);
//-----------------------------------------------------------------------------
static void init_first_polynomial(siqs_context_t* const context);
//-----------------------------------------------------------------------------
static void init_next_polynomial(siqs_context_t* const context);
//-----------------------------------------------------------------------------
static uint32_t get_sieve_threshold(siqs_context_t* const context);
//-----------------------------------------------------------------------------
static uint32_t best_multiplier(const mpz_t n);
//-----------------------------------------------------------------------------
static void collect_relations(siqs_context_t* const context);
//-----------------------------------------------------------------------------
static void update_threshold(siqs_context_t* const context);
//-----------------------------------------------------------------------------
static void update_current_polynomial(siqs_context_t* const context);
//-----------------------------------------------------------------------------
inline static void fill_sieve_chunk(siqs_context_t* const context);
//-----------------------------------------------------------------------------
inline static void scan_sieve_chunk(siqs_context_t* const context);
//-----------------------------------------------------------------------------
inline static void keep_relations_with_smooth_gx(siqs_context_t* const context);
//-----------------------------------------------------------------------------
static void prefill_matrix_with_a_decomps(siqs_context_t* const context);
//-----------------------------------------------------------------------------
static void multiply_by_smooth_a(siqs_context_t* const context);
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//                 "PUBLIC" FUNCTION(S) IMPLEMENTATION
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
ecode_t siqs(mpz_array_t* const factors, uint32_array_t* const multis,
             const mpz_t n, const siqs_params_t* const params,
             const factoring_mode_t mode) {

    PRINT_INTRO_MSG("self-initializing quadratic sieve");
    INIT_TIMER;
    START_TIMER;

    ecode_t ecode = NO_FACTOR_FOUND;

    factoring_machine_t machine;

    mpz_init_set(machine.n, n);
    machine.mode                = mode;
    machine.params              = (void*) params;
    machine.init_context_func   = init_siqs_context;
    machine.perform_algo_func   = perform_siqs;
    machine.update_context_func = update_siqs_context;
    machine.clear_context_func  = clear_siqs_context;
    machine.recurse_func        = recurse;
    machine.factors             = factors;
    machine.multis              = multis;

    ecode = run_machine(&machine);

    mpz_clear(machine.n);

    STOP_TIMER;
    PRINT_STATUS(machine.success, ecode);

    return ecode;
}
//------------------------------------------------------------------------------
void set_siqs_params_to_default(siqs_params_t* const params) {
    //
    // _NOTE_: For the time being we do not have enough empirical data to
    //         find good default values depending on the size of the number
    //         to factor... This should definitely be fixed before any serious
    //         use of the SIQS algorithm!
    //
    // _TO_DO_: Gather enough data to provide relevant and more or less optimal
    //          parameter values depending on the size of the number to factor.
    //
    params->sieve_half_width = SIQS_DFLT_SIEVE_HALF_WIDTH;
    params->nprimes_in_base  = SIQS_DFLT_NPRIMES_IN_BASE;
    params->nprimes_tdiv     = SIQS_DFLT_NPRIMES_TDIV;
    params->nrelations       = SIQS_DFLT_NRELATIONS;
    params->lsr_method       = SIQS_DFLT_LSR_METHOD;
    params->use_large_primes = SIQS_DFLT_USE_LARGE_PRIMES;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                  "PRIVATE" FUNCTION(S) IMPLEMENTATION
//------------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                  Functions used by factoring machine
 *-----------------------------------------------------------------------------
 */

//------------------------------------------------------------------------------
static ecode_t init_siqs_context(factoring_machine_t* const machine) {
    //
    // Initialize the SIQS implementation specific variables...
    //
    PRINT_INIT_MSG;
    INIT_TIMER;
    START_TIMER;

    machine->context        = (void*) malloc(sizeof(siqs_context_t));
    siqs_context_t* context = (siqs_context_t*) machine->context;
    siqs_params_t* params   = (siqs_params_t*)  machine->params;

    context->params = params;

    mpz_init_set(context->n, machine->n);
    //
    // Compute the best multiplier. The strategy used is similar to the one
    // described by M. A. Morrison and J. Brillhart in the remark 5.3 of
    // the CFRAC-focused paper "A Method of Factoring and the Factorization of
    // F_7" (Mathematics of Computation, vol 29, #129, Jan 1975, pages 183-205).
    //
    context->multiplier = best_multiplier(context->n);

    mpz_init(context->kn);
    mpz_mul_ui(context->kn, context->n, context->multiplier);

    context->factor_base     = alloc_uint32_array(params->nprimes_in_base);
    context->log_factor_base = alloc_uint32_array(
                                   context->factor_base->alloced
                               );
    context->sqrtm_pi        = alloc_uint32_array(
                                   context->factor_base->alloced
                               );
    //
    // Initialize all the data common to all the polynomials we are likely to
    // use for the sieving.
    //
    init_startup_data(context);
    context->ptree  = prod_tree_ui(context->factor_base);
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
    uint32_t e             = most_significant_bit(context->npairs_wanted);
    //
    // We'll collect more than the number of relations needed before attempting
    // a batch smoothness test a la Bernstein. See remark in the init_qs_context
    // function of QS in the file qs.c.
    //
    context->to_collect = 1 << (e + 1);

    context->xpool = alloc_int32_array(context->to_collect);

    context->size_sieve = (2 * params->sieve_half_width + 1);
    if (context->size_sieve > SIEVE_CHUNK_MAX_SIZE) {
        //
        // We limit the size of the sieve to make sure it can fit in the
        // cache memory of the processor, otherwise performance will suffer.
        // This is processor dependant so you may need to fine tune the value
        // of SIEVE_CHUNK_MAX_SIZE depending on the target processor.
        //
        context->size_sieve = SIEVE_CHUNK_MAX_SIZE;
    }
    context->sieve         = alloc_uint32_array(context->size_sieve);
    context->sieve->length = context->sieve->alloced;
    context->cand_redgx    = alloc_mpz_array(context->to_collect);
    context->cand_u        = alloc_mpz_array(context->to_collect);

    //
    context->cand_a_array = alloc_mpz_array(context->to_collect);

    for (uint32_t i = 0; i < context->to_collect; i++) {
        //
        // We fully initialize these arrays since they will be used over
        // and over again. Of course, we should remember to reset their
        // lengths to their alloced values before clering them!
        //
        mpz_init(context->cand_redgx->data[i]);
        mpz_init(context->cand_u->data[i]);
        //
        mpz_init(context->cand_a_array->data[i]);

    }
    context->smooth_redgx       = alloc_mpz_array(context->matrix_nrows);
    context->a_for_smooth_redgx = alloc_mpz_array(context->matrix_nrows);
    context->u                  = alloc_mpz_array(context->matrix_nrows);

    context->htable = NULL;

    if (params->use_large_primes) {
        uint32_t size = HTABLE_SIZE_MULTIPLIER * context->smooth_redgx->alloced;
        context->htable = alloc_init_hashtable(
                              size, uint32_cmp_func, hash_rj_32
                          );
    }

    mpz_init(context->a);
    mpz_init(context->b);
    mpz_init(context->c);

    context->imin = 0;
    context->imax = 0;

    //
    // Set 'to_approx' to sqrt(2 * 'kn')/'sieve_half_width'
    //
    mpz_init(context->to_approx);
    mpz_mul_ui(context->to_approx, context->kn, 2);
    mpz_sqrt(context->to_approx, context->to_approx);
    mpz_fdiv_q_ui(
        context->to_approx,
        context->to_approx,
        params->sieve_half_width
    );

    context->sieve_begin = -(int32_t)params->sieve_half_width;
    context->sieve_end   =   context->sieve_begin
                           + (int32_t)context->sieve->alloced - 1;

    if (context->sieve_end > (int32_t)params->sieve_half_width) {
        context->sieve_end = (int32_t)params->sieve_half_width;
    }
    context->scan_begin = context->sieve_begin;
    context->sol1       = alloc_uint32_array(context->factor_base->length);
    context->sol2       = alloc_uint32_array(context->factor_base->length);

    //
    // Compute the 'a' param of the first polynomial we'll use to sieve.
    //
    determine_first_a(context);

    context->polynomial_number = 1;
    context->nprimes_in_a      = context->imax - context->imin + 1;
    context->nb_poly_same_a    = 1 << (context->nprimes_in_a - 1);

    context->Bl = alloc_mpz_array(context->nprimes_in_a);

    for (uint32_t i = 0; i < context->nprimes_in_a; i++) {
        //
        // We fully initialize Bl since this array will be used over and over
        // again. No need to worry when clearing it since its length will always
        // equals its alloced value.
        //
        mpz_init(context->Bl->data[i]);
    }

    context->threshold = 0;
    context->loga      = 0;
    context->av_x      = 0;

    context->log_sqrtkn = mpz_sizeinbase(context->kn, 2) / 2;

    context->Bainv2 = malloc(context->factor_base->length * sizeof(uint32_t*));
    for (uint32_t i = 0; i < context->factor_base->length; i++) {
        context->Bainv2[i] = malloc(context->nprimes_in_a * sizeof(uint32_t));
    }

    init_first_polynomial(context);

    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//------------------------------------------------------------------------------
static ecode_t clear_siqs_context(factoring_machine_t* const machine) {

    PRINT_CLEAN_MSG;
    INIT_TIMER;
    START_TIMER;

    siqs_context_t* context = (siqs_context_t*) machine->context;

    if (context != NULL) {
        for (uint32_t i = 0; i < context->factor_base->length; i++) {
            free(context->Bainv2[i]);
        }
        free(context->Bainv2);
        //
        // Since the folowing arrays were fully initialized, don't forget to
        // set their current length to their alloced length before clearing
        // them...
        //
        context->xpool->length      = context->xpool->alloced;
        context->cand_redgx->length = context->cand_redgx->alloced;
        context->cand_u->length     = context->cand_u->alloced;

        context->cand_a_array->length     = context->cand_a_array->alloced;
        clear_mpz_array(context->cand_a_array);
        free(context->cand_a_array);


        clear_mpz_array(context->cand_u);
        free(context->cand_u);

        clear_mpz_array(context->u);
        free(context->u);

        clear_mpz_array(context->cand_redgx);
        free(context->cand_redgx);

        clear_mpz_array(context->smooth_redgx);
        free(context->smooth_redgx);

        clear_mpz_array(context->a_for_smooth_redgx);
        free(context->a_for_smooth_redgx);

        clear_mpz_array(context->Bl);
        free(context->Bl);

        clear_mpz_tree(context->ptree);
        free(context->ptree);

        clear_int32_array(context->xpool);
        free(context->xpool);

        clear_uint32_array(context->factor_base);
        free(context->factor_base);

        clear_uint32_array(context->log_factor_base);
        free(context->log_factor_base);

        clear_uint32_array(context->sol1);
        free(context->sol1);

        clear_uint32_array(context->sol2);
        free(context->sol2);

        clear_uint32_array(context->sieve);
        free(context->sieve);
        
        clear_uint32_array(context->sqrtm_pi);
        free(context->sqrtm_pi);

        clear_binary_matrix(context->matrix);
        free(context->matrix);

        mpz_clear(context->a);
        mpz_clear(context->b);
        mpz_clear(context->c);
        mpz_clear(context->n);
        mpz_clear(context->kn);
        mpz_clear(context->to_approx);

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
//------------------------------------------------------------------------------
static ecode_t update_siqs_context(factoring_machine_t* const machine) {
    //
    // _TO_DO_: No fallback strategy is implemented yet. We could for example
    //          try to gather more relations...
    //
    siqs_context_t* context = (siqs_context_t*) machine->context;
    context = NULL;
    return GIVING_UP;
}
//------------------------------------------------------------------------------
static ecode_t perform_siqs(factoring_machine_t* const machine) {

    INIT_TIMER;

    siqs_context_t* context = (siqs_context_t*) machine->context;
    siqs_params_t*  params  = context->params;
    
    PRINT_COLLECT_RELS_MSG;
    START_TIMER;

    collect_relations(context);

    STOP_TIMER;
    PRINT_COLLECT_RELS_DONE_MSG;
    PRINT_TIMING;

    //
    // We have now all the relations we need to fill the matrix and solve the
    // resulting linear system. The rest of the algorithm is now completely
    // standard to all congruence of squares methods.
    //

    PRINT_FACTOR_RES_MSG;
    RESET_TIMER;
    START_TIMER;

    //
    // Cofactors of generated smooth g_{a,b}(xi) after trial divisions.
    //
    mpz_array_t* partial_gx_array = alloc_mpz_array(context->matrix_nrows);
    //
    // Restrict the trial division to only params->nprimes_tdiv primes
    //
    uint32_array_t primes_array;
    primes_array.alloced = 0;
    primes_array.length  = params->nprimes_tdiv;
    primes_array.data    = (uint32_t*)(context->factor_base->data);

    //
    // Partly fill the matrix via trial divisions of the 'a's...
    //
    prefill_matrix_with_a_decomps(context);

    //
    // Partly fill the matrix via trial divisions of the 'g_{a,b}/a'...
    //
    fill_matrix_trial_div(
        context->matrix,
        partial_gx_array,
        context->smooth_redgx,
        &primes_array
    );

    if (params->nprimes_tdiv < params->nprimes_in_base) {
        //
        // The g_{a,b}/a are not completely factored on our factor base...
        //
        uint32_array_list_t*
        decomp_list = alloc_uint32_array_list(context->smooth_redgx->length);
        //
        // The primes in our factor base that was _not_ used in the
        // previous trial divisions are now pointed by primes_array
        //
        primes_array.length  =  params->nprimes_in_base - params->nprimes_tdiv;
        primes_array.data = (uint32_t*)
                            &(context->factor_base->data[params->nprimes_tdiv]);
        //
        // Use Bernstein's batch algo for the complete factorization of the
        // residues...
        //
        bern_71(decomp_list, partial_gx_array, &primes_array);
        //
        // ... and finish to fill the matrix...
        //
        fill_matrix_from_list(context->matrix, partial_gx_array, decomp_list,
                              context->factor_base);

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
    // relations of the form: A^2 = Q^2 mod n where  A^2 = prod(ui)^2 and
    // Q^2  = prod(g_{a,b}(xi)).
    //
    uint32_array_list_t *relations;
    relations = find_dependencies(context->matrix, params->lsr_method);

    STOP_TIMER;
    PRINT_TIMING;
    PRINT_DED_FACTORS_MSG;
    RESET_TIMER;
    START_TIMER;

    multiply_by_smooth_a(context);

    //
    // Deduce factors of n using the found relations.
    //
    ecode_t ecode = find_factors(
                        machine->factors,
                        context->n,
                        context->u,
                        context->smooth_redgx,
                        relations
                    );

    STOP_TIMER;
    PRINT_TIMING;

    clear_uint32_array_list(relations);
    free(relations);

    clear_mpz_array(partial_gx_array);
    free(partial_gx_array);

    return ecode;
}
//------------------------------------------------------------------------------
static ecode_t recurse(mpz_array_t* const factors, uint32_array_t* const multis,
                       const mpz_t n, factoring_mode_t mode) {

    siqs_params_t params;
    //
    // _NOTE_: For the time being, setting the parameters to their default
    //         values is pretty bad since this does not take into account
    //         the size of the number to factor. This should be changed.
    //
    //         See note in the set_siqs_params_to_default implementation.
    //
    set_siqs_params_to_default(&params);

    return siqs(factors, multis, n, &params, mode);
}
//------------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                     Sieving utilitarian functions
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
static void fill_sieve_array(
                        uint32_array_t* const sieve,
                        const int32_t sieve_begin,
                        const int32_t sieve_end,
                        const uint32_array_t* const factor_base,
                        const uint32_array_t* const log_factor_base,
                        const uint32_array_t* const sol1,
                        const uint32_array_t* const sol2,
                        const uint32_t a_pmin,
                        const uint32_t a_pmax) {
    //
    // Fill the sieve array 'sieve'.
    //
    // The sieve is computed by chunks of size sieve->alloced to avoid
    // a performance degradation due to the array not fitting in the
    // cache memory of the processor. This is completely architecture
    // dependant, so the value of sieve->alloced should theoretically
    // been adapted for each processor.
    //
    // The current chunk is given by sieve_begin and sieve_end (included)
    // such that sieve->data[0] really correspond to the (sieve_begin)-th
    // position of the real sieve.
    //
    // _NOTE_: Contrary to the fill_sieve_array function for the quadratic
    //         "mono-polynomial" sieve, sieve_begin and sieve_end can be
    //         negative since we're sieving on the interval -sieve_half_width
    //         to +sieve_half_width.
    //
    uint32_t logp     = 0;
    int32_t  imax     = 0;
    int32_t  imin     = 0;
    uint32_t curprime = 0;
    uint32_t cursol   = 0;
    uint32_t sindex   = 0;

    //
    // _WARNING_: This function uses explicit casts extensively.
    //            Indeed automatic type coercions _will_ give unwanted results!
    //

    //
    // (Re-)initialize the sieve array with zeroes...
    //
#if TIFA_USE_CALLOC_MEMSET
    memset(sieve->data, 0x00, sieve->length * sizeof(uint32_t));
#else
    for (uint32_t i = 0; i < sieve->length; i++) {
        sieve->data[i] = 0;
    }
#endif

    //
    // Sieve with the prime 2 only with sol1.
    //
    // imax (respectively imin) is the maximum (resp. minimum) value of i
    // for which (sol1->data[0] + 2*i) >= sieve_begin and
    //           (sol1->data[0] + 2*i) <= sieve_end
    //
    cursol = sol1->data[0];

    imax = (sieve_end - (int32_t)cursol);
    imax = imax >> 1;

    imin = (sieve_begin - (int32_t)cursol);
    imin = (imin >> 1) + (imin & 1);

    sindex = cursol + (2 * imin) - sieve_begin;

    for (int32_t i = imin; i <= imax; i++) {
       sieve->data[sindex] += 1;
       sindex += 2;
    }

    //
    // Sieve with primes > 2 with sol1 and sol2 but skip the primes that
    // divide the coefficient 'a'
    //
    for (uint32_t p = 1; p < factor_base->length; p++) {

        if (p == a_pmin) {
            //
            // Skip the primes that divide the coefficient 'a'
            //
            p = a_pmax + 1;
            if (p >= factor_base->length) {
               break;
            }
        }

        curprime = factor_base->data[p];
        cursol   = sol1->data[p];
        logp     = log_factor_base->data[p];
        //
        // Sieve with sol1
        //
        // imax (respectively imin) is the maximum (resp. minimum) value of i
        // for which:
        //   (sol1->data[p] + factor_base->data[p]*i) >= sieve_begin
        // and:
        //   (sol1->data[p] + factor_base->data[p]*i) <= sieve_end
        //
        imax = (sieve_end - (int32_t)cursol);
        if (imax >= 0) {
            imax = imax / (int32_t)curprime;
        } else {
            imax = ((imax + 1) / (int32_t)curprime) - 1;
        }

        imin = (sieve_begin - (int32_t)cursol);
        if (imin <= 0) {
            imin = (imin / (int32_t)curprime);
        } else {
            imin = ((imin - 1) / (int32_t)curprime) + 1;
        }

        sindex = cursol + (curprime * imin) - sieve_begin;

        for (int32_t i = imin; i <= imax; i++) {
            sieve->data[sindex] += logp;
            sindex += curprime;
        }
        //
        // Sieve with sol2
        //
        // imax (respectively imin) is the maximum (resp. minimum) value of i
        // for which:
        //   (sol2->data[p] + factor_base->data[p]*i) >= sieve_begin
        // and
        //   (sol2->data[p] + factor_base->data[p]*i) <= sieve_end
        //
        cursol = sol2->data[p];

        imax = (sieve_end - (int32_t)cursol);
        if (imax >= 0) {
            imax = imax / (int32_t)curprime;
        } else {
            imax = ((imax + 1) / (int32_t)curprime) - 1;
        }

        imin = (sieve_begin - (int32_t)cursol);
        if (imin <= 0) {
            imin = (imin / (int32_t)curprime);
        } else {
            imin = ((imin - 1) / (int32_t)curprime) + 1;
        }

        sindex = cursol + (curprime * imin) - sieve_begin;

        for (int32_t i = imin; i <= imax; i++) {
            sieve->data[sindex] += logp;
            sindex += curprime;
        }
    }
}
//-----------------------------------------------------------------------------
static void scan_sieve_array(
                        int32_array_t* const xpool,
                        const uint32_array_t* const sieve,
                        const int32_t sieve_begin,
                        const uint32_t sieve_half_width,
                        int32_t* scan_begin,
                        const uint32_t threshold) {

    //
    // Scans the sieve array sieve from the position scan_begin (included)
    // to the position min(sieve_begin + sieve->length - 1, sieve_half_width)
    // (also included) and stores in xpool only the locations x for which:
    //     sieve->data[x] >= threshold
    //
    // sieve_begin and scan_begin can be negative in contrast to what happens
    // in the standard QS version of this function.
    //
    // Here, the only contraints are:
    //     sieve_begin <= sieve_half_width
    // and:
    //     sieve_begin <= *scan_begin <= sieve_half_width
    //
    int32_t  sieve_end = 0;

    if (   (sieve_begin + (int32_t)sieve->length - 1)
         < (int32_t)sieve_half_width) {
        sieve_end = sieve_begin + (int32_t)sieve->length - 1;
    } else {
        sieve_end = (int32_t)sieve_half_width;
    }

    //
    // _NOTE_: Use local variables to avoid numerous readings/writings
    //         of pointed data...
    //
    uint32_t tmp_scan_begin = *scan_begin;
    uint32_t tmp_length     = xpool->length;

    for (int32_t xindex = *scan_begin; xindex <= sieve_end; xindex++) {

        tmp_scan_begin++;
        //
        // _REMINDER_: scan_begin is given as a value in the _global_ conceptual
        //             sieve, and _not_ as a _relative_ value in the given chunk
        //             as that was the case with the QS version of the
        //             scan_sieve function.
        //
        if (sieve->data[xindex - sieve_begin] >= threshold) {
            //
            // Check passed: store this x in the xpool array...
            //
            xpool->data[tmp_length] = (int32_t)xindex;
            tmp_length++;
            if (tmp_length == xpool->alloced) {
                break;
            }
        }
    }
    *scan_begin   = tmp_scan_begin;
    xpool->length = tmp_length;
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                Determination of factor base and misc. data
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
static void init_startup_data(siqs_context_t* const context) {

    uint32_array_t* const factor_base     = context->factor_base;
    uint32_array_t* const log_factor_base = context->log_factor_base;
    uint32_array_t* const sqrtm_pi        = context->sqrtm_pi;
    const mpz_srcptr kn                   = context->kn;

    //
    // Initializes the factor base factor_base by keeping only those primes
    // pi such that the number to factor n has a square root mod pi, and keeps
    // one of these square roots in the array sqrtm_pi.
    //
    // Also precomputes the logarithms of the primes in the factor base and
    // stores them in the array log_factor_base.
    //

    //
    // Always put 2 in the factor base...
    //
    factor_base->data[0] = 2;
    factor_base->length  = 1;

    log_factor_base->data[0] = 1;
    log_factor_base->length  = 1;

    if (0 == mpz_tstbit(kn, 0)) {
        sqrtm_pi->data[0] = 0;
    } else {
        sqrtm_pi->data[0] = 1;
    }
    sqrtm_pi->length  = 1;

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

            sqrtm_pi->data[sqrtm_pi->length] = t;
            sqrtm_pi->length++;
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
}
//-----------------------------------------------------------------------------
static void determine_first_a(siqs_context_t* const context) {

    mpz_ptr approxed                 = context->a;
    uint32_t* const imin             = &(context->imin);
    uint32_t* const imax             = &(context->imax);
    const mpz_ptr to_approx          = context->to_approx;
    const uint32_array_t* const base = context->factor_base;

    //
    // Determines 'a', i.e the leading coefficient of the first serie of
    // polynoms g_{a,b}(x) = (ax + b)^2 - n  used by SIQS.
    //
    // More precisely, given the integer 'to_approx', determines 'approxed'
    // (i.e. the coefficient 'a' in our case), 'imin' and 'imax' such that:
    //      'approxed' = base->data[imin] x ... x base->data[imax]
    // and:
    //      'approxed' <= 'to_approx'
    // and:
    //      'approxed' x base->data[imax+1] > 'to_approx'.
    // and:
    //      base->data[imin] >= min(MIN_PRIME_DECOMP,
    //                              'to_approx'/DIVISOR_DECOMP)
    //
    // Of course, 'approxed' is NOT the best possible approximation of
    // 'to_approx' using the integers from the 'base' array. Actually, for
    // this SIQS implementation, we want 'approxed' to be composed by as many
    // factors as possible without being wildly smaller or larger that
    // 'to_approx'.
    //
    // Also note that the present implementation requires the 'base' array to
    // be sorted.
    //

    *imin = 0;
    *imax = 0;

    //
    // In his master's thesis "Factoring Integers with the Self-initializing
    // Quadratic Sieve", S.P. Contini recommends to exclude small factors in
    // the product of primes. To account for possible small numbers, we require
    // the smallest prime allowed to be larger than MIN_PRIME_DECOMP or
    // than the NTH_ROOT_MIN_PRIME_DECOMP-th root of 'to_approx' , whichever is
    // smaller.
    //
    mpz_t threshold;
    mpz_init(threshold);
    mpz_root(threshold, to_approx, NTH_ROOT_MIN_PRIME_DECOMP);

    if (mpz_cmp_ui(threshold, MIN_PRIME_DECOMP) >= 0) {
        mpz_set_ui(threshold, MIN_PRIME_DECOMP);
    }

    while (mpz_cmp_ui(threshold, base->data[*imin]) >= 0) {
        (*imin)++;
    }
    *imax = *imin;
    //
    // Finding 'imin' and 'imax' in the simplest fashion: just multiply the
    // smallest numbers from the base (without taking twice the same factor)
    // until 'approxed' is larger than 'to_approx'...
    //
    mpz_set_ui(approxed, base->data[*imax]);

    while (mpz_cmp(approxed, to_approx) <= 0) {
        (*imax)++;
        mpz_mul_ui(approxed, approxed, base->data[*imax]);
    }
    //
    // ...and then, if possible, backtrap one step so that
    // 'approxed' <= 'to_approx'
    //
    if (*imax > *imin) {
        mpz_divexact_ui(approxed, approxed, base->data[*imax]);
        (*imax)--;
    }
    mpz_clear(threshold);

    context->loga = mpz_sizeinbase(context->a, 2);
}
//-----------------------------------------------------------------------------
static void determine_next_a(siqs_context_t* const context) {

    mpz_ptr a                        = context->a;
    uint32_t* const imin             = &(context->imin);
    uint32_t* const imax             = &(context->imax);
    const uint32_array_t* const base = context->factor_base;

    //
    // In the quite frequent occasion where all the polynomials having 'a'
    // as their leading coefficient have been used but not enough congruences
    // of square have been collected, it is neccesary to generate a new "family"
    // of polynomials by using a different leading coefficient 'a'.
    //
    // To keep the coefficients at more or less the same order of magnitude,
    // the new 'a' is derived from the previous 'a' simply by incrementing the
    // 'imin' and 'imax' indexes by 2 (and updating 'a' accordingly). Why not
    // incrementing these indexes by 1? Well, it seems that we get a lot of
    // redundant relations if several values of 'a' differ only by a single
    // prime. This is explained in S.P. Contini's thesis "Factoring
    // integers with the Self-initializing Quadratic Sieve" and was indeed
    // confirmed by our own experimentations.
    //
    mpz_divexact_ui(a, a, base->data[*imin]);
    (*imin)++;
    mpz_divexact_ui(a, a, base->data[*imin]);
    (*imin)++;

    if ((*imax) >= (base->length + 1)) {
        fprintf(stderr, "%s: ERROR: imax >= factor_base->length!\n",
                __func__);
        exit(-1);
    }

    (*imax)++;
    mpz_mul_ui(a, a, base->data[*imax]);
    (*imax)++;
    mpz_mul_ui(a, a, base->data[*imax]);

    context->loga = mpz_sizeinbase(context->a, 2);
}
//------------------------------------------------------------------------------
static void init_first_polynomial(siqs_context_t* const context) {

    const uint32_array_t* factor_base = context->factor_base;
    const uint32_array_t* sqrtm_pi    = context->sqrtm_pi;
    const uint32_t imin               = context->imin;
    const uint32_t imax               = context->imax;
    const mpz_ptr a                   = context->a;
    mpz_ptr b                         = context->b;
    mpz_ptr c                         = context->c;
    mpz_array_t* const Bl             = context->Bl;
    uint32_t** Bainv2                 = context->Bainv2;
    uint32_array_t* const sol1        = context->sol1;
    uint32_array_t* const sol2        = context->sol2;

    //
    // Determines the first polynomial having 'a' as its leading coefficient
    // and performs all needed precomputations. Notations are mainly that of
    // S.P. Contini's thesis "Factoring integers with the Self-initializing
    // Quadratic Sieve".
    //
    mpz_t inv;
    mpz_init(inv);
    mpz_t mpzp;
    mpz_init(mpzp);

    uint64_t gamma    = 0;
    int64_t  tmp64    = 0;
    uint32_t curprime = 0;
    uint32_t i        = 0;

    mpz_set_ui(b, 0);
    //
    // Computes the {Bl}, 1 <= l <= s such that:
    //     Bl^2 = N mod ql
    //     Bl = 0 mod ql forall j != l
    // where:
    //     qi are the primes such that a = q1 x q2 x ... x qs
    //
    for (uint32_t l = 0; l <= (imax - imin); l++) {

        i = l + imin;
        if (i >= factor_base->length) {
            break;
        }
        curprime = factor_base->data[i];
        mpz_set_ui(mpzp, curprime);
        //
        // Don't forget that Bl has already been completely initialized, so
        // no mpz_init calls are necessary
        //
        mpz_divexact_ui(Bl->data[l], a, curprime);

        if (0 == mpz_invert(inv, Bl->data[l], mpzp)) {
            gmp_fprintf(stderr,
                        "%s 1: ERROR: %Zd has no inverse mod %Zd!\n",
                        __func__, Bl->data[l], mpzp);
            fprintf(stderr, "This should not happen - Program aborted!\n");
            fflush(stderr);
            exit(-1);
        };
        gamma  = mpz_get_ui(inv) * sqrtm_pi->data[i];
        //
        // _NOTE_: One has to be careful when dealing with a mix of several
        //         integer types and pay attention to the possible automatic
        //         type coercion. For example a int64_t % a uint32_t will give
        //         the expected result, but not a int32_t % a uint32_t...
        //
        gamma %= factor_base->data[i];
        if (gamma > (curprime >> 1)) {
            gamma = curprime - gamma;
        }
        mpz_mul_ui(Bl->data[l], Bl->data[l], (unsigned long int)gamma);
        mpz_add(b, b, Bl->data[l]);
    }
    Bl->length = Bl->alloced;

    mpz_mul(c, b, b);
    mpz_sub(c, c, context->kn);
    mpz_divexact(c, c, a);

    context->logb = mpz_sizeinbase(b, 2);
    context->logc = mpz_sizeinbase(c, 2);

    //
    // Computes the two solutions to g_{a,b}(x) = 0 mod p for each prime
    // and stores them in the arrays sol1 and sol2.
    //
    // Also, for each prime p that does not divide the polynomial coefficient
    // 'a', computes and store the values 2 x Bj x inv(a).mod(p) mod p foreach
    // prime pj dividing 'a'. These values will be kept in the array
    // Bainv2 indexed first by the primes p, then by the indexes j.
    //
    for (uint32_t i = 0; i < factor_base->length; i++) {

        if (i == imin) {
            //
            // Skip the primes that divide the coefficient 'a'...
            //
            i = imax + 1;
            if (i >= factor_base->length) {
                break;
            }
        }
        curprime = factor_base->data[i];
        mpz_set_ui(mpzp, curprime);

        if (0 == mpz_invert(inv, a, mpzp)) {
            gmp_fprintf(stderr,
                        "%s 2: ERROR: %Zd has no inverse mod %"PRIu32"!\n",
                        __func__, a, curprime);
            fprintf(stderr, "This should not happen - Program aborted!\n");
            fflush(stderr);
            exit(-1);
        };

        for (uint32_t j = 0; j <= (imax - imin); j++) {
            gamma  = mpz_fdiv_ui(Bl->data[j], curprime);
            gamma *= (2 * mpz_get_ui(inv));
            gamma %= factor_base->data[i];
            Bainv2[i][j] = (uint32_t) gamma;
        }
        //
        // Computes the first solution to the g_{a,b}(x) = 0 mod p congruence
        //
        tmp64  = 0;
        tmp64  = (int64_t)(sqrtm_pi->data[i]);
        tmp64 -= (int64_t)(mpz_fdiv_ui(b, curprime));
        tmp64 %= curprime;
        tmp64 *= mpz_get_ui(inv);
        tmp64 %= curprime;
        if (tmp64 < 0) {
            tmp64 += curprime;
        }
        sol1->data[i] = (uint32_t) tmp64;
        //
        // Computes the second solution to the g_{a,b}(x) = 0 mod p congruence
        //
        tmp64  = 0;
        tmp64 -= (int64_t)(sqrtm_pi->data[i]);
        tmp64 -= (int64_t)(mpz_fdiv_ui(b, curprime));
        tmp64 %= curprime;
        tmp64 *= mpz_get_ui(inv);
        tmp64 %= curprime;
        if (tmp64 < 0) {
            tmp64 += curprime;
        }
        sol2->data[i] = (uint32_t) tmp64;
    }
    //
    // Setting the length are not strictly needed (as we don't use them), but
    // let's try to maintain a bit of consistency...
    //
    sol1->length   = sol1->alloced;
    sol2->length   = sol2->alloced;

    mpz_clear(inv);
    mpz_clear(mpzp);
}
//------------------------------------------------------------------------------
static void init_next_polynomial(siqs_context_t* const context) {

    const uint32_t ipol              = context->polynomial_number;
    const uint32_array_t* const base = context->factor_base;
    const uint32_t imin              = context->imin;
    const uint32_t imax              = context->imax;
    mpz_ptr b                        = context->b;
    mpz_ptr c                        = context->c;
    const mpz_array_t* const Bl      = context->Bl;
    uint32_t** const Bainv2          = context->Bainv2;
    uint32_array_t* const sol1       = context->sol1;
    uint32_array_t* const sol2       = context->sol2;

    //
    // Determines the next polynomial having 'a' as its leading coefficient and
    // updates the needed data. Notations are mainly that of S.P. Contini's
    // thesis "Factoring integers with the Self-initializing Quadratic Sieve".
    //
    // The current polynomial is indexed by ipol so we are now initializing
    // the polynomial indexed by (ipol + 1)
    //

    //
    // Compute the 'b' coefficient of the next polynomial to use
    //
    uint32_t nu = 0;
    while ( ((2 * ipol) & (1<<nu)) == 0) {
        nu++;
    }
    uint32_t io2nu = (ipol / (1<<nu)) + 1;

    if ((io2nu & 1) == 0) {
        mpz_add(b, b, Bl->data[nu-1]);
        mpz_add(b, b, Bl->data[nu-1]);
    } else {
        mpz_sub(b, b, Bl->data[nu-1]);
        mpz_sub(b, b, Bl->data[nu-1]);
    }
    mpz_mul(c, b, b);
    mpz_sub(c, c, context->kn);
    mpz_divexact(c, c, context->a);

    context->logb = mpz_sizeinbase(b, 2);
    context->logc = mpz_sizeinbase(c, 2);

    uint32_t curprime = 0;
    int32_t  tmp32    = 0;

    //
    // Compute the new sol1 and sol2 arrrays for sieving with the new polynomial
    //
    for (uint32_t i = 0; i < base->length; i++) {

        if (i == imin) {
            //
            // Skip the primes that divide the coefficient 'a'
            //
            i = imax + 1;
            if (i >= base->length) {
               break;
            }
        }
        curprime = base->data[i];

        if ((io2nu & 1) != 0) {
            sol1->data[i] += Bainv2[i][nu-1];
            sol1->data[i] %= curprime;
            sol2->data[i] += Bainv2[i][nu-1];
            sol2->data[i] %= curprime;
        } else {
            //
            // _NOTE_ : Be very careful with the possible negative case: use
            //          a temporary _signed_ variable and beware of
            //          automatic type coercion that can @#%!* results.
            //
            // _TODO_: Do we really need _unsigned_ variables anyway?
            //
            tmp32  = (int32_t)sol1->data[i] - (int32_t)Bainv2[i][nu-1];
            tmp32 %= (int32_t)curprime;
            if (tmp32 < 0) {
                tmp32 += (int32_t)curprime;
            }
            sol1->data[i] = (uint32_t)tmp32;

            tmp32  = (int32_t)sol2->data[i] - (int32_t)Bainv2[i][nu-1];
            tmp32 %= (int32_t) curprime;
            if (tmp32 < 0) {
                tmp32 += (int32_t)curprime;
            }
            sol2->data[i] = (uint32_t)tmp32;
        }
    }
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                  Miscellaneous utilitarian functions.
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
static void compute_reduced_polynomial_values(siqs_context_t* const context) {

    const mpz_ptr a = context->a;
    const mpz_ptr b = context->b;
    const mpz_ptr c = context->c;

    mpz_array_t* const cand_u       = context->cand_u;
    mpz_array_t* const cand_redgx   = context->cand_redgx;
    const int32_array_t* const x    = context->xpool;

    mpz_array_t* const cand_a_array = context->cand_a_array;

    //
    // Computes the "reduced" values g_{a,b}(xi)/a = a.xi^2 + 2.b.xi + c for
    // all xi values in the array x (with i >= cand_redgx->length) and stores
    // them in the array cand_redgx. Also keeps the (a.xi + b) in the array
    // cand_u.
    //
    mpz_t rgx;
    mpz_init(rgx);

    for (uint32_t i = cand_redgx->length; i < x->length; i++) {

        if ((x->data[i]) < 0) {
            mpz_mul_ui(rgx, a, (long unsigned int)(-(x->data[i])));
            mpz_neg(rgx, rgx);
        } else {
            mpz_mul_ui(rgx, a, x->data[i]);
        }
        mpz_add(rgx, rgx, b);
        mpz_set(cand_u->data[i], rgx);

        mpz_add(rgx, rgx, b);

        if ((x->data[i]) < 0) {
            mpz_mul_ui(rgx, rgx, (long unsigned int)(-(x->data[i])));
            mpz_neg(rgx, rgx);
        } else {
            mpz_mul_ui(rgx, rgx, x->data[i]);
        }
        mpz_add(rgx, rgx, c);

        mpz_set(cand_redgx->data[i], rgx);
        mpz_set(cand_a_array->data[i], a);
    }
    cand_redgx->length   = x->length;
    cand_u->length       = x->length;
    cand_a_array->length = x->length;

    mpz_clear(rgx);
}
//-----------------------------------------------------------------------------
static uint32_t get_sieve_threshold(siqs_context_t* const context) {
    //
    // Get the 'theoretical' sieve threshold: this is an approximate value of
    // the logarithm of g_{a,b}/a in base 2.
    //
    if (context->av_x == 0) {
        return context->logc;
    }
    uint32_t logx      = floor_log2(ABS(context->av_x));
    uint32_t threshold = context->loga + 2 * logx;

    threshold = MAX(threshold, 1 + context->logb + logx) + 1;
    threshold = MAX(threshold, context->logc);

    return threshold;
};
//-----------------------------------------------------------------------------
static void multiply_by_smooth_a(siqs_context_t* const context) {
    //
    // Multiply each of the smooth value of g_{a,b}/a in the smooth_redgx array
    // by their corresponding 'a'. This is needed at the very end of the SIQS
    // algorithm when we actually compute the (hopefully non trivial) factors...
    //
    for (uint32_t i = 0; i < context->smooth_redgx->length; i++) {
        mpz_mul(
            context->smooth_redgx->data[i],
            context->smooth_redgx->data[i],
            context->a_for_smooth_redgx->data[i]
        );
    }
}
//-----------------------------------------------------------------------------
static void prefill_matrix_with_a_decomps(siqs_context_t* const context) {
    //
    // 'Prefill' the binary matrix using the decomposition of the values of
    // 'a' associated to each smooth values of g_{a,b}/a in the array
    // smooth_redgx...
    //
    const uint32_array_t* const factor_base = context->factor_base;
    binary_matrix_t* const matrix  = context->matrix;
    mpz_array_t*     const a_array = context->a_for_smooth_redgx;

    uint32_t nrows_to_add = a_array->length;
    const uint32_t imax   = context->imax + 1;

    mpz_t tmp;
    mpz_init(tmp);

    for (uint32_t i = 0; i < nrows_to_add; i++) {

        mpz_set(tmp, a_array->data[i]);

        for (uint32_t j = 0; j < imax; j++) {

            if  (0 != mpz_divisible_ui_p(tmp, factor_base->data[j]) ) {
                mpz_divexact_ui(tmp, tmp, factor_base->data[j]);
                flip_matrix_bit(i, j + 1, matrix);
            }
            if (mpz_cmp_ui(tmp, 1) == 0) {
                break;
            }
        }
    }
    mpz_clear(tmp);
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                  Determination of best multiplier
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
static uint32_t best_multiplier(const mpz_t n) {
    //
    // Chooses the best multiplier to use for the SIQS algorithm.
    //
    // Chooses the best multiplier to use for the SIQS algorithm, following the
    // suggestion of M. A. Morrison and J. Brillhart in the remark 5.3 of
    // the CFRAC-focused paper "A Method of Factoring and the Factorization of
    // F_7" (Mathematics of Computation, vol 29, #129, Jan 1975, pages 183-205).
    //
    // \param[in] n The integer to factor.
    //
    mult_data_t* mdata = malloc(LARGEST_MULTIPLIER*sizeof(mult_data_t));

    mpz_t kn;
    mpz_init(kn);

    uint32_t ret_val;

    for (uint32_t k = 0; k < LARGEST_MULTIPLIER; k++) {

        mdata[k].multiplier = k+1;
        mdata[k].count      = 0;
        mdata[k].sum_inv_pi = 0.0;

        mpz_mul_ui(kn, n, mdata[k].multiplier);

        for (uint32_t j = 0; j < MAX_IPRIME_IN_MULT_CALC; j++) {

            if (    (-1 != mpz_kronecker_ui(kn, 3))
                 || (-1 != mpz_kronecker_ui(kn, 5))
               ) {

                if (-1 != mpz_kronecker_ui(kn, first_primes[j])) {
                    mdata[k].count++;
                    mdata[k].sum_inv_pi += 1/(double)first_primes[j];
                }
            }
        }
    }
    //
    // The cmp_mult_data function gives the meaning of what a "better"
    // multiplier is.
    //
    qsort(mdata, LARGEST_MULTIPLIER, sizeof(mult_data_t), cmp_mult_data);

    ret_val = mdata[LARGEST_MULTIPLIER-1].multiplier;

    mpz_clear(kn);
    free(mdata);

    return ret_val;
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *          Collection of relations (a.xi + b)^2 = g_{a,b}(xi) mod (kn)
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
static void collect_relations(siqs_context_t* const context) {

    //
    // Collect the relations (a.xi + b)^2 = g_{a,b}(xi) mod (kn). This is
    // of course the 'core' of the SIQS algorithm...
    //

    siqs_params_t*  params  = context->params;

    bool chunk_partially_scanned = false;
    bool sieve_fully_scanned     = false;

    //
    // Collect relations that will ultimately lead to congruences of squares,
    // i.e. keep the (a.x + b) and their associated g_{a,b}(x) values iif
    // g_{a,b}(x) is smooth on the factor base.
    //
    while (context->smooth_redgx->length != context->smooth_redgx->alloced) {
        //
        // Note that we may collect more than the remaining number of relations
        // needed before attempting a batch smoothness test a la Bernstein.
        //
        while (context->xpool->length != context->xpool->alloced) {

            if (! chunk_partially_scanned) {
                //
                // The sieve slice has either been completely scanned for
                // candidate xi values or not been computed at all. In
                // either case, we need to fill a (new) sieve slice.
                //
                fill_sieve_chunk(context);
            }

            //
            // Scan the yet unscanned part of the slice to keep xi values
            // for which g_{a,b}(xi) _may_ be smooth...
            //
            update_threshold(context);
            scan_sieve_chunk(context);

            if (context->scan_begin == (int32_t)params->sieve_half_width + 1) {
                //
                // We have completely sieved the interval
                // [-params->sieve_half_width, params->sieve_half_width].
                //
                // _NOTE_: For the time being we perform a smoothness batch
                //         whenever we have to change to another polynomial
                //         even if the ideal batch size is not reached. This
                //         makes for a simpler logic (we don't have to keep
                //         track of all the polynomials used) but may be less
                //         efficient.
                //
                sieve_fully_scanned     = true;
                chunk_partially_scanned = false;

            } else {

                if (context->scan_begin == context->sieve_end + 1) {
                    //
                    // We just finished to scan a whole chunk of the sieve,
                    // so let's prepare for the next one...
                    //
                    context->sieve_begin = context->sieve_end + 1;
                    context->sieve_end  += (int32_t)context->sieve->alloced;

                    if (context->sieve_end>(int32_t)params->sieve_half_width) {
                        //
                        // Let's stay within the sieving interval boundaries...
                        //
                        context->sieve_end = (int32_t)params->sieve_half_width;
                    }
                    context->scan_begin = context->sieve_begin;
                    chunk_partially_scanned = false;
                } else {
                    //
                    // The scan stopped in the middle of the sieve slice
                    // because we have filled the xpool array. The next
                    // scan (if any is needed) will start from where we
                    // stopped.
                    //
                    chunk_partially_scanned = true;
                }
            }
            if (sieve_fully_scanned) {
                //
                // We should compute the "reduced" g_{a,b}(xi)/a BEFORE updating
                // the current polynomial...
                //
                compute_reduced_polynomial_values(context);

                update_current_polynomial(context);

                context->scan_begin  = -(int32_t)params->sieve_half_width;
                context->sieve_begin = -(int32_t)params->sieve_half_width;
                context->sieve_end   = context->sieve_begin;
                context->sieve_end  += (int32_t)context->sieve->alloced - 1;

                if (context->sieve_end > (int32_t)params->sieve_half_width) {
                    context->sieve_end = (int32_t)params->sieve_half_width;
                }
                sieve_fully_scanned     = false;
                chunk_partially_scanned = false;
            } else {
                //
                // Compute the "reduced" g_{a,b}(xi)/a...
                //
                compute_reduced_polynomial_values(context);
            }
        }
        //
        // After computing the "reduced" g_{a,b}(xi)/a, keep the ones that are
        // smooth in the smooth_redgx array together with the associated
        // u = (a.x + b) values...
        //
        keep_relations_with_smooth_gx(context);

        context->xpool->length        = 0;
        context->cand_redgx->length   = 0;
        context->cand_u->length       = 0;
        context->cand_a_array->length = 0;

        PRINT_NRELS_FOUND(
            context->smooth_redgx->length,
            context->smooth_redgx->alloced
        );
    }
}
//-----------------------------------------------------------------------------
static void update_threshold(siqs_context_t* const context) {
    //
    // Update the SIQS context to change the sieving criterion.
    //
    context->av_x = context->sieve_begin + (int32_t)(context->sieve->length/2);

    //
    // Get an approximationf of log_2(g_{a,b}(avx)) for avx the absolute value
    // of the average of x in the scanning interval.
    //
    context->threshold = get_sieve_threshold(context);

    //
    // Correct this 'theoretical' threshold by substracting the value of the
    // logarithm (in base 2) of the smallest prime greater than every primes
    // in the factor base...
    //
    if (context->loglast < context->threshold) {
        context->threshold -= context->loglast;
        //
        // Also multiply this 'corrected' threshold by a final multiplicative
        // correction factor to account for the fact that we do not sieve with
        // prime powers...
        //
        context->threshold  = (uint32_t)(
                                  (float)context->threshold
                                * SIEVE_THRESHOLD_MULTIPLIER
                              );
    } else {
        context->threshold = 1;
    }
}
//-----------------------------------------------------------------------------
static void update_current_polynomial(siqs_context_t* const context) {
    //
    // Update the SIQS context to use the next polynomial.
    //
    if (context->polynomial_number != context->nb_poly_same_a) {
        //
        // We can keep the coefficient 'a' of the current polynomial and
        // just use another 'b': this is quite fast.
        //
        init_next_polynomial(context);
        context->polynomial_number++;

    } else {
        //
        // All of the suitable 'b' values have been used to construct
        // polynomials with the current 'a'. Determine another value for
        // the coefficient 'a' and proceed to a complete polynomial
        // initialization. This will obviously be costlier than the previous
        // case...
        //
        determine_next_a(context);
        init_first_polynomial(context);
        context->polynomial_number = 1;
    }
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                        Syntaxic sugar functions
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
inline static void fill_sieve_chunk(siqs_context_t* const context) {
    fill_sieve_array(
        context->sieve,
        context->sieve_begin,
        context->sieve_end,
        context->factor_base,
        context->log_factor_base,
        context->sol1,
        context->sol2,
        context->imin,
        context->imax
    );
}
//-----------------------------------------------------------------------------
inline static void scan_sieve_chunk(siqs_context_t* const context) {

    scan_sieve_array(
        context->xpool,
        context->sieve,
        context->sieve_begin,
        context->params->sieve_half_width,
        &context->scan_begin,
        context->threshold
    );
}
//-----------------------------------------------------------------------------
inline static void
keep_relations_with_smooth_gx(siqs_context_t* const context) {

    if (context->params->use_large_primes) {
        bern_21_rt_pairs_lp_siqs(
            context->n,
            context->htable,
            context->u,
            context->smooth_redgx,
            context->a_for_smooth_redgx,
            context->cand_u,
            context->cand_redgx,
            context->cand_a_array,
            context->ptree->data[0]
        );
    } else {
        bern_21_rt_pairs_siqs(
            context->u,
            context->smooth_redgx,
            context->a_for_smooth_redgx,
            context->cand_u,
            context->cand_redgx,
            context->cand_a_array,
            context->ptree->data[0]
        );
    }
}
//-----------------------------------------------------------------------------

