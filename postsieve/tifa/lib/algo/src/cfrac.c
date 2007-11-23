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
 * \file    cfrac.c
 * \author  Jerome Milan
 * \date    Early November 2007
 * \version 1.2
 */

/*
 * History:
 *
 * 1.2: Early November 2007 by JM
 *      - Now uses a smooth_filter_t object to perform smooth residue
 *        selection with the choice of trial division, trial division +
 *        early abort or batch.
 * 1.1: Early March 2007 by JM
 *      - Completely refactored to use a factoring_machine_t.
 * 1.0: Tue Mar 14 2006 by JM
 *      - Initial version.
 */

#include <stdlib.h>
#include <inttypes.h>

#include "tifa_config.h"
#include "gmp_utils.h"
#include "first_primes.h"
#include "matrix.h"
#include "x_array_list.h"
#include "funcs.h"
#include "hashtable.h"
#include "sqrt_cont_frac.h"
#include "bernsteinisms.h"
#include "x_tree.h"
#include "macros.h"
#include "smooth_filter.h"
#include "cfrac.h"

#define __PREFIX__      "cfrac: "
#define __VERBOSE__     TIFA_VERBOSE_CFRAC
#define __TIMING__      TIFA_TIMING_CFRAC

#include "messages.h"

//------------------------------------------------------------------------------
//                         NON PUBLIC DEFINE(S)
//------------------------------------------------------------------------------

//
// Before collecting residues from the continued fraction expansion, computes
// MIN_NSTEP_CFRAC terms (Warning: should be higher than or equal to 2).
//
#define MIN_NSTEP_CFRAC 4
//
// In order to not degrade too much the performance of the hashtable used
// for the large prime variation (due to too many collisions), we set its size
// to the size of the factor base multiplied by HTABLE_SIZE_MULTIPLIER.
//
#define HTABLE_SIZE_MULTIPLIER 8
//
// Maximum number of times to try to get more relations to find a factor
// before giving up and starting again with another multiplier.
//
#define MAX_NFIND_MORE_RELS 3
//
// Maximum number of multiplier to use sequentially. If no factor is found
// after using that many multipliers, just give up the factorization and admit
// defeat.
//
#define MAX_NMULT_USED 1
//
// Number of ranges for the size of the number to factor. Each range is
// define by the optimal value of the factor base's size to use.
//
#define NRANGES 33
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                   NON PUBLIC STRUCTURE(S) AND TYPEDEF(S)
//------------------------------------------------------------------------------

//
// An ad-hoc structure holding the data variables used by the CFRAC
// implementation.
//
struct struct_cfrac_context_t {
    //
    // Number to factor.
    //
    mpz_t n;
    //
    // Data used in the determination of the multiplier to use.
    //
    mult_data_t* mult_data;
    //
    // Multiplier to use.
    //
    uint32_t multiplier;
    //
    // Index of multiplier to use in the mult_data array.
    //
    uint32_t multiplier_index;
    //
    // Number to factor x multiplier.
    //
    mpz_t kn;
    //
    // State of continued fraction generator.
    //
    cont_frac_state_t* cfstate;
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
    binary_matrix_t* matrix;
    //
    // Current number of matrix rows to generate.
    //
    uint32_t nrows_to_collect;
    //
    // Cofactors of generated residues after trial divisions.
    //
    mpz_array_t* partial_yi_array;
    //
    // Relations to find: xi^2 = yi (mod kn) with yi smooth.
    // These are the yi to check for smoothness together with their
    // corresponding xi.
    //
    mpz_array_t* candidate_yi;
    mpz_array_t* candidate_xi;
    //
    // All accepted relations, i.e. the yi and their corresponding xi, such that
    // xi^2 = yi (mod kn) with yi smooth.
    //
    mpz_array_t* all_accepted_yi;
    mpz_array_t* all_accepted_xi;
    //
    // Last accepted relations collected. Actually used as pointers to
    // sub arrays of all_accepted_yi and all_accepted_xi.
    //
    mpz_array_t new_accepted_yi;
    mpz_array_t new_accepted_xi;
    //
    // Factor base used.
    //
    uint32_array_t* factor_base;
    //
    // Number of time additionnal relations have been seached for.
    //
    uint32_t nfind_more_rels;
    //
    // Number of multiplier used.
    //
    uint32_t nmult_used;
    //
    // Hashtable used in large primes variation. It will hold mpz_pair_t's
    // such as (xi, yi) with an index in first_primes_array as key.
    //
    hashtable_t* htable;
    //
    // Smoothness filter with multi-step early abort.
    //
    smooth_filter_t filter;
    //
    // Product tree of primes in factor base.
    //
    mpz_tree_t* ptree;
    //
    // true if large prime variation is used.
    //
    bool use_large_primes;
};
//------------------------------------------------------------------------------
typedef struct struct_cfrac_context_t cfrac_context_t;
//------------------------------------------------------------------------------
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
   // These optimal values were computed by:
   //    - experimental determination for some given size of numbers
   //    - linear interpolation between experimental values (a linear
   //      interpolation is actually not the best choice since it tends to over
   //      estimate the optimal size but it's good enough for the time being)
   //
   {.min_size_n = 0	 , .size_base = 64   },
   {.min_size_n = 15 , .size_base = 80   },
   {.min_size_n = 30 , .size_base = 96   },
   {.min_size_n = 45 , .size_base = 112  },
   {.min_size_n = 60 , .size_base = 128  }, // Empirical determination
   {.min_size_n = 68 , .size_base = 160  },
   {.min_size_n = 76 , .size_base = 192  },
   {.min_size_n = 84 , .size_base = 224  },
   {.min_size_n = 92 , .size_base = 256  }, // Empirical determination
   {.min_size_n = 96 , .size_base = 320  },
   {.min_size_n = 101, .size_base = 384  },
   {.min_size_n = 105, .size_base = 448  },
   {.min_size_n = 110, .size_base = 512  }, // Empirical determination
   {.min_size_n = 114, .size_base = 640  },
   {.min_size_n = 118, .size_base = 768  },
   {.min_size_n = 121, .size_base = 896  },
   {.min_size_n = 125, .size_base = 1024 }, // Empirical determination
   {.min_size_n = 130, .size_base = 1280 },
   {.min_size_n = 135, .size_base = 1536 },
   {.min_size_n = 140, .size_base = 1792 },
   {.min_size_n = 145, .size_base = 2048 }, // Empirical determination
   {.min_size_n = 152, .size_base = 2559 },
   {.min_size_n = 159, .size_base = 3070 },
   {.min_size_n = 165, .size_base = 3581 },
   {.min_size_n = 172, .size_base = 4092 }, // Empirical determination
   {.min_size_n = 181, .size_base = 5117 },
   {.min_size_n = 191, .size_base = 6142 },
   {.min_size_n = 200, .size_base = 7167 },
   {.min_size_n = 210, .size_base = 8192 }, // Empirical determination
   {.min_size_n = 222, .size_base = 10240},
   {.min_size_n = 235, .size_base = 12288},
   {.min_size_n = 247, .size_base = 14336},
   {.min_size_n = 260, .size_base = 16384}
};
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                 PROTOTYPES OF NON PUBLIC FUNCTION(S)
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
static ecode_t init_cfrac_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t clear_cfrac_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t update_cfrac_context(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t perform_cfrac(factoring_machine_t* const machine);
//------------------------------------------------------------------------------
static ecode_t recurse(
    mpz_array_t* const,
    uint32_array_t* const,
    const mpz_t,
    factoring_mode_t
);
//------------------------------------------------------------------------------
static void collect_xi_yi_pairs(cfrac_context_t* context);
//------------------------------------------------------------------------------
static void generate_xi_yi_pairs(
    cont_frac_state_t* const cfstate,
    mpz_array_t* const xi_array,
    mpz_array_t* const yi_array,
    uint32_t npairs
);
//------------------------------------------------------------------------------
inline static uint32_t select_xi_yi_pairs(cfrac_context_t* context);
//------------------------------------------------------------------------------
static void compute_factor_base(uint32_array_t* factor_base, const mpz_t kn);
//------------------------------------------------------------------------------
static void sort_multipliers(cfrac_context_t* const context);
//------------------------------------------------------------------------------
static void choose_multiplier(cfrac_context_t* const context);
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//                 "PUBLIC" FUNCTION(S) IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
ecode_t cfrac(mpz_array_t* const factors, uint32_array_t* const multis,
              const mpz_t n, const cfrac_params_t* const params,
              const factoring_mode_t mode) {

    PRINT_INTRO_MSG("continued fraction");
    INIT_TIMER;
    START_TIMER;

    ecode_t ecode = NO_FACTOR_FOUND;

    factoring_machine_t machine;

    mpz_init_set(machine.n, n);
    machine.mode                = mode;
    machine.params              = (void*) params;
    machine.init_context_func   = init_cfrac_context;
    machine.perform_algo_func   = perform_cfrac;
    machine.update_context_func = update_cfrac_context;
    machine.clear_context_func  = clear_cfrac_context;
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
void set_cfrac_params_to_default(const mpz_t n, cfrac_params_t* const params) {
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
    params->nrelations         = CFRAC_DFLT_NRELATIONS;
    params->linalg_method      = CFRAC_DFLT_LINALG_METHOD;
    params->use_large_primes   = CFRAC_DFLT_USE_LARGE_PRIMES;
    params->filter_method      = DJB_BATCH;
    params->nsteps_early_abort = 0;
    
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
static ecode_t init_cfrac_context(factoring_machine_t* const machine) {
    //
    // Initialize the CFRAC implementation specific variables...
    //
    PRINT_INIT_MSG;
    INIT_TIMER;
    START_TIMER;

    machine->context         = (void*) malloc(sizeof(cfrac_context_t));
    cfrac_context_t* context = (cfrac_context_t*) machine->context;
    cfrac_params_t*  params  = (cfrac_params_t*)  machine->params;

    mpz_init_set(context->n, machine->n);
    mpz_init2(
        context->kn,
        mpz_sizeinbase(context->n, 2) + BITSIZE_LARGEST_MULTIPLIER
    );

    context->mult_data     = malloc(LARGEST_MULTIPLIER * sizeof(mult_data_t));
    context->factor_base   = alloc_uint32_array(params->nprimes_in_base);
    //
    // Number of columns is indeed params->nprimes_in_base + 1 since the first
    // row is used to keep track of the sign of the residue.
    //
    context->matrix_ncols  = params->nprimes_in_base + 1;
    //
    // Number of rows is context->matrix_ncols + MAX_NFIND_MORE_RELS *
    // params->nrelations. In other words, allocate _now_ extra memory space
    // in case the function fails and we want to try to find more relations
    // after updating the context.
    //
    context->matrix_nrows  = context->matrix_ncols;
    context->matrix_nrows += MAX_NFIND_MORE_RELS * params->nrelations;
    context->matrix        = alloc_binary_matrix(
                                 context->matrix_nrows,
                                 context->matrix_ncols
                             );
    context->matrix->nrows = 0;
    //
    // For the first run, we only want to find context->matrix_ncols +
    // params->nrelations relations xi^2 = yi (mod kn).
    //
    context->nrows_to_collect = context->matrix_ncols + params->nrelations;
    context->all_accepted_xi  = alloc_mpz_array(context->matrix_nrows);
    context->all_accepted_yi  = alloc_mpz_array(context->matrix_nrows);
    context->partial_yi_array = alloc_mpz_array(context->matrix_nrows);
    context->multiplier       = 0;
    context->multiplier_index = 0;
    context->nfind_more_rels  = 0;
    context->nmult_used       = 0;
    context->use_large_primes = params->use_large_primes;
    //
    // Sort the multipliers and select the better one.
    //
    sort_multipliers(context);
    choose_multiplier(context);
    //
    // Compute the factor base using the found multiplier.
    //
    compute_factor_base(context->factor_base, context->kn);
    //
    // Compute the product tree of all primes in the base once and for all...
    //
    context->ptree = prod_tree_ui(context->factor_base);
    //
    // Init state of continued fraction generator.
    //
    context->cfstate = malloc(sizeof(cont_frac_state_t));

	init_cont_frac_state(context->cfstate, context->kn);
	step_cont_frac_state(context->cfstate, MIN_NSTEP_CFRAC - 2);

    context->htable = NULL;
    if (params->use_large_primes) {
        uint32_t size = context->all_accepted_yi->alloced;
        size *= HTABLE_SIZE_MULTIPLIER;
        context->htable = alloc_init_hashtable(
                              size, uint32_cmp_func, hash_rj_32
                          );
    }
    // const uint32_t e = ceil_log2(context->nrows_to_collect);
    //
    // Uses smaller batches for the time being until we have better studied
    // the impact of the batch sizes.
    //
    const uint32_t e = 6;

    context->candidate_yi = alloc_mpz_array(1 << e);
    context->candidate_xi = alloc_mpz_array(1 << e);
    //
    // _TO_DO_: Since the maximum size of the yi is known we can actually
    //          allocate all the space needed by hand with only one call to
    //          malloc. Fix this is next code revision...
    //
    for (uint32_t i = 0; i < context->candidate_yi->alloced; i++) {
        mpz_init(context->candidate_yi->data[i]);
        mpz_init(context->candidate_xi->data[i]);
    }
    //
    // Init filter used for smoothness detection.
    //
    context->filter.n             = context->n;
    context->filter.kn            = context->kn;
    context->filter.batch_size    = 1 << e;
    context->filter.nsteps        = params->nsteps_early_abort;
    context->filter.htable        = context->htable;
    context->filter.base_size     = context->factor_base->length;
    context->filter.candidate_xi  = context->candidate_xi;
    context->filter.candidate_yi  = context->candidate_yi;
    context->filter.accepted_xi   = &(context->new_accepted_xi);
    context->filter.accepted_yi   = &(context->new_accepted_yi);
    context->filter.complete_base = context->factor_base;
    
    context->filter.method = params->filter_method;
    context->filter.nsteps = params->nsteps_early_abort;
    
    //
    // _WARNING_: Large prime variation is mandatory for the time being...
    //
    //context->filter.use_large_primes = params->use_large_primes;
    context->filter.use_large_primes = true;
    context->filter.use_siqs_variant = false;

    complete_filter_init(&(context->filter), context->factor_base);

    STOP_TIMER;
    PRINT_TIMING;

    return SUCCESS;
}
//------------------------------------------------------------------------------
static ecode_t clear_cfrac_context(factoring_machine_t* const machine) {

    PRINT_CLEAN_MSG;
    INIT_TIMER;
    START_TIMER;

    cfrac_context_t* context = (cfrac_context_t*) machine->context;

    if (context != NULL) {
        mpz_clear(context->n);
        mpz_clear(context->kn);
        //
        // Since the candidate_* arrays were fully initialized, don't forget to
        // set their current length to their alloced length before clearing
        // them...
        //
        context->candidate_yi->length = context->candidate_yi->alloced;
        context->candidate_xi->length = context->candidate_xi->alloced;

        clear_mpz_array(context->all_accepted_yi);
        free(context->all_accepted_yi);

        clear_mpz_array(context->all_accepted_xi);
        free(context->all_accepted_xi);

        clear_mpz_array(context->candidate_yi);
        free(context->candidate_yi);

        clear_mpz_array(context->candidate_xi);
        free(context->candidate_xi);

        clear_uint32_array(context->factor_base);
        free(context->factor_base);

        clear_binary_matrix(context->matrix);
        free(context->matrix);

        clear_mpz_array(context->partial_yi_array);
        free(context->partial_yi_array);

        clear_cont_frac_state(context->cfstate);
        free(context->cfstate);

        clear_mpz_tree(context->ptree);
        free(context->ptree);

        free(context->mult_data);
        
        clear_smooth_filter(&(context->filter));

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
static ecode_t update_cfrac_context(factoring_machine_t* const machine) {

    //
    // If no factors (or not enough factors) are found, update the CFRAC
    // context according to the following strategy:
    //
    //    - First, try to find more relations by generating new residues.
    //    - If, after MAX_NFIND_MORE_RELS such tries, no extra factors are
    //      found, then, change the multiplier used and restart the whole
    //      process.
    //    - If after MAX_NMULT_USED used multipliers, no factor is found, just
    //      give up. At this stage it is unlikely that we'll find a factor
    //      anyway, so no need to spend too much time...
    //

    INIT_TIMER;
    START_TIMER;

    cfrac_context_t* context = (cfrac_context_t*) machine->context;
    cfrac_params_t*  params  = (cfrac_params_t*)  machine->params;

    ecode_t ecode = SUCCESS;

    if (context->nfind_more_rels == MAX_NFIND_MORE_RELS) {
        if (context->nmult_used == MAX_NMULT_USED) {
            PRINT_UPDATE_GIVEUP_MSG;
            //
            // At that point, it seems unlikely that we'll succeed finding
            // a factor, so let's just give up the factorization...
            //
            ecode = GIVING_UP;
        } else {
            PRINT_UPDATE_NEW_MULT_MSG;
            //
            // We'll try with another multiplier.
            //
            // _KLUDGE_: Right now, we just forget everything we've done so far
            //           and prepare to do everything from scratch with a new
            //           multiplier. Actually we could also keep the relations
            //           found but we'll have to change the base (various
            //           strategies are possible) and update the matrix
            //           accordingly. This is quite cumbersome and would make
            //           for a rather complex code.
            //
            // _TO_DO_: Find about the aforementionned problem, and try to fix
            //          it in a not-that-ugly way.
            //
            context->nfind_more_rels  = 0;

            reset_binary_matrix(context->matrix);
            context->matrix->nrows    = 0;
            context->nrows_to_collect =   context->matrix_ncols
                                        + params->nrelations;

            //
            // The candidate multipliers were already sorted. Just use the next
            // best one.
            //
            choose_multiplier(context);
            compute_factor_base(context->factor_base, context->kn);

            clear_mpz_tree(context->ptree);
            free(context->ptree);
            context->ptree = prod_tree_ui(context->factor_base);

            clear_cont_frac_state(context->cfstate);
        	init_cont_frac_state(context->cfstate, context->kn);
        	step_cont_frac_state(context->cfstate, MIN_NSTEP_CFRAC - 2);
            //
            // We forget everything done so far. A shame, really...
            //
            uint32_t len = context->all_accepted_xi->length;

            for (uint32_t i = 0U; i != len; i++) {
                if (context->all_accepted_xi->data[i] != NULL) {
                    mpz_clear(context->all_accepted_xi->data[i]);
                }
                if (context->all_accepted_yi->data[i] != NULL) {
                    mpz_clear(context->all_accepted_yi->data[i]);
                }
            }
            context->all_accepted_yi->length = 0;
            context->all_accepted_xi->length = 0;
        }
    } else {
        PRINT_UPDATE_MORE_RELS_MSG;
        //
        // We'll try to find params->nrelations more relations using the same
        // multiplier
        //
        context->nfind_more_rels++;
        context->nrows_to_collect = params->nrelations;
    }

    STOP_TIMER;
    PRINT_TIMING;

    return ecode;
}
//------------------------------------------------------------------------------
static ecode_t perform_cfrac(factoring_machine_t* const machine) {
    //
    // Perform the CFRAC algorithm using the current CFRAC context.
    //
    cfrac_context_t* context = (cfrac_context_t*) machine->context;
    cfrac_params_t*  params  = (cfrac_params_t*)  machine->params;

    INIT_TIMER;
    PRINT_COLLECT_RELS_MSG;
    START_TIMER;
    //
    // Collect xi^2 = yi (mod kn) relations with yi smooth
    //
    collect_xi_yi_pairs(context);

    STOP_TIMER;
    PRINT_COLLECT_RELS_DONE_MSG;
    PRINT_TIMING;
    PRINT_FACTOR_RES_MSG;
    RESET_TIMER;
    START_TIMER;
    
    uint32_array_t primes_array;
    primes_array.alloced = 0;
    primes_array.length  = params->nprimes_tdiv;
    primes_array.data    = (uint32_t*)(context->factor_base->data);
    
    //
    // Fill the matrix: first by trial divisions...
    //
    fill_matrix_trial_div(
        context->matrix,
        context->partial_yi_array,
        &(context->new_accepted_yi),
        &primes_array
    );

    if (params->nprimes_tdiv < params->nprimes_in_base) {
        //
        // Finish to fill the matrix with one of Bernstein's batch algo...
        //
        uint32_array_list_t*
        decomp_list = alloc_uint32_array_list(context->all_accepted_yi->length);

        primes_array.length = params->nprimes_in_base - params->nprimes_tdiv;
        primes_array.data   = (uint32_t*) &(
                                context->factor_base->data[params->nprimes_tdiv]
                              );
        
        bern_71(decomp_list, context->partial_yi_array, &primes_array);
        
        fill_matrix_from_list(
            context->matrix,
            context->partial_yi_array,
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
    uint32_array_list_t* relations;
    
    if (machine->mode == SINGLE_RUN) {
        relations = find_dependencies(context->matrix, params->linalg_method);
    } else {
        //
        // We have to clone the matrix and perform the resolution on the copy
        // since we might need the matrix in its original form should we have
        // to add more relations if no factors are found...
        //
        binary_matrix_t* clone = clone_binary_matrix(context->matrix);
        
        relations = find_dependencies(clone, params->linalg_method);
        
        clear_binary_matrix(clone);
        free(clone);
    }
    STOP_TIMER;
    PRINT_TIMING;
    PRINT_DED_FACTORS_MSG;
    RESET_TIMER;
    START_TIMER;

    //
    // Deduce factors of n using the found relations.
    //
    ecode_t ecode = find_factors(
                        machine->factors,
                        context->n,
                        context->all_accepted_xi,
                        context->all_accepted_yi,
                        relations
                    );
    clear_uint32_array_list(relations);
    free(relations);

    STOP_TIMER;
    PRINT_TIMING;

    return ecode;
}
//------------------------------------------------------------------------------
static ecode_t recurse(mpz_array_t* const factors, uint32_array_t* const multis,
                       const mpz_t n, factoring_mode_t mode) {

    cfrac_params_t params;
    set_cfrac_params_to_default(n, &params);

    return cfrac(factors, multis, n, &params, mode);
}
//------------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                  Residues generation and selection
 *-----------------------------------------------------------------------------
 */

//------------------------------------------------------------------------------
static void collect_xi_yi_pairs(cfrac_context_t* context) {    
    //
    // Collect relations xi^2 = yi (mod kn) with yi smooth.
    //
    INIT_NAMED_TIMER(generation);
    INIT_NAMED_TIMER(selection);
    
    //
    // We have already accepted 'accepted' relations. (Used if the CFRAC
    // context was updated after a failure to try to find more relations)
    //
    uint32_t accepted = context->all_accepted_yi->length;
    //
    // Make new_accepted_yi and new_accepted_xi point to subarrays of
    // all_accepted_yi and all_accepted_xi.
    //
    context->new_accepted_yi.length = 0;
    context->new_accepted_xi.length = 0;
    //
    // _NOTE_: context->nrows_to_collect is updated by the
    //         update_cfrac_context() function. Consequently, it is not
    //         necessarily the number of rows of the matrix.
    //
    context->new_accepted_yi.alloced = context->nrows_to_collect;
    context->new_accepted_xi.alloced = context->nrows_to_collect;

    context->new_accepted_yi.data = &(context->all_accepted_yi->data[accepted]);
    context->new_accepted_xi.data = &(context->all_accepted_xi->data[accepted]);

    while (context->new_accepted_yi.length != context->nrows_to_collect) {
        //
        // If context->candidate_yi->length != 0, it means that not all
        // relations were used to reach the amount needed to go to the
        // linear algebra phase.
        //
        // However if the context needs to be updated later on (because none of 
        // the found relations yields a non trivial factor) we'll come back 
        // in this loop to collect more residues.
        //
        START_NAMED_TIMER(generation);
        //
        // Compute candidate yi and xi such that yi = xi^2 (mod kn)...
        //
        generate_xi_yi_pairs(
            context->cfstate,
            context->candidate_xi,
            context->candidate_yi,
            1 << ceil_log2(context->nrows_to_collect)
        );

        STOP_NAMED_TIMER(generation);
        START_NAMED_TIMER(selection);

        //
        // ... and keep only pairs for which yi is smooth and pairs found
        // via the large prime variation (if this variation is used).
        //        
        filter_new_relations(&(context->filter));

        STOP_NAMED_TIMER(selection);

        PRINT_NRELS_FOUND(
            context->new_accepted_yi.length,
            context->nrows_to_collect
        );
    }
    
    PRINT_RES_GENERATED_MSG(GET_NAMED_TIMING(generation));
    PRINT_RES_SELECTED_MSG(GET_NAMED_TIMING(selection));

    //
    // We now have all the pairs (xi, yi) we need...
    //
    context->all_accepted_xi->length += context->nrows_to_collect;
    context->all_accepted_yi->length += context->nrows_to_collect;
}
//------------------------------------------------------------------------------
static void generate_xi_yi_pairs(
                        cont_frac_state_t* const cfstate,
                        mpz_array_t* const xi_array,
                        mpz_array_t* const yi_array,
                        uint32_t npairs) {

    //
    // Generate npairs relations xi^2 = yi (mod kn) and append them in
    // the xi_array and yi_array.
    //
    uint32_t const free_space = xi_array->alloced - xi_array->length;

    npairs = MIN(npairs, free_space);

    uint32_t const from = xi_array->length;
    uint32_t const to   = xi_array->length + npairs;
    uint32_t       i    = from;

    while (i < to) {

		step_cont_frac_state(cfstate, 1);
		//
        // Now, (-1)^m * cfstate.q = cfstate.a^2 - k*N*cfstate.b^2
        // i.e. cfstate.a^2 = (-1)^m cfstate.q (mod N)
        //
        // with m = 1 if (cfstate.nsteps_performed) is odd
        //          0 otherwise
        //
        if (mpz_cmpabs_ui(cfstate->q, first_primes[0]) >= 0) {

        	if (IS_ODD(cfstate->nsteps_performed)) {
        	    mpz_neg(yi_array->data[i], cfstate->q);
            } else {
                mpz_set(yi_array->data[i], cfstate->q);
            }
            //
            // Uncomment this 'if' loop to keep only positive residues...
            //
            //if (mpz_sgn(yi_array->data[i]) == 1) {

        	mpz_set(xi_array->data[i], cfstate->a);
        	i++;

            //}
		}
	}
	xi_array->length += npairs;
	yi_array->length += npairs;
}
//------------------------------------------------------------------------------
inline static uint32_t select_xi_yi_pairs(cfrac_context_t* context) {
    //
    // Select relations xi^2 = yi (mod kn) from the candidate_xi and
    // candidate_yi arrays so that yi is smooth.
    //
    if (context->use_large_primes) {
        return bern_21_rt_pairs_lp(
            context->n,
            context->htable,
            &(context->new_accepted_xi),
            &(context->new_accepted_yi),
            context->candidate_xi,
            context->candidate_yi,
            context->ptree->data[0]
        );
    } else {
        return bern_21_rt_pairs(
            &(context->new_accepted_xi),
            &(context->new_accepted_yi),
            context->candidate_xi,
            context->candidate_yi,
            context->ptree->data[0]
        );
    }
}
//------------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *             Determination of best multiplier and factor base
 *-----------------------------------------------------------------------------
 */

//------------------------------------------------------------------------------
static void sort_multipliers(cfrac_context_t* const context) {
    //
    // Sorts the multipliers to use for the CFRAC algorithm.
    //
    // _NOTE_: We follow the suggestion of  M. A. Morrison and J. Brillhart
    //         in the remark 5.3 of the paper "A Method of Factoring and the
    //         Factorization of F_7" (Mathematics of Computation, vol 29,
    //         #129, Jan 1975, pages 183-205)
    //
    mult_data_t* mdata = context->mult_data;

    mpz_t kn;
    mpz_init2(kn, mpz_sizeinbase(context->n, 2) + BITSIZE_LARGEST_MULTIPLIER);

    for (uint32_t k = 0; k < LARGEST_MULTIPLIER; k++) {

        mdata[k].multiplier = k + 1;
        mdata[k].count      = 0;
        mdata[k].sum_inv_pi = 0.0;

        mpz_mul_ui(kn, context->n, mdata[k].multiplier);

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
    mpz_clear(kn);
    //
    // The cmp_mult_data function gives the meaning of what a "better"
    // multiplier is.
    //
    qsort(mdata, LARGEST_MULTIPLIER, sizeof(mult_data_t), cmp_mult_data);
}
//------------------------------------------------------------------------------
static void choose_multiplier(cfrac_context_t* const context) {
    //
    // Chooses the multiplier to use for the CFRAC algorithm.
    //
    // _NOTE_: We should also check that k*n is NOT a square, otherwise
    //         the CFRAC algorithm will fail. Such a check may not be needed
    //         for rather "large" number since the probability for k*n to be
    //         square is quite small. However for really tiny integers, such a
    //         case can occur more often than desired...
    //
    uint32_t i = LARGEST_MULTIPLIER - 1;
    if (context->multiplier != 0) {
        //
        // context->multiplier != 0 implies that this is not the first
        // multiplier we're using. So use the best multiplier coming after the
        // ones already used.
        //
        i = context->multiplier_index - 1;
    }

    mpz_mul_ui(context->kn, context->n, context->mult_data[i].multiplier);

    while ((0 != mpz_perfect_square_p(context->kn)) && (i != 0)) {
        i--;
        mpz_mul_ui(context->kn, context->n, context->mult_data[i].multiplier);
    }
    context->multiplier       = context->mult_data[i].multiplier;
    context->multiplier_index = i;
    context->nmult_used++;
}
//------------------------------------------------------------------------------
static void compute_factor_base(uint32_array_t* factor_base, const mpz_t kn) {
    //
    // Initializes the factor base by keeping only those primes pi such that
    // the legendre symbol (kn/pi) is either 0 or 1. Note that 2 is always
    // included in the factor base.
    //
    factor_base->data[0] = 2;
    factor_base->length  = 1;

    uint32_t pindex = 1;

    while (   (factor_base->length != factor_base->alloced)
           && (pindex != NFIRST_PRIMES)
          ) {
        if (-1 != mpz_kronecker_ui(kn, first_primes[pindex])) {
            factor_base->data[factor_base->length] = first_primes[pindex];
            factor_base->length++;
        }
        pindex++;
    }
}
//------------------------------------------------------------------------------
