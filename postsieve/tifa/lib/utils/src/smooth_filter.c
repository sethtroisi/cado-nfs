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
 * \file    smooth_filter.c
 * \author  Jerome Milan
 * \date    Fri Oct 12 2007
 * \version 1.0
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "tifa_config.h"
#include "x_tree.h"
#include "bernsteinisms.h"
#include "smooth_filter.h"
#include "res_tdiv.h"
#include "gmp_utils.h"
#include "first_primes.h"

//------------------------------------------------------------------------------
uint32_t select_res(smooth_filter_t* const filter, unsigned long int step);
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void complete_filter_init(smooth_filter_t* const filter,
                          uint32_array_t*  const base) {   
    
    switch (filter->method) {
        
    case TDIV:
        filter->nsteps = 0;
        filter->factor_base[0] = base;
        break;
        
    case TDIV_EARLY_ABORT:
        filter->nsteps = 1;
        //
        // Allocate the memory space needed by the arrays.
        //
        for (unsigned int i = 0; i < filter->nsteps; i++) {
            filter->filtered_xi[i] = alloc_mpz_array(filter->batch_size);
            filter->filtered_yi[i] = alloc_mpz_array(filter->batch_size);
            filter->cofactors[i]   = alloc_mpz_array(filter->batch_size);
        }
        //
        // Compute the various bounds and partial factor base needed for the
        // multi-step early abort strategy.
        //
        unsigned long int nprimes = base->length;
        
        uint32_t lp     = base->data[base->length - 1];
        uint32_t rootlp = sqrt(lp);
        
        uint32_t ip  = 0;
        uint32_t eap = base->data[ip];
        
        while (eap < rootlp) {
            ip++;
            eap = base->data[ip];
        }
        filter->factor_base[0]          = malloc(sizeof(uint32_array_t));
        filter->factor_base[0]->data    = base->data;
        filter->factor_base[0]->length  = ip + 1;
        filter->factor_base[0]->alloced = 0;
        
        filter->factor_base[1]          = malloc(sizeof(uint32_array_t));
        filter->factor_base[1]->data    = &(base->data[ip + 1]);
        filter->factor_base[1]->length  = nprimes - ip - 2;
        filter->factor_base[1]->alloced = 0;
        
        mpz_init_set(filter->bounds[0], filter->kn);
        mpz_root(filter->bounds[0], filter->bounds[0], 7);
        mpz_pow_ui(filter->bounds[0], filter->bounds[0], 3);
        
        break;

    case DJB_BATCH:
    default:
        filter->method = DJB_BATCH;
        filter->nsteps = 0;
        mpz_tree_t* ptree = prod_tree_ui(base);
        mpz_init_set(filter->prod_pj[0], ptree->data[0]);
        clear_mpz_tree(ptree);
        break;
    }
}
//------------------------------------------------------------------------------
void clear_smooth_filter(smooth_filter_t* const filter) {
    
    switch (filter->method) {
        
    case TDIV_EARLY_ABORT:
    
        for (unsigned int i = 0; i < filter->nsteps; i++) {
            mpz_t* const xidat = filter->filtered_xi[i]->data;
            mpz_t* const yidat = filter->filtered_yi[i]->data;
            mpz_t* const codat = filter->cofactors  [i]->data;
            
            for (unsigned int j = 0; j < filter->batch_size; j++) {
                mpz_clear(xidat[j]);
                mpz_clear(yidat[j]);
                mpz_clear(codat[j]);
            }
            free(filter->filtered_xi[i]);
            free(filter->filtered_yi[i]);
            free(filter->cofactors[i]);
        }
        if (filter->use_siqs_variant) {
            for (unsigned int i = 0; i < filter->nsteps; i++) {
                mpz_t* const xidat = filter->filtered_ai[i]->data;
                for (unsigned int j = 0; j < filter->batch_size; j++) {
                    mpz_clear(xidat[j]);
                }
                free(filter->filtered_ai[i]);
            }
        }
        for (unsigned int i = 0; i < filter->nsteps; i++) {
            mpz_clear(filter->bounds[i]);
        }
        free(filter->factor_base[0]);
        free(filter->factor_base[1]);
        break;
        
    case DJB_BATCH:
        mpz_clear(filter->prod_pj[0]);
        break;
    
    default:
        break;
    }
}
//------------------------------------------------------------------------------
void filter_new_relations(smooth_filter_t* const filter) {

    uint32_t used = 0;
    unsigned int step = 1;
        
    switch (filter->method) {
    
    case TDIV:
    case DJB_BATCH:
        
        used = select_res(filter, 0);
        
        filter->candidate_xi->length -= used;
        filter->candidate_yi->length -= used;
        
        break;
        
    case TDIV_EARLY_ABORT:    

        used = select_res(filter, 0);
    
        filter->candidate_xi->length -= used;
        filter->candidate_yi->length -= used;

        while (step <= filter->nsteps) {

            if (! ARRAY_IS_FULL(filter->filtered_yi[step - 1])) {
                break;
            }
            if (ARRAY_IS_FULL(filter->accepted_yi)) {
                break;
            }
            used = select_res(filter, step);

            filter->filtered_yi[step - 1]->length -= used;
            filter->filtered_xi[step - 1]->length -= used;
            filter->cofactors  [step - 1]->length -= used;

            step++;
        }
        break;
    
    default:
        break;
    }
}
//------------------------------------------------------------------------------
uint32_t select_res(smooth_filter_t* const filter, unsigned long int step) {
    //
    // Performs a round of smoothness detection with 'filter'.
    //
    // The parameter 'step' is only used if filter->method = TDIV_EARLY_ABORT,
    // in which case we have the following behaviour:
    //
    //   step == 0:
    //     First early abort step. Performs a test where the 'good'
    //     relations given by the filter->candidate_* arrays are stored in
    //     filter->filtered_*[0] and/or filter->accepted_*.
    //
    //   1 <= step < filter->nsteps:
    //     Intermediate steps corresponding to early abort. Performs a test
    //     where the 'good' relations given by the filter->filtered_*[step - 1] 
    //     arrays are stored in filter->filtered_*[step] and/or 
    //     filter->accepted_*.
    //
    //   step == filter->nsteps:
    //     Last step. Performs a test where the 'good' relations from the 
    //     last filter->filtered_* arrays are stored in filter->accepted_*.
    //    
    switch (filter->method) {
    case TDIV:
    case TDIV_EARLY_ABORT: 
        return res_tdiv(filter, step);
        break;
    
    case DJB_BATCH:
        return djb_batch_rt(filter, step);
        break;
       
    default:
        break;
    }
    return 0;
}
//------------------------------------------------------------------------------
void print_filter_status(smooth_filter_t* const filter) {
    //
    // Just for debug purposes...
    //
    printf("  ------------------------------\n");
    printf("     accepted: %6"PRIu32" / %"PRIu32"\n",
        filter->accepted_yi->length, 
        filter->accepted_yi->alloced
    );
    printf("  ------------------------------\n");

    for (short i = filter->nsteps - 1; i >= 0; i--) {
        printf("     filter %1i: %6"PRIu32" / %"PRIu32"\n",
            i, 
            filter->filtered_yi[i]->length, 
            filter->filtered_yi[i]->alloced
        );
        printf("  ------------------------------\n");
    }
    printf("    candidate: %6"PRIu32" / %"PRIu32"\n",
        filter->candidate_yi->length, 
        filter->candidate_yi->alloced
    );
    printf("  ------------------------------\n");
}
//------------------------------------------------------------------------------

