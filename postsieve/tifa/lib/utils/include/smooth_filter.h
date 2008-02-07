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
 * \file    smooth_filter.h
 * \author  Jerome Milan
 * \date    Fri Oct 12 2007
 * \version 1.0
 *
 * \brief Smooth integer filter.
 *
 * The \c smooth_filter_t structure and its associated functions implement the
 * multi-step early abort strategy in a way reminiscent of Pomerance's 
 * suggestion in "Analysis and Comparison of Some Integer Factoring Algorithm" 
 * with the exception that the smoothness tests are performed by batch 
 * (see bernsteinisms.h) instead of trial division.
 * 
 * <h3>How to use a \c smooth_filter_t structure?</h3>
 *
 * The following code snippet, while incomplete, illustrates the way a
 * \c smooth_filter_t should be used.
 *
 * \verbatim
  //
  // Fill the with the smooth_filter_t with our parameters...
  //
  filter.n          = n;    // number to factor
  filter.kn         = kn;   // number to factor x multiplier
  filter.batch_size = 1024; // number of relations to perform a batch
  filter.methid     = TDIV_EARLY_ABORT; // use the early abort strategy
  filter.nsteps     = 2;    // use a 2-step early abort strategy
  
  filter.htable                 = htable;
  filter.use_large_primes       = true;
  filter.use_siqs_batch_variant = false;
  
  filter.base_size    = factor_base->length; // size of factor base
  filter.candidate_xi = candidate_xi; // candidate relations stored
  filter.candidate_yi = candidate_yi; // in candidate_* arrays
  filter.accepted_xi  = accepted_xi;  // 'good' relations stored
  filter.accepted_yi  = accepted_yi;  // in accepted_* arrays
  //
  // Complete the filter initialization by allocating its internal
  // buffers, computing the early abort bounds and the intermediate
  // factor bases...
  //
  complete_filter_init(&filter, factor_base);
  
  while (accepted_yi->length != nrels_to_collect) {
      //
      // While we don't have enough relations, create new ones and
      // stores them in the candidate_* arrays such that
      // yi = xi^2 (mod kn). (The generate_relations function here
      // is completely fictitious)
      //
      generate_relations(candidate_xi, candidate_yi);
      
      //
      // Select 'good' relations such that yi = xi^2 (mod kn) with
      // yi smooth. Note that pointers to the candidate_* and
      // accepted_* arrays were given in the filter structure.
      //
      filter_new_relations(&filter);
  }
 \endverbatim
 *
 */

#if !defined(_TIFA_SMOOTH_FILTER_H_)
   /**
    * \def _TIFA_SMOOTH_FILTER_H_
    * Standard include guard.
    */
#define _TIFA_SMOOTH_FILTER_H_

#include <inttypes.h>
#include <array.h>
#include <hashtable.h>
#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

   /**
    * \def MAX_NSTEPS
    * Maximum number of steps used in the multi-step early abort strategy.
    */
#define MAX_NSTEPS 4

   /**
    * \enum smooth_filter_method_enum
    *
    * An enumeration of the possible methods used to test residue smoothness.
    */
enum smooth_filter_method_enum {
        /**
         * Simple trial division.
         */
    TDIV = 0,
        /**
         * Simple trial division with (multi-step) early abort.
         */
    TDIV_EARLY_ABORT,
        /**
         * D. Bernstein's batch method described in "How to find smooth parts
         * of integers", http://cr.yp.to/factorization/smoothparts-20040510.pdf.
         */
    DJB_BATCH,
};

   /**
    * \typedef smooth_filter_t
    * \brief Equivalent to <tt>struct struct_smooth_filter_t</tt>.
    */
typedef struct struct_smooth_filter_t smooth_filter_t;

   /**
    * Global constant array mapping filter methods to their string 
    * representations.
    */
static const char* const filter_method_to_str[3] = {
    "trial division",
    "trial division + early abort",
    "batch",
};

   /**
    * \typedef smooth_filter_method_t
    * \brief Equivalent to <tt>enum smooth_filter_method_enum</tt>.
    */
typedef enum smooth_filter_method_enum smooth_filter_method_t;

   /**
    * \struct struct_smooth_filter_t smooth_filter.h
    *                                lib/utils/include/smooth_filter.h
    *
    * \brief Structure grouping variables needed for multi-step early abort
    * strategy.
    *
    * This structure and its associated functions implement the multi-step
    * early abort strategy in a way reminiscent of Pomerance's suggestion
    * in "Analysis and Comparison of Some Integer Factoring Algorithm" with
    * the exception that the smoothness tests are performed by batch instead
    * of trial division.
    *
    * \See C. Pomerance, <i>Analysis and Comparison of Some Integer Factoring
    * Algorithm</i>, in Mathematical Centre Tracts 154.
    */
struct struct_smooth_filter_t {
       /**
        * The number to factor.
        */
    mpz_ptr n;
       /**
        * The number to factor multiplied by a multiplier.
        */
    mpz_ptr kn;
       /**
        * The method to use for smooth residue detection.
        */
    smooth_filter_method_t method;
       /**
        * The number of steps in the early abort strategy. If <tt>nsteps ==
        * 0</tt> no early abort is performed.
        *
        * \note \c nsteps should be less than or equal to \c MAX_NSTEPS.
        */
    unsigned short int nsteps;
       /**
        * Number of relations to accumulate before testing for smoothness.
        */
    unsigned long int batch_size;
       /**
        * Size of the complete factor base.
        */
    unsigned long int base_size;
       /**
        * Pointer to the complete factor base.
        */
    uint32_array_t* complete_base;
       /**
        * Array giving the factor base to use at each step if we use the
        * early-abort strategy.
        */
    uint32_array_t* factor_base[MAX_NSTEPS];
       /**
        * The candidate x's. Together with \c candidate_yi, stores candidate
        * relations verifying \c x^2 (mod \c kn) == \c y (mod \c kn)
        *
        * \note See \c candidate_ai if \c use_siqs_batch_variant is true.
        * In that case the relations become \c x^2 (mod \c kn) == \c y * \c a 
        * (mod \c kn).
        */
    mpz_array_t* candidate_xi;
       /**
        * The candidate y's. Together with \c candidate_xi, stores candidate
        * relations verifying \c x^2 (mod \c kn) == \c y (mod \c kn)
        *
        * \note See \c candidate_ai if \c use_siqs_batch_variant is true.
        * In that case the relations become \c x^2 (mod \c kn) == \c y * \c a 
        * (mod \c kn).
        */
    mpz_array_t* candidate_yi;
       /**
        * Used only if \c use_siqs_batch_variant is true.
        *
        * \c candidate_ai stores the a's (i.e. the value of the first 
        * parameter of the SIQS polynomial used). Together with \c candidate_xi 
        * and \c candidate_yi, stores candidate relations verifying
        * \c x^2 (mod \c kn) == \c y * \c a (mod \c kn).
        */
    mpz_array_t* candidate_ai;
       /**
        * The accepted x's. Together with \c accepted_yi, stores 'good'
        * relations verifying \c x^2 (mod \c kn) == \c y (mod \c kn) with
        * \c y smooth over the factor base.
        *
        * \note See \c accepted_ai if \c use_siqs_batch_variant is true.
        * In that case the relations become \c x^2 (mod \c kn) == \c y * \c a 
        * (mod \c kn).
        */
    mpz_array_t* accepted_xi;
       /**
        * The accepted y's. Together with \c accepted_xi, stores 'good'
        * relations verifying \c x^2 (mod \c kn) == \c y (mod \c kn) with
        * \c y smooth over the factor base.
        *
        * \note See \c accepted_ai if \c use_siqs_batch_variant is true.
        * In that case the relations become \c x^2 (mod \c kn) == \c y * \c a 
        * (mod \c kn).
        */
    mpz_array_t* accepted_yi;
       /**
        * Used only if \c use_siqs_batch_variant is true.
        *
        * The accepted a's (i.e. the value of the first 
        * parameter of the SIQS polynomial used). Together with \c accepted_xi, 
        * and \c accepted_yi, stores candidate relations verifying
        * \c x^2 (mod \c kn) == \c y * \c a (mod \c kn) with \c y smooth over 
        * the factor base.
        */
    mpz_array_t* accepted_ai; 
       /**
        * The \c filtered_xi[s] array gives the filtered x's after \c s 
        * early abort steps. Together with \c filtered_yi[s] and
        * \c cofactors[s], stores relations verifying \c x^2 (mod \c kn)
        * == \c y * \c cofactor (mod \c kn) with \c cofactor smooth over
        * the base composed by the partial_factor bases used in the
        * early abort steps up to (and including) the \c s-th one, and with
        * \c y less than \c bounds[s].
        *
        * \note See \c filtered_ai if \c use_siqs_batch_variant is true.
        * In that case the relations become \c x^2 (mod \c kn) == \c y * \c a *
        * \c cofactor (mod \c kn).
        */
    mpz_array_t* filtered_xi[MAX_NSTEPS];
       /**
        * The \c filtered_yi[s] array gives the filtered x's after \c s 
        * early abort steps. See filtered_xi.
        */
    mpz_array_t* filtered_yi[MAX_NSTEPS];
       /**
        * The \c filtered_ai[s] array gives the filtered a's after \c s 
        * early abort steps. See filtered_xi.
        */
    mpz_array_t* filtered_ai[MAX_NSTEPS];
       /**
        * The \c cofactor[s] array gives the cofactors of the y's in
        * \c filtered_yi[s] after \c s early abort steps. See filtered_xi.
        */
    mpz_array_t* cofactors  [MAX_NSTEPS];
       /**
        * \c bounds[s] is the upper bound of the \c s-th early abort step.
        * A relation will not pass this step if it has a 'y' greater than this
        * bound.
        */
    mpz_t bounds[MAX_NSTEPS];
       /**
        * \c prod_pj[s] is the product of all the elements in the partial
        * factor base at the \c s-th early abort step.
        */
    mpz_t prod_pj[MAX_NSTEPS + 1];
       /**
        * Hashtable used if the large prime variation is used. Can be \c NULL
        * if the variation is not used.
        */
    hashtable_t* htable;
       /**
        * true if and only if we are using the large prime variation.
        */
    bool use_large_primes;
       /**
        * true if and only if we are using this \c smooth_filter_t with
        * SIQS.
        */
    bool use_siqs_variant;
};

   /**
    * \brief Complete initialization of a \c smooth_filter_t.
    *
    * Complete the initialization of a \c smooth_filter_t by allocating
    * needed memory space.
    *
    * \warning It is the responsability of the client code to set the
    * following structure variables <em>before</em> calling this function:
    *
    *   \li \c n
    *   \li \c kn
    *   \li \c nsteps
    *   \li \c batch_size
    *   \li \c base_size
    *   \li \c htable
    *   \li \c candidate_xi
    *   \li \c candidate_yi
    *   \li \c accepted_xi
    *   \li \c accepted_yi
    *   \li \c use_large_primes
    *   \li \c use_siqs_batch_variant
    *
    * No pointer ownership is transfered. For example, it is still the
    * responsability of the client code to properly clears the \c candidate_*
    * arrays since the structure just <em>refers</em> to them, but does
    * not <em>own</em> them.
    *
    * \note If <tt>nsteps > MAX_NSTEPS</tt> then \c complete_filter_init
    * will set it to \c MAX_NSTEPS.
    *
    * \warning If \c method == \c DJB_BATCH then \c nsteps will be set to 0
    * and no early abort will be performed.
    *
    * \param filter a pointer to the \c smooth_filter_t to initialize
    * \param base   a pointer to the complete factor base used
    */
void complete_filter_init(
    smooth_filter_t* const filter,
    uint32_array_t*  const base
);

   /**
    * \brief Clears a <tt>smooth_filter_t</tt>.
    *
    * Clears the memory of a \c smooth_filter_t that was allocated by
    * \c complete_filter_init().
    *
    * \warning This clears <em>only</em> the internal buffers allocated by
    * by complete_filter_init(), and not the whole structure. For example, it
    * is still the responsability of the client code to properly clears the \c
    * candidate_* arrays or the \c htable hashtable.
    */
void clear_smooth_filter(smooth_filter_t* const filter);

   /**
    * \brief Filters new relations to keep 'good' ones.
    *
    * Filters the relations given by <tt>filter->candidate_*</tt> via
    * a smoothness detecting batch using a <tt>filter->nsteps</tt> steps
    * early abort strategy and stores the 'good' relations in
    * <tt>filter->accepted_*</tt>. Has no effect if the
    * <tt>filter->candidate_*</tt> are not full since we need
    * <tt>filter->batch_size</tt> relations to perform a batch.
    *
    * \param filter a pointer to the \c smooth_filter_t used
    */
void filter_new_relations(smooth_filter_t* const filter);

   /**
    * \brief Prints the status of the buffers of a \c smooth_filter_t.
    *
    * Prints a status summary of the internal buffers of a \c smooth_filter_t
    * on the standard output.
    *
    * \note This function is mostly intended for debugging purposes as the
    * output is not particularly well structured.
    *
    * \param filter a pointer to the \c smooth_filter_t to inspect
    */
void print_filter_status(smooth_filter_t* const filter);

#ifdef __cplusplus
}
#endif

#endif
