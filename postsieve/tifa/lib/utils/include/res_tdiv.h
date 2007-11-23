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
 * \file    res_tdiv.h
 * \author  Jerome Milan
 * \date    Thu Nov 15 2007
 * \version 1.0
 *
 * \brief Trial division of residues with optional early abort.
 *
 * This file defines functions used to trial divide residues on a factor
 * base using optional multi-step early abort.
 *
 * \See C. Pomerance, <i>Analysis and Comparison of Some Integer Factoring
 * Algorithm</i>, in Mathematical Centre Tracts 154.
 */

 /*
  *  History:
  *    1.0: Thu Nov 15 2007 by JM:
  *         - Initial version.
  */

#if !defined(_TIFA_RES_TDIV_H_)
   /**
    * \def _TIFA_RES_TDIV_H_
    * Standard include guard.
    */
#define _TIFA_RES_TDIV_H_

#include <inttypes.h>
#include "smooth_filter.h"

#ifdef __cplusplus
extern "C" {
#endif

   /**
    * \brief Trial divide residues using data from a \c smooth_filter_t.
    *
    * Filters the relations given by <tt>filter->candidate_*</tt> via
    * trial division at the <tt>step</tt>-th early abort step and stores the
    * 'good' relations in <tt>filter->accepted_*</tt>. The \c step parameter
    * has no effect if <tt>filter->method == TDIV</tt> (i.e. no early abort
    * variation).
    *
    * \warning This function is only meant to be used if <tt>filter->method == 
    * TDIV</tt> is either \c TDIV or \c TDIV_EARLY_ABORT.
    *
    * \param[in] filter a pointer to the \c smooth_filter_t to use.
    * \param[in] step the step in the early abort strategy to perform.
    */
inline uint32_t res_tdiv(smooth_filter_t* const filter, unsigned long int step);

#ifdef __cplusplus
}
#endif

#endif
