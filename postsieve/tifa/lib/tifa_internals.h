//
// Copyright (C) 2006, 2007, 2008 INRIA (French National Institute for Research
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
 * \file    tifa_internals.h
 * \author  Jerome Milan
 * \date    Thu Jan 10 2008
 * \version 1.0
 *
 * \brief Library wide include file (complete with internal structures /
 * functions).
 *
 * Includes all TIFA's structures and functions.
 *
 * \warning Usually, only the tifa.h include file is needed. tifa_internals.h
 * should only be included to access some internal structures or functions.
 * Be warned that conflicts with client code or external libraries are then
 * more likely to occur.
 */

#if !defined(_TIFA_TIFA_H_)
#define _TIFA_TIFA_H_

//
// The following symbols need to be defined before including messages.h
// and timer.h.
//
#if !defined(__VERBOSE__)
    #define __VERBOSE__ 0
#endif
#if !defined(__TIMING__)
    #define __TIMING__ 0
#endif
#if !defined(__PREFIX__)
    #define __PREFIX__ ""
#endif

//
// Configuration file
//
#include "tifa_config.h"
//
// Includes from lib/algo
//
#include "cfrac.h"
#include "ecm.h"
#include "fermat.h"
#include "qs.h"
#include "siqs.h"
#include "squfof.h"
#include "tdiv.h"
#include "tifa_factor.h"
//
// Includes from lib/data
//
#include "first_primes.h"
//
// Includes from lib/utils
//
#include "array.h"
#include "bernsteinisms.h"
#include "bitstring_t.h"
#include "exit_codes.h"
#include "factoring_machine.h"
#include "funcs.h"
#include "gauss_elim.h"
#include "gmp_utils.h"
#include "hashtable.h"
#include "lindep.h"
#include "linked_list.h"
      //
      // _NOTE_: we also include macros.h but conflicts are highly likely!
      //
#include "macros.h"
#include "matrix.h"
#include "messages.h"
#include "print_error.h"
#include "res_tdiv.h"
#include "smooth_filter.h"
#include "sqrt_cont_frac.h"
#include "timer.h"
#include "x_array_list.h"
#include "x_tree.h"

#endif
