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
 * \file    tifa.h
 * \author  Jerome Milan
 * \date    Thu Jan 25 2007
 * \version 1.0
 *
 * \brief Library wide include file.
 *
 * Includes all TIFA's structures and functions.
 */

 /*
  *  Copyright (C) 2006, 2007  INRIA
  *  Licence: GNU Lesser General Public License (LGPL)
  *  History:
  *
  *  Mon Sep 4 2006: Initial version by JM
  */

#if !defined(_TIFA_TIFA_H_)
#define _TIFA_TIFA_H_

//
// _TO_DO_: Not everything needs to be visible from outside of the TIFA
//          library, so we should certainly sort this out and #include
//          in this header only the absolutely necessary structures and
//          functions.
//
#include "tifa_config.h"

#include "cfrac.h"
#include "qs.h"
#include "siqs.h"
#include "squfof.h"
#include "tdiv.h"

#include "first_primes.h"

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
#include "matrix.h"
#include "messages.h"
#include "print_error.h"
#include "sqrt_cont_frac.h"
#include "timer.h"
#include "x_array_list.h"
#include "x_tree.h"

#endif
