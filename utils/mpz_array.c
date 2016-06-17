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
 * \file    array.c
 * \author  Jerome Milan
 * extracted by Paul Zimmermann from the TIFA library, December 8, 2008
 * \date    Circa Fri Mar 29 2008
 * \version 1.2.3
 */

/*
 * History:
 * --------
 * 1.2.3: Fri Mar 29 (?) 2008 by JM:
 *        - Inlined functions pertaining to binary_array_t's in array.h and
 *          removed them from this file.
 * 1.2.2: Mon Mar 10 2008 by JM:
 *        - Inlined is_in_sorted_*_array(...) functions.
 * 1.2.1: Fri Mar 7 2008 by JM:
 *        - WARNING: Changed memory managment. clear_*_array functions now
 *                   frees the pointer to the array (that's more logical).
 *   1.2: Fri Feb 29 2008 by JM:
 *        - WARNING: Semantic of mpz_array_t's length field has changed!.
 *                   Now length is only the "useful part" of the array, but
 *                   the data field is completely mpz_init'ed (i.e. up to
 *                   the alloced field, not longer up to length).
 *   1.1: Wed Oct 17 2007 by JM:
 *        - Added byte_array_t with corresponding functions.
 *   1.0: Wed Mar 1 2006 by JM:
 *        - Initial version.
 */

#include "cado.h"
#include <stdlib.h>

#include "mpz_array.h"
#include "macros.h"

/* Note: this routine is only used in append_mpz_to_array, with
   alloced = array->alloced + ELONGATION, thus alloced >= array->alloced.
   Since the initial allocation has always alloced > 0, and alloced cannot
   decrease, we have array->alloced > 0. */
//-----------------------------------------------------------------------------
static void resize_mpz_array (mpz_array_t* const array, uint32_t alloced) {
    ASSERT(array->alloced > 0);
    array->data = realloc (array->data, alloced * sizeof(mpz_t));
    //
    // _WARNING_: Since version 1.2, the mpz_t array given by array->data is
    //            fully mpz_init'ed once and for all!
    //
    for (uint32_t i = array->alloced; i < alloced; i++)
      mpz_init (array->data[i]);
    array->alloced = alloced;
}

//-----------------------------------------------------------------------------
void append_mpz_to_array (mpz_array_t* array, const mpz_t to_append) {
    if (array->length >= array->alloced)
      /* necessarily array->length = array->alloced in this case */
      resize_mpz_array (array, array->length + ELONGATION);
    //
    // _WARNING_: Since version 1.2, the mpz_t array given by array->data is
    //            fully mpz_init'ed once and for all!
    //
    mpz_set (array->data[array->length], to_append);
    array->length++;
}

/* this routine is only called from append_uint32_to_array, with
   alloced = array->length + ELONGATION = array->alloced + ELONGATION.
   Since the initial allocation has always alloced > 0, and alloced cannot
   decrease, we have array->alloced > 0. */
static void resize_uint32_array (uint32_array_t* array, uint32_t alloced) {

    ASSERT(array->alloced > 0);
    array->data = realloc (array->data, alloced * sizeof(uint32_t));
    array->alloced = alloced;
}

//-----------------------------------------------------------------------------
void append_uint32_to_array (uint32_array_t* array, const uint32_t to_append) {
    if (array->length >= array->alloced)
      /* necessarily array->length = array->alloced here */
      resize_uint32_array (array, array->length + ELONGATION);
    array->data[array->length] = to_append;
    array->length++;
}

//-----------------------------------------------------------------------------
mpz_array_t* alloc_mpz_array (uint32_t length) {

    mpz_array_t* array = malloc(sizeof(mpz_array_t));

    ASSERT(length > 0);

    array->alloced = length;
    array->length  = 0U;
    array->data    = malloc(array->alloced * sizeof(mpz_t));

    //
    // _WARNING_: Since version 1.2, the mpz_t array given by array->data is
    //            fully mpz_init'ed once and for all to avoid possible
    //            memory deallocations and reallocations.
    //
    for (uint32_t i = 0U; i < array->alloced; i++) {
        mpz_init(array->data[i]);
    }

    return array;
}

//-----------------------------------------------------------------------------
uint32_array_t* alloc_uint32_array (uint32_t length) {

    uint32_array_t* array = malloc(sizeof(uint32_array_t));

    ASSERT(length > 0);

    array->alloced = length;
    array->length  = 0U;

#if USE_CALLOC_MEMSET
    array->data = calloc (array->alloced, sizeof(uint32_t));
#else
    array->data = malloc (array->alloced * sizeof(uint32_t));
    for (uint32_t i = 0U; i < array->alloced; i++)
      array->data[i] = 0U;
#endif

    return array;
}

//-----------------------------------------------------------------------------
void clear_mpz_array(mpz_array_t* array) {
    if (array->alloced != 0U) {
        //
        // _WARNING_: Since version 1.2, the mpz_t array given by array->data
        //            is fully mpz_init'ed once and for all. Consequently
        //            array->data has to be completely mpz_clear'ed (and not
        //            just up to array->length as in older versions).
        //
        for (uint32_t i = 0U; i < array->alloced; i++) {
            mpz_clear(array->data[i]);
        }
        free(array->data);
        free(array);
    }
}

//-----------------------------------------------------------------------------
void clear_uint32_array(uint32_array_t* array) {
    if (array->alloced != 0U) {
        free(array->data);
        free(array);
    }
}

