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
 * \date    Mon Mar 10 2008
 * \version 1.2.2
 */

 /*
  * History:
  * --------
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "tifa_config.h"

#if USE_CALLOC_MEMSET
    #include <string.h>
#endif

#include "array.h"
#include "funcs.h"
#include "macros.h"
#include "bitstring_t.h"

/*
 *-----------------------------------------------------------------------------
 *              byte_array_t and associated functions
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
byte_array_t* alloc_byte_array(uint32_t length) {

    byte_array_t* array = malloc(sizeof(byte_array_t));

    array->alloced = length;
    array->length  = 0U;

#if USE_CALLOC_MEMSET
    array->data = calloc(array->alloced, sizeof(unsigned char));
#else
    array->data = malloc(array->alloced * sizeof(unsigned char));
    for (uint32_t i = 0; i < array->alloced; i++) {
        array->data[i] = 0;
    }
#endif

    return array;
}
//-----------------------------------------------------------------------------
void clear_byte_array(byte_array_t* array) {
    if (array->alloced != 0U) {
        free(array->data);
        free(array);
    }
}
//-----------------------------------------------------------------------------
void resize_byte_array(byte_array_t* array, uint32_t alloced) {

    if (alloced < array->length) {
        array->length = alloced;
    }
    if (array->alloced > 0) {
        array->data = realloc(array->data, alloced * sizeof(unsigned char));
    } else {
        array->data = malloc(alloced * sizeof(unsigned char));
    }
    array->alloced = alloced;
}
//-----------------------------------------------------------------------------
void append_byte_to_array(byte_array_t* array, const unsigned char to_append) {
    if (array->length >= array->alloced) {
        resize_byte_array(array, array->length + ELONGATION);
    }
    array->data[array->length] = to_append;
    array->length++;
}
//-----------------------------------------------------------------------------
void append_byte_array(byte_array_t* const array,
                       const byte_array_t* const to_append) {

    if ((array->length + to_append->length) > array->alloced) {
        resize_byte_array(
            array,
            array->length + to_append->length + ELONGATION
        );
    }
    uint32_t current = array->length;

    for (uint32_t i = 0; i < to_append->length; i++ ) {
        array->data[current] = to_append->data[i];
        current++;
    }
    array->length += to_append->length;
}
//-----------------------------------------------------------------------------
void swap_byte_array(byte_array_t* const a, byte_array_t* const b) {

    uint32_t tmpalloced = a->alloced;
    uint32_t tmplength  = a->length;
    unsigned char* tmpdata = a->data;

    a->alloced = b->alloced;
    a->length  = b->length;
    a->data    = b->data;

    b->alloced = tmpalloced;
    b->length  = tmplength;
    b->data    = tmpdata;
}
//-----------------------------------------------------------------------------
void print_byte_array(const byte_array_t* const array) {
    //
    // Mostly for debugging purposes as the output is not particularly
    // well structured...
    //
    for (uint32_t i = 0U; i < array->length; i++) {
        printf("%u ", (unsigned int)array->data[i]);
    }
    printf("\n");
}
//-----------------------------------------------------------------------------
void ins_sort_byte_array(byte_array_t* const array) {
    //
    // A very basic insertion sort...
    //
    unsigned char to_insert;
    uint32_t j;
    for (uint32_t i = 1; i < array->length; i++) {
        to_insert = array->data[i];
        j = i;
        while ((j > 0) && (array->data[j-1] > to_insert)) {
            array->data[j] = array->data[j-1];
            j--;
        }
        array->data[j] = to_insert;
    }
}
//-----------------------------------------------------------------------------
void qsort_byte_array(byte_array_t* const array) {
    //
    // Quick sort using the qsort function from the standard C library...
    //
    qsort(array->data, array->length, sizeof(unsigned char), uint32_cmp_func);
}
//-----------------------------------------------------------------------------
uint32_t index_in_byte_array(unsigned char to_find,
                             const byte_array_t* const array) {

    for (uint32_t i = 0; i < array->length; i++) {
        if (to_find == array->data[i]) {
            return i;
        }
    }
    //
    // _KLUDGE_: If the integer to_find is not found, we'd like to return -1.
    //           However, the return type being unsigned, returning -1 would
    //           not be very elegant, so we return UINT32_MAX (i.e. the value
    //           of the NOT_IN_ARRAY symbol) which, in a signed context, is 
    //           indeed -1. Of course, the other possibility would be to return 
    //           a int32_t but that would limit the length of the uint32_array_t 
    //           to INT32_MAX instead of UINT32_MAX.
    //
    return NOT_IN_ARRAY;
}
//-----------------------------------------------------------------------------
uint32_t index_in_sorted_byte_array(unsigned char to_find,
                                    const byte_array_t* const sorted_array,
                                    uint32_t min_index, uint32_t max_index) {
    //
    // Basic binary search...
    //
    uint32_t min = min_index - 1;
    //
    // _KLUDGE_: What if min_index = 0 ? Substracting one will yield
    //           -1, or rather, UINT32_MAX since min is unsigned. This is
    //           not very elegant but it works...
    //
    uint32_t max = max_index;
    uint32_t index;

    while ((max - min) > 1) {

        index = (max + min)/2;

        if (to_find <= sorted_array->data[index]) {
            max = index;
        } else {
            min = index;
        }
    }
    if (to_find == sorted_array->data[max]) {
        return max;
    } else {
        //
        // _KLUDGE_: Return type is unsigned so -1 is the same as UINT32_MAX 
        //           (i.e. the value of the NOT_IN_ARRAY symbol).
        //
        return NOT_IN_ARRAY;
    }
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *              uint32_array_t and associated functions
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
uint32_array_t* alloc_uint32_array(uint32_t length) {

    uint32_array_t* array = malloc(sizeof(uint32_array_t));

    array->alloced = length;
    array->length  = 0U;

#if USE_CALLOC_MEMSET
    array->data = calloc(array->alloced, sizeof(uint32_t));
#else
    array->data = malloc(array->alloced * sizeof(uint32_t));
    for (uint32_t i = 0U; i < array->alloced; i++) {
        array->data[i] = 0U;
    }
#endif

    return array;
}
//-----------------------------------------------------------------------------
void clear_uint32_array(uint32_array_t* array) {
    if (array->alloced != 0U) {
        free(array->data);
        free(array);
    }
}
//-----------------------------------------------------------------------------
void resize_uint32_array(uint32_array_t* array, uint32_t alloced) {

    if (alloced < array->length) {
        array->length = alloced;
    }
    if (array->alloced > 0) {
        array->data = realloc(array->data, alloced * sizeof(uint32_t));
    } else {
        array->data = malloc(alloced * sizeof(uint32_t));
    }
    array->alloced = alloced;
}
//-----------------------------------------------------------------------------
void append_uint32_to_array(uint32_array_t* array, const uint32_t to_append) {
    if (array->length >= array->alloced) {
        resize_uint32_array(array, array->length + ELONGATION);
    }
    array->data[array->length] = to_append;
    array->length++;
}
//-----------------------------------------------------------------------------
void append_uint32_array(uint32_array_t* const array,
                         const uint32_array_t* const to_append) {

    if ((array->length + to_append->length) > array->alloced) {
        resize_uint32_array(
            array,
            array->length + to_append->length + ELONGATION
        );
    }
    uint32_t current = array->length;

    for (uint32_t i = 0; i < to_append->length; i++ ) {
        array->data[current] = to_append->data[i];
        current++;
    }
    array->length += to_append->length;
}
//-----------------------------------------------------------------------------
void swap_uint32_array(uint32_array_t* const a, uint32_array_t* const b) {

    uint32_t  tmpalloced = a->alloced;
    uint32_t  tmplength  = a->length;
    uint32_t* tmpdata    = a->data;

    a->alloced = b->alloced;
    a->length  = b->length;
    a->data    = b->data;

    b->alloced = tmpalloced;
    b->length  = tmplength;
    b->data    = tmpdata;
}
//-----------------------------------------------------------------------------
void print_uint32_array(const uint32_array_t* const array) {
    //
    // Mostly for debugging purposes as the output is not particularly
    // well structured...
    //
    for (uint32_t i = 0U; i < array->length; i++) {
        printf("%"PRIu32" ", array->data[i]);
    }
    printf("\n");
}
//-----------------------------------------------------------------------------
void ins_sort_uint32_array(uint32_array_t* const array) {
    //
    // A very basic insertion sort...
    //
    uint32_t to_insert;
    uint32_t j;
    for (uint32_t i = 1; i < array->length; i++) {
        to_insert = array->data[i];
        j = i;
        while ((j > 0) && (array->data[j-1] > to_insert)) {
            array->data[j] = array->data[j-1];
            j--;
        }
        array->data[j] = to_insert;
    }
}
//-----------------------------------------------------------------------------
void qsort_uint32_array(uint32_array_t* const array) {
    //
    // Quick sort using the qsort function from the standard C library...
    //
    qsort(array->data, array->length, sizeof(uint32_t), uint32_cmp_func);
}
//-----------------------------------------------------------------------------
uint32_t index_in_uint32_array(uint32_t to_find,
                               const uint32_array_t* const array) {

    for (uint32_t i = 0; i < array->length; i++) {
        if (to_find == array->data[i]) {
            return i;
        }
    }
    //
    // _KLUDGE_: If the integer to_find is not found, we'd like to return -1.
    //           However, the return type being unsigned, returning -1 would
    //           not be very elegant, so we return UINT32_MAX (i.e. the value
    //           of the NOT_IN_ARRAY symbol) which, in a signed context, is 
    //           indeed -1. Of course, the other possibility would be to return 
    //           a int32_t but that would limit the length of the uint32_array_t 
    //           to INT32_MAX instead of UINT32_MAX.
    //
    return NOT_IN_ARRAY;
}
//-----------------------------------------------------------------------------
uint32_t index_in_sorted_uint32_array(uint32_t to_find,
                                      const uint32_array_t* const sorted_array,
                                      uint32_t min_index, uint32_t max_index) {
    //
    // Basic binary search...
    //
    uint32_t min = min_index - 1;
    //
    // _KLUDGE_: What if min_index = 0 ? Substracting one will yield
    //           -1, or rather, UINT32_MAX since min is unsigned. This is
    //           not very elegant but it works...
    //
    uint32_t max = max_index;
    uint32_t index;

    while ((max - min) > 1) {

        index = (max + min)/2;

        if (to_find <= sorted_array->data[index]) {
            max = index;
        } else {
            min = index;
        }
    }
    if (to_find == sorted_array->data[max]) {
        return max;
    } else {
        //
        // _KLUDGE_: Return type is unsigned so -1 is the same as UINT32_MAX 
        //           (i.e. the value of the NOT_IN_ARRAY symbol).
        //
        return NOT_IN_ARRAY;
    }
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *              int32_array_t and associated functions
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
int32_array_t* alloc_int32_array(uint32_t length) {

    int32_array_t* array = malloc(sizeof(int32_array_t));

    array->alloced = length;
    array->length  = 0U;

#if USE_CALLOC_MEMSET
    array->data = calloc(array->alloced, sizeof(int32_t));
#else
    array->data = malloc(array->alloced * sizeof(int32_t));
    for (uint32_t i = 0U; i < array->alloced; i++) {
        array->data[i] = 0;
    }
#endif

    return array;
}
//-----------------------------------------------------------------------------
void clear_int32_array(int32_array_t* array) {
    if (array->alloced != 0U) {
        free(array->data);
        free(array);
    }
}
//-----------------------------------------------------------------------------
void resize_int32_array(int32_array_t* array, uint32_t alloced) {

    if (alloced < array->length) {
        array->length = alloced;
    }
    if (array->alloced > 0) {
        array->data = realloc(array->data, alloced * sizeof(int32_t));
    } else {
        array->data = malloc(alloced * sizeof(int32_t));
    }
    array->alloced = alloced;
}
//-----------------------------------------------------------------------------
void append_int32_to_array(int32_array_t* array, const int32_t to_append) {
    if (array->length >= array->alloced) {
        resize_int32_array(array, array->length + ELONGATION);
    }
    array->data[array->length] = to_append;
    array->length++;
}
//-----------------------------------------------------------------------------
void append_int32_array(int32_array_t* const array,
                        const int32_array_t* const to_append) {

    if ((array->length + to_append->length) > array->alloced) {
        resize_int32_array(
            array,
            array->length + to_append->length + ELONGATION
        );
    }
    uint32_t current = array->length;

    for (uint32_t i = 0; i < to_append->length; i++ ) {
        array->data[current] = to_append->data[i];
        current++;
    }
    array->length += to_append->length;
}
//-----------------------------------------------------------------------------
void print_int32_array(const int32_array_t* const array) {
    //
    // Mostly for debugging purposes as the output is not particularly
    // well structured...
    //
    for (uint32_t i = 0U; i < array->length; i++) {
        printf("%"PRId32" ", array->data[i]);
    }
    printf("\n");
}
//-----------------------------------------------------------------------------
void swap_int32_array(int32_array_t* const a, int32_array_t* const b) {

    uint32_t tmpalloced = a->alloced;
    uint32_t tmplength  = a->length;
    int32_t* tmpdata    = a->data;

    a->alloced = b->alloced;
    a->length  = b->length;
    a->data    = b->data;

    b->alloced = tmpalloced;
    b->length  = tmplength;
    b->data    = tmpdata;
}
//-----------------------------------------------------------------------------
uint32_t index_in_int32_array(int32_t to_find,
                              const int32_array_t* const array) {

    for (uint32_t i = 0; i < array->length; i++) {
        if (to_find == array->data[i]) {
            return i;
        }
    }
    //
    // _KLUDGE_: If the integer to_find is not found, we'd like to return -1.
    //           However, the return type being unsigned, returning -1 would
    //           not be very elegant, so we return UINT32_MAX (i.e. the value
    //           of the NOT_IN_ARRAY symbol) which, in a signed context, is 
    //           indeed -1. Of course, the other possibility would be to return 
    //           a int32_t but that would limit the length of the uint32_array_t 
    //           to INT32_MAX instead of UINT32_MAX.
    //
    return NOT_IN_ARRAY;
}
//-----------------------------------------------------------------------------
uint32_t index_in_sorted_int32_array(int32_t to_find,
                                     const int32_array_t* const sorted_array,
                                     uint32_t min_index, uint32_t max_index) {
    //
    // Basic binary search...
    //
    uint32_t min = min_index - 1;
    //
    // _KLUDGE_: What if min_index = 0 ? Substracting one will yield
    //           -1, or rather, UINT32_MAX since min is unsigned. This is
    //           not very elegant but it works...
    //
    uint32_t max = max_index;
    uint32_t index;

    while ((max - min) > 1) {

        index = (max + min)/2;

        if (to_find <= sorted_array->data[index]) {
            max = index;
        } else {
            min = index;
        }
    }
    if (to_find == sorted_array->data[max]) {
        return max;
    } else {
        //
        // _KLUDGE_: Return type is unsigned so -1 is the same as UINT32_MAX
        //           (i.e. the value of the NOT_IN_ARRAY symbol).
        //
        return NOT_IN_ARRAY;
    }
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                  mpz_array_t and associated functions
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
mpz_array_t* alloc_mpz_array(uint32_t length) {

    mpz_array_t* array = malloc(sizeof(mpz_array_t));

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
void resize_mpz_array(mpz_array_t* const array, uint32_t alloced) {
    if (alloced < array->alloced) {
        //
        // Don't forget to free the memory used by the mpz_t's that will be
        // discarded in case the new size is less than the original!
        //
        // _WARNING_: Since version 1.2, the mpz_t array given by array->data
        //            is fully mpz_init'ed once and for all. Consequently
        //            array->data has to be mpz_clear'ed from alloced to
        //            array->alloced (and not just up to array->length as in 
        //            older versions).
        //
        for (uint32_t i = alloced; i < array->alloced; i++) {
            mpz_clear(array->data[i]);
        }
        array->length = alloced;
    }
    if (array->alloced > 0) {
        array->data = realloc(array->data, alloced * sizeof(mpz_t));
    } else {
        array->data = malloc(alloced * sizeof(mpz_t));
    }
    //
    // _WARNING_: Since version 1.2, the mpz_t array given by array->data is
    //            fully mpz_init'ed once and for all!
    //
    for (uint32_t i = array->alloced; i < alloced; i++) {
        mpz_init(array->data[i]);
    }
    array->alloced = alloced;
}
//-----------------------------------------------------------------------------
void append_mpz_to_array(mpz_array_t* array, const mpz_t to_append) {
    if (array->length >= array->alloced) {
        resize_mpz_array(array, array->length + ELONGATION);
    }
    //
    // _WARNING_: Since version 1.2, the mpz_t array given by array->data is
    //            fully mpz_init'ed once and for all!
    //
    mpz_set(array->data[array->length], to_append);
    array->length++;
}
//-----------------------------------------------------------------------------
void append_mpz_array(mpz_array_t* const array,
                      const mpz_array_t* const to_append) {

    if ((array->length + to_append->length) > array->alloced) {
        resize_mpz_array(array, array->length + to_append->length + ELONGATION);
    }
    uint32_t current = array->length;

    for (uint32_t i = 0; i < to_append->length; i++ ) {
        //
        // _WARNING_: Since version 1.2, the mpz_t array given by array->data
        //            is fully mpz_init'ed once and for all!
        //
        mpz_set(array->data[current], to_append->data[i]);
        current++;
    }
    array->length += to_append->length;
}
//-----------------------------------------------------------------------------
void swap_mpz_array(mpz_array_t* const a, mpz_array_t* const b) {

    uint32_t tmpalloced = a->alloced;
    uint32_t tmplength  = a->length;
    mpz_t*   tmpdata    = a->data;

    a->alloced = b->alloced;
    a->length  = b->length;
    a->data    = b->data;

    b->alloced = tmpalloced;
    b->length  = tmplength;
    b->data    = tmpdata;
}
//-----------------------------------------------------------------------------
void print_mpz_array(const mpz_array_t* const array) {
    //
    // Mostly for debugging purposes as the output is not particularly
    // well structured...
    //
    for (uint32_t i = 0U; i != array->length; i++) {
        gmp_printf("%Zd\n", array->data[i]);
    }
}
//-----------------------------------------------------------------------------
uint32_t index_in_mpz_array(const mpz_t to_find,
                            const mpz_array_t* const array) {

    for (uint32_t i = 0; i < array->length; i++) {
        if (0 == mpz_cmp(to_find, array->data[i])) {
            return i;
        }
    }
    //
    // _KLUDGE_: Return type is unsigned so -1 is the same as UINT32_MAX
    //           (i.e. the value of the NOT_IN_ARRAY symbol).
    //
    return NOT_IN_ARRAY;
}
//-----------------------------------------------------------------------------
uint32_t index_in_sorted_mpz_array(const mpz_t to_find,
                                   const mpz_array_t* const sorted_array,
                                   uint32_t min_index, uint32_t max_index) {
    //
    // Basic binary search...
    //
    uint32_t min = min_index - 1;
    //
    // _KLUDGE_: What if min_index = 0 ? Substracting one will yield
    //           -1, or rather, UINT32_MAX since min is unsigned. This is
    //           not very elegant but it works...
    //
    uint32_t max = max_index;
    uint32_t index;

    while ((max - min) > 1) {

        index = (max + min)/2;

        if (mpz_cmp(to_find, sorted_array->data[index]) <= 0) {
            max = index;
        } else {
            min = index;
        }
    }
    if (mpz_cmp(to_find, sorted_array->data[max]) == 0) {
        return max;
    } else {
        //
        // _KLUDGE_: Return type is unsigned so -1 is the same as UINT32_MAX
        //           (i.e. the value of the NOT_IN_ARRAY symbol).
        //
        return NOT_IN_ARRAY;
    }
}
//-----------------------------------------------------------------------------
void ins_sort_mpz_array(mpz_array_t* const array) {
    //
    // A very basic insertion sort...
    //
    mpz_t to_insert;
    mpz_init(to_insert);
    uint32_t j;
    for (uint32_t i = 1; i < array->length; i++) {
        mpz_set(to_insert, array->data[i]);
        j = i;
        while ((j > 0) && (mpz_cmp(array->data[j - 1], to_insert) > 0)) {
            mpz_set(array->data[j], array->data[j - 1]);
            j--;
        }
        mpz_set(array->data[j], to_insert);
    }
    mpz_clear(to_insert);
}
//-----------------------------------------------------------------------------
void qsort_mpz_array(mpz_array_t* const array) {
    //
    // Quick sort using the qsort function from the standard C library...
    //
    qsort(array->data, array->length, sizeof(mpz_t), mpz_cmp_func);
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *                  binary_array_t and associated functions
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
binary_array_t* alloc_binary_array(uint32_t length) {

    binary_array_t* array = malloc(sizeof(binary_array_t));
    //
    // Remember now that array->alloced is not the number of bits but the
    // smallest number of TIFA_BITSTRING_T needed to store length bits.
    //
    array->alloced = ceil(((float)length)/BITSTRING_T_BITSIZE);
    array->length  = 0U;

#if USE_CALLOC_MEMSET
    array->data = calloc(array->alloced, sizeof(TIFA_BITSTRING_T));
#else
    array->data = malloc(array->alloced * sizeof(TIFA_BITSTRING_T));
    for (uint32_t i = 0U; i < array->alloced; i++) {
        array->data[i]  = 0;
        //
        // Use an xor instead of simply 0 to make sure we have a bitstring
        // made of 0 bits only. Granted, TIFA is not likely to compile on
        // the exotic architectures where the value 0 is not made of 0 bits
        // only anyways, but still...
        //
        array->data[i] ^= 0;
    }
#endif

    return array;
}
//-----------------------------------------------------------------------------
void clear_binary_array(binary_array_t* array) {
    if (array->alloced != 0U) {
        free(array->data);
        free(array);
    }
}
//-----------------------------------------------------------------------------
void resize_binary_array(binary_array_t* array, uint32_t alloced) {

    uint32_t nalloc = ceil(((float)alloced)/BITSTRING_T_BITSIZE);

    if (alloced < array->length) {
        array->length = alloced;
    }
    if (array->alloced > 0) {
        array->data = realloc(array->data, nalloc * sizeof(TIFA_BITSTRING_T));
    } else {
        array->data = malloc(nalloc * sizeof(TIFA_BITSTRING_T));
    }
    if (array->alloced < nalloc) {
        //
        // For binary array, proper initialization with 0 bits is important.
        //
        for (uint32_t i = array->alloced; i < nalloc; i++ ) {
            array->data[i] = 0;
            //
            // Use an xor instead of simply 0 to make sure we have a bitstring
            // made of 0 bits only. Granted, TIFA is not likely to compile on
            // the exotic architectures where the value 0 is not made of 0 bits
            // only anyways, but still...
            //
            array->data[i] ^= 0;
        }
    }
    array->alloced = nalloc;
}
//-----------------------------------------------------------------------------
void append_bit_to_array(binary_array_t* array, const unsigned int to_append) {

    if (array->length >= (array->alloced * BITSTRING_T_BITSIZE)) {
        resize_binary_array(
            array,
            array->length + ELONGATION * BITSTRING_T_BITSIZE
        );
    }
    if (to_append != 0) {
        set_array_bit_to_one(array->length, array);
    } else {
        set_array_bit_to_zero(array->length, array);
    }
    array->length++;
}
//-----------------------------------------------------------------------------
void print_binary_array(const binary_array_t* const array) {
    if (array->length == 0U) {
        printf("<empty binary_array>");
    }
    for (uint32_t i = 0U; i != array->length; i++) {
        printf("%d", get_array_bit(i, array));
    }
    printf("\n");
}
//-----------------------------------------------------------------------------
uint8_t get_array_bit(uint32_t index, const binary_array_t* const array) {

#if BITSTRING_T_SIZE_IS_POW_OF_TWO
    uint32_t offset = index & (BITSTRING_T_BITSIZE - 1);
    index >>= POW_TWO_BITSTRING_T_SIZE;
#else
    uint32_t offset = index % BITSTRING_T_BITSIZE;
    index /= BITSTRING_T_BITSIZE;
#endif

    offset = BITSTRING_T_BITSIZE - 1 - offset;

    if (0 == ((((TIFA_BITSTRING_T)1)<<offset) & array->data[index])) {
        return 0;
    } else {
        return 1;
    }
}
//-----------------------------------------------------------------------------
void set_array_bit_to_one(uint32_t index, binary_array_t* const array) {

#if BITSTRING_T_SIZE_IS_POW_OF_TWO
    uint32_t offset = index & (BITSTRING_T_BITSIZE - 1);
    index >>= POW_TWO_BITSTRING_T_SIZE;
#else
    uint32_t offset = index % BITSTRING_T_BITSIZE;
    index /= BITSTRING_T_BITSIZE;
#endif

    offset = BITSTRING_T_BITSIZE - 1 - offset;

    array->data[index] |= (((TIFA_BITSTRING_T)1)<<offset);
}
//-----------------------------------------------------------------------------
void set_array_bit_to_zero(uint32_t index, binary_array_t* const array) {

#if BITSTRING_T_SIZE_IS_POW_OF_TWO
    uint32_t offset = index & (BITSTRING_T_BITSIZE - 1);
    index >>= POW_TWO_BITSTRING_T_SIZE;
#else
    uint32_t offset = index % BITSTRING_T_BITSIZE;
    index /= BITSTRING_T_BITSIZE;
#endif

    offset = BITSTRING_T_BITSIZE - 1 - offset;

    array->data[index] &= !(((TIFA_BITSTRING_T)1)<<offset);
}
//-----------------------------------------------------------------------------
void flip_array_bit(uint32_t index, binary_array_t* const array) {

#if BITSTRING_T_SIZE_IS_POW_OF_TWO
    uint32_t offset = index & (BITSTRING_T_BITSIZE - 1);
    index >>= POW_TWO_BITSTRING_T_SIZE;
#else
    uint32_t offset = index % BITSTRING_T_BITSIZE;
    index /= BITSTRING_T_BITSIZE;
#endif

    offset = BITSTRING_T_BITSIZE - 1 - offset;

    array->data[index] ^= (((TIFA_BITSTRING_T)1)<<offset);
}
//-----------------------------------------------------------------------------
