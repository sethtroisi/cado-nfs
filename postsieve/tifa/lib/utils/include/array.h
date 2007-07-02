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
 * \file    array.h
 * \author  Jerome Milan
 * \date    Wed Mar 28 2007
 * \version 1.1
 *
 * \brief Higher level arrays and associated functions.
 *
 * Defines higher level arrays and their associated functions.
 */

 /*
  *  License: GNU Lesser General Public License (LGPL)
  *  History:
  *    1.1: Wed Mar 28 2007 by JM:
  *           - Added factor_array_t and associated functions.
  *    1.0: Wed Mar 1 2006 by JM:
  *           - Initial version.
  */

#if !defined(_TIFA_ARRAY_H_)
   /**
    * \def _TIFA_ARRAY_H_
    * Standard include guard.
    */
#define _TIFA_ARRAY_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "tifa_config.h"

#include <inttypes.h>
#include <stdbool.h>

#include <gmp.h>

   /**
    * \def ELONGATION
    * Incremental size used when automatically expanding the capacity of a
    * <tt>*_array_t<\tt>.
    *
    * \note This is, of course, an hint to GMP's limbs and nails :-)
    */
#define ELONGATION 16

/*
 *-----------------------------------------------------------------------------
 *              uint32_array_t and associated functions
 *-----------------------------------------------------------------------------
 */

   /**
    * \struct struct_uint32_array_t array.h lib/utils/include/array.h
    * \brief  Defines an array of <tt>uint32</tt>.
    *
    * This structure defines a special kind of \c uint32 array which knows
    * its current length and its allocated memory space.
    */
struct struct_uint32_array_t {
       /**
        * Memory space allocated for this array's \c data field, given as
        * a multiple of <tt>sizeof(uint32_t)</tt>. This is the maximum
        * number of \c uint32_t that the array can accommodate.
        */
    uint32_t alloced;
       /**
        * Current number of \c uint32_t hold in the array pointed by the
        * structure's \c data field.
        */
    uint32_t length;
       /**
        * Array of \c uint32_t whose size is given by the \c alloced field.
        */
    uint32_t* data;
};

   /**
    * \typedef uint32_array_t
    * \brief Equivalent to <tt>struct struct_uint32_array_t</tt>.
    */
typedef struct struct_uint32_array_t uint32_array_t;

   /**
    * \brief Allocates and returns a new <tt>uint32_array_t</tt>.
    *
    * Allocates and returns a new <tt>uint32_array_t</tt> such that:
    * \li its \c alloced field is set to the parameter length.
    * \li its \c length field is set to zero.
    * \li its \c data array is completely filled with zeroes.
    *
    * \param[in] length The maximum length of the \c uint32_array_t to allocate.
    * \return A pointer to the newly allocated \c uint32_array_t structure.
    */
uint32_array_t* alloc_uint32_array(uint32_t length);

   /**
    * \brief Clears a <tt>uint32_array_t</tt>.
    *
    * Clears a <tt>uint32_array_t</tt>, or, more precisely, clears the memory
    * space used by the array pointed by the \c data field of a
    * <tt>uint32_array_t</tt>. Also set its \c alloced and \c length fields
    * to zero.
    *
    * \param[in] array A pointer to the <tt>uint32_array_t</tt> to clear.
    */
void clear_uint32_array(uint32_array_t* const array);

   /**
    * \brief Resizes the allocated memory of an <tt>uint32_array_t</tt>.
    *
    * Resizes the storage available to an <tt>uint32_array_t</tt> to make room
    * for \c alloced integers, while preserving its content. If \c alloced is
    * less than the length of the array, then obviously some of its content
    * will be lost.
    *
    * \param[in] alloced The new maximum length of the \c uint32_array_t to
    *                    resize.
    * \return array A pointer to the \c uint32_array_t to resize.
    */
void resize_uint32_array(uint32_array_t* const array, uint32_t alloced);

   /**
    * \brief Appends a \c uint32_t to an <tt>uint32_array_t</tt>.
    *
    * Appends the \c uint32_t integer \c to_append to <tt>array</tt>.
    * If \c array has not enough capacity to accommodate this extra element it
    * will be resized via a call to \c resize_uint32_array adding \c ELONGATION
    * \c uint32_t slots to avoid too frequent resizes.
    *
    * \param[in] array A pointer to the recipient <tt>uint32_array_t</tt>.
    * \param[in] to_append The integer to append.
    */
void append_uint32_to_array(uint32_array_t* array, uint32_t to_append);

   /**
    * \brief Appends the content of a <tt>uint32_array_t</tt> to another one.
    *
    * Appends the content of the \c to_append array to the \c uint32_array_t
    * named \c array. If \c array has not enough capacity to accommodate all
    * elements from \c to_append, it will be resized via a call to
    * \c resize_uint32_array with extra room for \c ELONGATION unused
    * \c uint32_t slots to avoid too frequent resizes.
    *
    * \param[in] array A pointer to the recipient <tt>uint32_array_t</tt>.
    * \param[in] to_append A pointer to the <tt>uint32_array_t</tt> to append.
    */
void append_uint32_array(uint32_array_t* const array,
                         const uint32_array_t* const to_append);

   /**
    * \brief Swaps two <tt>uint32_array_t</tt>'s contents.
    *
    * Swaps the contents of \c a and <tt>b</tt>, two <tt>uint32_array_t</tt>'s.
    *
    * \note In some case, pointer swapping is inappropriate (for example, if
    * the pointers are passed as function arguments!), hence the need for such
    * a swapping function.
    *
    * \param[in] a A pointer to the first <tt>uint32_array_t</tt> to swap.
    * \param[in] b A pointer to the second <tt>uint32_array_t</tt> to swap.
    */
void swap_uint32_array(uint32_array_t* const a, uint32_array_t* const b);

   /**
    * \brief Prints a <tt>uint32_array_t</tt>.
    *
    * Prints a <tt>uint32_array_t</tt>'s \c data elements on the standard
    * output.
    *
    * \note This function is mostly intended for debugging purposes as the
    * output is not particularly well structured.
    *
    * \param[in] array A pointer to the <tt>uint32_array_t</tt> to print.
    */
void print_uint32_array(const uint32_array_t* const array);

   /**
    * \brief Sorts the \c uint32_t elements of a <tt>uint32_array_t</tt>.
    *
    * Sorts the uint32_t elements of a <tt>uint32_array_t</tt> in natural order
    * using a basic insertion sort.
    *
    * \param[in] array A pointer to the <tt>uint32_array_t</tt> to sort.
    */
void ins_sort_uint32_array(uint32_array_t* const array);

   /**
    * \brief Sorts the uint32_t elements of a <tt>uint32_array_t</tt> with a
    * quick sort.
    *
    * Sorts the uint32_t elements of a <tt>uint32_array_t</tt> in natural order
    * using the quick sort algorithm.
    *
    * \note This function relies on the C library implementation of the
    * quick sort provided by the function <tt>qsort</tt>.
    *
    * \param[in] array A pointer to the <tt>uint32_array_t</tt> to sort.
    */
void qsort_uint32_array(uint32_array_t* const array);


   /**
    * \brief Returns the position of an integer in a <tt>uint32_array_t</tt>.
    *
    * Returns the position of the integer \c to_find in the
    * <tt>uint32_array_t</tt> pointed to by \c array. If the integer \c to_find
    * is not found in the <tt>uint32_array_t</tt>, returns <tt>UINT32_MAX</tt>.
    *
    * \note The <tt>UINT32_MAX</tt> value is actually -1 if interpreted as a
    * signed <tt>int32_t</tt>.
    * \note The fact that \c UINT32_MAX is returned when the integer is not
    * found in the array means that the practical, usable maximum length of
    * a <tt>uint32_array_t</tt> is <tt>UINT32_MAX - 1</tt> and not
    * <tt>UINT32_MAX</tt>.
    * \note If the array is already sorted, the more efficient function
    * <tt>index_in_sorted_uint32_array</tt> can be used as it uses a basic
    * binary search instead of a complete scanning of the array.
    *
    * \param[in] to_find The integer to find in the <tt>uint32_array_t</tt>.
    * \param[in] array   A pointer to the <tt>uint32_array_t</tt>.
    * \returns The index of \c to_find in the array if \c to_find is found.
    * \returns <tt>UINT32_MAX</tt> otherwise.
    */
uint32_t index_in_uint32_array(uint32_t to_find,
                               const uint32_array_t* const array);

   /**
    * \brief Returns true if a given integer is in a given array.
    *
    * Returns true if the integer \c to_find is in the <tt>uint32_array_t</tt>
    * pointed to by \c array. Returns false otherwise.
    *
    * \note If the array is already sorted, the more efficient function
    * <tt>is_in_sorted_uint32_array</tt> can be used as it uses a basic
    * binary search instead of a complete scanning of the array.
    *
    * \param[in] to_find The integer to find in the <tt>uint32_array_t</tt>.
    * \param[in] array   A pointer to the <tt>uint32_array_t</tt>.
    * \returns true if \c to_find is in the array \c array.
    * \returns false otherwise.
    */
bool is_in_uint32_array(uint32_t to_find, const uint32_array_t* const array);

   /**
    * \brief Returns the position of an integer in a sorted portion of a
    * <tt>uint32_array_t</tt>.
    *
    * Returns the position of the integer \c to_find in a \e sorted portion
    * of the <tt>uint32_array_t</tt> pointed to by \c array. If the integer
    * \c to_find is not found in the portion delimited by \c min_index and
    * \c max_index, returns <tt>UINT32_MAX</tt>.
    *
    * \note The <tt>UINT32_MAX</tt> value is actually -1 if interpreted as a
    * signed <tt>int32_t</tt>.
    * \note The fact that \c UINT32_MAX is returned when the integer is not
    * found in the array means that the practical, usable maximum length of
    * a <tt>uint32_array_t</tt> is <tt>UINT32_MAX - 1</tt> and not
    * <tt>UINT32_MAX</tt>.
    *
    * \param[in] to_find The integer to find in the <tt>uint32_array_t</tt>.
    * \param[in] sorted_array A pointer to the <tt>uint32_array_t</tt>.
    * \param[in] min_index The beginning of the sorted array portion to
    *                      search in.
    * \param[in] max_index The end of the sorted array portion to search in.
    * \returns The index of \c to_find in the array if \c to_find is found in
    *          the sorted array portion.
    * \returns <tt>UINT32_MAX</tt> otherwise.
    */
uint32_t index_in_sorted_uint32_array(uint32_t to_find,
                                      const uint32_array_t* const sorted_array,
                                      uint32_t min_index, uint32_t max_index);

   /**
    * \brief Returns true if a given integer is in a given array.
    *
    * Returns true if the integer \c to_find is in the (\e already \e sorted)
    * <tt>uint32_array_t</tt> pointed to by \c array. Returns false otherwise.
    *
    * \param[in] to_find The integer to find in the <tt>uint32_array_t</tt>.
    * \param[in] array   A pointer to the \e sorted <tt>uint32_array_t</tt>.
    * \returns true if \c to_find is in the array \c array.
    * \returns false otherwise.
    */
bool is_in_sorted_uint32_array(uint32_t to_find,
                               const uint32_array_t* const array);

/*
 *-----------------------------------------------------------------------------
 *              int32_array_t and associated functions
 *-----------------------------------------------------------------------------
 */

   /**
    * \struct struct_int32_array_t array.h lib/utils/include/array.h
    * \brief  Defines an array of <tt>int32</tt>.
    *
    * This structure defines a special kind of \c int32 array which knows
    * its current length and its allocated memory space.
    */
struct struct_int32_array_t {
       /**
        * Memory space allocated for this array's \c data field, given as
        * a multiple of <tt>sizeof(int32_t)</tt>. This is the maximum
        * number of \c int32_t that the array can accommodate.
        */
    uint32_t alloced;
       /**
        * Current number of \c int32_t hold in the array pointed by the
        * structure's \c data field.
        */
    uint32_t length;
       /**
        * Array of \c int32_t whose size is given by the \c alloced field.
        */
    int32_t* data;
};

   /**
    * \typedef int32_array_t
    * \brief Equivalent to <tt>struct struct_int32_array_t</tt>.
    */
typedef struct struct_int32_array_t int32_array_t;

   /**
    * \brief Allocates and returns a new <tt>int32_array_t</tt>.
    *
    * Allocates and returns a new <tt>int32_array_t</tt> such that:
    * \li its \c alloced field is set to the parameter length.
    * \li its \c length field is set to zero.
    * \li its \c data array is completely filled with zeroes.
    *
    * \param[in] length The maximum length of the \c int32_array_t to allocate.
    * \return A pointer to the newly allocated \c int32_array_t structure.
    */
int32_array_t* alloc_int32_array(uint32_t length);

   /**
    * \brief Clears a <tt>int32_array_t</tt>.
    *
    * Clears a <tt>int32_array_t</tt>, or, more precisely, clears the memory
    * space used by the array pointed by the \c data field of a
    * <tt>int32_array_t</tt>. Also set its \c alloced and \c length fields
    * to zero.
    *
    * \param[in] array A pointer to the <tt>int32_array_t</tt> to clear.
    */
void clear_int32_array(int32_array_t* const array);

   /**
    * \brief Resizes the allocated memory of an <tt>int32_array_t</tt>.
    *
    * Resizes the storage available to an <tt>int32_array_t</tt> to make room
    * for \c alloced integers, while preserving its content. If \c alloced is
    * less than the length of the array, then obviously some of its content
    * will be lost.
    *
    * \param[in] alloced The new maximum length of the \c int32_array_t to
    *                    resize.
    * \return array A pointer to the \c int32_array_t to resize.
    */
void resize_int32_array(int32_array_t* const array, uint32_t alloced);

   /**
    * \brief Appends a \c int32_t to an <tt>int32_array_t</tt>.
    *
    * Appends the \c int32_t integer \c to_append to <tt>array</tt>.
    * If \c array has not enough capacity to accommodate this extra element it
    * will be resized via a call to \c resize_int32_array adding \c ELONGATION
    * \c int32_t slots to avoid too frequent resizes.
    *
    * \param[in] array A pointer to the recipient <tt>int32_array_t</tt>.
    * \param[in] to_append The integer to append.
    */
void append_int32_to_array(int32_array_t* array, int32_t to_append);

   /**
    * \brief Appends the content of an <tt>int32_array_t</tt> to another one.
    *
    * Appends the content of the \c to_append array to the \c int32_array_t
    * named \c array. If \c array has not enough capacity to accommodate all
    * elements from \c to_append, it will be resized via a call to
    * \c resize_int32_array with extra room for \c ELONGATION unused \c int32_t
    * slots to avoid too frequent resizes.
    *
    * \param[in] array A pointer to the recipient <tt>int32_array_t</tt>.
    * \param[in] to_append A pointer to the <tt>int32_array_t</tt> to append.
    */
void append_int32_array(int32_array_t* const array,
                        const int32_array_t* const to_append);

   /**
    * \brief Swaps two <tt>int32_array_t</tt>'s contents.
    *
    * Swaps the contents of \c a and <tt>b</tt>, two <tt>int32_array_t</tt>'s.
    *
    * \note In some case, pointer swapping is inappropriate (for example, if
    * the pointers are passed as function arguments!), hence the need for such
    * a swapping function.
    *
    * \param[in] a A pointer to the first <tt>int32_array_t</tt> to swap.
    * \param[in] b A pointer to the second <tt>int32_array_t</tt> to swap.
    */
void swap_int32_array(int32_array_t* const a, int32_array_t* const b);

   /**
    * \brief Prints a <tt>int32_array_t</tt>.
    *
    * Prints a <tt>int32_array_t</tt>'s \c data elements on the standard
    * output.
    *
    * \note This function is mostly intended for debugging purposes as the
    * output is not particularly well structured.
    *
    * \param[in] array A pointer to the <tt>int32_array_t</tt> to print.
    */
void print_int32_array(const int32_array_t* const array);

   /**
    * \brief Returns the position of an integer in a <tt>int32_array_t</tt>.
    *
    * Returns the position of the integer \c to_find in the
    * <tt>int32_array_t</tt> pointed to by \c array. If the integer \c to_find
    * is not found in the <tt>int32_array_t</tt>, returns <tt>UINT32_MAX</tt>.
    *
    * \note The <tt>UINT32_MAX</tt> value is actually -1 if interpreted as a
    * signed <tt>int32_t</tt>.
    * \note The fact that \c UINT32_MAX is returned when the integer is not
    * found in the array means that the practical, usable maximum length of
    * a <tt>int32_array_t</tt> is <tt>UINT32_MAX - 1</tt> and not
    * <tt>UINT32_MAX</tt>.
    * \note If the array is already sorted, the more efficient function
    * <tt>index_in_sorted_int32_array</tt> can be used as it uses a basic
    * binary search instead of a complete scanning of the array.
    *
    * \param[in] to_find The integer to find in the <tt>int32_array_t</tt>.
    * \param[in] array   A pointer to the <tt>int32_array_t</tt>.
    * \returns The index of \c to_find in the array if \c to_find is found.
    * \returns <tt>UINT32_MAX</tt> otherwise.
    */
uint32_t index_in_int32_array(int32_t to_find,
                              const int32_array_t* const array);

   /**
    * \brief Returns true if a given integer is in a given array.
    *
    * Returns true if the integer \c to_find is in the <tt>int32_array_t</tt>
    * pointed to by \c array. Returns false otherwise.
    *
    * \note If the array is already sorted, the more efficient function
    * <tt>is_in_sorted_int32_array</tt> can be used as it uses a basic
    * binary search instead of a complete scanning of the array.
    *
    * \param[in] to_find The integer to find in the <tt>int32_array_t</tt>.
    * \param[in] array   A pointer to the <tt>int32_array_t</tt>.
    * \returns true if \c to_find is in the array \c array.
    * \returns false otherwise.
    */
bool is_in_int32_array(int32_t to_find, const int32_array_t* const array);

   /**
    * \brief Returns the position of an integer in a sorted portion of a
    * <tt>int32_array_t</tt>.
    *
    * Returns the position of the integer \c to_find in a \e sorted portion
    * of the <tt>int32_array_t</tt> pointed to by \c array. If the integer
    * \c to_find is not found in the portion delimited by \c min_index and
    * \c max_index, returns <tt>UINT32_MAX</tt>.
    *
    * \note The <tt>UINT32_MAX</tt> value is actually -1 if interpreted as a
    * signed <tt>int32_t</tt>.
    * \note The fact that \c UINT32_MAX is returned when the integer is not
    * found in the array means that the practical, usable maximum length of
    * a <tt>int32_array_t</tt> is <tt>UINT32_MAX - 1</tt> and not
    * <tt>UINT32_MAX</tt>.
    *
    * \param[in] to_find The integer to find in the <tt>int32_array_t</tt>.
    * \param[in] sorted_array A pointer to the <tt>int32_array_t</tt>.
    * \param[in] min_index The beginning of the sorted array portion to
    *                      search in.
    * \param[in] max_index The end of the sorted array portion to search in.
    * \returns The index of \c to_find in the array if \c to_find is found in
    *          the sorted array portion.
    * \returns <tt>UINT32_MAX</tt> otherwise.
    */
uint32_t index_in_sorted_int32_array(int32_t to_find,
                                     const int32_array_t* const sorted_array,
                                     uint32_t min_index, uint32_t max_index);

   /**
    * \brief Returns true if a given integer is in a given array.
    *
    * Returns true if the integer \c to_find is in the (\e already \e sorted)
    * <tt>int32_array_t</tt> pointed to by \c array. Returns false otherwise.
    *
    * \param[in] to_find The integer to find in the <tt>int32_array_t</tt>.
    * \param[in] array   A pointer to the \e sorted <tt>int32_array_t</tt>.
    * \returns true if \c to_find is in the array \c array.
    * \returns false otherwise.
    */
bool is_in_sorted_int32_array(int32_t to_find,
                              const int32_array_t* const array);

/*
 *-----------------------------------------------------------------------------
 *                  mpz_array_t and associated functions
 *-----------------------------------------------------------------------------
 */

   /**
    * \struct struct_mpz_array_t array.h lib/utils/include/array.h
    * \brief  Defines an array of <tt>mpz_t</tt> elements from the GMP
    * library.
    *
    * This structure defines a special kind of \c mpz array which knows
    * its current length and its allocated memory space.
    */
struct struct_mpz_array_t {
       /**
        * Memory space allocated for this array's \c data field, given as
        * a multiple of <tt>sizeof(mpz_t)</tt>. This is the maximum
        * number of \c mpz_t elements that the array can accommodate.
        */
    uint32_t alloced;
       /**
        * Current number of \c mpz_t elements hold in the array pointed by the
        * structure's \c data field.
        */
    uint32_t length;
       /**
        * Array of \c mpz_t elements whose size is given by the \c alloced
        * field.
        */
    mpz_t* data;
};
   /**
    * \typedef mpz_array_t
    * \brief Equivalent to <tt>struct struct_mpz_array_t</tt>.
    */
typedef struct struct_mpz_array_t mpz_array_t;

   /**
    * \brief Allocates and returns a new <tt>mpz_array_t</tt>.
    *
    * Allocates and returns a new <tt>mpz_array_t</tt> such that:
    * \li its \c alloced field is set to the parameter length.
    * \li its \c length field is set to zero.
    * \li its \c data array is left \e uninitialized.
    *
    * \param[in] length The maximum length of the \c mpz_array_t to allocate.
    * \return A pointer to the newly allocated \c mpz_array_t structure.
    */
mpz_array_t* alloc_mpz_array(uint32_t length);

   /**
    * \brief Clears a <tt>mpz_array_t</tt>.
    *
    * Clears a <tt>mpz_array_t</tt>, or, more precisely, clears the memory
    * space used by the array pointed by the \c data field of a
    * <tt>mpz_array_t</tt>. Also set its \c alloced and \c length fields
    * to zero.
    *
    * \param[in] array A pointer to the <tt>mpz_array_t</tt> to clear.
    */
void clear_mpz_array(mpz_array_t* const array);

   /**
    * \brief Resizes the allocated memory of an <tt>mpz_array_t</tt>.
    *
    * Resizes the storage available to an <tt>mpz_array_t</tt> to make room
    * for \c alloced integers, while preserving its content. If \c alloced is
    * less than the length of the array, then obviously some of its content
    * will be freed and lost.
    *
    * \param[in] alloced The new maximum length of the \c mpz_array_t to
    *                    resize.
    * \return array A pointer to the \c mpz_array_t to resize.
    */
void resize_mpz_array(mpz_array_t* const array, uint32_t alloced);

   /**
    * \brief Resets an <tt>mpz_array_t</tt>.
    *
    * Resets an <tt>mpz_array_t</tt>, i.e.:
    *
    *   \li Clears all of its \c mpz_t elements via \c mpz_clear;
    *   \li Sets its \c length field to zero.
    *
    * Note that its \c alloced field is left unchanged and that memory for
    * \c alloced \c mpz_t elements is still allocated (even tough the elements
    * still have to be initialized).
    *
    * \param[in] array A pointer to the <tt>mpz_array_t</tt> to clear.
    */
void reset_mpz_array(mpz_array_t* const array);

   /**
    * \brief Swaps two <tt>mpz_array_t</tt>'s contents.
    *
    * Swaps the contents of \c a and <tt>b</tt>, two <tt>mpz_array_t</tt>'s.
    *
    * \note In some case, pointer swapping is inappropriate (for example, if
    * the pointers are passed as function arguments!), hence the need for such
    * a swapping function.
    *
    * \param[in] a A pointer to the first <tt>mpz_array_t</tt> to swap.
    * \param[in] b A pointer to the second <tt>mpz_array_t</tt> to swap.
    */
void swap_mpz_array(mpz_array_t* const a, mpz_array_t* const b);

   /**
    * \brief Appends an \c mpz_t to an <tt>mpz_array_t</tt>.
    *
    * Appends the \c mpz_t integer \c to_append to <tt>array</tt>. If \c array
    * has not enough capacity to accommodate this extra element it will be
    * resized via a call to \c resize_mpz_array adding \c ELONGATION
    * \c mpz_t slots to avoid too frequent resizes.
    *
    * \param[in] array A pointer to the recipient <tt>mpz_array_t</tt>.
    * \param[in] to_append The <tt>mpz_t</tt> to append.
    */
void append_mpz_to_array(mpz_array_t* array, mpz_t to_append);

   /**
    * \brief Appends the content of an <tt>mpz_array_t</tt> to another one.
    *
    * Appends the content of the \c to_append array to the \c mpz_array_t
    * named \c array. If \c array has not enough capacity to accommodate all
    * elements from \c to_append, it will be resized via a call to
    * \c resize_mpz_array with extra room for \c ELONGATION unused \c mpz_t
    * slots to avoid too frequent resizes.
    *
    * \param[in] array A pointer to the recipient <tt>mpz_array_t</tt>.
    * \param[in] to_append A pointer to the <tt>mpz_array_t</tt> to append.
    */
void append_mpz_array(mpz_array_t* const array,
                      const mpz_array_t* const to_append);

   /**
    * \brief Prints a <tt>mpz_array_t</tt>.
    *
    * Prints a <tt>mpz_array_t</tt>'s \c data elements on the standard
    * output.
    *
    * \note This function is mostly intended for debugging purposes as the
    * output is not particularly well structured.
    *
    * \param[in] array A pointer to the <tt>mpz_array_t</tt> to print.
    */
void print_mpz_array(const mpz_array_t* const array);

   /**
    * \brief Returns the position of a \c mpz_t in a <tt>mpz_array_t</tt>.
    *
    * Returns the position of the \c mpz_t \c to_find in the
    * <tt>mpz_array_t</tt> pointed to by \c array. If the integer \c to_find
    * is not found in the <tt>mpz_array_t</tt>, returns <tt>UINT32_MAX</tt>.
    *
    * \note The <tt>UINT32_MAX</tt> value is actually -1 if interpreted as a
    * signed <tt>int32_t</tt>.
    * \note The fact that \c UINT32_MAX is returned when the integer is not
    * found in the array means that the practical, usable maximum length of
    * a <tt>mpz_array_t</tt> is <tt>UINT32_MAX - 1</tt> and not
    * <tt>UINT32_MAX</tt>.
    *
    * \param[in] to_find The \c mpz_t integer to find in the
    *                    <tt>mpz_array_t</tt>.
    * \param[in] array   A pointer to the <tt>mpz_array_t</tt>.
    * \returns The index of \c to_find in the array if \c to_find is found.
    * \returns <tt>UINT32_MAX</tt> otherwise.
    */
uint32_t index_in_mpz_array(const mpz_t to_find,
                            const mpz_array_t* const array);
   /**
    * \brief Returns the position of an \c mpz_t in a \e sorted portion of an
    * <tt>mpz_array_t</tt>.
    *
    * Returns the position of the \c mpz_t \c to_find in the \e sorted portion
    * of the <tt>mpz_array_t</tt> pointed to by \c array. If the integer \c
    * to_find is not found in this portion, returns <tt>UINT32_MAX</tt>.
    *
    * \note The <tt>UINT32_MAX</tt> value is actually -1 if interpreted as a
    * signed <tt>int32_t</tt>.
    * \note The fact that \c UINT32_MAX is returned when the integer is not
    * found in the array means that the practical, usable maximum length of
    * a <tt>mpz_array_t</tt> is <tt>UINT32_MAX - 1</tt> and not
    * <tt>UINT32_MAX</tt>.
    *
    * \param[in] to_find The \c mpz_t integer to find in the \c mpz_array_t.
    * \param[in] array   A pointer to the \e sorted <tt>mpz_array_t</tt>.
    * \param[in] min_index The beginning of the sorted array portion to
    *                      search in.
    * \param[in] max_index The end of the sorted array portion to search in.
    * \returns The index of \c to_find in the array if \c to_find is found.
    * \returns <tt>UINT32_MAX</tt> otherwise.
    */
uint32_t index_in_sorted_mpz_array(const mpz_t to_find,
                                   const mpz_array_t* const sorted_array,
                                   uint32_t min_index, uint32_t max_index);

   /**
    * \brief Returns true if a given integer is in a given array.
    *
    * Returns true if the \c mpz_t integer \c to_find is in the
    * <tt>mpz_array_t</tt> pointed to by \c array. Returns false otherwise.
    *
    * \note If the array is already sorted, the more efficient function
    * <tt>is_in_sorted_mpz_array</tt> can be used as it uses a basic
    * binary search instead of a complete scanning of the array.
    *
    * \param[in] to_find The integer to find in the <tt>mpz_array_t</tt>.
    * \param[in] array   A pointer to the <tt>mpz_array_t</tt>.
    * \returns true if \c to_find is in the array \c array.
    * \returns false otherwise.
    */
bool is_in_mpz_array(const mpz_t to_find, const mpz_array_t* const array);

   /**
    * \brief Sorts the mpz_t elements of a <tt>mpz_array_t</tt>.
    *
    * Sorts the mpz_t elements of a <tt>mpz_array_t</tt> in natural order using
    * a basic insertion sort.
    *
    * \param[in] array A pointer to the <tt>mpz_array_t</tt> to sort.
    */
void ins_sort_mpz_array(mpz_array_t* const array);

   /**
    * \brief Sorts the mpz_t elements of a <tt>mpz_array_t</tt> with a quick
    * sort.
    *
    * Sorts the mpz_t elements of a <tt>mpz_array_t</tt> in natural order using
    * the quick sort algorithm.
    *
    * \note This function relies on the C library implementation of the
    * quick sort provided by the function <tt>qsort</tt>.
    *
    * \param[in] array A pointer to the <tt>mpz_array_t</tt> to sort.
    */
void qsort_mpz_array(mpz_array_t* const array);

   /**
    * \brief Returns true if a given integer is in a given array.
    *
    * Returns true if the \c mpz_t integer \c to_find is in the
    * <tt>mpz_array_t</tt> pointed to by \c array. Returns false otherwise.
    *
    * \param[in] to_find The integer to find in the <tt>mpz_array_t</tt>.
    * \param[in] array   A pointer to the \e sorted <tt>mpz_array_t</tt>.
    * \returns true if \c to_find is in the array \c array.
    * \returns false otherwise.
    */
bool is_in_sorted_mpz_array(const mpz_t to_find,
                            const mpz_array_t* const array);

/*
 *-----------------------------------------------------------------------------
 *                  binary_array_t and associated functions
 *-----------------------------------------------------------------------------
 */

   /**
    * \struct struct_binary_array_t array.h lib/utils/include/array.h
    * \brief  Defines an array of bits.
    *
    * This structure defines an array of bits which knows
    * its current length and its allocated memory space.
    *
    * \note Internally, bits are packed in a \c TIFA_BITSTRING_T array.
    */
struct struct_binary_array_t {
      /**
       * Memory space allocated for this array's \c data field, given as
       * a multiple of <tt>sizeof(TIFA_BITSTRING_T)</tt>. This is the maximum
       * number of \c TIFA_BITSTRING_T that the array can accommodate.
       * The number of bits that the array can hold is hence
       * <tt>CHAR_BIT * sizeof(TIFA_BITSTRING_T)</tt> times this value
       * (\c CHAR_BIT being the number of bits used to represent a
       * <tt>char</tt>, usually 8 on most current architectures).
       */
   uint32_t alloced;
      /**
       * Current number of bits hold in the array pointed by the
       * structure's \c data field.
       */
   uint32_t length;
      /**
       * Array of \c TIFA_BITSTRING_T whose size is given by the
       * \c alloced field.
       */
   TIFA_BITSTRING_T* data;
};
   /**
    * \typedef binary_array_t
    * \brief Equivalent to <tt>struct struct_binary_array_t</tt>.
    */
typedef struct struct_binary_array_t binary_array_t;

   /**
    * \brief Allocates and returns a new <tt>binary_array_t</tt>.
    *
    * Allocates and returns a new <tt>binary_array_t</tt> such that:
    * \li its \c alloced field is set to the minimum number of
    * \c TIFA_BITSTRING_T variables needed to store length bits.
    * \li its \c length field is set to zero.
    * \li its \c data array is completely filled with zeroes.
    *
    * \param[in] length The maximum bitlength of the \c uint32_array_t to
    * allocate.
    * \return A pointer to the newly allocated \c uint32_array_t structure.
    * Note that this array may hold more that length bits if length is not a
    * multiple of 8 * <tt>sizeof(TIFA_BITSTRING_T)</tt>.
    */
binary_array_t* alloc_binary_array(uint32_t length);

   /**
    * \brief Clears a <tt>binary_array_t</tt>.
    *
    * Clears a <tt>binary_array_t</tt>, or, more precisely, clears the memory
    * space used by the array pointed by the \c data field of a
    * <tt>binary_array_t</tt>. Also set its \c alloced and \c length fields
    * to zero.
    *
    * \param[in] array A pointer to the <tt>binary_array_t</tt> to clear.
    */
void clear_binary_array(binary_array_t* const array);

   /**
    * \brief Resizes the allocated memory of a <tt>binary_array_t</tt>.
    *
    * Resizes the storage available to an <tt>binary_array_t</tt> to make room
    * for \c alloced integers, while preserving its content. If \c alloced is
    * less than the length of the array, then obviously some of its content
    * will be lost.
    *
    * \param[in] alloced The new maximum length of the \c binary_array_t to
    *                    resize.
    * \return array A pointer to the \c binary_array_t to resize.
    */
void resize_binary_array(binary_array_t* const array, uint32_t alloced);

   /**
    * \brief Appends a bit to a <tt>binary_array_t</tt>.
    *
    * Appends a bit (set to one if <tt>to_append != 0</tt>, set to zero
    * otherwise) to <tt>array</tt>. If \c array has not enough capacity to
    * accommodate this extra element it will be resized via a call to
    * \c resize_binary_array adding <tt>ELONGATION * BITSTRING_T_BITSIZE</tt>
    * bit slots to avoid too frequent resizes.
    *
    * \param[in] array A pointer to the recipient <tt>binary_array_t</tt>.
    * \param[in] to_append The bit to append (1 if <tt>to_append != 0</tt>,
    *                      0 otherwise).
    */
void append_binary_to_array(binary_array_t* array, unsigned int to_append);


   /**
    * \brief Prints a <tt>binary_array_t</tt>.
    *
    * Prints a <tt>binary_array_t</tt>'s on the standard output.
    *
    * \param[in] array A pointer to the <tt>binary_array_t</tt> to print.
    */
void print_binary_array(const binary_array_t* const array);

   /**
    * \brief Returns the value of a given bit in a <tt>binary_array_t</tt>.
    *
    * Returns the value of the <tt>index</tt>-th bit of the
    * <tt>binary_array_t</tt> pointed to by \c array, as either 0 or 1.
    *
    * \param[in] index The position of the bit to read.
    * \param[in] array A pointer to the <tt>binary_array_t</tt>.
    * \return The value of the <tt>index</tt>-th bit: either 0 or 1.
    */
uint8_t get_array_bit(uint32_t index, const binary_array_t* const array);

   /**
    * \brief Sets a given bit to one in a <tt>binary_array_t</tt>.
    *
    * Sets the <tt>index</tt>-th bit of the
    * <tt>binary_array_t</tt> pointed to by \c array to 1.
    *
    * \param[in] index The position of the bit to set.
    * \param[in] array A pointer to the <tt>binary_array_t</tt>.
    */
void set_array_bit_to_one(uint32_t index, binary_array_t* const array);

   /**
    * \brief Sets a given bit to zero in a <tt>binary_array_t</tt>.
    *
    * Sets the <tt>index</tt>-th bit of the
    * <tt>binary_array_t</tt> pointed to by \c array to 0.
    *
    * \param[in] index The position of the bit to set.
    * \param[in] array A pointer to the <tt>binary_array_t</tt>.
    */
void set_array_bit_to_zero(uint32_t index, binary_array_t* const array);

   /**
    * \brief Flips a given bit to zero in a <tt>binary_array_t</tt>.
    *
    * Flips the <tt>index</tt>-th bit of the
    * <tt>binary_array_t</tt> pointed to by \c array.
    *
    * \param[in] index The position of the bit to flip.
    * \param[in] array A pointer to the <tt>binary_array_t</tt>.
    */
void flip_array_bit(uint32_t index, binary_array_t* const array);

#ifdef __cplusplus
}
#endif

#endif
