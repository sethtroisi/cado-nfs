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
 * \date    Circa Fri Mar 29 2008
 * \version 1.2.3
 *
 * \brief Higher level arrays and associated functions.
 *
 * This file defines higher level arrays together with some associated
 * functions.
 *
 * The <tt>*_array_t</tt> types and their associated functions are quite
 * similar, the only differences being the type of the elements these arrays
 * hold. Each <tt>*_array_t</tt> type is a structure composed of three fields:
 *
 * \li \c alloced - The maximum number of element the array can accomodate
 * \li \c length - The current number of element in the array
 * \li \c data - A pointer to the allocated memory space of \c alloced elements
 *
 * \warning Since version 1.2.1 memory management changed. See the
 * \c alloc_*_array and \c clear_*_array functions for more information.
 */

 /*
  * History:
  * --------
  * 1.2.3: Fri Mar 29 (?) 2008 by JM:
  *        - Inlined functions pertaining to binary_array_t's.
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

#if !defined(_TIFA_ARRAY_H_)
   /**
    * \def _TIFA_ARRAY_H_
    * Standard include guard.
    */
#define _TIFA_ARRAY_H_

#include <inttypes.h>
#include <stdbool.h>
#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

   /**
    * \def ELONGATION
    * Incremental size used when automatically expanding the capacity of a
    * \c *_array_t.
    *
    * \note This is, of course, an hint to GMP's limbs and nails :-)
    */
#define ELONGATION 16

   /**
    * \def NOT_IN_ARRAY
    * Value returned by the <tt>index_in_*_array(x, array, ...)</tt>
    * functions if the element \c x is not in the array <tt>array</tt>.
    */
#define NOT_IN_ARRAY UINT32_MAX

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
    * Clears the <tt>uint32_array_t</tt> pointed to by \c array, \e i.e.
    * clears the memory space used by the C-style array pointed by
    * <tt>array->data</tt> and frees the \c array pointer.
    *
    * \warning Before version 1.2.1, the \c array pointer was not freed which
    *          required explicit calls to \c free(...) in client code.
    *
    * \param[in] array A pointer to the <tt>uint32_array_t</tt> to clear.
    */
void clear_uint32_array(uint32_array_t* array);

   /**
    * \def reset_uint32_array(ARRAY)
    * \brief Resets a <tt>uint32_array_t</tt>.
    *
    * Resets the \c length field of \c array to zero.
    *
    * Note that its \c alloced field is left unchanged and that memory for
    * \c alloced \c uint32_t elements is still allocated.
    *
    * \param[in] array A pointer to the <tt>uint32_array_t</tt> to reset.
    */
#define reset_uint32_array(ARRAY) do {(ARRAY)->length = 0;} while (0)

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
void append_uint32_to_array(uint32_array_t* array, const uint32_t to_append);

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
    * is not found in the <tt>uint32_array_t</tt>, returns
    * <tt>NOT_IN_ARRAY</tt>.
    *
    * \note The <tt>NOT_IN_ARRAY</tt> value is actually -1 if interpreted as a
    * signed <tt>int32_t</tt>.
    *
    * \note If the array is already sorted, the more efficient function
    * <tt>index_in_sorted_uint32_array</tt> can be used as it uses a basic
    * binary search instead of a complete scanning of the array.
    *
    * \param[in] to_find The integer to find in the <tt>uint32_array_t</tt>.
    * \param[in] array   A pointer to the <tt>uint32_array_t</tt>.
    * \returns The index of \c to_find in the array if \c to_find is found.
    * \returns <tt>NOT_IN_ARRAY</tt> otherwise.
    */
uint32_t index_in_uint32_array(uint32_t to_find,
                               const uint32_array_t* const array);

   /**
    * \brief Returns the position of an integer in a sorted portion of a
    * <tt>uint32_array_t</tt>.
    *
    * Returns the position of the integer \c to_find in a \e sorted portion
    * of the <tt>uint32_array_t</tt> pointed to by \c array. If the integer
    * \c to_find is not found in the portion delimited by \c min_index and
    * \c max_index, returns <tt>NOT_IN_ARRAY</tt>.
    *
    * \note The <tt>NOT_IN_ARRAY</tt> value is actually -1 if interpreted as a
    * signed <tt>int32_t</tt>.
    *
    * \param[in] to_find The integer to find in the <tt>uint32_array_t</tt>.
    * \param[in] sorted_array A pointer to the <tt>uint32_array_t</tt>.
    * \param[in] min_index The beginning of the sorted array portion to
    *                      search in.
    * \param[in] max_index The end of the sorted array portion to search in.
    * \returns The index of \c to_find in the array if \c to_find is found in
    *          the sorted array portion.
    * \returns <tt>NOT_IN_ARRAY</tt> otherwise.
    */
uint32_t index_in_sorted_uint32_array(uint32_t to_find,
                                      const uint32_array_t* const sorted_array,
                                      uint32_t min_index, uint32_t max_index);

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
        * Current number of "useful" \c mpz_t elements hold in the array
        * pointed by the structure's \c data field.
        *
        * \warning Prior to version 1.2, the \c length field also indicated
        * which positions had been \c mpz_init'ed in the \c data field. Since
        * version 1.2 this is no longer true. Now all positions in the \c data
        * array are \c mpz_init'ed and \c length only gives which part of the
        * array is useful from the client standpoint.
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
    * \li its \c data array is fully \c mpz_init'ed.
    *
    * \param[in] length The maximum length of the \c mpz_array_t to allocate.
    * \return A pointer to the newly allocated \c mpz_array_t structure.
    *
    * \warning Since version 1.2, the \c data field is completely
    * \c mpz_init'ed (from \c data[0] to \c data[\c alloced -1]) whereas
    * older versions did not \c mpz_init anything. This change in behaviour
    * was prompted by the need to avoid multiple memory deallocations and
    * reallocations when using the same \c mpz_array_t repeatedly.
    */
mpz_array_t* alloc_mpz_array(uint32_t length);

   /**
    * \brief Clears a <tt>mpz_array_t</tt>.
    *
    * Clears the <tt>mpz_array_t</tt> pointed to by \c array, \e i.e.
    * clears the memory space used by the C-style array pointed by
    * <tt>array->data</tt> and frees the \c array pointer.
    *
    * \warning Before version 1.2.1, the \c array pointer was not freed which
    *          required explicit calls to \c free(...) in client code.
    *
    * \param[in] array A pointer to the <tt>mpz_array_t</tt> to clear.
    */
void clear_mpz_array(mpz_array_t* array);

   /**
    * \brief Resets an <tt>mpz_array_t</tt>.
    *
    * Resets the \c length field of \c array to zero.
    *
    * Note that its \c alloced field is left unchanged and that memory for
    * \c alloced \c mpz_t elements is still allocated (all the elements
    * remaining fully \c mpz_init'ed).
    *
    * \warning Prior to 1.2 when the semantic was different, this function used
    * to \c mpz_clear all positions in \c array->data. This is no longer true.
    *
    * \param[in] array A pointer to the <tt>mpz_array_t</tt> to clear.
    */
#define reset_mpz_array(ARRAY) do {(ARRAY)->length = 0;} while (0)

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
void append_mpz_to_array(mpz_array_t* array, const mpz_t to_append);

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
    * is not found in the <tt>mpz_array_t</tt>, returns <tt>NOT_IN_ARRAY</tt>.
    *
    * \note The <tt>NOT_IN_ARRAY</tt> value is actually -1 if interpreted as a
    * signed <tt>int32_t</tt>.
    *
    * \param[in] to_find The \c mpz_t integer to find in the
    *                    <tt>mpz_array_t</tt>.
    * \param[in] array   A pointer to the <tt>mpz_array_t</tt>.
    * \returns The index of \c to_find in the array if \c to_find is found.
    * \returns <tt>NOT_IN_ARRAY</tt> otherwise.
    */
uint32_t index_in_mpz_array(const mpz_t to_find,
                            const mpz_array_t* const array);
   /**
    * \brief Returns the position of an \c mpz_t in a \e sorted portion of an
    * <tt>mpz_array_t</tt>.
    *
    * Returns the position of the \c mpz_t \c to_find in the \e sorted portion
    * of the <tt>mpz_array_t</tt> pointed to by \c array. If the integer \c
    * to_find is not found in this portion, returns <tt>NOT_IN_ARRAY</tt>.
    *
    * \note The <tt>NOT_IN_ARRAY</tt> value is actually -1 if interpreted as a
    * signed <tt>int32_t</tt>.
    *
    * \param[in] to_find The \c mpz_t integer to find in the \c mpz_array_t.
    * \param[in] sorted_array   A pointer to the \e sorted <tt>mpz_array_t</tt>.
    * \param[in] min_index The beginning of the sorted array portion to
    *                      search in.
    * \param[in] max_index The end of the sorted array portion to search in.
    * \returns The index of \c to_find in the array if \c to_find is found.
    * \returns <tt>NOT_IN_ARRAY</tt> otherwise.
    */
uint32_t index_in_sorted_mpz_array(const mpz_t to_find,
                                   const mpz_array_t* const sorted_array,
                                   uint32_t min_index, uint32_t max_index);

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

#ifdef __cplusplus
}
#endif

#endif
