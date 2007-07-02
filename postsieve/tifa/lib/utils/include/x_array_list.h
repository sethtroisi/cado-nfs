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
 * \file    x_array_list.h
 * \author  Jerome Milan
 * \date    Thu Mar 2 2006
 * \version 1.0
 *
 * \brief Higher level lists of arrays and associated functions.
 *
 * Defines higher level lists of arrays and their associated functions. This
 * terminology is actually a bit confusing since they are actually more
 * arrays of arrays rather than lists strictly speaking.
 */

 /*
  *  License: GNU Lesser General Public License (LGPL)
  *  History:
  */

#if !defined(_TIFA_X_ARRAY_LIST_H_)
   /**
    * \def _TIFA_X_ARRAY_LIST_H_
    * Standard include guard.
    */
#define _TIFA_X_ARRAY_LIST_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <inttypes.h>
#include "array.h"

/*
 *-----------------------------------------------------------------------------
 *              uint32_array_list_t and associated functions
 *-----------------------------------------------------------------------------
 */

   /**
    * \struct struct_uint32_array_list_t x_array_list.h
    *                                    lib/utils/include/x_array.h
    *
    * \brief  Defines a list of <tt>uint32_array_t</tt>.
    *
    * This structure defines an array of pointers to <tt>uint32_array_t</tt>
    * elements. Its name is a bit confusing since it is actually more of an
    * array than a list strictly speaking.
    */
struct struct_uint32_array_list_t {
    /**
     * This is the maximum number of \c uint32_array_t pointers that the array
     * can accommodate.
     */
  uint32_t alloced;
    /**
     * Current number of \c uint32_array_t pointers hold in the array pointed
     * by the structure's \c data field.
     */
  uint32_t length;
    /**
     * Array of pointers to \c uint32_array_t whose size is given by the
     * \c alloced field.
     */
  uint32_array_t** data;
};

   /**
    * \typedef uint32_array_list_t
    * \brief Equivalent to <tt>struct struct_uint32_array_list_t</tt>.
    */
typedef struct struct_uint32_array_list_t uint32_array_list_t;

   /**
    * \brief Allocates and returns a new <tt>uint32_array_list_t</tt>.
    *
    * Allocates and returns a new <tt>uint32_array_list_t</tt> such that:
    * \li its \c alloced field is set to the parameter alloced.
    * \li its \c length field set to zero.
    * \li its \c data array is left \e uninitialized.
    *
    * \param[in] alloced The allocated length of the \c uint32_array_list_t
    *                    to allocate.
    * \return A pointer to the newly allocated \c uint32_array_list_t structure.
    */
uint32_array_list_t* alloc_uint32_array_list(uint32_t alloced);

   /**
    * \brief Adds an entry to a <tt>uint32_array_list_t</tt>.
    *
    * Adds the \c entry pointer to \c list and increments its
    * \c length field.
    *
    * \warning This function tranfers the ownership of the \c uint32_array_t
    * pointed to by \c entry to <tt>list</tt>. This means that any client code
    * \e should \e not clear any \c uint32_array_t that has been added to a
    * <tt>uint32_array_list_t</tt> since this is the exclusive responsability of
    * the <tt>uint32_array_list_t</tt>.
    *
    * \param[in] entry A pointer to the array to add.
    * \param[in] list A pointer to the <tt>uint32_array_list_t</tt>.
    * \return A pointer to the newly allocated \c uint32_array_list_t structure.
    */
inline void add_entry_in_uint32_array_list(uint32_array_t* const entry,
                                           uint32_array_list_t* const list);

   /**
    * \brief Clears a <tt>uint32_array_list_t</tt>.
    *
    * Clears a <tt>uint32_array_list_t</tt>, or, more precisely, clears
    * the memory space used by the array pointed by the \c data field of a
    * <tt>uint32_array_list_t</tt>. Also set its \c alloced and \c length
    * fields to zero.
    *
    * \param[in] list A pointer to the <tt>uint32_array_list_t</tt> to clear.
    */
void clear_uint32_array_list(uint32_array_list_t* const list);

   /**
    * \brief Prints a <tt>uint32_array_list_t</tt>.
    *
    * Prints a <tt>uint32_array_list_t</tt> on the standard output.
    *
    * \note This function is mostly intended for debugging purposes as the
    * output is not particularly well structured
    *
    * \param[in] list A pointer to the <tt>uint32_array_list_t</tt> to print.
    */
void print_uint32_array_list(const uint32_array_list_t* const list);


/*
 *-----------------------------------------------------------------------------
 *              mpz_array_list_t and associated functions
 *-----------------------------------------------------------------------------
 */

   /**
    * \struct struct_mpz_array_list_t x_array_list.h
    *                                 lib/utils/include/x_array.h
    *
    * \brief  Defines a list of <tt>mpz_array_t</tt>.
    *
    * This structure defines an array of pointers to <tt>mpz_array_t</tt>
    * elements. Its name is a bit confusing since it is actually more of an
    * array than a list strictly speaking.
    */
struct struct_mpz_array_list_t {
      /**
       * This is the maximum number of \c mpz_array_t pointers that the array
       * can accommodate.
       */
    uint32_t alloced;
      /**
       * Current number of \c mpz_array_t pointers hold in the array pointed
       * by the structure's \c data field.
       */
    uint32_t length;
      /**
       * Array of pointers to \c mpz_array_t whose size is given by the
       * \c alloced field.
       */
    mpz_array_t** data;
};

   /**
    * \typedef mpz_array_list_t
    * \brief Equivalent to <tt>struct struct_mpz_array_list_t</tt>.
    */
typedef struct struct_mpz_array_list_t mpz_array_list_t;

   /**
    * \brief Allocates and returns a new <tt>mpz_array_list_t</tt>.
    *
    * Allocates and returns a new <tt>mpz_array_list_t</tt> such that:
    * \li its \c alloced field is set to the parameter alloced.
    * \li its \c length field set to zero.
    * \li its \c data array is left \e uninitialized.
    *
    * \param[in] alloced The allocated length of the \c mpz_array_list_t
    *                    to allocate.
    * \return A pointer to the newly allocated \c mpz_array_list_t structure.
    */
mpz_array_list_t* alloc_mpz_array_list(uint32_t alloced);

   /**
    * \brief Adds an entry to a <tt>mpz_array_list_t</tt>.
    *
    * Adds the \c entry pointer to \c list and increments its
    * \c length field.
    *
    * \warning This function tranfers the ownership of the \c mpz_array_t
    * pointed to by \c entry to <tt>list</tt>. This means that any client code
    * \e should \e not clear any \c mpz_array_t that has been added to a
    * <tt>mpz_array_list_t</tt> since this is the exclusive responsability of
    * the <tt>mpz_array_list_t</tt>.
    *
    * \param[in] entry A pointer to the \c mpz_array_t to add in the list.
    * \param[in] list A pointer to the <tt>mpz_array_list_t</tt>.
    * \return A pointer to the newly allocated \c mpz_array_list_t structure.
    */
inline void add_entry_in_mpz_array_list(mpz_array_t* const entry,
                                        mpz_array_list_t* const list);

   /**
    * \brief Clears a <tt>mpz_array_list_t</tt>.
    *
    * Clears a <tt>mpz_array_list_t</tt>, or, more precisely, clears
    * the memory space used by the array pointed by the \c data field of a
    * <tt>mpz_array_list_t</tt>. Also set its \c alloced and \c length
    * fields to zero.
    *
    * \param[in] list A pointer to the <tt>mpz_array_list_t</tt> to clear.
    */
void clear_mpz_array_list(mpz_array_list_t* const list);

   /**
    * \brief Prints a <tt>mpz_array_list_t</tt>.
    *
    * Prints a <tt>mpz_array_list_t</tt> on the standard output.
    *
    * \note This function is mostly intended for debugging purposes as the
    * output is not particularly well structured
    *
    * \param[in] list A pointer to the <tt>mpz_array_list_t</tt> to print.
    */
void print_mpz_array_list(const mpz_array_list_t* const list);

#ifdef __cplusplus
}
#endif

#endif
