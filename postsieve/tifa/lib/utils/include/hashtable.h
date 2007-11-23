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
 * \file    hashtable.h
 * \author  Jerome Milan
 * \date    Thu Mar 30 2006
 * \version 1.0
 *
 * \brief Generic hashtable.
 *
 * Yet another implementation of a generic hashtable.
 */

 /*
  *  License: GNU Lesser General Public License (LGPL)
  *  History:
  */

#if !defined(_TIFA_HASHTABLE_H_)
   /**
    * \def _TIFA_HASHTABLE_H_
    * Standard include guard.
    */
#define _TIFA_HASHTABLE_H_

#include <inttypes.h>
#include "linked_list.h"

#ifdef __cplusplus
extern "C" {
#endif

   /**
    * \struct struct_hashtable_t hashtable.h lib/utils/include/hashtable.h
    * \brief  A basic implementation of a hashtable.
    *
    * This structure defines a simple generic hashtable. It can store any
    * type of elements, provided that suitable comparison and hash functions
    * exist for the type of the keys used.
    *
    * This hashtable implementation uses a simple sequential search in a linked
    * list to solve the collisions.
    *
    * \warning Due to limitations in the current implementation, it is
    * strongly advised to use \e only pointers to integers or strings as
    * keys.
    */
struct struct_hashtable_t {
       /**
        * Number of allocated buckets (always a power of two).
        */
    uint32_t alloced;
       /**
        * Current number of entries in the hashtable.
        */
    uint32_t nb_entries;
       /**
        * Array of \c linked_list_t of size <tt>alloced</tt> used to store the
        * hashtable's entries.
        */
    linked_list_t* buckets;
       /**
        * Pointer to a comparison function for the keys stored
        * in the hashtable. The function's signature should be:
        * <tt>int cmp_func(const void* const, const void* const)</tt>. The
        * function should take two keys of \e identical type passed as
        * pointers to <tt>void</tt>.
        *
        * For now, the only requirement if that the pointed function should
        * return 0 if two identical keys are compared.
        */
    int      (*cmp_func)  (const void* const key_a, const void* const key_b);
       /**
        * Pointer to the hash function used by the hashtable. The hash
        * function's signature should be:
        * <tt>uint32_t hash_func(const void* const)</tt>. The
        * function should take a key passed as a pointer to <tt>void</tt>.
        *
        * Needless to say, the real type of the key handled by this function
        * should be the same as the one handled by the comparison function
        * pointed by <tt>cmp_int</tt>.
        *
        */
    uint32_t (*hash_func) (const void* const key);
};
   /**
    * \typedef hashtable_t
    * \brief Equivalent to <tt>struct struct_hashtable_t</tt>.
    */
typedef struct struct_hashtable_t hashtable_t;

   /**
    * \struct struct_hashtable_entry_t lib/utils/include/hashtable.h
    * \brief  The structure of a hashtable's entry.
    *
    * This structure defines a hashtable entry, i.e. some data and
    * its associated key.
    */
struct struct_hashtable_entry_t {
    /**
     * Key associated to this entry, as a pointer to <tt>void</tt>.
     */
    void* key;
    /**
     * Data associated to this entry, as a pointer to <tt>void</tt>.
     */
    void* data;
};
   /**
    * \typedef hashtable_entry_t
    * \brief Equivalent to <tt>struct struct_hashtable_entry_t</tt>.
    */
typedef struct struct_hashtable_entry_t hashtable_entry_t;

   /**
    * \brief Allocates and returns a new <tt>hashtable_t</tt>.
    *
    * Allocates and returns a pointer to a new <tt>hashtable_t</tt> with
    * \c size allocated buckets, using the comparison function pointed by
    * \c cmp_func and the hash function given by \c hash_func.
    *
    * \note If \c size is not a power of two, the lowest power of two greater
    * than \c size is used instead.
    *
    * \param[in] size The number of allocated buckets.
    * \param[in] cmp_func A pointer to the comparison function.
    * \param[in] hash_func A pointer to the hash function.
    */
hashtable_t* alloc_init_hashtable(
    uint32_t size,
    int      (*cmp_func)  (const void* const key_a, const void* const key_b),
    uint32_t (*hash_func) (const void* const key)
);

   /**
    * \brief Clears a <tt>hashtable_t</tt>.
    *
    * Clears the <tt>hashtable_t</tt> pointed by \c htable.
    *
    * \warning This function merely clears the memory used by the linked
    * lists but \e does \e not frees the memory space used by the \c key
    * and \c data pointers of the <tt>hashtable_entry_t</tt>.
    *
    * \param[in] htable A pointer to the <tt>hashtable_t</tt> to clear.
    */
void clear_hashtable(hashtable_t* const htable);

   /**
    * \brief Adds an entry in a <tt>hashtable_t</tt>.
    *
    * Creates an entry with its \c data field given by the pointer \c data
    * and its \c key field given by the pointer \c key and adds it in a
    * <tt>hashtable_t</tt>.
    *
    * \warning This function does not copy the actual content of the variable
    * referenced by \c key and \c data but merely copies these references in
    * the \c hashtable_entry_t structure.
    *
    * \param[in] htable A pointer to the <tt>hashtable_t</tt>.
    * \param[in] key A pointer to the key of the new entry.
    * \param[in] data A pointer to the data of the new entry.
    */
void add_entry_in_hashtable(hashtable_t* const htable,
                            const void* const key, const void* const data);
   /**
    * \brief Gets an entry's data from a <tt>hashtable_t</tt>.
    *
    * Gets the \c data field of the entry from the hashtable \c htable whose
    * key is given by <tt>key</tt>.
    *
    * \param[in] htable A pointer to the <tt>hashtable_t</tt>.
    * \param[in] key A pointer to the key of the entry to read.
    * \return The \c data field of the entry whose key is given by <tt>key</tt>,
    *         if any.
    * \return NULL if no such entry is found in the hashtable.
    */
void* get_entry_in_hashtable(hashtable_t* const htable,
                             const void* const key);

   /**
    * \brief Removes an entry from a <tt>hashtable_t</tt>.
    *
    * Removes the entry from the hashtable \c htable whose
    * key is given by <tt>key</tt> and returns its \c data field.
    *
    * \warning The \c key field of the removed entry is freed if a matching
    * entry is found. The means that if \c key is a pointer to a structure
    * containing some other pointers, all the memory may not be freed since
    * the real type of \c key is not known by the hashtable. This
    * is why it is strongly recommended to use \e only pointers to integers or
    * strings as keys.
    *
    * \param[in] htable A pointer to the <tt>hashtable_t</tt>.
    * \param[in] key A pointer to the key of the entry to remove.
    * \return The \c data field of the removed entry.
    * \return NULL if no matching entry is found in the hashtable.
    */
void* remove_entry_in_hashtable(hashtable_t* const htable,
                                const void* const key);

#ifdef __cplusplus
}
#endif

#endif
