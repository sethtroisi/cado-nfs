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
 * \file    hashtable.c
 * \author  Jerome Milan
 * \date    Thu Mar 30 2006
 * \version 1.0
 */

#include <stdlib.h>
#include <inttypes.h>

#include "hashtable.h"
#include "macros.h"
#include "funcs.h"

//
// _TO_DO_: Consider adding in the hashtable_t structure (and hence in the
//          linked_list_t structure) a pointer to a function clearing the
//          memory space used by the entries?
//

//-----------------------------------------------------------------------------
hashtable_t* alloc_init_hashtable(
                        uint32_t size,
                        int (*cmp_func) (const void* const key_a,
                                         const void* const key_b),
                        uint32_t (*hash_func) (const void* const key) ) {

    hashtable_t *htable = malloc(sizeof(hashtable_t));
    //
    // Make sure that htable->alloced is a power of 2 so that we can just use
    // the binary & operator instead of computing a modulo...
    //
    if (!IS_POWER_OF_2_UI(size)) {
        size = ceil_log2(size);
        size = 1 << size;
    }
    htable->alloced  = size;
    htable->nentries = 0;
    htable->buckets  = malloc(size * sizeof(linked_list_t));

    htable->hash_func = hash_func;
    htable->cmp_func  = cmp_func;

    for (uint32_t i = 0; i < size; i++) {
        init_linked_list(&(htable->buckets[i]), cmp_func);
    }
    return htable;
}
//-----------------------------------------------------------------------------
void clear_hashtable(hashtable_t* const htable) {
    //
    // _WARNING_: For the time being, clear_linked_list may not clear
    //            completely the memory used by the list's node. See
    //            comments at the beginning of this file and in the file
    //            linked_list.c
    //
    for (uint32_t i = 0; i < htable->alloced; i++) {
        clear_linked_list(&(htable->buckets[i]));
    }
    free(htable->buckets);
    htable->alloced    = 0;
    htable->nentries = 0;
}
//-----------------------------------------------------------------------------
void add_entry_in_hashtable(hashtable_t* const htable,
                            const void* const key,
                            const void* const data) {

    hashtable_entry_t* entry = malloc(sizeof(hashtable_entry_t));

    entry->key  = (void*)key;
    entry->data = (void*)data;

    uint32_t ind = htable->hash_func(key) & (htable->alloced - 1);

    append_to_linked_list(&(htable->buckets[ind]), (void*)entry);
}
//-----------------------------------------------------------------------------
void* get_entry_in_hashtable(hashtable_t* const htable,
                             const void* const key) {

    uint32_t ind = htable->hash_func(key) & (htable->alloced - 1);

    if (htable->buckets[ind].length == 0) {
        return NULL;
    }
    linked_list_node_t* node = htable->buckets[ind].head;
    hashtable_entry_t* entry = NULL;

    while (node != NULL) {
        entry = ((hashtable_entry_t*)node->data);
        if (htable->cmp_func(entry->key, key) == 0) {
            return entry->data;
        }
        node = node->next;
    }
    return NULL;
}
//-----------------------------------------------------------------------------
void* remove_entry_in_hashtable(hashtable_t* const htable,
                                const void* const  key) {

    uint32_t ind = htable->hash_func(key) & (htable->alloced - 1);

    if (htable->buckets[ind].length == 0) {
        return NULL;
    }
    linked_list_node_t* node = htable->buckets[ind].head;
    hashtable_entry_t* entry = NULL;

    while (node != NULL) {

        entry = ((hashtable_entry_t*)node->data);

        if (htable->cmp_func(entry->key, key) == 0) {
            remove_node_from_linked_list(&(htable->buckets[ind]), node);
            void* data = entry->data;
            free(entry->key);
            free(entry);
            free(node);
            htable->nentries--;
            return data;
        }
        node = node->next;
    }
    return NULL;
}
//-----------------------------------------------------------------------------
