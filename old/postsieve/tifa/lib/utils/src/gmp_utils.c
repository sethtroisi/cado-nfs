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
 * \file    gmp_utils.c
 * \author  Jerome Milan
 * \date    Mon Mar 10 2008
 * \version 1.0.2
 */

/*
 * History:
 * --------
 *    1.0.2: Mon Mar 10 2008 by JM
 *        - Inlined some functions.
 *    1.0.1: Wed Mar 7 2007 by JM
 *        - Added *_mpzpair_htable functions.
 *    1.0: Wed Mar 1 2006 by JM
 *        - Initial version.
 */

#include <stdlib.h>
#include <gmp.h>

#include "gmp_utils.h"

//------------------------------------------------------------------------------
void empty_mpzpair_htable(hashtable_t* const htable) {
    //
    // Empties the hashtable and all of its entries... It's a bit messy
    // so should we...
    //
    // _NOTE_: Consider adding in the hashtable_t structure (and hence
    //         in the linked_list_t structure) a pointer to a function
    //         clearing the memory space used by the entries?
    //
    for (uint32_t i = 0; i < htable->alloced; i++) {
        linked_list_node_t *node = htable->buckets[i].head;
        linked_list_node_t *next = NULL;
        while (NULL != node) {
            hashtable_entry_t* entry = (hashtable_entry_t*)node->data;
            clear_mpz_pair((mpz_pair_t*)entry->data);
            free(entry->data);
            free(entry->key);
            free(entry);
            next = node->next;
            free(node);
            node = next;
        }
        htable->buckets[i].head   = NULL;
        htable->buckets[i].tail   = NULL;
        htable->buckets[i].length = 0;
    }
}
//------------------------------------------------------------------------------
