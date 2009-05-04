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
 * \file    linked_list.c
 * \author  Jerome Milan
 * \date    Thu Mar 30 2006
 * \version 1.0
 */

#include <stdlib.h>
#include <inttypes.h>

#include "linked_list.h"

//
// _TO_DO_: Consider adding in the linked_list_t structure a pointer to a
//          function clearing completely the memory space referenced by
//          linked_list_node_t's data pointer ?
//

//-----------------------------------------------------------------------------
void init_linked_list(linked_list_t* const list,
                      int (*cmp_func) (const void* const data_a,
                                       const void* const data_b) ) {
    list->head     = NULL;
    list->tail     = NULL;
    list->cmp_func = cmp_func;
    list->length   = 0;
}
//-----------------------------------------------------------------------------
void clear_linked_list(linked_list_t* const list) {
    linked_list_node_t *node = list->head;
    linked_list_node_t *next = NULL;
    while (NULL != node) {
        //
        // _WARNING_: For the time being, we merely free node->data.
        //            However it is possible that the object referenced
        //            by this pointer contains other pointers that needs to
        //            be freed. See _TO_DO_ at the beginning of this file.
        //
        free(node->data);
        next = node->next;
        free(node);
        node = next;
    }
    list->head   = NULL;
    list->tail   = NULL;
    list->length = 0;
}
//-----------------------------------------------------------------------------
void append_to_linked_list(linked_list_t* const list, void* const data) {
    linked_list_node_t* to_append = malloc(sizeof(linked_list_node_t));

    to_append->data = data;
    to_append->next = NULL;

    if (0 == list->length) {
        list->head = to_append;
        list->tail = to_append;
    } else {
        list->tail->next = to_append;
        list->tail       = to_append;
    }
    list->length++;
}
//-----------------------------------------------------------------------------
void prepend_to_linked_list(linked_list_t* const list, void* const data) {
    linked_list_node_t* to_prepend = malloc(sizeof(linked_list_node_t));

    to_prepend->data = data;
    to_prepend->next = list->head;

    if (0 == list->length) {
        list->head = to_prepend;
        list->tail = to_prepend;
    } else {
        list->head = to_prepend;
    }
    list->length++;
}
//-----------------------------------------------------------------------------
void* pop_linked_list(linked_list_t* const list) {

    switch (list->length) {
        //
        // _NOTE_: The brackets are actually necessary to be able to
        //         declare new local variables without conflicts between
        //         the different cases.
        //
        case 0: {
            return NULL;
            break;
        }
        case 1: {
            void* const data = list->tail->data;
            free(list->tail);

            list->head   = NULL;
            list->tail   = NULL;
            list->length = 0;

            return data;

            break;
        }
        default: {
            linked_list_node_t* last = list->tail;
            void* const         data = list->tail->data;

            list->tail = list->head;
            while (list->tail->next != last) {
                list->tail = list->tail->next;
            }
            list->tail->next = NULL;

            free(last);

            list->length--;

            return data;
        }
    }
}
//-----------------------------------------------------------------------------
void* push_linked_list(linked_list_t* const list) {

    switch (list->length) {
        //
        // _NOTE_: The brackets are actually necessary to be able to
        //         declare new local variables without conflicts between
        //         the different cases.
        //
        case 0: {
            return NULL;
            break;
        }
        case 1: {
            void* const data = list->head->data;
            free(list->head);

            list->head   = NULL;
            list->tail   = NULL;
            list->length = 0;

            return data;

            break;
        }
        default: {
            linked_list_node_t* first = list->head;
            void* const         data  = list->head->data;

            list->head = list->head->next;
            free(first);

            list->length--;

            return data;
        }
    }
}
//-----------------------------------------------------------------------------
void insert_in_linked_list(linked_list_t* const list, void* const data) {

    linked_list_node_t* to_insert = malloc(sizeof(linked_list_node_t));

    to_insert->data = data;

    if (0 == list->length) {
        list->head = to_insert;
        list->tail = to_insert;
        to_insert->next = NULL;
        list->length++;
        return;
    }
    //
    // As is usual with linked list, handle head's case separately...
    //
    if (list->cmp_func(data, list->head->data) < 0) {
        to_insert->next = list->head;
        list->head      = to_insert;
        list->length++;
        return;
    }
    linked_list_node_t* current = list->head;
    linked_list_node_t* next    = current->next;
    while (NULL != current->next) {
        if (list->cmp_func(data, current->next->data) > 0) {
            current = current->next;
            next    = current->next;
        } else {
            to_insert->next = current->next;
            current->next   = to_insert;
            list->length++;
            return;
        }
    }
    list->tail->next = to_insert;
    list->tail       = to_insert;
    list->length++;
}
//-----------------------------------------------------------------------------
linked_list_node_t*
get_node_in_linked_list(linked_list_t* const list, void* const data) {

    if (0 == list->length) {
        return NULL;
    }
    linked_list_node_t* node = list->head;
    while (NULL != node) {
        if (list->cmp_func(data, node->data) == 0) {
            return node;
        } else {
            node = node->next;
        }
    }
    return NULL;
}
//-----------------------------------------------------------------------------
linked_list_node_t*
remove_from_linked_list(linked_list_t* const list, void* const data) {

    if (0 == list->length) {
        return NULL;
    }
    //
    // Handles head's case separately...
    //
    linked_list_node_t* node = list->head;
    if (list->cmp_func(node->data, data) == 0) {
        list->head = node->next;
        list->length--;
        if (0 == list->length) {
            list->tail = NULL;
        }
        return node;
    }
    //
    // ... from the other cases...
    //
    linked_list_node_t* next = NULL;
    while (NULL != node->next) {
        if (list->cmp_func(node->next->data, data) == 0) {
            if (node->next == list->tail) {
                list->tail = node;
            }
            next = node->next;
            node->next = node->next->next;
            list->length--;
            return next;
        }
        node = node->next;
    }
    return NULL;
}
//-----------------------------------------------------------------------------
void remove_node_from_linked_list(linked_list_t* const list,
                                  linked_list_node_t* const node) {

    if (0 == list->length) {
        return;
    }
    //
    // Handles head's case separately...
    //
    if (node == list->head) {
        list->head = node->next;
        list->length--;
        if (0 == list->length) {
            list->tail = NULL;
        }
        return;
    }
    //
    // ... from the other cases...
    //
    linked_list_node_t *curnode = list->head;
    while (NULL != curnode->next) {
        if (curnode->next == node) {
            if (node == list->tail) {
                list->tail = curnode;
            }
            curnode->next = node->next;
            list->length--;
            return;
        }
        curnode = curnode->next;
    }
}
//-----------------------------------------------------------------------------
void delete_in_linked_list(linked_list_t* const list, void* const data) {

    if (0 == list->length) {
        return;
    }
    //
    // Handles head's case separately...
    //
    linked_list_node_t *curnode = list->head;
    if (0 == list->cmp_func(data, curnode->data)) {
        list->head = curnode->next;
        //
        // _WARNING_: For the time being, we merely free curnode->data.
        //            However it is possible that the object referenced
        //            by this pointer contains other pointers that needs to
        //            be freed. See _TO_DO_ at the beginning of this file.
        //
        free(curnode->data);
        free(curnode);
        list->length--;
        if (0 == list->length) {
            list->tail = NULL;
        }
        return;
    }
    //
    // ... from the other cases...
    //
    linked_list_node_t *next = NULL;
    while (NULL != curnode->next) {
        if (0 == list->cmp_func(curnode->next->data, data)) {
            if (curnode->next == list->tail) {
                list->tail = curnode;
            }
            next = curnode->next->next;
            //
            // _WARNING_: For the time being, we merely free curnode->data.
            //            However it is possible that the object referenced
            //            by this pointer contains other pointers that needs to
            //            be freed. See _TO_DO_ at the beginning of this file.
            //
            free(curnode->next->data);
            free(curnode->next);
            curnode->next = next;
            list->length--;
            return;
        }
        curnode = curnode->next;
    }
}
//-----------------------------------------------------------------------------
