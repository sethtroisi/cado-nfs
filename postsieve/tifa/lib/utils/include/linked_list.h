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
 * \file    linked_list.h
 * \author  Jerome Milan
 * \date    Thu Mar 30 2006
 * \version 1.0
 *
 * \brief Standard singly-linked list.
 *
 * Defines generic singly-linked lists and their associated functions.
 */

 /*
  *  License: GNU Lesser General Public License (LGPL)
  *  History:
  */

#if !defined(_TIFA_LINKED_LIST_H_)
   /**
    * \def _TIFA_LINKED_LIST_H_
    * Standard include guard.
    */
#define _TIFA_LINKED_LIST_H_

#include <inttypes.h>

#ifdef __cplusplus
extern "C" {
#endif

   /**
    * \struct struct_linked_list_node_t linked_list.h
    *                                   lib/utils/include/linked_list.h
    * \brief  A basic implementation of a linked list node.
    *
    * This structure defines a simple linked list node, i.e. some data together
    * with a pointer to the next node.
    */
struct struct_linked_list_node_t {
       /**
        * Pointer to the next node of the linked list.
        */
    struct struct_linked_list_node_t* next;
       /**
        * Pointer to the node data.
        */
    void*  data;
};

   /**
    * \struct struct_linked_list_t linked_list.h lib/utils/include/linked_list.h
    * \brief  A basic implementation of a linked list.
    *
    * This structure defines a simple linked list. It can store any
    * type of elements, provided that there is a suitable comparison function
    * for the type of the data used.
    */
struct struct_linked_list_t {
       /**
        * Pointer to the head of the linked list.
        */
    struct struct_linked_list_node_t* head;
       /**
        * Pointer to the tail of the linked list.
        */
    struct struct_linked_list_node_t* tail;
       /**
        * Pointer to the comparison function used to compare the node data.
        * This function, which take two pointers to the data to compare,
        * should return:
        * \li 0 if the datas pointed by \c data_a and \c data_b are the same.
        * \li A positive number if the data pointed by \c data_a is greater
        *     than the data pointed by <tt>data_b</tt>.
        * \li A negative number if the data pointed by \c data_a is less
        *     than the data pointed by <tt>data_b</tt>.
        */
    int    (*cmp_func) (const void* const data_a, const void* const data_b);
       /**
        * Number of nodes composing the linked list.
        */
    uint32_t length;
};

   /**
    * \typedef linked_list_t
    * \brief Equivalent to <tt>struct struct_linked_list_t</tt>.
    */
typedef struct struct_linked_list_t linked_list_t;

   /**
    * \typedef linked_list_node_t
    * \brief Equivalent to <tt>struct struct_linked_list_node_t</tt>.
    */
typedef struct struct_linked_list_node_t linked_list_node_t;

   /**
    * \brief Initializes a <tt>linked_list_t</tt>.
    *
    * Initializes a linked list <tt>list</tt>:
    * \li Sets its \c head to <tt>NULL</tt>.
    * \li Sets its \c tail to <tt>NULL</tt>.
    * \li Sets its \c cmp_func to the <tt>cmp_func</tt> argument.
    * \li Sets its \c length to 0.
    *
    * \param[in] list A pointer to the \c linked_list_t to initialize.
    * \param[in] cmp_func A pointer to the comparison function.
    */
void init_linked_list(linked_list_t* const list,
                      int (*cmp_func) (const void* const data_a,
                                       const void* const data_b) );

   /**
    * \brief Clears a <tt>linked_list_t</tt>.
    *
    * Clears a linked list <tt>list</tt> by freeing its nodes.
    *
    * \warning Each node's \c data field is freed. However, if \c data is a
    * pointer to a structure containing some other pointers, all the memory
    * may not be freed.
    *
    * \param[in] list A pointer to the \c linked_list_t to clear.
    */
void clear_linked_list(linked_list_t* const list);

   /**
    * \brief Appends a node in a <tt>linked_list_t</tt>.
    *
    * Appends a node (whose data is given by the pointer \c data) in the
    * linked list <tt>list</tt>.
    *
    * \param[in] list A pointer to the <tt>linked_list_t</tt>.
    * \param[in] data A pointer to the new node's data.
    */
void append_to_linked_list(linked_list_t* const list, void* const data);

   /**
    * \brief Prepends a node in a <tt>linked_list_t</tt>.
    *
    * Prepends a node (whose data is given by the pointer \c data) in the
    * linked list <tt>list</tt>.
    *
    * \param[in] list A pointer to the <tt>linked_list_t</tt>.
    * \param[in] data A pointer to the new node's data.
    */
void prepend_to_linked_list(linked_list_t* const list, void* const data);


   /**
    * \brief Deletes the last node of a <tt>linked_list_t</tt>.
    *
    * Returns the data of the last node of the linked list <tt>list</tt>
    * and deletes this node, similar to Perl's pop function.
    *
    * \param[in] list A pointer to the <tt>linked_list_t</tt>.
    */
void* pop_linked_list(linked_list_t* const list);

   /**
    * \brief Deletes the first node of a <tt>linked_list_t</tt>.
    *
    * Returns the data of the first node of the linked list <tt>list</tt>
    * and deletes this node, similar to Perl's push function.
    *
    * \param[in] list A pointer to the <tt>linked_list_t</tt>.
    */
void* push_linked_list(linked_list_t* const list);

   /**
    * \brief Inserts a node in a <tt>linked_list_t</tt>.
    *
    * Inserts a node (whose data is given by the pointer \c data) in the
    * linked list <tt>list</tt>, so that all the previous nodes have \c data
    * fields pointing to datas less than the new node's data.
    *
    * \param[in] list A pointer to the <tt>linked_list_t</tt>.
    * \param[in] data A pointer to the new node's data.
    */
void insert_in_linked_list(linked_list_t* const list, void* const data);

   /**
    * \brief Gets a node in a <tt>linked_list_t</tt>.
    *
    * Returns a pointer to a node (whose data is given by the pointer \c data)
    * from the linked list <tt>list</tt>.
    *
    * \param[in] list A pointer to the <tt>linked_list_t</tt>.
    * \param[in] data A pointer to the seeked node data.
    * \return A pointer to the node with data pointed by <tt>data</tt>, if such
    *         a node exists.
    * \return NULL if no matching node is found in the linked list.
    */
linked_list_node_t* get_node_in_linked_list(linked_list_t* const list,
                                            void* const data);

   /**
    * \brief Gets a node in a <tt>linked_list_t</tt> and removes it from
    * the list.
    *
    * Returns a pointer to a node (whose data is given by the pointer \c data)
    * from the linked list <tt>list</tt> and removes this node from the list.
    *
    * \param[in] list A pointer to the <tt>linked_list_t</tt>.
    * \param[in] data A pointer to the seeked node data.
    * \return A pointer to the node with data pointed by <tt>data</tt>, if such
    *         a node exists.
    * \return NULL if no matching node is found in the linked list.
    */
linked_list_node_t* remove_from_linked_list(linked_list_t* const list,
                                            void* const data);

   /**
    * \brief Removes a given node from a <tt>linked_list_t</tt>.
    *
    * Removes the node whose data is given by the pointer \c data
    * from the linked list <tt>list</tt>. If the node is not found, the linked
    * list is left unchanged.
    *
    * \param[in] list A pointer to the <tt>linked_list_t</tt>.
    * \param[in] node A pointer to the node to remove from the list.
    */
void remove_node_from_linked_list(linked_list_t* const list,
                                  linked_list_node_t* const node);

   /**
    * \brief Finds and deletes a node from a <tt>linked_list_t</tt>.
    *
    * Deletes the node whose data is given by the pointer \c data
    * from the linked list <tt>list</tt>. If the node is not found, the linked
    * list is left unchanged.
    *
    * \warning If a matching node is found, its \c data field is freed.
    * However, if \c data is a pointer to a structure containing some
    * other pointers, all the memory may not be freed.
    *
    * \param[in] list A pointer to the <tt>linked_list_t</tt>.
    * \param[in] data A pointer to the data of the node to delete.
    */
void delete_in_linked_list(linked_list_t* const list, void* const data);

#ifdef __cplusplus
}
#endif

#endif
