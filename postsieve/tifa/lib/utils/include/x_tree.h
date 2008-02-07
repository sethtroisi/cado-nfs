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
 * \file    x_tree.h
 * \author  Jerome Milan
 * \date    Wed Nov 28 2007
 * \version 1.0.1
 *
 * \brief Product and remainder trees.
 *
 * Implementation of the product and remainder trees used in D. J. Bernstein's
 * algorithms.
 */

 /*
  *  History:
  *
  *  1.0.1: Wed Nov 28 2007 by JM:
  *         - Added the (currently unused) 'prod_tree_mod' function protoype.
  *
  *  1.0.0: Mon Mar 6 2006 by JM:
  *         - Initial version.
  */

#if !defined(_TIFA_X_TREE_H_)
   /**
    * \def _TIFA_X_TREE_H_
    * Standard include guard.
    */
#define _TIFA_X_TREE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <gmp.h>
#include "array.h"

   /**
    * \typedef mpz_tree_t
    * \brief Equivalent to <tt>mpz_array_t</tt>.
    *
    * While an \c mpz_tree_t is just an <tt>mpz_array_t</tt>, its memory is
    * allocated in a different manner than in the \c mpz_array_t case. Indeed,
    * the elements of an \c mpz_tree_t array should NOT be modified later on
    * since the memory used is allocated in one huge block to prevent overhead
    * from multiple malloc calls. So the allocated memory of the mpz_t's in the
    * tree can NOT be increased.
    *
    * The \c mpz_tree_t typedef is introduced only as a reminder of this
    * different underlying memory allocation. \c clear_mpz_tree should be
    * used to clear the memory space occupied by an <tt>mpz_tree_t</tt>. Do
    * NOT call \c clear_mpz_array on an <tt>mpz_tree_t</tt>!
    */
typedef mpz_array_t mpz_tree_t;

   /**
    * \brief Computes the product tree of some \c mpz_t integers.
    *
    * Computes the product tree of the \c mpz_t integers given by \c array
    * and returns it as a newly allocated <tt>mpz_tree_t</tt>.
    *
    * \note The product tree is implemented as a single \c mpz_array_t \c tree
    * with the usual compact representation: \c tree->data[2i+1] and
    * \c tree->data[2i+2] are the children of the node \c tree->data[i].
    *
    * Hence, in order to avoid useless nodes (i.e nodes with value 1),
    * it is recommended to have \c array->length equals to a power of 2.
    * If this is not the case, the product tree will be computed as if
    * array was completed by as many useless nodes as necessary until
    * a power of 2 is reached.
    *
    * This choice was made to keep a space efficient representation
    * and to avoid dynamic allocations of nodes.
    *
    * \warning
    * Although the product tree returned is actually a pointer to an
    * \c mpz_tree_t structure (i.e. an \c mpz_array_t), the elements of the
    * array should NOT be modified later on since the memory used is allocated
    * in one huge block to prevent overhead from multiple malloc calls. So the
    * allocated memory of the mpz_t's in the array can NOT be increased...
    *
    * \param[in] array Pointer to the \c mpz_array_t containing the
    *                  \c mpz_t integers to multiply.
    * \return A pointer to a newly allocated \c mpz_tree_t holding the
    * computed product tree.
    */
mpz_tree_t* prod_tree(const mpz_array_t* const array);

   /**
    * \brief Computes the product tree of some \c mpz_t integers modulo
    * a positive integer.
    *
    * Similar to \c prod_tree but each node is reduced mod \c n.
    *
    * \warning n should be strictly positive or results will be unpredictable.
    *
    * \see The function <tt>prod_tree(const mpz_array_t* const array)</tt>.
    *
    * \param[in] array Pointer to the \c mpz_array_t containing the
    *                  \c mpz_t integers to multiply.
    * \return A pointer to a newly allocated \c mpz_tree_t holding the
    * computed product tree.
    */
mpz_tree_t* prod_tree_mod(const mpz_array_t* const array, const mpz_t n);

   /**
    * \brief Computes the product tree of some \c uint32_t integers.
    *
    * Computes the product tree of the \c uint32_t integers given by \c array
    * and returns it as a newly allocated <tt>mpz_tree_t</tt>.
    *
    * \note The product tree is implemented as a single \c mpz_array_t \c tree
    * with the usual compact representation: \c tree->data[2i+1] and
    * \c tree->data[2i+2] are the children of the node \c tree->data[i].
    *
    * Hence, in order to avoid useless nodes (i.e nodes with value 1),
    * it is recommended to have \c array->length equals to a power of 2.
    * If this is not the case, the product tree will be computed as if
    * array was completed by as many useless nodes as necessary until
    * a power of 2 is reached.
    *
    * This choice was made to keep a space efficient representation
    * and to avoid dynamic allocations of nodes.
    *
    * \warning
    * Although the product tree returned is actually a pointer to an
    * \c mpz_tree_t structure (i.e. an \c mpz_array_t), the elements of the
    * array should NOT be modified later on since the memory used is allocated
    * in one huge block to prevent overhead from multiple malloc calls. So the
    * allocated memory of the mpz_t's in the array can NOT be increased...
    *
    * \param[in] array Pointer to the \c uint32_array_t containing the
    *                  \c mpz_t integers to multiply.
    * \return A pointer to a newly allocated \c mpz_tree_t holding the
    * computed product tree.
    */
mpz_tree_t* prod_tree_ui(const uint32_array_t* const array);

   /**
    * \brief Computes a remainder tree.
    *
    * Computes the remainder tree of \c z by the \c mpz_t integers whose
    * product tree is given by \c ptree and returns it as a newly allocated
    * <tt>mpz_tree_t</tt>. If \c rtree is the returned remainder tree, then
    * one has: \c rtree->data[i] = \c z mod <tt>ptree->data[i]</tt>.
    *
    * \note The remainder tree is implemented as a single \c mpz_array_t \c tree
    * with the usual compact representation: \c tree->data[2i+1] and
    * \c tree->data[2i+2] are the children of the node \c tree->data[i].
    *
    * \warning
    * Although the remainder tree returned is actually a pointer to an
    * \c mpz_tree_t structure (i.e. an \c mpz_array_t), the elements of the
    * array should NOT be modified later on since the memory used is allocated
    * in one huge block to prevent overhead from multiple malloc calls. So the
    * allocated memory of the mpz_t's in the array can NOT be increased...
    *
    * \param[in] z The integer to divide.
    * \param[in] ptree Pointer to the \c mpz_tree_t containing the product
    *                  tree of the \c mpz_t integers to divide \c z by.
    * \return A pointer to a newly allocated \c mpz_tree_t holding the
    * computed remainder tree.
    */
mpz_tree_t* rem_tree(const mpz_t z, const mpz_tree_t* const ptree);

   /**
    * \brief Clears a tree of \c mpz_t integers.
    *
    * Clears a tree of \c mpz_t integers returned by <tt>prod_tree</tt>,
    * <tt>prod_tree_ui</tt> or <tt>rem_tree</tt>.
    *
    * \note This function is actually different from <tt>clear_mpz_array</tt>.
    * Indeed, even if the \c mpz_tree_t type is merely a typedef of
    * <tt>mpz_array_t</tt>, the memory used by the \c mpz_t elements is
    * allocated in a completely different manner, hence the need for a
    * different function.
    *
    * \param[in] tree Pointer to the \c mpz_tree_t to clear.
    */
void clear_mpz_tree(mpz_tree_t* tree);

   /**
    * \brief Prints a tree of \c mpz_t integers.
    *
    * Prints a tree of \c mpz_t integers on the standard output. Useful
    * only for debugging purposes and for relatively small trees.
    *
    * \param[in] tree Pointer to the \c mpz_array_t containing the tree
    *                 to print.
    */
void print_mpz_tree(const mpz_tree_t* const tree);

#ifdef __cplusplus
}
#endif

#endif
