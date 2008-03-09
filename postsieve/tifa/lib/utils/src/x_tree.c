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
 * \file    x_tree.c
 * \author  Jerome Milan
 * \date    Wed Nov 28 2007
 * \version 1.1.1
 */

 /*
  *  History:
  *
  *  1.1.1: Wed Nov 28 2007 by JM:
  *         - Added the (currently unused) 'prod_tree_mod' function.
  *
  *  1.1.0: Tue Dec 12 2006 by JM:
  *         - new typedef: mpz_tree_t.
  *         - prod_tree, prod_tree_ui and rem_tree modified to avoid multiple
  *           small memory allocations. Product tree functions are a few
  *           percents faster as a result. However, rem_tree is now more than
  *           twice as fast as a straitforward mpz_t-based implementation.
  *         - removed unused print_uint32_tree function.
  *
  *  1.0.2: Tue Nov 21 2006 by JM:
  *         - rem_tree rewritten to perform the division with the mpn
  *           layer. Speed-up of the order of 30%.
  *
  *  1.0.1: Wed Nov 15 2006 by JM:
  *         - prod_tree rewritten to perform the multiplication with the mpn
  *           layer. Amazing speed-up of the order of... a couple of percents!
  *
  *  1.0.0: Mon Mar 6 2006 by JM:
  *         - Initial version.
  */

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>

#include "funcs.h"
#include "macros.h"
#include "gmp_utils.h"
#include "x_tree.h"

//-----------------------------------------------------------------------------
void print_mpz_tree_rec(const mpz_tree_t* const, uint32_t, uint32_t);
//-----------------------------------------------------------------------------
mpz_tree_t* prod_tree(const mpz_array_t* const array) {
    //
    // _NOTE_: The product tree is implemented as a single mpz_array_t tree
    //         with the usual compact representation: tree->data[2i+1]
    //         and tree->data[2i+2] are the children of the node tree->data[i].
    //
    //         Hence, in order to avoid useless nodes (i.e nodes with value 1),
    //         it is recommended to have array->length equals to a power of 2.
    //         If this is not the case, the product tree will be computed as if
    //         array was completed by as many useless nodes as necessary until
    //         a power of 2 is reached.
    //
    //         This choice was made to keep a space efficient representation
    //         and to avoid dynamic allocation of nodes.
    //
    // _WARNING_:
    //
    //         Although the product tree returned is actually a pointer
    //         to an mpz_array_t structure, the array should NOT be modified
    //         later on since the memory it uses is allocated in one huge
    //         block to prevent overhead from multiple malloc calls. So the
    //         allocated memory of the mpz_t's in the array can NOT be
    //         increased...
    //
    // _NOTE_: The differences between this function and the prod_tree_ui
    //         function are quite minimal. In other words, there is a lot of
    //         code reuse here which, from a purely software engineering point
    //         of view, is pretty ugly. It may be worthwhile to rewrite these
    //         two functions in a clever way to avoid such bad code reuse...
    //
    // _UPDATE_: Wed Nov 15 2006 by JM
    //
    //         This function has been modified to perform the
    //         multiplication directly with the mpn layer. This can
    //         speed-up the function up to... a couple of percents! Wow!
    //
    // _UPDATE_: Mon Dec 11 2006 by JM
    //
    //         This function has been modified to allocate all the needed
    //         memory with a single call to malloc (See the warning). Together
    //         with the rewrite using the mpn layer, we gain about 10 percents.
    //         Not much, but still better than nothing...
    //
    mpz_tree_t* const tree = malloc(sizeof(mpz_tree_t));
    //
    // Handle the special cases separately...
    //
    if (0 == array->length) {
        tree->alloced = 0;
        tree->length = 0;
        tree->data = NULL;
        return tree;
    }
    if (1 == array->length) {
        tree->alloced = 1;
        tree->length = 1;
        tree->data = malloc(sizeof(mpz_t));
        mpz_init_set(tree->data[0], array->data[0]);
        return tree;
    }
    //
    // General case...
    //
    // minpow is the least integer such that 2^minpow >= array->length
    //
    uint32_t minpow = most_significant_bit(array->length);
    if (((1U<<minpow) ^ array->length) != 0) {
        minpow++;
    }
    //
    // The number of leaves is (1<<minpow), which is different from
    // array->length if array->length is not a power of 2.
    //
    tree->alloced = 2*(1<<minpow) - 1;
    tree->length = tree->alloced;
    tree->data = malloc(tree->alloced*sizeof(mpz_t));
    //
    // Use a local constant pointer to the data to bypass the reading of a
    // non constant pointer. This should hopefully speed-up things a bit...
    // A (very) little bit...
    //
    mpz_t* const treedata = tree->data;
    //
    // i >= offset => tree[i] is a leaf
    // i <  offset => tree[i] is a node
    //
    uint32_t const offset = tree->length/2;
    //
    // Computes the memory needed to hold the product tree and allocate the
    // whole memory needed in only one call to malloc
    //
    uint32_t size_tree   = 0;
    uint32_t size_leaves = 0;
    for (uint32_t inode = 0U; inode < array->length; inode++) {
        size_leaves += ABSIZ(array->data[inode]);
    }
    for (uint32_t inode = tree->length - 1;
         inode >= offset + array->length;
         inode--) {
        size_leaves++;
    }
    size_tree = size_leaves * (minpow + 1);

    PTR(treedata[0]) = malloc(size_tree * sizeof(mp_limb_t));
    mp_limb_t* const treelimbs = PTR(treedata[0]);

    uint32_t next_free_space_index = size_tree;
    //
    // If necessary (i.e. array->length is not a power of 2), sets the
    // remaining tree leaves to 1. Note that in that case the allocated memory
    // will be strictly greater than what's really needed, but any
    // "optimization" here is likely to be more of a waste that just living
    // with it...
    //
    for (uint32_t inode = tree->length - 1;
         inode >= offset + array->length;
         inode--) {

        next_free_space_index--;

        treelimbs[next_free_space_index] = (mp_limb_t)1;

        PTR(treedata[inode])   = &(treelimbs[next_free_space_index]);
        SIZ(treedata[inode])   = 1;
        ALLOC(treedata[inode]) = 1;
    }
    uint32_t iarray = array->length - 1;
    //
    // Sets the value of the leaves...
    //
    for (uint32_t inode = offset + array->length - 1;
         inode >= offset; inode--) {

        next_free_space_index -= ABSIZ(array->data[iarray]);

        memcpy(
            &(treelimbs[next_free_space_index]),
            PTR(array->data[iarray]),
            ABSIZ(array->data[iarray]) * sizeof(mp_limb_t)
        );
        PTR(treedata[inode])   = &(treelimbs[next_free_space_index]);
        SIZ(treedata[inode])   = SIZ(array->data[iarray]);
        ALLOC(treedata[inode]) = ABSIZ(array->data[iarray]);

        iarray--;
    }
    uint32_t ichild_1 = 0;
    uint32_t ichild_2 = 0;
    //
    // Computes the products for each node...
    //
    mp_size_t size_1   = 0;
    mp_size_t size_2   = 0;
    mp_size_t size_res = 0;

    for (uint32_t inode = offset - 1; inode != 0U; inode--) {

        ichild_1 = (inode<<1) | 1;
        ichild_2 = (inode<<1) + 2;
        //
        // Set the mp_limb_t pointers to their "dedicated" memory spaces...
        //
        size_1   = ABSIZ(treedata[ichild_1]);
        size_2   = ABSIZ(treedata[ichild_2]);
        size_res = size_1 + size_2;

        next_free_space_index -= size_res;

        PTR(treedata[inode]) = &(treelimbs[next_free_space_index]);

        ALLOC(treedata[inode]) = size_res;

        if (size_1 > size_2) {
            mpn_mul(PTR(treedata[inode]),
                    PTR(treedata[ichild_1]), size_1,
                    PTR(treedata[ichild_2]), size_2);
        } else {
            mpn_mul(PTR(treedata[inode]),
                    PTR(treedata[ichild_2]), size_2,
                    PTR(treedata[ichild_1]), size_1);
        }
        //
        // Normalize the size of the result...
        //
        MPN_NORMALIZE(PTR(treedata[inode]), size_res);

        if ((SIZ(treedata[ichild_1]) ^ SIZ(treedata[ichild_2])) >= 0) {
            SIZ(treedata[inode]) = size_res;
        } else {
            SIZ(treedata[inode]) = -size_res;
        }
    }
    //
    // Root node is computed separately...
    //
    size_1   = ABSIZ(treedata[1]);
    size_2   = ABSIZ(treedata[2]);
    size_res = size_1 + size_2;

    treedata[0][0]._mp_alloc = size_res;

    if (size_1 > size_2) {
        mpn_mul(PTR(treedata[0]),
                PTR(treedata[1]), size_1,
                PTR(treedata[2]), size_2);
    } else {
        mpn_mul(PTR(treedata[0]),
                PTR(treedata[2]), size_2,
                PTR(treedata[1]), size_1);
    }
    MPN_NORMALIZE(PTR(treedata[0]), size_res);

    if ((SIZ(treedata[1]) ^ SIZ(treedata[2])) >= 0) {
        SIZ(treedata[0]) = size_res;
    } else {
        SIZ(treedata[0]) = -size_res;
    }
    return tree;
}
//-----------------------------------------------------------------------------
mpz_tree_t* prod_tree_mod(const mpz_array_t* const array, const mpz_t n) {
    //
    // _NOTE_: See introductory comments in the body of the prod_tree function.
    //
    // _WARNING_: n should be strictly positive.
    //
    // _NOTE_: The differences between this function and the prod_tree
    //         function are quite minimal. In other words, there is a lot of
    //         code reuse here which, from a purely software engineering point
    //         of view, is pretty ugly. It may be worthwhile to rewrite these
    //         two functions in a clever way to avoid such bad code reuse...
    //
    mpz_tree_t* const tree = malloc(sizeof(mpz_tree_t));
    //
    // Handle the special cases separately...
    //
    if (0 == array->length) {
        tree->alloced = 0;
        tree->length = 0;
        tree->data = NULL;
        return tree;
    }
    if (1 == array->length) {
        tree->alloced = 1;
        tree->length = 1;
        tree->data = malloc(sizeof(mpz_t));
        mpz_init_set(tree->data[0], array->data[0]);
        return tree;
    }
    //
    // General case...
    //
    // minpow is the least integer such that 2^minpow >= array->length
    //
    uint32_t minpow = most_significant_bit(array->length);
    if (((1U<<minpow) ^ array->length) != 0) {
        minpow++;
    }
    //
    // The number of leaves is (1<<minpow), which is different from
    // array->length if array->length is not a power of 2.
    //
    tree->alloced = 2*(1<<minpow) - 1;
    tree->length = tree->alloced;
    tree->data = malloc(tree->alloced*sizeof(mpz_t));
    //
    // Use a local constant pointer to the data to bypass the reading of a
    // non constant pointer. This should hopefully speed-up things a bit...
    // A (very) little bit...
    //
    mpz_t* const treedata = tree->data;
    //
    // i >= offset => tree[i] is a leaf
    // i <  offset => tree[i] is a node
    //
    uint32_t const offset = tree->length/2;
    //
    // Computes the memory needed to hold the product tree and allocate the
    // whole memory needed in only one call to malloc
    //
    uint32_t size_tree   = 0;
    uint32_t size_leaves = 0;
    for (uint32_t inode = 0U; inode < array->length; inode++) {
        size_leaves += ABSIZ(array->data[inode]);
    }
    for (uint32_t inode = tree->length - 1;
         inode >= offset + array->length;
         inode--) {
        size_leaves++;
    }
    size_tree = size_leaves * (minpow + 1);

    PTR(treedata[0]) = malloc(size_tree * sizeof(mp_limb_t));
    mp_limb_t* const treelimbs = PTR(treedata[0]);

    uint32_t next_free_space_index = size_tree;
    //
    // If necessary (i.e. array->length is not a power of 2), sets the
    // remaining tree leaves to 1. Note that in that case the allocated memory
    // will be strictly greater than what's really needed, but any
    // "optimization" here is likely to be more of a waste that just living
    // with it...
    //
    for (uint32_t inode = tree->length - 1;
         inode >= offset + array->length;
         inode--) {

        next_free_space_index--;

        treelimbs[next_free_space_index] = (mp_limb_t)1;

        PTR(treedata[inode])   = &(treelimbs[next_free_space_index]);
        SIZ(treedata[inode])   = 1;
        ALLOC(treedata[inode]) = 1;
    }
    uint32_t iarray = array->length - 1;
    //
    // Sets the value of the leaves...
    //
    for (uint32_t inode = offset + array->length - 1;
         inode >= offset; inode--) {

        next_free_space_index -= ABSIZ(array->data[iarray]);

        memcpy(
            &(treelimbs[next_free_space_index]),
            PTR(array->data[iarray]),
            ABSIZ(array->data[iarray]) * sizeof(mp_limb_t)
        );
        PTR(treedata[inode])   = &(treelimbs[next_free_space_index]);
        SIZ(treedata[inode])   = SIZ(array->data[iarray]);
        ALLOC(treedata[inode]) = ABSIZ(array->data[iarray]);

        iarray--;
    }
    
    mp_limb_t* quotient = malloc(
                              (size_leaves - SIZ(n) + 1) * sizeof(mp_limb_t)
                          );
    mp_limb_t* remainder   = malloc(SIZ(n) * sizeof(mp_limb_t));

    uint32_t ichild_1 = 0;
    uint32_t ichild_2 = 0;
    //
    // Computes the products for each node...
    //
    mp_size_t size_1   = 0;
    mp_size_t size_2   = 0;
    mp_size_t size_res = 0;

    for (uint32_t inode = offset - 1; inode != 0U; inode--) {

        ichild_1 = (inode<<1) | 1;
        ichild_2 = (inode<<1) + 2;
        //
        // Set the mp_limb_t pointers to their "dedicated" memory spaces...
        //
        size_1   = ABSIZ(treedata[ichild_1]);
        size_2   = ABSIZ(treedata[ichild_2]);
        size_res = size_1 + size_2;

        next_free_space_index -= size_res;

        PTR(treedata[inode]) = &(treelimbs[next_free_space_index]);

        treedata[inode][0]._mp_alloc = size_res;

        if (size_1 > size_2) {
            mpn_mul(PTR(treedata[inode]),
                    PTR(treedata[ichild_1]), size_1,
                    PTR(treedata[ichild_2]), size_2);
        } else {
            mpn_mul(PTR(treedata[inode]),
                    PTR(treedata[ichild_2]), size_2,
                    PTR(treedata[ichild_1]), size_1);
        }
        //
        // Normalize the size of the result...
        //
        MPN_NORMALIZE(PTR(treedata[inode]), size_res);

        if ((SIZ(treedata[ichild_1]) ^ SIZ(treedata[ichild_2])) >= 0) {
            SIZ(treedata[inode]) = size_res;
        } else {
            SIZ(treedata[inode]) = -size_res;
        }   
        //
        // Reduce modulo n
        //
        if (ABSIZ(treedata[inode]) >= SIZ(n)) {
            mpn_tdiv_qr(
                quotient, remainder,
                0,
                PTR(treedata[inode]), ABSIZ(treedata[inode]),
                PTR(n), SIZ(n)
            );        
            for (int i = 0; i < SIZ(n); i++) {
                treedata[inode][0]._mp_d[i] = remainder[i];
            }
            SIZ(treedata[inode]) = SIZ(n);
            
            MPN_NORMALIZE(PTR(treedata[inode]), SIZ(treedata[inode]));
            
            if ((SIZ(treedata[ichild_1]) ^ SIZ(treedata[ichild_2])) < 0) {
                SIZ(treedata[inode]) = -SIZ(treedata[inode]);
            }
        }
    }
    //
    // Root node is computed separately...
    //
    size_1   = ABSIZ(treedata[1]);
    size_2   = ABSIZ(treedata[2]);
    size_res = size_1 + size_2;

    ALLOC(treedata[0]) = size_res;

    if (size_1 > size_2) {
        mpn_mul(PTR(treedata[0]),
                PTR(treedata[1]), size_1,
                PTR(treedata[2]), size_2);
    } else {
        mpn_mul(PTR(treedata[0]),
                PTR(treedata[2]), size_2,
                PTR(treedata[1]), size_1);
    }
    MPN_NORMALIZE(PTR(treedata[0]), size_res);

    if ((SIZ(treedata[1]) ^ SIZ(treedata[2])) >= 0) {
        SIZ(treedata[0]) = size_res;
    } else {
        SIZ(treedata[0]) = -size_res;
    }
    //
    // Reduce modulo n
    //
    if (ABSIZ(treedata[0]) >= SIZ(n)) {
        mpn_tdiv_qr(
            quotient, remainder,
            0,
            PTR(treedata[0]), ABSIZ(treedata[0]),
            PTR(n), SIZ(n)
        );   
        for (int i = 0; i < SIZ(n); i++) {
            treedata[0][0]._mp_d[i] = remainder[i];
        }
        SIZ(treedata[0]) = SIZ(n);
        
        MPN_NORMALIZE(PTR(treedata[0]), SIZ(treedata[0]));
        
        if ((SIZ(treedata[1]) ^ SIZ(treedata[2])) < 0) {
            SIZ(treedata[0]) = -SIZ(treedata[0]);
        }
    }
    free(quotient);
    free(remainder);
    
    return tree;
}
//-----------------------------------------------------------------------------
mpz_tree_t* prod_tree_ui(const uint32_array_t* const array) {
    //
    // _NOTE_: The product tree is implemented as a single mpz_array_t tree
    //         with the usual compact representation: tree->data[2i+1]
    //         and tree->data[2i+2] are the children of the node tree->data[i].
    //
    //         Hence, in order to avoid useless nodes (i.e nodes with value 1),
    //         it is recommended to have array->length equals to a power of 2.
    //         If this is not the case, the product tree will be computed as if
    //         array was completed by as many useless nodes as necessary until
    //         a power of 2 is reached.
    //
    //         This choice was made to keep a space efficient representation
    //         and to avoid dynamic allocation of nodes.
    //
    // _WARNING_:
    //
    //         Although the product tree returned is actually a pointer
    //         to an mpz_array_t structure, the array should NOT be modified
    //         later on since the memory it uses is allocated in one huge
    //         block to prevent overhead from multiple malloc calls. So the
    //         allocated memory of the mpz_t's in the array can NOT be
    //         increased...
    //
    // _NOTE_: The differences between this function and the prod_tree function
    //         are quite minimal. In other words, there is a lot of code reuse
    //         here which, from a purely software engineering point of view, is
    //         pretty ugly. It may be worthwhile to rewrite these two functions
    //         in a clever way to avoid such bad code reuse...
    //
    // _UPDATE_: Wed Nov 15 2006 by JM
    //
    //         This function has been rewritten to perform the
    //         multiplication directly with the mpn layer. This can
    //         speed-up the function up to... a couple of percents! Wow!
    //
    // _UPDATE_: Mon Dec 11 2006 by JM
    //
    //         This function has been modified to allocate all the needed
    //         memory with a single call to malloc (See the warning). Together
    //         with the rewrite using the mpn layer, we gain about 10 percents.
    //         Not much, but still better than nothing...
    //
    mpz_tree_t* const tree = malloc(sizeof(mpz_tree_t));
    //
    // Handle the special cases separately...
    //
    if (0 == array->length) {
        tree->alloced = 0;
        tree->length = 0;
        tree->data = NULL;
        return tree;
    }
    if (1 == array->length) {
        tree->alloced = 1;
        tree->length = 1;
        tree->data = malloc(sizeof(mpz_t));
        mpz_init_set_ui(tree->data[0], array->data[0]);
        return tree;
    }
    //
    // General case...
    //
    // minpow is the least integer such that 2^minpow >= array->length
    //
    uint32_t minpow = most_significant_bit(array->length);
    if (((1U<<minpow) ^ array->length) != 0) {
        minpow++;
    }
    //
    // The number of leaves is (1<<minpow), which is different from
    // array->length if array->length is not a power of 2.
    //
    tree->alloced = 2*(1<<minpow) - 1;
    tree->length = tree->alloced;
    tree->data = malloc(tree->alloced*sizeof(mpz_t));
    //
    // Use a local constant pointer to the data to bypass the reading of a
    // non constant pointer. This should hopefully speed-up things a bit...
    // A (very) little bit...
    //
    mpz_t* const treedata = tree->data;
    //
    // i >= offset => tree[i] is a leaf
    // i <  offset => tree[i] is a node
    //
    uint32_t const offset = tree->length/2;
    //
    // Computes the memory needed to hold the product tree and allocate the
    // whole memory needed in only one call to malloc
    //
    uint32_t size_tree   = 0;
    uint32_t size_leaves = 0;
    for (uint32_t inode = 0U; inode < array->length; inode++) {
        size_leaves ++;
    }
    for (uint32_t inode = tree->length - 1;
         inode >= offset + array->length;
         inode--) {
        size_leaves++;
    }
    size_tree = size_leaves * (minpow + 1);

    PTR(treedata[0]) = malloc(size_tree * sizeof(mp_limb_t));
    mp_limb_t* const treelimbs = PTR(treedata[0]);

    uint32_t next_free_space_index = size_tree;
    //
    // If necessary (i.e. array->length is not a power of 2), sets the
    // remaining tree leaves to 1.
    //
    for (uint32_t inode = tree->length - 1;
         inode >= offset + array->length;
         inode--) {

        next_free_space_index--;

        treelimbs[next_free_space_index] = (mp_limb_t)1;

        PTR(treedata[inode])   = &(treelimbs[next_free_space_index]);
        SIZ(treedata[inode])   = 1;
        ALLOC(treedata[inode]) = 1;
    }
    uint32_t iarray = array->length - 1;
    //
    // Sets the value of the leaves... Note that in that case the allocated
    // memory will be strictly greater than what's really needed, but any
    // "optimizations" here are likely to be more of a waste that just living
    // with it...
    //
    for (uint32_t inode = offset + array->length - 1;
         inode >= offset; inode--) {
        //
        // The following loop is (together with the code in the
        // 1==array->length special case) is really the only difference
        // with the prod_tree function. This means a lot of copy'n'paste so it
        // would actually be worthwhile to rewrite these two functions in a
        // clever way to avoid such bad code reuse... Maybe the preprocessor
        // could be used, even if this is ugly...
        //
        next_free_space_index--;

        treelimbs[next_free_space_index] = array->data[iarray];

        PTR(treedata[inode]) = &(treelimbs[next_free_space_index]);
        if (array->data[iarray] > 0) {
            SIZ(treedata[inode]) = 1;
        } else {
            SIZ(treedata[inode]) = -1;
        }
        ALLOC(treedata[inode]) = 1;

        iarray--;
    }
    uint32_t ichild_1 = 0;
    uint32_t ichild_2 = 0;

    next_free_space_index = size_tree - 1 - size_leaves;
    //
    // Computes the products for each node...
    //
    mp_size_t size_1   = 0;
    mp_size_t size_2   = 0;
    mp_size_t size_res = 0;

    for (uint32_t inode = offset - 1; inode != 0U; inode--) {

        ichild_1 = (inode<<1) | 1;
        ichild_2 = (inode<<1) + 2;
        //
        // Set the mp_limb_t pointers to their "dedicated" memory spaces...
        //
        size_1   = ABSIZ(treedata[ichild_1]);
        size_2   = ABSIZ(treedata[ichild_2]);
        size_res = size_1 + size_2;

        next_free_space_index -= size_res;

        PTR(treedata[inode]) = &(treelimbs[next_free_space_index]);

        ALLOC(treedata[inode]) = size_res;

        if (size_1 > size_2) {
            mpn_mul(PTR(treedata[inode]),
                    PTR(treedata[ichild_1]), size_1,
                    PTR(treedata[ichild_2]), size_2);
        } else {
            mpn_mul(PTR(treedata[inode]),
                    PTR(treedata[ichild_2]), size_2,
                    PTR(treedata[ichild_1]), size_1);
        }
        //
        // Normalize the size of the result...
        //
        MPN_NORMALIZE(PTR(treedata[inode]), size_res);

        if ((SIZ(treedata[ichild_1]) ^ SIZ(treedata[ichild_2])) >= 0) {
            SIZ(treedata[inode]) = size_res;
        } else {
            SIZ(treedata[inode]) = -size_res;
        }
    }
    //
    // Root node is computed separately...
    //
    size_1   = ABSIZ(treedata[1]);
    size_2   = ABSIZ(treedata[2]);
    size_res = size_1 + size_2;

    ALLOC(treedata[0]) = size_res;

    if (size_1 > size_2) {
        mpn_mul(PTR(treedata[0]),
                PTR(treedata[1]), size_1,
                PTR(treedata[2]), size_2);
    } else {
        mpn_mul(PTR(treedata[0]),
                PTR(treedata[2]), size_2,
                PTR(treedata[1]), size_1);
    }
    MPN_NORMALIZE(PTR(treedata[0]), size_res);
    if ((SIZ(treedata[1]) ^ SIZ(treedata[2])) >= 0) {
        SIZ(treedata[0]) = size_res;
    } else {
        SIZ(treedata[0]) = -size_res;
    }
    return tree;
}
//-----------------------------------------------------------------------------
mpz_tree_t* rem_tree(const mpz_t z, const mpz_tree_t* const ptree) {
    //
    // _NOTE_: The remainder tree is implemented as a single mpz_array_t tree
    //         with the usual compact representation: tree[2i+1] and tree[2i+2]
    //         are the children of the node tree[i].
    //
    //         See comments in body of function prod_tree for more information.
    //
    // _WARNING_:
    //
    //         Although the product tree returned is actually a pointer
    //         to an mpz_array_t structure, the array should NOT be modified
    //         later on since the memory it uses is allocated in one huge
    //         block to prevent overhead from multiple malloc calls. So the
    //         allocated memory of the mpz_t's in the array can NOT be
    //         increased...
    //
    // _UPDATE_: Mon Nov 20 2006 by JM
    //
    //         This function has been rewritten to perform the
    //         division directly with the mpn layer. Contrary to what
    //         happens with the prod_tree* functions, the switch from
    //         mpz to mpn leads to a much more noticeable speed-up...
    //         (Typically something like 30%)
    //
    // _UPDATE_: Tue Dec 12 2006 by JM
    //
    //         This function has been modified to allocate all the needed
    //         memory with a single call to malloc (see the warning). This
    //         is more than twice as fast as a straitforward mpz implementation.
    //
    mpz_tree_t* const res = malloc(sizeof(mpz_tree_t));

    if (0 == ptree->length) {
        res->alloced = 0;
        res->length  = 0;
        res->data    = NULL;
        return res;
    }
    res->alloced = ptree->length;
    res->length  = res->alloced;
    res->data    = malloc(res->alloced * sizeof(mpz_t));
    //
    // Use local constant pointers to the data to by pass the reading of
    // non constant pointers. This should hopefully speed-up things a bit...
    // A (very) little bit...
    //
    mpz_t* const resdata   = res->data;
    mpz_t* const ptreedata = ptree->data;

    uint32_t ichild_1    = 0;
    uint32_t ichild_2    = 0;
    uint32_t const max_i = res->length >> 1;

    //
    // _KLUDGE_: We use our knowledge of the internal representation of the
    //           mpz_t type to switch from the mpz to the mpn layer to speed-up
    //           computation. Maintainers do not need to panic: the mpz_t type
    //           is completely stable and has no reason to change in the
    //           foreseeable future.
    //
    mp_limb_t* numdata  = NULL;
    mp_size_t  numsize  = 0;
    mp_size_t  numalloc = 0;

    mp_limb_t* dendata  = NULL;
    mp_size_t  denalloc = 0;

    mp_size_t quotalloc = 0;
    mp_size_t remalloc  = 0;
    //
    // We don't need to keep the value of the quotients. So allocate enough
    // memory for the largest possible quotient and reuse this array...
    //
    mp_limb_t* quotdata = malloc(ABSIZ(z) * sizeof(mp_limb_t));

    //
    // Computes the memory needed to hold the remainder tree and allocate the
    // whole memory needed in only one call to malloc.
    //
    uint32_t size_tree   = 0;
    for (uint32_t inode = 0U; inode < ptree->length; inode++) {
        size_tree += ABSIZ(ptree->data[inode]);
    }
    PTR(resdata[0]) = malloc(size_tree * sizeof(mp_limb_t));
    mp_limb_t* const reslimbs = PTR(resdata[0]);
    //
    // Sets each mpz_t's data pointer to its "dedicated" memory location
    //
    uint32_t  next_index = 0;
    for (uint32_t inode = 0U; inode < ptree->length; inode++) {
        PTR(resdata[inode])   = &(reslimbs[next_index]);
        ALLOC(resdata[inode]) = ABSIZ(ptree->data[inode]);

        next_index += ABSIZ(ptree->data[inode]);
    }
    //
    // Compute z mod ptreedata[0] and store the result in resdata[0]
    //
    numdata  = PTR(z);
    numsize  = SIZ(z);
    numalloc = ABS(numsize);

    dendata  = PTR(ptreedata[0]);
    denalloc = ABSIZ(ptreedata[0]);

    quotalloc = numalloc - denalloc + 1;
    remalloc  = denalloc;

    mp_limb_t* remdata = PTR(resdata[0]);

    if (quotalloc <= 0) {
        if (numsize >= 0) {
            //
            // The result is the same as the numerator, so just copy it. memcpy
            // works here because the data we copy does not contain pointers to
            // other data locations...
            //
            memcpy(remdata, numdata, numalloc * sizeof(mp_limb_t));
            resdata[0]->_mp_size  = numsize;

        } else {
            mpn_sub(remdata, dendata, denalloc, numdata, numalloc);

            MPN_NORMALIZE(remdata, denalloc);

            resdata[0]->_mp_size = denalloc;
        }
    } else {

        mpn_tdiv_qr(quotdata, remdata, 0, numdata, numalloc, dendata, denalloc);

        if (numsize < 0) {
            mpn_sub(remdata, dendata, denalloc, remdata, denalloc);
        }

        MPN_NORMALIZE(remdata, remalloc);

        resdata[0]->_mp_size  = remalloc;
    }
    for (uint32_t i = 0U; i < max_i; i++) {

        ichild_1 = (i<<1) | 1;
        ichild_2 = (i<<1) + 2;
        //
        // Compute resdata[i] mod ptreedata[2*i+1] and store the result in
        // resdata[2*i+1]
        //
        numdata  = PTR(resdata[i]);
        numsize  = SIZ(resdata[i]);
        numalloc = ABS(numsize);

        dendata  = PTR(ptreedata[ichild_1]);
        denalloc = ABSIZ(ptreedata[ichild_1]);

        quotalloc = numalloc - denalloc + 1;
        remalloc  = denalloc;

        remdata = PTR(resdata[ichild_1]);

        if (quotalloc <= 0) {
            if (numsize >= 0) {
                memcpy(remdata, numdata, numalloc * sizeof(mp_limb_t));
                resdata[ichild_1]->_mp_size  = numsize;
            } else {
                mpn_sub(remdata, dendata, denalloc, numdata, numalloc);

                MPN_NORMALIZE(remdata, denalloc);

                resdata[ichild_1]->_mp_size = denalloc;
            }
        } else {
            mpn_tdiv_qr(quotdata, remdata, 0,
                        numdata, numalloc, dendata, denalloc);

            if (numsize < 0) {
                mpn_sub(remdata, dendata, denalloc, remdata, denalloc);
            }

            MPN_NORMALIZE(remdata, remalloc);

            resdata[ichild_1]->_mp_size  = remalloc;
        }
        //
        // Compute resdata[i] mod ptreedata[2*i+2] and store the result in
        // resdata[2*i+2]
        //
        denalloc = ABSIZ(ptreedata[ichild_2]);
        dendata  = PTR(ptreedata[ichild_2]);

        quotalloc = numalloc - denalloc + 1;
        remalloc  = denalloc;

        remdata = PTR(resdata[ichild_2]);

        if (quotalloc <= 0) {
            if (numsize >= 0) {
                memcpy(remdata, numdata, numalloc * sizeof(mp_limb_t));
                resdata[ichild_2]->_mp_size  = numsize;
            } else {
                mpn_sub(remdata, dendata, denalloc, numdata, numalloc);

                MPN_NORMALIZE(remdata, denalloc);

                resdata[ichild_2]->_mp_size = denalloc;
            }

        } else {
            mpn_tdiv_qr(quotdata, remdata, 0,
                        numdata, numalloc, dendata, denalloc);

            if (numsize < 0) {
                mpn_sub(remdata, dendata, denalloc, remdata, denalloc);
            }
            MPN_NORMALIZE(remdata, remalloc);

            resdata[ichild_2]->_mp_size  = remalloc;
        }
    }
    free(quotdata);

    return res;
}
//-----------------------------------------------------------------------------
void clear_mpz_tree(mpz_tree_t* tree) {
    free(PTR(tree->data[0]));
    free(tree->data);
    free(tree);
}
//-----------------------------------------------------------------------------
void print_mpz_tree(const mpz_tree_t* const tree) {
    //
    // Mostly for debugging purposes as the output is not particularly
    // well structured...
    //
    print_mpz_tree_rec(tree, 0U, 0U);
}
//-----------------------------------------------------------------------------
void print_mpz_tree_rec(const mpz_tree_t* const tree,
                        uint32_t curpos, uint32_t indent) {

    if (curpos < tree->length) {
        for (uint32_t i = 0U; i < indent; i++) {
            printf("    ");
        }
        gmp_printf("> %Zd\n", tree->data[curpos]);
    }
    if (curpos < (tree->length-1)/2) {
        indent++;
        print_mpz_tree_rec(tree, 2*curpos+1, indent);
        print_mpz_tree_rec(tree, 2*curpos+2, indent);
    }
}
//-----------------------------------------------------------------------------
