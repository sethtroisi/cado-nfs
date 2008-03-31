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
 * \file    matrix.c
 * \author  Jerome Milan
 * \date    Circa Fri Mar 29 2008
 * \version 1.1.1
 */

/*
 *  History:
 *  1.1.1: Fri Mar 29 (?) 2008 by JM:
 *        - Inlined functions pertaining to binary_matrix_t's in matrix.h and
 *          removed them from this file.
 *    1.1: Tue Jun 19 2007 by JM:
 *          - Added clone_binary_matrix function.
 *    1.0: Wed Mar 1 2006 by JM:
 *          - Initial version.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "matrix.h"

#if TIFA_USE_CALLOC_MEMSET
    #include <string.h>
#endif

#include "bitstring_t.h"

/*
 *-----------------------------------------------------------------------------
 *              binary_matrix_t and associated functions
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
void reset_binary_matrix(binary_matrix_t* const matrix) {
    for (uint32_t i = 0U; i < matrix->nrows_alloced; i++) {
#if TIFA_USE_CALLOC_MEMSET
        memset(matrix->data[i], 0x00,
               matrix->ncols_alloced * SIZEOF_BITSTRING_T);
#else
        for (uint32_t j = 0U; j < matrix->ncols_alloced; j++) {
            matrix->data[i][j] = 0;
            //
            // Use an xor instead of simply 0 to make sure we have a bitstring
            // made of 0 bits only. Granted, TIFA is not likely to compile on
            // the exotic architectures where the value 0 is not made of 0 bits
            // only anyways, but still...
            //
            matrix->data[i][j] ^= 0;
        }
#endif
    }
}
//-----------------------------------------------------------------------------
binary_matrix_t* alloc_binary_matrix(uint32_t nrows, uint32_t ncols) {

    binary_matrix_t* matrix = malloc(sizeof(binary_matrix_t));

    matrix->nrows_alloced = nrows;

#if BITSTRING_T_SIZE_IS_POW_OF_TWO
    matrix->ncols_alloced = ((ncols - 1) >> POW_TWO_BITSTRING_T_SIZE) + 1;
#else
    matrix->ncols_alloced = (ncols - 1)/BITSTRING_T_BITSIZE + 1;
#endif

    matrix->nrows = nrows;
    matrix->ncols = ncols;
    matrix->data  = malloc(matrix->nrows_alloced * sizeof(TIFA_BITSTRING_T*));

    for (uint32_t i = 0U; i < matrix->nrows_alloced; i++) {
        matrix->data[i] = malloc(matrix->ncols_alloced * SIZEOF_BITSTRING_T);
    }

    reset_binary_matrix(matrix);

    return matrix;
}
//-----------------------------------------------------------------------------
binary_matrix_t* clone_binary_matrix(const binary_matrix_t * const matrix) {

    binary_matrix_t* clone = malloc(sizeof(binary_matrix_t));

    clone->nrows_alloced = matrix->nrows_alloced;
    clone->ncols_alloced = matrix->ncols_alloced;
    clone->nrows = matrix->nrows;
    clone->ncols = matrix->ncols;
    clone->data  = malloc(clone->nrows_alloced * sizeof(TIFA_BITSTRING_T*));

    size_t rowsize = clone->ncols_alloced * SIZEOF_BITSTRING_T;

    for (uint32_t i = 0U; i < clone->nrows_alloced; i++) {
        clone->data[i] = malloc(rowsize);
        memcpy(clone->data[i], matrix->data[i], rowsize);
    }
    return clone;
}
//-----------------------------------------------------------------------------
void clear_binary_matrix(binary_matrix_t* const matrix) {
    for (uint32_t i = 0U; i < matrix->nrows_alloced; i++) {
        free(matrix->data[i]);
    }
    free(matrix->data);

    matrix->nrows_alloced = 0U;
    matrix->ncols_alloced = 0U;
    matrix->nrows = 0U;
    matrix->ncols = 0U;
}
//-----------------------------------------------------------------------------
void print_binary_matrix(const binary_matrix_t* const matrix) {
    //
    // Mostly for debugging purposes...
    //
    for (uint32_t i = 0U; i < matrix->nrows; i++) {
        printf("row %4u : ", i);
        for (uint32_t j = 0U; j < matrix->ncols; j++) {
            if (0 != get_matrix_bit(i, j, matrix)) {
                printf("1");
            } else {
                printf(".");
            }
        }
        printf("\n");
    }
}
//-----------------------------------------------------------------------------

/*
 *-----------------------------------------------------------------------------
 *              byte_matrix_t and associated functions
 *-----------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
void reset_byte_matrix(byte_matrix_t* const matrix) {
    for (uint32_t i = 0U; i < matrix->nrows_alloced; i++) {
#if TIFA_USE_CALLOC_MEMSET
        memset(
            matrix->data[i],
            0x00,
            matrix->ncols_alloced * sizeof(unsigned char)
        );
#else
        for (uint32_t j = 0U; j < matrix->ncols_alloced; j++) {
            matrix->data[i][j] = 0;
        }
#endif
    }
}
//-----------------------------------------------------------------------------
byte_matrix_t* alloc_byte_matrix(uint32_t nrows, uint32_t ncols) {

    byte_matrix_t* matrix = malloc(sizeof(byte_matrix_t));

    matrix->nrows_alloced = nrows;
    matrix->ncols_alloced = ncols;

    matrix->nrows = nrows;
    matrix->ncols = ncols;

    matrix->data  = malloc(matrix->nrows_alloced * sizeof(unsigned char*));

    for (uint32_t i = 0U; i < matrix->nrows_alloced; i++) {
#if TIFA_USE_CALLOC_MEMSET
        matrix->data[i] = calloc(matrix->ncols_alloced, sizeof(unsigned char));
#else
        matrix->data[i] = malloc(matrix->ncols_alloced * sizeof(unsigned char));
#endif
    }

#if ! TIFA_USE_CALLOC_MEMSET
    reset_byte_matrix(matrix);
#endif

    reset_byte_matrix(matrix);

    return matrix;
}
//-----------------------------------------------------------------------------
byte_matrix_t* clone_byte_matrix(const byte_matrix_t * const matrix) {

    byte_matrix_t* clone = malloc(sizeof(byte_matrix_t));

    clone->nrows_alloced = matrix->nrows_alloced;
    clone->ncols_alloced = matrix->ncols_alloced;
    clone->nrows = matrix->nrows;
    clone->ncols = matrix->ncols;
    clone->data  = malloc(clone->nrows_alloced * sizeof(unsigned char));

    size_t rowsize = clone->ncols_alloced * sizeof(unsigned char);

    for (uint32_t i = 0U; i < clone->nrows_alloced; i++) {
        clone->data[i] = malloc(rowsize);
        memcpy(clone->data[i], matrix->data[i], rowsize);
    }
    return clone;
}
//-----------------------------------------------------------------------------
void clear_byte_matrix(byte_matrix_t* const matrix) {
    for (uint32_t i = 0U; i < matrix->nrows_alloced; i++) {
        free(matrix->data[i]);
    }
    free(matrix->data);

    matrix->nrows_alloced = 0U;
    matrix->ncols_alloced = 0U;
    matrix->nrows = 0U;
    matrix->ncols = 0U;
}
//-----------------------------------------------------------------------------
void print_byte_matrix(const byte_matrix_t* const matrix) {
    //
    // Mostly for debugging purposes...
    //
    for (uint32_t i = 0U; i < matrix->nrows; i++) {
        printf("row %4u : ", i);
        for (uint32_t j = 0U; j < matrix->ncols; j++) {
            printf("%4u", (unsigned int)matrix->data[i][j]);
        }
        printf("\n");
    }
}
//-----------------------------------------------------------------------------
