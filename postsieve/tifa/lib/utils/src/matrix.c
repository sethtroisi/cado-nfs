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
 * \file    matrix.c
 * \author  Jerome Milan
 * \date    Tue Jun 19 2007
 * \version 1.1
 */

 /*
  *  History:
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
uint8_t get_matrix_bit(uint32_t row, uint32_t col,
                       const binary_matrix_t* const matrix) {

#if BITSTRING_T_SIZE_IS_POW_OF_TWO
    uint32_t col_offset = col & (BITSTRING_T_BITSIZE - 1);

    col_offset = BITSTRING_T_BITSIZE - 1 - col_offset;
    col        = col >> POW_TWO_BITSTRING_T_SIZE;
#else
    uint32_t col_offset = col % BITSTRING_T_BITSIZE;

    col_offset = BITSTRING_T_BITSIZE - 1 - col_offset;
    col       /= BITSTRING_T_BITSIZE;
#endif

    if (0 == ((((TIFA_BITSTRING_T)1)<<col_offset) & matrix->data[row][col])) {
        return 0;
    } else {
        return 1;
    }
}
//-----------------------------------------------------------------------------
void set_matrix_bit_to_one(uint32_t row, uint32_t col,
                           binary_matrix_t* const matrix) {

#if BITSTRING_T_SIZE_IS_POW_OF_TWO
    uint32_t col_offset = col & (BITSTRING_T_BITSIZE - 1);

    col_offset = BITSTRING_T_BITSIZE - 1 - col_offset;
    col        = col >> POW_TWO_BITSTRING_T_SIZE;
#else
    uint32_t col_offset = col % BITSTRING_T_BITSIZE;

    col_offset = BITSTRING_T_BITSIZE - 1 - col_offset;
    col       /= BITSTRING_T_BITSIZE;
#endif

    matrix->data[row][col] |= ((TIFA_BITSTRING_T)1<<col_offset);
}
//-----------------------------------------------------------------------------
void set_matrix_bit_to_zero(uint32_t row, uint32_t col,
                            binary_matrix_t* const matrix) {

#if BITSTRING_T_SIZE_IS_POW_OF_TWO
    uint32_t col_offset = col & (BITSTRING_T_BITSIZE - 1);

    col_offset = BITSTRING_T_BITSIZE - 1 - col_offset;
    col        = col >> POW_TWO_BITSTRING_T_SIZE;
#else
    uint32_t col_offset = col % BITSTRING_T_BITSIZE;

    col_offset = BITSTRING_T_BITSIZE - 1 - col_offset;
    col       /= BITSTRING_T_BITSIZE;
#endif

    matrix->data[row][col] &= !(((TIFA_BITSTRING_T)1)<<col_offset);
}
//-----------------------------------------------------------------------------
void flip_matrix_bit(uint32_t row, uint32_t col,
                     binary_matrix_t* const matrix){

#if BITSTRING_T_SIZE_IS_POW_OF_TWO
    uint32_t col_offset = col & (BITSTRING_T_BITSIZE - 1);

    col_offset = BITSTRING_T_BITSIZE - 1 - col_offset;
    col        = col >> POW_TWO_BITSTRING_T_SIZE;
#else
    uint32_t col_offset = col % BITSTRING_T_BITSIZE;

    col_offset = BITSTRING_T_BITSIZE - 1 - col_offset;
    col       /= BITSTRING_T_BITSIZE;
#endif

    matrix->data[row][col] ^= (((TIFA_BITSTRING_T)1)<<col_offset);
}
//-----------------------------------------------------------------------------
uint32_t first_row_with_one_on_col(uint32_t col,
                                   const binary_matrix_t* const matrix) {
    //
    // This function returns the index of the first row which has a 1 in
    // its col-th column. It returns NOT_IN_ARRAY if no such row is found.
    //
#if BITSTRING_T_SIZE_IS_POW_OF_TWO
    uint32_t col_offset = col & (BITSTRING_T_BITSIZE - 1);

    col_offset = BITSTRING_T_BITSIZE - 1 - col_offset;
    col        = col >> POW_TWO_BITSTRING_T_SIZE;
#else
    uint32_t col_offset = col % BITSTRING_T_BITSIZE;

    col_offset = BITSTRING_T_BITSIZE - 1 - col_offset;
    col       /= BITSTRING_T_BITSIZE;
#endif

    for (uint32_t irow = 0; irow < matrix->nrows; irow++) {
        if (0 != ((((TIFA_BITSTRING_T)1)<<col_offset) & matrix->data[irow][col])
           ) {
            return irow;
        }
    }
    //
    // _KLUDGE_: If no row with a one on column col is found, returns -1 or
    //           rather NOT_IN_ARRAY (i.e. UINT32_MAX) which, in a signed 
    //           context, is indeed -1.
    //
    return NO_SUCH_ROW;
}
//-----------------------------------------------------------------------------
