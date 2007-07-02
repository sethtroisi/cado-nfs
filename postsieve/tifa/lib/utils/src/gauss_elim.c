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
 * \file    gauss_elim.c
 * \author  Jerome Milan
 * \date    Wed Mar 1 2006
 * \version 1.0
 */

#include <stdio.h>
#include <stdlib.h>

#include "gauss_elim.h"
#include "array.h"

//-----------------------------------------------------------------------------
// Ref. "A compact algorithm for Gaussian elimination over GF(2)
//       implemented on highly parallel computers",
//       Dennis Parkinson, Marvin Wunderlich,
//       Parallel Computing 1 (1984) 65-73
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void gaussian_elim(uint32_array_list_t* relations,
                   binary_matrix_t* const matrix) {
    //
    // Parkinson and Wunderlich's algorithm is straitforwardly implemented
    // without any kind of enhancement.
    //
    // _WARNING_: The aformentionned paper contains an omission in the
    //            description of the algorithm. See comments for futher
    //            explanation...
    //
    uint32_t nrow = matrix->nrows;
    uint32_t ncol = matrix->ncols;

    uint32_t row = 0;

    binary_array_t* used  = alloc_binary_array(nrow);
    uint32_array_t* pivot = alloc_uint32_array(ncol);
    used->length  = nrow;
    //
    // Stage 1: Process the matrix
    //
    // _NOTE_: We process the matrix from right to left since the rightmost
    //         part of the matrix is sparser than the leftmost part
    //
    for (uint32_t icol = ncol - 1; icol != UINT32_MAX; icol--) {

        row = first_row_with_one_on_col(icol, matrix);
        if (row == UINT32_MAX) {
            //
            // No row found, so skip to the next column...
            //
            continue;
        }
        pivot->data[icol] = row;

        for (uint32_t irow = row+1; irow < nrow; irow++) {

            if (1 == get_matrix_bit(irow, icol, matrix)) {

                for (uint32_t index = 0U; index < matrix->ncols_alloced;
                     index++) {
                    matrix->data[irow][index] ^= matrix->data[row][index];
                }
                set_matrix_bit_to_one(irow, icol, matrix);
            }
        }
        //
        // _WARNING_: The following is completely omitted in the algorithm
        //            description of P & W's paper. We need to keep the
        //            icol-th bit of the row "row" to one but to zero
        //            all the other bits in this row.
        //
        for (uint32_t index = 0U; index < matrix->ncols_alloced;
             index++) {
            matrix->data[row][index] = 0U;
        }
        set_matrix_bit_to_one(row, icol, matrix);
        //
        // Back to the algorithm as described in the paper...
        //
        set_array_bit_to_one(row, used);
    }
    pivot->length = ncol;
    //
    // Stage 2: Find dependencies
    //
    for (uint32_t irow = 0U; irow < matrix->nrows; irow++) {
        if (relations->length == relations->alloced) {
            break;
        }
        if (0U == get_array_bit(irow, used)) {
            //
            // Computes the size of the uint32_array_t to hold the relation
            //
            uint32_t nrows_in_rel = 0U;
            for (uint32_t icol = 0U; icol < matrix->ncols; icol++) {
                if (1 == get_matrix_bit(irow, icol, matrix)) {
                    nrows_in_rel++;
                }
            }
            //
            // rel's size is indeed (nrows_in_rel+1) as we also have to
            // add the current irow...
            //
            uint32_array_t* rel = alloc_uint32_array(nrows_in_rel+1);
            rel->length = nrows_in_rel+1;

            nrows_in_rel = 0U;
            for (uint32_t icol = 0U; icol < matrix->ncols; icol++) {
                if (1 == get_matrix_bit(irow, icol, matrix)) {
                    rel->data[nrows_in_rel] = pivot->data[icol];
                    nrows_in_rel++;
                }
            }
            rel->data[nrows_in_rel] = irow;
            add_entry_in_uint32_array_list(rel, relations);
        }
    }
    clear_binary_array(used);
    clear_uint32_array(pivot);

    free(used);
    free(pivot);
}
//-----------------------------------------------------------------------------
