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
 * \file    test_gauss_elim.c
 * \author  Jerome Milan
 * \date    Fri Mar 2 2006
 * \version 1.0
 */

 /*
  *  Copyright (C) 2006, 2007 INRIA
  *  License: GNU Lesser General Public License (LGPL)
  *  History:
  */

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>

#include "array.h"
#include "x_array_list.h"
#include "matrix.h"
#include "gauss_elim.h"

//------------------------------------------------------------------------------
void do_test() {

    uint32_t nrows = 9;
    uint32_t ncols = 7;

    binary_matrix_t* matrix = alloc_binary_matrix(nrows, ncols);
/*
    set_matrix_bit_to_one(0, 2, matrix);
    set_matrix_bit_to_one(0, 5, matrix);

    set_matrix_bit_to_one(1, 1, matrix);
    set_matrix_bit_to_one(1, 5, matrix);
    set_matrix_bit_to_one(1, 6, matrix);

    set_matrix_bit_to_one(2, 0, matrix);
    set_matrix_bit_to_one(2, 4, matrix);

    set_matrix_bit_to_one(3, 2, matrix);
    set_matrix_bit_to_one(3, 3, matrix);
    set_matrix_bit_to_one(3, 6, matrix);

    set_matrix_bit_to_one(4, 0, matrix);
    set_matrix_bit_to_one(4, 1, matrix);
    set_matrix_bit_to_one(4, 4, matrix);

    set_matrix_bit_to_one(5, 3, matrix);
    set_matrix_bit_to_one(5, 6, matrix);

    set_matrix_bit_to_one(6, 3, matrix);
    set_matrix_bit_to_one(6, 4, matrix);

    set_matrix_bit_to_one(7, 0, matrix);
    set_matrix_bit_to_one(7, 2, matrix);
    set_matrix_bit_to_one(7, 5, matrix);

    set_matrix_bit_to_one(8, 1, matrix);
    set_matrix_bit_to_one(8, 3, matrix);
    set_matrix_bit_to_one(8, 6, matrix);
*/
    set_matrix_bit_to_one(0, 3, matrix);
    set_matrix_bit_to_one(0, 6, matrix);

    set_matrix_bit_to_one(1, 1, matrix);
    set_matrix_bit_to_one(1, 2, matrix);
    set_matrix_bit_to_one(1, 4, matrix);

    set_matrix_bit_to_one(2, 1, matrix);
    set_matrix_bit_to_one(2, 4, matrix);

    set_matrix_bit_to_one(3, 0, matrix);
    set_matrix_bit_to_one(3, 6, matrix);

    set_matrix_bit_to_one(4, 4, matrix);
    set_matrix_bit_to_one(4, 6, matrix);

    set_matrix_bit_to_one(5, 1, matrix);
    set_matrix_bit_to_one(5, 2, matrix);

    set_matrix_bit_to_one(6, 0, matrix);
    set_matrix_bit_to_one(6, 5, matrix);

    set_matrix_bit_to_one(7, 3, matrix);
    set_matrix_bit_to_one(7, 1, matrix);

    set_matrix_bit_to_one(8, 1, matrix);
    set_matrix_bit_to_one(8, 3, matrix);
    set_matrix_bit_to_one(8, 5, matrix);

    printf("---------------------------\n");
    printf("Original matrix\n");
    printf("---------------------------\n");
    print_binary_matrix(matrix);

    uint32_array_list_t* relations;
    relations = alloc_uint32_array_list(matrix->nrows - matrix->ncols);
    gaussian_elim(relations, matrix);

    printf("---------------------------\n");
    printf("Matrix completely processed\n");
    printf("---------------------------\n");
    print_binary_matrix(matrix);

    printf("---------------------------\n");
    printf("Relations obtained\n");
    printf("---------------------------\n");

    print_uint32_array_list(relations);

    clear_binary_matrix(matrix);
    clear_uint32_array_list(relations);
}
//------------------------------------------------------------------------------
int main() {
    do_test();
    return 0;
}
//------------------------------------------------------------------------------
