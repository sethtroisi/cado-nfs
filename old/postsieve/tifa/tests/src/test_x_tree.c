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
 * \file    test_x_tree.c
 * \author  Jerome Milan
 * \date    Mon Mar 6 2006
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
#include "x_tree.h"
#include "first_primes.h"

//------------------------------------------------------------------------------
void do_test() {

    uint32_t length = 16;

    mpz_array_t array;
    array.alloced = length;
    array.length  = length;
    array.data    = malloc(array.alloced*sizeof(mpz_t));

    for (uint32_t i = 0U; i < array.length; i++) {
        mpz_init_set_ui(array.data[i], first_primes[i]);
    }
    printf("------------ array ------------\n");
    print_mpz_array(&array);

    mpz_array_t* ptree = prod_tree(&array);
    printf("--------- product tree --------\n");
    print_mpz_tree(ptree);

    mpz_t z;
    mpz_init_set_ui(z, 17534);
    gmp_printf("\nz = %Zd\n\n", z);

    mpz_array_t* rtree = rem_tree(z, ptree);
    printf("-------- remainder tree -------\n");
    print_mpz_tree(rtree);

    clear_mpz_array(ptree);
    clear_mpz_array(rtree);
    mpz_clear(z);
}
//------------------------------------------------------------------------------
int main() {
    do_test();
    return 0;
}
//------------------------------------------------------------------------------
