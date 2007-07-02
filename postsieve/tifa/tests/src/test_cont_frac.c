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
 * \file    test_cont_frac.c
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
#include "gmp.h"
#include "sqrt_cont_frac.h"

//------------------------------------------------------------------------------
#define NB_STEP 32
#define MAX_NB_DIGITS 32
//------------------------------------------------------------------------------
void do_test() {

    char str_buffer[MAX_NB_DIGITS];
    printf("> Enter an integer: ");
    char* str_n = fgets(str_buffer, MAX_NB_DIGITS, stdin);

    mpz_t n;
    mpz_init_set_str(n, str_n, 10);

    cont_frac_state_t state;

    mpf_t sqrtn;
    mpf_init(sqrtn);
    mpf_set_z(sqrtn, n);
    mpf_sqrt(sqrtn, sqrtn);

    gmp_printf("\nn       = %Zd\n", n);
    gmp_printf("sqrt(n) = %Ff\n\n", sqrtn);

    init_cont_frac_state(&state, n);
    step_cont_frac_state(&state, NB_STEP-1);

    //mpf_t fa;
    //mpf_t fb;

    //mpf_init(fa);
    //mpf_init(fb);

    //mpf_set_z(fa, state.a);
    //mpf_set_z(fb, state.b);

    //mpf_div(fa, fa, fb);

    printf("nb_step_performed = %"PRIu32"\n\n", state.nb_step_performed);

    gmp_printf("a     = %Zd\n", state.a);
    //gmp_printf("b     = %Zd\n", state.b);
    //gmp_printf("a/b   = %Ff\n\n", fa);

    //mpz_t a2mb2;
    //mpz_init(a2mb2);
    //mpz_mul(a2mb2, state.b, state.b);
    //mpz_mul(a2mb2, a2mb2, n);
    //mpz_neg(a2mb2, a2mb2);
    //mpz_addmul(a2mb2, state.a, state.a);

    if (state.nb_step_performed & 0x1) {
        printf("q     = - a2mb2\n");
    } else {
        printf("q     = + a2mb2\n");
    }
    gmp_printf("q     = %Zd\n", state.q);
    //gmp_printf("a2mb2 = %Zd\n", a2mb2);

    mpz_clear(n);
    mpf_clear(sqrtn);
    //mpf_clear(fa);
    //mpf_clear(fb);
    //mpz_clear(a2mb2);

    clear_cont_frac_state(&state);
}
//------------------------------------------------------------------------------
int main() {
    do_test();
    return 0;
}
//------------------------------------------------------------------------------
