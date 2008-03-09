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
 * \file    test_bernsteinisms.c
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
#include <string.h>
#include <inttypes.h>

#include "array.h"
#include "x_tree.h"
#include "bernsteinisms.h"
#include "test_bernsteinisms.h"

#include "first_primes.h"

//------------------------------------------------------------------------------
int main(int argc, char** argv) {

    if (2 != argc) {
        printf("Usage: %s (21|51|53|63|71)\n", argv[0]);
        return -1;
    }
    if (0 == strcmp(argv[1], "21")) {
        do_test_21();
        return 0;
    }
    if (0 == strcmp(argv[1], "51")) {
        do_test_51();
        return 0;
    }
    if (0 == strcmp(argv[1], "53")) {
        do_test_53();
        return 0;
    }
    if (0 == strcmp(argv[1], "63")) {
        do_test_63();
        return 0;
    }
    if (0 == strcmp(argv[1], "71")) {
        do_test_71();
        return 0;
    }
    printf("%s : invalid test number\n", argv[1]);
    return -1;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void do_test_71() {

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);

    FILE *devurand;
    uint32_t seed;
    if ((devurand = fopen("/dev/urandom","r")) == 0) {
        printf("Cannot open /dev/urandom");
        exit(1);
    }
    fread(&seed, sizeof(seed), 1, devurand);
    fclose(devurand);
    printf("/dev/urandom gives : %u\n", seed);

    gmp_randseed_ui(rstate, seed);

    mpz_array_t int_to_fact;
    int_to_fact.alloced = 8;
    int_to_fact.length = int_to_fact.alloced;
    int_to_fact.data = malloc(int_to_fact.alloced*sizeof(mpz_t));
    for (uint32_t i = 0U; i < int_to_fact.length; i++) {
        mpz_init(int_to_fact.data[i]);
        mpz_urandomb(int_to_fact.data[i], rstate, 12);
        if (0 == mpz_tstbit(int_to_fact.data[i], 0)) {
            mpz_add_ui(int_to_fact.data[i], int_to_fact.data[i], 1);
        }
        while (0 != mpz_probab_prime_p(int_to_fact.data[i], 10)) {
            mpz_init(int_to_fact.data[i]);
            mpz_urandomb(int_to_fact.data[i], rstate, 12);
            if (0 == mpz_tstbit(int_to_fact.data[i], 0)) {
                mpz_add_ui(int_to_fact.data[i], int_to_fact.data[i], 1);
            }
        }
    }

    uint32_array_t primes;
    primes.alloced = 16;
    primes.length  = primes.alloced;
    primes.data    = malloc(primes.alloced*sizeof(uint32_t));
    for (uint32_t i = 0; i < primes.length; i++) {
        primes.data[i] = first_primes[i+1];
    }
    printf("Prime factor base:\n");
    printf("------------------------------\n");
    print_uint32_array((const uint32_array_t*)&primes);
    printf("------------------------------\n");
    printf("Integers to factor :\n");
    printf("------------------------------\n");
    print_mpz_array((const mpz_array_t*)&int_to_fact);
    printf("------------------------------\n");

    uint32_array_list_t* list = alloc_uint32_array_list(int_to_fact.length);

    algo_71(list, (const mpz_array_t*)&int_to_fact,
            (const uint32_array_t*)&primes);

    printf("Smooth factors of the integers to factor :\n");
    printf("------------------------------------------\n");
    print_uint32_array_list(list);

    gmp_randclear(rstate);
    clear_uint32_array_list(list);
    clear_uint32_array(&primes);
    clear_mpz_array(&int_to_fact);
}
//------------------------------------------------------------------------------
void do_test_63() {

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);

    FILE *devurand;
    uint32_t seed;
    if ((devurand = fopen("/dev/urandom","r")) == 0) {
        printf("Cannot open /dev/urandom");
        exit(1);
    }
    fread(&seed, sizeof(seed), 1, devurand);
    fclose(devurand);
    printf("/dev/urandom gives : %u\n", seed);

    gmp_randseed_ui(rstate, seed);
    mpz_t x;
    mpz_init(x);
    mpz_urandomb(x, rstate, 9);
    if (0 == mpz_tstbit(x, 0)) {
        mpz_add_ui(x, x, 1);
    }

    uint32_t length = 12;
    uint32_array_t primes;
    primes.alloced = length;
    primes.length  = length;
    primes.data    = malloc(primes.alloced*sizeof(uint32_t));
    for (uint32_t i = 0; i < primes.length; i++) {
        primes.data[i] = first_primes[i+1];
    }
    printf("---- Used primes ------\n");
    print_uint32_array(&primes);
    printf("-----------------------\n");

    gmp_printf("\n Integer to factor : %Zd\n\n", x);


    mpz_tree_t* tree = prod_tree_ui(&primes);

    uint32_array_t* divs = algo_63(x, tree);
    print_uint32_array(divs);

    gmp_randclear(rstate);

    clear_mpz_tree(tree);
    clear_uint32_array(divs);
    clear_uint32_array(&primes);
}
//------------------------------------------------------------------------------
void do_test_53() {

    mpz_t u;
    mpz_init_set_ui(u, 271);

    uint32_t b = 7;
    mpz_t x;
    mpz_init_set_ui(x, 103);

    mpz_t *r = algo_53(b, u, x);

    gmp_printf("b = %d\n", b);
    gmp_printf("u = %Zd\n", u);
    gmp_printf("x = %Zd\n", x);
    gmp_printf("r = %Zd\n\n", *r);

    mpz_mul_2exp(*r, *r, b);
    gmp_printf("r.2^b        = %Zd\n", *r);
    mpz_mod(*r, *r, u);
    gmp_printf("r.2^b mod u  = %Zd\n", *r);
    mpz_mod(x, x, u);
    gmp_printf("x mod u      = %Zd\n", x);

    mpz_clear(*r);
    free(r);
    mpz_clear(u);
    mpz_clear(x);

}
//------------------------------------------------------------------------------
void do_test_51() {

    FILE *devurand;
    uint32_t seed;
    if ((devurand = fopen("/dev/urandom","r")) == 0) {
        printf("Cannot open /dev/urandom");
        exit(1);
    }
    fread(&seed, sizeof(seed), 1, devurand);
    printf("/dev/urandom gives : %u\n", seed);


    uint32_t b = (seed & 0xF) + 1;

    fread(&seed, sizeof(seed), 1, devurand);
    fclose(devurand);
    printf("/dev/urandom gives : %u\n", seed);
    fflush(stdout);


    mpz_t u;
    mpz_init_set_ui(u, (seed & 0xFFFF) | 1);

    gmp_printf("b = %d\n", b);
    gmp_printf("u = %Zd\n", u);
    fflush(stdout);

    mpz_t *v = algo_51(b, u);

    gmp_printf("v = %Zd\n", *v);

    mpz_mul(u, u, *v);
    mpz_add_ui(u, u, 1);
    gmp_printf("1+uv         = %Zd\n", u);

    mpz_fdiv_r_2exp(u, u, b);
    gmp_printf("1+uv mod 2^b = %Zd\n", u);

    mpz_clear(*v);
    free(v);
    mpz_clear(u);
}
//------------------------------------------------------------------------------
void do_test_21() {
    uint32_t length_to_fact = 12;
    uint32_t length_primes  = 5;

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);

    FILE *devurand;
    uint32_t seed;
    if ((devurand = fopen("/dev/urandom","r")) == 0) {
        printf("Cannot open /dev/urandom");
        exit(1);
    }
    fread(&seed, sizeof(seed), 1, devurand);
    fclose(devurand);
    printf("/dev/urandom gives : %u\n", seed);

    gmp_randseed_ui(rstate, seed);

    mpz_array_t int_to_fact;
    int_to_fact.alloced = length_to_fact;
    int_to_fact.length = int_to_fact.alloced;
    int_to_fact.data = malloc(int_to_fact.alloced*sizeof(mpz_t));
    for (uint32_t i = 0U; i < int_to_fact.length; i++) {
        mpz_init(int_to_fact.data[i]);
        mpz_urandomb(int_to_fact.data[i], rstate, 8);
    }

    uint32_array_t primes;
    primes.alloced = length_primes;
    primes.length  = primes.alloced;
    primes.data    = malloc(primes.alloced*sizeof(mpz_t));

    for (uint32_t i = 0U; i < primes.length; i++) {
        primes.data[i] = first_primes[i];
    }
    printf("--- Primes used ----\n");
    print_uint32_array(&primes);
    printf("---- Numbers -------\n");
    print_mpz_array(&int_to_fact);

    mpz_array_t* smooth = alloc_mpz_array(int_to_fact.length);
    mpz_tree_t*  ptree  = prod_tree_ui(&primes);

    algo_21_rtree(smooth, &int_to_fact, ptree->data[0]);
    printf("-- Smooth Numbers --\n");
    print_mpz_array(smooth);

    gmp_randclear(rstate);
    clear_mpz_array(&int_to_fact);
    clear_uint32_array(&primes);
    clear_mpz_array(smooth);
    clear_mpz_tree(ptree);
}
//------------------------------------------------------------------------------
