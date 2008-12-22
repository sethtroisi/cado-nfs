#define _BSD_SOURCE

/* This program is the simplest interface to the bare matrix
 * multiplication routine. It's meant to provide an easy way of benching,
 * and comparing, different matrix product implementations.
 *
 * It must be used with the INNER matrix, not the .info one. So either use
 * the vanilla cado format (before balance), or one of the .h<i>.v<j>
 * matrices as output by balance.
 */

#include <stdio.h>
#include <time.h>

#include "matmul.h"
#include "abase.h"
#include "macros.h"

int main(int argc, char * argv[])
{
    abobj_t xx MAYBE_UNUSED;
    abobj_init(xx);

    matmul_t mm;

    if (argc != 2) {
        fprintf(stderr, "Usage: ./bench <file>\n");
        exit(1);
    }

    mm = matmul_reload_cache(xx, argv[1]);
    if (mm) {
        fprintf(stderr, "Reusing cache file for %s\n", argv[1]);
    } else {
        fprintf(stderr, "Building cache file for %s\n", argv[1]);
        mm = matmul_build(xx, argv[1]);
        matmul_save_cache(mm, argv[1]);
    }

    fprintf(stderr, "%s: %u rows %u cols %lu coeffs\n",
            argv[1], mm->nrows, mm->ncols, mm->ncoeffs);

    abt * src = abinit(xx, mm->ncols);
    abt * dst = abinit(xx, mm->nrows);

    clock_t t0 = clock();

    for(int n = 0 ;  ; n++ ) {
        matmul(mm, dst, src);
        double dt = clock() - t0;
        dt /= CLOCKS_PER_SEC;
        printf("%d iterations in %2.fs, %.2f/1, %.2f ns/coeff        \r",
                n, dt, dt/n, 1.0e9 * dt/n/mm->ncoeffs);
        if (dt > 100)
            break;
    }
    printf("\n");

    matmul_clear(mm);

    return 0;
}
