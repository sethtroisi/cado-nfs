#define _BSD_SOURCE

/* This program is the simplest interface to the bare matrix
 * multiplication routine. It's meant to provide an easy way of benching,
 * and comparing, different matrix product implementations.
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

    mm = matmul_build(xx, argv[1]);
    matmul_save_cache(mm, argv[1]);

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
