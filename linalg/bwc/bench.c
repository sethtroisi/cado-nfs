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
#include <errno.h>
#include <limits.h>
#include <inttypes.h>

#include "matmul.h"
#include "abase.h"
#include "macros.h"
#include "params.h"
#include "debug.h"

void usage()
{
    fprintf(stderr,
            "Usage: ./bench [--impl <implementation>] [--transpose] <file>\n");
    exit(1);
}

int main(int argc, char * argv[])
{
    abobj_t xx MAYBE_UNUSED;
    abobj_init(xx);

    matmul_t mm;

    const char * impl = "sliced";
    const char * file = NULL;
    int transpose = 0;
    int nocheck = 0;
    int nchecks = 4;

    param_list pl;

    param_list_init(pl);

    argv++,argc--;
    param_list_configure_knob(pl, "--transpose", &transpose);
    param_list_configure_knob(pl, "--nocheck", &nocheck);
    param_list_configure_alias(pl, "--transpose", "-t");

    int wild = 0;
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        if (wild == 0) {
            file = argv[0];
            argv++,argc--;
            wild++;
            continue;
        }
        fprintf (stderr, "Unknown option: %s\n", argv[0]);
        usage();
    }

    unsigned int nmax = UINT_MAX;
    param_list_parse_uint(pl, "nmax", &nmax);

    double tmax = 100.0;
    param_list_parse_double(pl, "tmax", &tmax);
    param_list_parse_int(pl, "nchecks", &nchecks);

    const char * tmp;
    if ((tmp = param_list_lookup_string(pl, "impl")) != NULL) {
        impl = tmp;
    }

    unsigned int nbys;

    if (param_list_parse_uint(pl, "nbys", &nbys)) {
        /* The api mandates that we set the desired value for nbys. Here,
         * that's not really our intent, since we really want to bench
         * the layer in its favorite working context. Most of the time,
         * setting nbys is pointless.
         */
        abobj_set_nbys(xx,nbys);
    }

    if (param_list_warn_unused(pl)) {
        usage();
    }

    fprintf(stderr, "Using implementation \"%s\"\n", impl);

    clock_t t0 = clock();

    mm = matmul_reload_cache(xx, file, impl, pl, !transpose);
    if (mm) {
        fprintf(stderr, "Reusing cache file for %s\n", file);
        fprintf(stderr, "Cache load time %.2fs\n",
                (double) (clock()-t0) / CLOCKS_PER_SEC);
    } else {
        fprintf(stderr, "Building cache file for %s\n", file);
        mm = matmul_build(xx, file, impl, pl, !transpose);
        matmul_save_cache(mm, file);
        fprintf(stderr, "Cache build time %.2fs\n",
                (double) (clock()-t0) / CLOCKS_PER_SEC);
    }

    fprintf (stderr, "%s: %u rows %u cols %" PRIu64 " coeffs\n",
            file, mm->dim[0], mm->dim[1], mm->ncoeffs);

    unsigned int nc = mm->dim[1];
    matmul_aux(mm, MATMUL_AUX_GET_READAHEAD, &nc);
    abt * src = abinit(xx, nc);
    abzero(xx, src, nc);

    unsigned int nr = mm->dim[0];
    matmul_aux(mm, MATMUL_AUX_GET_READAHEAD, &nr);
    abt * dst = abinit(xx, nr);
    abzero(xx, dst, nr);

    setup_seeding(1);

    abrandom(xx, src, mm->dim[1]);

    if (!nocheck) {
        abt * src2 = abinit(xx, nc);
        abt * dst2 = abinit(xx, nr);
        abzero(xx, src2, nc);
        abzero(xx, dst2, nr);

        abt * checkA = abinit(xx, abnbits(xx));
        abt * checkB = abinit(xx, abnbits(xx));

        for(int t = 0; t < nchecks ; t++) {
            abrandom(xx, src, mm->dim[1]);
            abrandom(xx, dst, mm->dim[0]);
            abcopy(xx, src2, src, mm->dim[1]);
            abcopy(xx, dst2, dst, mm->dim[0]);
            matmul_mul(mm, dst, src, 1);
            debug_write(dst, abbytes(xx, mm->dim[0]), "/tmp/Lmul.%d", t);

            abdotprod(xx, checkA, dst, dst2, mm->dim[0]);
            matmul_mul(mm, src2, dst2, 0);
            debug_write(src2, abbytes(xx, mm->dim[1]), "/tmp/Rmul.%d", t);

            abdotprod(xx, checkB, src, src2, mm->dim[1]);

            if (memcmp(checkA, checkB, aboffset(xx, abnbits(xx)) * sizeof(abt)) != 0) {
                fprintf(stderr, "Check %d failed\n", t);
                abort();
            }
        }
        printf("All %d checks passed\n", nchecks);

        abclear(xx, checkA, abnbits(xx));
        abclear(xx, checkB, abnbits(xx));

        abclear(xx, src2, nc);
        abclear(xx, dst2, nr);
    }

    t0 = clock();
    double next = 0.25 * CLOCKS_PER_SEC;
    if (next > tmax * CLOCKS_PER_SEC) { next = tmax * CLOCKS_PER_SEC; }
    for(unsigned int n = 0 ; n < nmax ; n++ ) {
        matmul_mul(mm, transpose ? src : dst, transpose ? dst : src, !transpose);
        double dt = clock() - t0;
        if (dt > next) {
            do { next += 0.25 * CLOCKS_PER_SEC; } while (dt > next);
            if (next > tmax * CLOCKS_PER_SEC) { next = tmax * CLOCKS_PER_SEC; }
            dt /= CLOCKS_PER_SEC;
            printf("%d iterations in %2.fs, %.2f/1, %.2f ns/coeff        \r",
                    n, dt, dt/n, 1.0e9 * dt/n/mm->ncoeffs);
            fflush(stdout);
            if (dt > tmax)
                break;
        }
    }
    printf("\n");

    matmul_clear(mm);

    param_list_clear(pl);

    return 0;
}
