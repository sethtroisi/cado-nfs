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

#include "matmul.h"
#include "abase.h"
#include "macros.h"
#include "params.h"

/* Include all matmul implementations here */
#include "matmul-basic.h"
#include "matmul-sliced.h"

struct matmul_bindings_s {
    matmul_ptr (*build)(abobj_ptr, const char * filename);
    matmul_ptr (*reload_cache)(abobj_ptr, const char * filename);
    void (*save_cache)(matmul_ptr, const char * filename);
    void (*mul)(matmul_ptr, abt *, abt const *, int);
    void (*report)(matmul_ptr);
    void (*clear)(matmul_ptr mm);
};

struct matmul_bindings_s bind[1];


#define REBIND_F(kind,func) bind->func = & MATMUL_NAME(kind, func)
#define REBIND_ALL(kind) do {						\
        REBIND_F(kind, build);						\
        REBIND_F(kind, reload_cache);					\
        REBIND_F(kind, save_cache);					\
        REBIND_F(kind, mul);						\
        REBIND_F(kind, report);						\
        REBIND_F(kind, clear);						\
    } while (0)

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

    param_list pl;

    param_list_init(pl);

    argv++,argc--;
    param_list_configure_knob(pl, "--transpose", &transpose);
    param_list_configure_alias(pl, "--transpose", "-t");

    int wild = 0;
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        if (wild == 0) {
            file = argv[0];
            argc--;
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

    const char * tmp;
    if ((tmp = param_list_lookup_string(pl, "impl")) != NULL) {
        impl = tmp;
    }

    if (param_list_warn_unused(pl)) {
        usage();
    }

#define CHECK_REBIND(K) \
    if (strcmp(impl, #K) == 0) { REBIND_ALL(K); } else

    CHECK_REBIND(sliced)        // no semicolon !
    CHECK_REBIND(basic)
    { fprintf(stderr, "no implementation %s known\n", impl); exit(1); }

    fprintf(stderr, "Using implementation \"%s\"\n", impl);

    clock_t t0 = clock();

    mm = (*bind->reload_cache)(xx, file);
    if (mm) {
        fprintf(stderr, "Reusing cache file for %s\n", file);
        fprintf(stderr, "Cache load time %.2fs\n",
                (double) (clock()-t0) / CLOCKS_PER_SEC);
    } else {
        fprintf(stderr, "Building cache file for %s\n", file);
        mm = (*bind->build)(xx, file);
        (*bind->save_cache)(mm, file);
        fprintf(stderr, "Cache build time %.2fs\n",
                (double) (clock()-t0) / CLOCKS_PER_SEC);
    }

    fprintf(stderr, "%s: %u rows %u cols %lu coeffs\n",
            file, mm->nrows, mm->ncols, mm->ncoeffs);

    abt * src = abinit(xx, mm->ncols);
    abt * dst = abinit(xx, mm->nrows);

    t0 = clock();
    double next = 0.25 * CLOCKS_PER_SEC;
    if (next > tmax * CLOCKS_PER_SEC) { next = tmax * CLOCKS_PER_SEC; }
    for(unsigned int n = 0 ; n < nmax ; n++ ) {
        (*bind->mul)(mm, dst, src, 1);
        double dt = clock() - t0;
        if (dt > next) {
            do { next += 0.25 * CLOCKS_PER_SEC; } while (dt > next);
            if (next > tmax * CLOCKS_PER_SEC) { next = tmax * CLOCKS_PER_SEC; }
            dt /= CLOCKS_PER_SEC;
            printf("%d iterations in %2.fs, %.2f/1, %.2f ns/coeff        \r",
                    n, dt, dt/n, 1.0e9 * dt/n/mm->ncoeffs);
            if (dt > tmax)
                break;
        }
    }
    printf("\n");

    (*bind->clear)(mm);

    param_list_clear(pl);

    return 0;
}
