#define _DARWIN_C_SOURCE
#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>
#include "macros.h"
#include "utils.h"

// This program takes a matrix in the ascii format, ant produces:
//
// - the headerless matrix in binary format.
// - the .cw and .rw files that go with it.
//
// if th

int main(int argc, char * argv[])
{
    param_list pl;

    const char * rwfile = NULL;
    const char * cwfile = NULL;
    const char * mfile = NULL;
    const char * ofile = NULL;

    unsigned int wild =  0;
    int quiet =  0;

    param_list_init(pl);
    argv++,argc--;

    param_list_configure_knob(pl, "--quiet", &quiet);
    for(;argc;) {
        if (param_list_update_cmdline(pl, &argc, &argv)) continue;
        if (wild == 0) {
            mfile = argv[0];
            wild++;
            argv++,argc--;
            continue;
        }
        fprintf(stderr, "unknown option %s\n", argv[0]);
    }

    const char * tmp;
    if ((tmp = param_list_lookup_string(pl, "mfile")) != NULL) {
        mfile = tmp;
    }
    if ((tmp = param_list_lookup_string(pl, "rwfile")) != NULL) {
        rwfile = tmp;
    }
    if ((tmp = param_list_lookup_string(pl, "cwfile")) != NULL) {
        cwfile = tmp;
    }
    if ((tmp = param_list_lookup_string(pl, "ofile")) != NULL) {
        ofile = tmp;
    }

    if (!mfile) {
        fprintf(stderr, "usage ; tbd\n");
        exit(1);
    }

    FILE * mat;
    FILE * bmat;
    if (strcmp(mfile, "-") == 0) { mfile = NULL; }
    if (strcmp(ofile, "-") == 0) { ofile = NULL; }

    if (!ofile) {
        if (!rwfile) fprintf(stderr, "no rwfile given !\n");
        if (!cwfile) fprintf(stderr, "no cwfile given !\n");
    } else {
        char * leakme;
        if (!rwfile) { asprintf(&leakme, "%s.rw", ofile); rwfile = leakme; }
        if (!cwfile) { asprintf(&leakme, "%s.cw", ofile); cwfile = leakme; }
    }

    if (!mfile) {
        mat = stdin;
    } else {
        mat = fopen(mfile, "r");
        ASSERT_ALWAYS(mat != NULL);
    }
    if (!ofile) {
        bmat = stdout;
    } else {
        bmat = fopen(ofile, "w");
    }

    /* We're reading ascii, so the number of rows and cols is readily
     * available -- hurrah !
     */
    unsigned int nr, nc;
    fscanf(mat, "%u %u", &nr, &nc);
    uint32_t * rw, * cw;
    rw = (uint32_t *) malloc(nr * sizeof(uint32_t));
    cw = (uint32_t *) malloc(nc * sizeof(uint32_t));
    memset(rw, 0, nr * sizeof(uint32_t));
    memset(cw, 0, nc * sizeof(uint32_t));

    time_t t0 = time(NULL);

    for(unsigned int i = 0 ; i < nr ; i++) {
        uint32_t w;
        fscanf(mat, "%"SCNu32, &w);
        rw[i] = w;
        fwrite(&w, sizeof(uint32_t), 1, bmat);
        for( ; w-- ; ) {
            uint32_t c;
            fscanf(mat, "%" SCNu32, &c);
            fwrite(&c, sizeof(uint32_t), 1, bmat);
            cw[c]++;
        }
        if ((i+1) % 100000 == 0) {
            fprintf(stderr, "read %zu MB, wrote %zu MB in wct %d s\n",
                    ftell(mat) >> 20,
                    ftell(bmat) >> 20,
                    (int) (time(NULL) - t0));
        }
    }
    if (mfile) fclose(mat);
    if (ofile) fclose(bmat);

    // writing rw and cw file

    if (rwfile) {
        FILE * f = fopen(rwfile, "w");
        fwrite(rw, sizeof(uint32_t), nr, f);
        fclose(f);
    }
    if (cwfile) {
        FILE * f = fopen(cwfile, "w");
        fwrite(cw, sizeof(uint32_t), nc, f);
        fclose(f);
    }

    param_list_clear(pl);
    return 0;
}
