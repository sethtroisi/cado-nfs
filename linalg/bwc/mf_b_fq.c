#define _DARWIN_C_SOURCE
#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>
#include "macros.h"
#include "utils.h"

// This program takes a matrix in the binary format, and keeps it. It
// only writes down the accompanying rw and cw files. No main output is
// produced -- most probably this is only intended as a recovery tool.

unsigned int nc = 0;
uint32_t * cw = NULL;

static inline void expand_cw(uint32_t nnc)
{
    cw = realloc(cw, nnc * sizeof(uint32_t));
    memset(cw + nc, 0, (nnc - nc) * sizeof(uint32_t));
    nc = nnc;
}


int main(int argc, char * argv[])
{
    param_list pl;

    const char * rwfile = NULL;
    const char * cwfile = NULL;
    const char * mfile = NULL;

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

    if (!mfile) {
        fprintf(stderr, "usage ; tbd\n");
        exit(1);
    }

    FILE * mat;
    if (strcmp(mfile, "-") == 0) { mfile = NULL; }

    if (!mfile) {
        if (!rwfile) fprintf(stderr, "no rwfile given !\n");
        if (!cwfile) fprintf(stderr, "no cwfile given !\n");
    } else {
        char * leakme;
        if (!rwfile) { asprintf(&leakme, "%s.rw", mfile); rwfile = leakme; }
        if (!cwfile) { asprintf(&leakme, "%s.cw", mfile); cwfile = leakme; }
    }

    if (!mfile) {
        mat = stdin;
    } else {
        mat = fopen(mfile, "r");
        ASSERT_ALWAYS(mat != NULL);
    }

    /* We're reading ascii, so the number of rows and cols is readily
     * available -- hurrah !
     */
    expand_cw(1000 * 1000);

    time_t t0 = time(NULL);
    unsigned int i = 0 ;

    FILE * frw = NULL;
    if (rwfile) frw = fopen(rwfile, "w");

    uint32_t lastcol = 0;
    for(; ; i++) {
        uint32_t w;
        fread(&w, sizeof(uint32_t), 1, mat);
        if (feof(mat)) break;
        if (frw) fwrite(&w, sizeof(uint32_t), 1, frw);
        for( ; w-- ; ) {
            uint32_t c;
            fread(&c, sizeof(uint32_t), 1, mat);
            if (c >= lastcol) lastcol = c;
            if (c >= nc) expand_cw(c + nc / 2);
            cw[c]++;
        }
        if ((i+1) % 100000 == 0) {
            fprintf(stderr, "read %zu MB in wct %d s\n",
                    ftell(mat) >> 20,
                    (int) (time(NULL) - t0));
        }
    }
    if (rwfile) fclose(frw);
    if (mfile) fclose(mat);
    nc = lastcol + 1;

    if (cwfile) {
        FILE * f = fopen(cwfile, "w");
        fwrite(cw, sizeof(uint32_t), nc, f);
        fclose(f);
    }

    param_list_clear(pl);
    return 0;
}
