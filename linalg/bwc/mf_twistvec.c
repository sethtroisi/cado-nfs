#include "cado.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include "balancing.h"
#include "utils.h"

void usage(const char * progname)
{
    fprintf(stderr, "Usage: %s [--twist|--untwist] [--nullspace left|right] [--expand|--truncate] [-x] <bfile> <vfile> [--out <vfile>]\n"
            "--twist\tCompute forward permutation of vfile from bfile (default)\n"
            "--untwist\tCompute backward permutation of vfile from bfile\n"
            "--expand\tWhen twisting, start from an input file with trailing padding block stripped out\n"
            "--truncate\tWhen untwisting, discard the trailing padding block of the result (output file is then shorter)\n"
            "--out\tSpecify output file (defaults to stdout)\n"
            "-x\tAssume files correspond to X vectors (different format)\n"
            "\n"
            "Note: specifying --nullspace does have an impact on the permutation used, in the case of the shuffled product\n", progname);
    exit(1);
}

int main(int argc, char * argv[])
{
    balancing bal;
    param_list pl;
    param_list_init(pl);
    const char * progname = argv[0];
    argv++, argc--;
    int nullspace_left = 0;
    int untwist = 0;
    int twist = 0;
    int truncate = 0;
    int expand = 0;
    int xmode = 0;
    param_list_configure_knob(pl, "--untwist", &untwist);
    param_list_configure_knob(pl, "--twist", &twist);
    param_list_configure_knob(pl, "--truncate", &truncate);
    param_list_configure_knob(pl, "--expand", &expand);
    param_list_configure_knob(pl, "-x", &xmode);
    unsigned int wild = 0;
    const char * bname = NULL;
    const char * vname = NULL;
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) continue;
        if (argv[0][0] != '-' && wild == 0) {
            bname = argv[0];
            wild++;
            argv++, argc--;
            continue;
        }
        if (argv[0][0] != '-' && wild == 1) {
            vname = argv[0];
            wild++;
            argv++, argc--;
            continue;
        }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        usage(progname);
    }
    ASSERT_ALWAYS(twist + untwist <= 1);
    if (!twist && !untwist) twist=1;
    // truncate only makes sense for untwisting
    ASSERT_ALWAYS(!truncate || untwist);
    // expand only makes sense for twisting
    ASSERT_ALWAYS(!expand || !untwist);
    {
        const char * tmp;
        if ((tmp = param_list_lookup_string(pl, "nullspace")) && strcmp(tmp, "left") == 0) {
            nullspace_left = 1;
        }
    }
    ASSERT_ALWAYS(bname);
    balancing_init(bal);
    balancing_read(bal, bname);
    unsigned long unpadded = nullspace_left ? bal->h->nrows : bal->h->ncols;
    unsigned long padded = nullspace_left ? bal->trows : bal->tcols;

    // in reality, both are the same permutation as long as only conjugated
    // permutations are implemented.
    uint32_t * xc = bal->colperm;
    uint32_t * xr = bal->rowperm;
    if (!(bal->h->flags & FLAG_REPLICATE)) {
        ASSERT_ALWAYS(bal->h->flags & FLAG_COLPERM);
        ASSERT_ALWAYS(bal->h->flags & FLAG_ROWPERM);
    } else {
        if (!xc) xc = xr;
        if (!xr) xr = xc;
    }
    ASSERT_ALWAYS(xc);
    ASSERT_ALWAYS(xr);
    uint32_t * perm = nullspace_left ? xr:xc;
    uint32_t * iperm = NULL;
    if (!untwist) {
        iperm = malloc(padded * sizeof(uint32_t));
        for(unsigned long k = 0 ; k < padded ; k++) {
            iperm[perm[k]]=k;
        }
    }

    uint32_t * pr_perm = NULL;
    if (nullspace_left && (bal->h->flags & FLAG_SHUFFLED_MUL)) {
        /* Presently, balancing_workhorse only permutes rows. This has the
         * consequence that only nullspace=left implies a permutation on the
         * vectors */
        pr_perm = malloc(padded * sizeof(uint32_t));

        uint32_t nh = bal->h->nh;
        uint32_t nv = bal->h->nv;
        ASSERT_ALWAYS(bal->trows % (nh * nv) == 0);
        uint32_t elem = bal->trows / (nh * nv);
        uint32_t ix = 0;
        uint32_t iy = 0;
        for(uint32_t i = 0 ; i < nh ; i++) {
            for(uint32_t j = 0 ; j < nv ; j++) {
                ix = (i * nv + j) * elem;
                iy = (j * nh + i) * elem;
                for(uint32_t k = 0 ; k < elem ; k++) {
                    if (untwist) {
                        pr_perm[ix+k]=iy+k;
                    } else {
                        pr_perm[iy+k]=ix+k;
                    }
                }
            }
        }
    }

    ASSERT_ALWAYS(vname);
    FILE * v;
    v = fopen(vname, "r");
    if (!v) {
        perror(vname);
        exit(1);
    }
    const char * outfilename = param_list_lookup_string(pl, "out");
    FILE * o = stdout;
    if (outfilename) {
        o = fopen(outfilename, "w");
        DIE_ERRNO_DIAG(!o, "fopen", outfilename);
    }
    if (!xmode) {
        struct stat sbuf[1];
        if (stat(vname, sbuf) < 0) {
            perror(vname);
            exit(1);
        }
        size_t input_items = expand ? unpadded : padded;
        size_t output_items = truncate ? unpadded : padded;
        ASSERT_ALWAYS(sbuf->st_size % input_items == 0);
        size_t chunk = sbuf->st_size / input_items;
        char * area = malloc(padded * chunk);
        memset(area, 0, padded * chunk);
        char * tmp = malloc(chunk);
        memset(tmp, 0, chunk);

        for(unsigned long k = 0 ; k < input_items ; k++) {
            if (fread(tmp, 1, chunk, v) < chunk) {
                fprintf(stderr, "short read\n");
                exit(1);
            }
            if (pr_perm) {
                if (untwist) {
                    memcpy(area + perm[pr_perm[k]] * chunk, tmp, chunk);
                } else {
                    memcpy(area + pr_perm[iperm[k]] * chunk, tmp, chunk);
                }
            } else {
                uint32_t * p = untwist ? perm : iperm;
                memcpy(area + p[k] * chunk, tmp, chunk);
            }
        }
        fwrite(area, chunk, output_items, o);
        free(tmp);
        free(area);
    } else {
        unsigned int nx;
        int rc;
        rc = fscanf(v, "%u", &nx); ASSERT_ALWAYS(rc == 1);
        // we have m*nx items. We don't want to bother reading m, so we just
        // apply the permutation to the complete file.
        fprintf(o, "%u\n", nx);
        unsigned int nr = 0;
        for( ; ; nr++) {
            unsigned int col;
            rc = fscanf(v,"%u", &col);
            if (rc <= 0) break;
            if (pr_perm) {
                if (untwist) {
                    col = perm[pr_perm[col]];
                } else {
                    col = pr_perm[iperm[col]];
                }
            } else {
                uint32_t * p = untwist ? perm : iperm;
                col = p[col];
            }
            fprintf(o, "%u\n", col);
        }
        ASSERT_ALWAYS(nr % nx == 0);
    }
    if (outfilename)
        fclose(o);
    fclose(v);



    free(iperm);
    free(pr_perm);
    return 0;
}

