#include "cado.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "balancing.h"
#include "utils.h"

int main(int argc, char * argv[])
{
    balancing bal;
    param_list pl;
    param_list_init(pl);
    argv++, argc--;
    // specifying nullspace-left or nullspace-right has no impact here,
    // as we have only one permutation to choose from (we're still
    // forcing the conjugated permutations case). However, future
    // evolutions might require this.
    int nullspace_left = 0;
    int untwist = 0;
    int truncate = 0;
    int expand = 0;
    param_list_configure_knob(pl, "--untwist", &untwist);
    param_list_configure_knob(pl, "--truncate", &truncate);
    param_list_configure_knob(pl, "--expand", &expand);
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
    }
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
    ASSERT_ALWAYS(vname);
    balancing_init(bal);
    balancing_read(bal, bname);
    FILE * v;
    v = fopen(vname, "r");
    if (!v) {
        perror(vname);
        exit(1);
    }
    struct stat sbuf[1];
    if (stat(vname, sbuf) < 0) {
        perror(vname);
        exit(1);
    }
    unsigned long unpadded = nullspace_left ? bal->h->nrows : bal->h->ncols;
    unsigned long padded = nullspace_left ? bal->trows : bal->tcols;
    size_t input_items = expand ? unpadded : padded;
    size_t output_items = truncate ? unpadded : padded;
    ASSERT_ALWAYS(sbuf->st_size % input_items == 0);
    size_t chunk = sbuf->st_size / input_items;
    char * area = malloc(padded * chunk);
    memset(area, 0, padded * chunk);
    char * tmp = malloc(chunk);
    memset(tmp, 0, chunk);

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
                    pr_perm[ix+k]=iy+k;
                }
            }
        }
    }
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

    uint32_t * p = untwist ? perm : iperm;
    for(unsigned long k = 0 ; k < input_items ; k++) {
        if (fread(tmp, 1, chunk, v) < chunk) {
            fprintf(stderr, "short read\n");
            exit(1);
        }
        if (pr_perm) {
            memcpy(area + p[pr_perm[k]] * chunk, tmp, chunk);
        } else {
            memcpy(area + p[k] * chunk, tmp, chunk);
        }
    }
    fclose(v);
    free(tmp);
    {
        const char * tmp;
        if ((tmp = param_list_lookup_string(pl, "out")) != NULL) {
            FILE * o;
            o= fopen(tmp, "w");
            fwrite(area, chunk, output_items, o);
            fclose(o);
        } else {
            fwrite(area, chunk, output_items, stdout);
        }
    }
    free(area);
    free(iperm);
    free(pr_perm);
    return 0;
}

