#include <cstdio>
#include <cerrno>

#include <string>

#include "info_file.h"
#include "matmul_top.h"
#include "intersections.hpp"

using namespace std;

#define ASSIGN_ONCE(data, value, counter) do {				\
        if ((counter) == 0) {						\
            (data) = (value);						\
        } else {							\
            ASSERT_ALWAYS((data) == (value));				\
        }								\
    } while (0)

#ifdef  CONJUGATED_PERMUTATIONS
void mmt_wiring_alloc(matmul_top_data_ptr mmt, unsigned int d)
{
    struct mmt_wiring * mw = mmt->wr[d];
    struct pi_wiring * pw = mmt->pi->wr[d];
    unsigned int n = mmt->wr[d]->n;

    unsigned int len = pw->ncores * pw->njobs;

    mw->fences = (unsigned int *) malloc((1+len) * sizeof(unsigned int));
    memset(mw->fences, 0, (1+len) * sizeof(unsigned int));
    mw->fences[pw->ncores * pw->njobs] = n;
}

void mmt_wiring_intersect(matmul_top_data_ptr mmt, unsigned int d)
{
    struct mmt_wiring * mw = mmt->wr[d];
    struct mmt_wiring * ow = mmt->wr[!d];
    unsigned int n = mmt->wr[d]->n;

    std::vector<isect_info> vv;
    vv = intersect(mw->fences, ow->i0, std::min(ow->i1, n));

    mw->xlen = vv.size();
    mw->x = (struct isect_info *) malloc((mw->xlen) * sizeof(struct isect_info));
    std::copy(vv.begin(), vv.end(), mw->x);
}
#endif

void read_info_file(matmul_top_data_ptr mmt, const char * filename)
{
    // we read the info file immediately in order to do some checking.
    string base(filename);
    base += ".info";

    FILE * f = fopen(base.c_str(), "r");
    char line[1024];
    char * rptr;

    if (f == NULL) {
        fprintf(stderr, "fopen(%s): %s\n", base.c_str(), strerror(errno));
        exit(1);
    }

    rptr = fgets(line, sizeof(line), f);
    ASSERT_ALWAYS(rptr != NULL);
    sscanf(line, "%u %u", & mmt->wr[0]->n, & mmt->wr[1]->n);

    rptr = fgets(line, sizeof(line), f);
    ASSERT_ALWAYS(rptr != NULL);
    sscanf(line, "%u %u", & mmt->pi->wr[0]->nslices, & mmt->pi->wr[1]->nslices);

    int ok=1;

    for(int d = 0 ; d < 2 ; d++)
        ok = ok && mmt->pi->wr[d]->nslices ==
            mmt->pi->wr[d]->njobs * mmt->pi->wr[d]->ncores;

    if (!ok) {
        if (mmt->pi->m->jrank == 0 && mmt->pi->m->trank == 0) {
            fprintf(stderr, "Configured split %ux%ux%ux%u does not match "
                    "on-disk data %ux%u\n",
                    mmt->pi->wr[1]->njobs,
                    mmt->pi->wr[0]->njobs,
                    mmt->pi->wr[1]->ncores,
                    mmt->pi->wr[0]->ncores,
                    mmt->pi->wr[1]->nslices,
                    mmt->pi->wr[0]->nslices);
        }
        exit(1);
    }

#ifdef  CONJUGATED_PERMUTATIONS
    for(int d = 0 ; d < 2 ; d++)
        mmt_wiring_alloc(mmt, d);
#endif

    mmt->ncoeffs_total = 0;
    for( ; !feof(f) ; ) {
        unsigned int i,j;
        unsigned int i0, i1, j0, j1, ncoeffs;
        char locfile[80];
        if (fgets(line, sizeof(line), f) == NULL) {
            break;
        }
        if (line[0] == '#') {
            continue;
        }
        int rc = sscanf(line, "%u %u %u %u %u %u %u %80s\n",
                &i, &j, &i0, &j0, &i1, &j1, &ncoeffs, locfile);
        if (rc != 8) {
            fprintf(stderr, "Could not parse input line:\n%s", line);
            continue;
        }
#ifdef  CONJUGATED_PERMUTATIONS
        ASSIGN_ONCE(mmt->wr[0]->fences[i], i0, j);
        ASSIGN_ONCE(mmt->wr[1]->fences[j], j0, i);
#endif
        if (i != mmt->pi->wr[0]->jrank * mmt->pi->wr[0]->ncores + mmt->pi->wr[0]->trank)
            continue;
        if (j != mmt->pi->wr[1]->jrank * mmt->pi->wr[1]->ncores + mmt->pi->wr[1]->trank)
            continue;
        mmt->wr[0]->i0 = i0;
        mmt->wr[0]->i1 = i1;
        mmt->wr[1]->i0 = j0;
        mmt->wr[1]->i1 = j1;
        mmt->ncoeffs = ncoeffs;
        mmt->ncoeffs_total += ncoeffs;
        mmt->filename = strdup(locfile);
    }
    fclose(f);

#ifdef  CONJUGATED_PERMUTATIONS
    for(int d = 0 ; d < 2 ; d++)
        mmt_wiring_intersect(mmt, d);
#endif

    if (!mmt->filename) {
        fprintf(stderr, "No file defined for chunk (%u,%u,%u,%u)\n",
                mmt->pi->wr[0]->jrank, mmt->pi->wr[1]->jrank,
                mmt->pi->wr[0]->trank, mmt->pi->wr[1]->trank);
        exit(1);
    }
}

