#define _GNU_SOURCE         /* asprintf */
#define _DARWIN_C_SOURCE    /* for asprintf. _ANSI_SOURCE must be undefined */
#define _POSIX_C_SOURCE 200112L

#include <stdio.h>
#include <errno.h>
#include <limits.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include "bwc_config.h"
#include "info_file.h"
#include "filenames.h"
#include "utils.h"      /* cado_strndup */

#if 0
void read_info_file_ng(matmul_top_data_ptr mmt, const char * filename)
{
    // this new function will be renamed at some point.  note that we
    // look for the balancing file in the current directory.
    DIR * dir;
    dir = opendir(".");
    struct dirent * de;
    for( ; (de = readdir(dir)) != NULL ; ) {
    }
}
// XXX TBC
#endif

void read_info_file(matmul_top_data_ptr mmt, const char * filename)
{
    FILE * f = fopen(filename, "r");
    char line[1024];
    char * rptr;

    char * last_slash = strrchr(filename, '/');

    DIE_ERRNO_DIAG(f == NULL, "fopen", filename);

    rptr = fgets(line, sizeof(line), f);
    ASSERT_ALWAYS(rptr != NULL);
    sscanf(line, "%u %u", & mmt->n[1], & mmt->n[0]);

    rptr = fgets(line, sizeof(line), f);
    ASSERT_ALWAYS(rptr != NULL);

    // nslices[0] is the number of slices in a row -- this is the number of
    // ``vertical slices'', in balance.c parlance.
    unsigned int nslices[2];
    sscanf(line, "%u %u", & nslices[1], & nslices[0]);

    int ok=1;

    for(int d = 0 ; d < 2 ; d++)
        ok = ok && mmt->pi->wr[d]->totalsize == nslices[d];

    if (!ok) {
        /* when displaying an error message which precedes exit(), it's
         * preferrable to have everybody shout. At the very very least all
         * threads must shout. Because as soon as one thread exits, all the
         * other threads are killed immediately. (an alternative could be to
         * serialize on pi->m prior to exiting. But that would be against the
         * rule of never trying to be smart upon fatal errors).
         */
        fprintf(stderr, "Configured split %ux%ux%ux%u does not match "
                "on-disk data %ux%u\n",
                mmt->pi->wr[1]->njobs,
                mmt->pi->wr[0]->njobs,
                mmt->pi->wr[1]->ncores,
                mmt->pi->wr[0]->ncores,
                nslices[1],
                nslices[0]);
        exit(1);
    }

    mmt->fences[0] = malloc((1 + nslices[1])*sizeof(unsigned int));
    mmt->fences[1] = malloc((1 + nslices[0])*sizeof(unsigned int));

    // we use UINT_MAX markers for sanity checking.
    memset(mmt->fences[0], -1, (1+nslices[1]) * sizeof(unsigned int));
    memset(mmt->fences[1], -1, (1+nslices[0]) * sizeof(unsigned int));

    mmt->fences[0][nslices[1]] = mmt->n[1];
    mmt->fences[1][nslices[0]] = mmt->n[0];

    mmt->ncoeffs_total = 0;
    for( ; !feof(f) ; ) {
        unsigned int i,j;
        unsigned int i0, i1, j0, j1, ncoeffs;
        char locfile[80] = { '\0', };
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
        ASSERT_ALWAYS(i < nslices[1]);
        ASSERT_ALWAYS(j < nslices[0]);
        ASSERT_ALWAYS(mmt->fences[0][i] == i0 || mmt->fences[0][i] == UINT_MAX);
        ASSERT_ALWAYS(mmt->fences[1][j] == j0 || mmt->fences[1][j] == UINT_MAX);
        mmt->fences[0][i] = i0;
        mmt->fences[1][j] = j0;
#endif
        if (i != mmt->pi->wr[1]->jrank * mmt->pi->wr[1]->ncores + mmt->pi->wr[1]->trank)
            continue;
        if (j != mmt->pi->wr[0]->jrank * mmt->pi->wr[0]->ncores + mmt->pi->wr[0]->trank)
            continue;

        /* In the info file, the number of rows/cols is given */
        i1 += i0;
        j1 += j0;

        mmt->wr[0]->i0 = i0;
        mmt->wr[0]->i1 = i1;
        mmt->wr[1]->i0 = j0;
        mmt->wr[1]->i1 = j1;
        mmt->ncoeffs = ncoeffs;
        mmt->ncoeffs_total += ncoeffs;
        
        if (strchr(locfile, '/') == NULL && last_slash != NULL) {
            // there is a dirname for the matrix info file, so honour it.
            char * dirname = cado_strndup(filename, last_slash+1-filename);
            int rc = asprintf(&(mmt->locfile), "%s%s", dirname, locfile);
            FATAL_ERROR_CHECK(rc < 0, "out of memory");
            free(dirname);
        } else {
            // keep it as it is.
            mmt->locfile = strdup(locfile);
        }
    }
    fclose(f);
}

