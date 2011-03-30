#define _POSIX_C_SOURCE 200112L
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include "matmul-mf.h"
#include "mf.h"

/* This is some scotch tape between the balancing layer and the
 * lower-level matmul layer. It also uses some code from mf. This, in
 * effect, bypasses matmul_common_read_stupid_data, which is no longer
 * used.
 */
void mf_prepare_matrix_u32(matmul_ptr mm, matrix_u32_ptr m, const char * file)
{
    struct mf_io_file mf[1];
    struct mf_io_file rw[1];
    struct mf_io_file cw[1];
    memset(mf, 0, sizeof(mf));
    memset(rw, 0, sizeof(rw));
    memset(cw, 0, sizeof(cw));
    struct stat sbuf[1];
    int rc = stat(file, sbuf);
    if (rc < 0) {
        fprintf(stderr, "stat(%s): %s\n", file, strerror(errno));
        exit(1);
    }
    mf->alloc = mf->size = sbuf->st_size / sizeof(uint32_t);
    mf->p = malloc(sbuf->st_size);
    ASSERT_ALWAYS(mf->p);
    FILE * f = fopen(file, "r");
    ASSERT_ALWAYS(f);
    int nread = fread(mf->p, sizeof(uint32_t), mf->size, f);
    if (nread < (int) mf->size) {
        fprintf(stderr, "%s: short read (%d < %zu)\n", file, nread, mf->size);
        exit(1);
    }
    fclose(f);
    matrix_read_pass(mf, NULL, rw, cw, 0, 0, 1);

    memset(m, 0, sizeof(matrix_u32));
    m->mfile = file;
    m->bfile = NULL;
    m->transpose = mm->store_transposed;
    m->size = sbuf->st_size / sizeof(uint32_t);
    m->ntwists = 0;
    m->twist = NULL;
    m->p = mf->p;

    mm->dim[m->transpose] = rw->size;
    mm->dim[!m->transpose] = cw->size;
}

