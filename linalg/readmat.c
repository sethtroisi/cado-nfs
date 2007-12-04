#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <limits.h>

#include "readmat.h"
#include "manu.h"

unsigned int next_alloc_size(unsigned int s, unsigned int needed)
{
    /* most probably only one spin */
    for( ; s < needed ; ) {
        s = MAX(s + 1024, s + s / 4);
    }
    return s;
}

unsigned int *read_matrix_row(FILE * file, sparse_mat_t mat,
			      unsigned int *dst, int compact)
{
    unsigned int nc;
    unsigned int j;
    int ret;


    long a;
    unsigned long b;

    if (compact == 0) {
	ret = fscanf(file, "%ld %lu", &a, &b);
	BUG_ON_MSG(ret != 2, "Could not read a and b");

        /* FIXME -- having a and b stored in the sparse mat struct simply
         * won't scale ! */

        BUG_ON_MSG(a > (long) INT_MAX, "Time to enlarge !");
        BUG_ON_MSG(a < (long) INT_MIN, "Time to enlarge !");
        BUG_ON_MSG(b > (unsigned long) INT_MAX, "Time to enlarge !");
    }

    ret = fscanf(file, "%d", &nc);
    BUG_ON_MSG(ret != 1, "Could not read ncoeffs");

    unsigned int offset = dst - mat->data;
    unsigned int needed = offset + nc + 1;

    if (compact == 0)
        needed += 2;

    if (needed > mat->alloc) {
	unsigned int next = next_alloc_size(mat->alloc, needed);
	mat->data = (unsigned int *) realloc(mat->data,
					     next * sizeof(unsigned int));
	mat->alloc = next;
	dst = mat->data + offset;
    }

    if (compact == 0) {
        *dst++ = (int) a;
        *dst++ = (int) b;
    }

    *dst++ = nc;
    for (j = 0; j < nc; ++j) {
	unsigned int x;
	ret = fscanf(file, "%u", &x);
	BUG_ON_MSG(ret != 1, "Could not read coefficient");
	*dst++ = x;
    }
    mat->wt += nc;
    return dst;
}

void read_matrix_header(FILE * file, sparse_mat_t mat)
{
    int ret;

    ret = fscanf(file, "%d %d", &(mat->nrows), &(mat->ncols));
    BUG_ON_MSG(ret != 2, "Could not read nrows/ncols");
    mat->alloc = mat->nrows * 10;
    mat->data = (unsigned int *) malloc(mat->alloc * sizeof(unsigned int));
    mat->wt = 0;
}

void readmat(FILE * file, sparse_mat_t mat, int compact)
{
    int i;
    unsigned int *dst;
    read_matrix_header(file, mat);
    dst = mat->data;
    for (i = 0; i < mat->nrows; ++i) {
	dst = read_matrix_row(file, mat, dst, compact);
    }
}

void sparse_mat_init(sparse_mat_t dst)
{
    memset(dst, 0, sizeof(sparse_mat_t));
}

void sparse_mat_clear(sparse_mat_t dst)
{
    free(dst->data);
    memset(dst, 0, sizeof(sparse_mat_t));
}

/* vim: set sw=4 sta et: */
