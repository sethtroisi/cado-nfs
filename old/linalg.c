#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "cado.h"
#include "readmat.h"

#include "gauss.h"

typedef struct {
    unsigned int nrows;
    unsigned int ncols;
    mp_limb_t *data;
    unsigned int limbs_per_row;
    unsigned int limbs_per_col;
} dense_mat_t[1];



void setcoeff(dense_mat_t Dmat, unsigned int i, unsigned int j)
{
    int index;
    int j0, j1;
    mp_limb_t mask;
    index = Dmat->limbs_per_row * i;
    j0 = j / GMP_NUMB_BITS;
    index += j0;
    j1 = j - (j0 * GMP_NUMB_BITS);
    mask = 1UL << j1;
    Dmat->data[index] |= mask;
}

void sparse2densetranspose(dense_mat_t Dmat, sparse_mat_t Smat)
{
    unsigned int i, j, *ptr;

    Dmat->nrows = Smat->ncols;
    Dmat->ncols = Smat->nrows;
    Dmat->limbs_per_row = (Dmat->ncols / GMP_NUMB_BITS) + 1;
    Dmat->limbs_per_col = (Dmat->nrows / GMP_NUMB_BITS) + 1;

    Dmat->data =
	(mp_limb_t *) malloc(Dmat->limbs_per_row * Dmat->nrows *
			     sizeof(mp_limb_t));
    for (i = 0; i < Dmat->limbs_per_row * Dmat->nrows; ++i)
	Dmat->data[i] = 0UL;

    ptr = Smat->data;
    for (i = 0; i < Smat->nrows; ++i) {
	unsigned int nc = *ptr++;
	for (j = 0; j < nc; ++j) {
	    unsigned int coeff = *ptr++;
	    setcoeff(Dmat, coeff, i);
	}
    }
}

void sparse2dense(dense_mat_t Dmat, sparse_mat_t Smat)
{
    unsigned int i, j, *ptr;

    Dmat->nrows = Smat->nrows;
    Dmat->ncols = Smat->ncols;
    Dmat->limbs_per_row = (Dmat->ncols / GMP_NUMB_BITS) + 1;
    Dmat->limbs_per_col = (Dmat->nrows / GMP_NUMB_BITS) + 1;

    Dmat->data =
	(mp_limb_t *) malloc(Dmat->limbs_per_row * Dmat->nrows *
			     sizeof(mp_limb_t));
    ASSERT(Dmat->data != NULL);
    for (i = 0; i < Dmat->limbs_per_row * Dmat->nrows; ++i)
	Dmat->data[i] = 0UL;

    ptr = Smat->data;
    for (i = 0; i < Smat->nrows; ++i) {
	unsigned int nc = *ptr++;
	for (j = 0; j < nc; ++j) {
	    unsigned int coeff = *ptr++;
	    setcoeff(Dmat, i, coeff);
	}
    }
}


int main(int argc, char **argv)
{
    FILE *file;
    sparse_mat_t mat;
    dense_mat_t dmat;
    mp_limb_t **ker;
    int compact = 0;
    unsigned int i, j, dim;

    if (argc != 3) {
	fprintf(stderr, "usage: %s [filename] [compact]\n", argv[0]);
	exit(1);
    }
    file = fopen(argv[1], "r");
    ASSERT(file != NULL);
    compact = atoi(argv[2]);

    sparse_mat_init(mat);
    readmat(file, mat, compact);
    fprintf(stderr, "have read matrix: nrows = %d, ncols = %d\n",
            mat->nrows, mat->ncols);
    fprintf(stderr, "average wt of rows is %f\n",
	    mat->wt / (double) mat->nrows);

    sparse2dense(dmat, mat);

    ker = (mp_limb_t **) malloc(dmat->nrows * sizeof(mp_limb_t *));
#if SAVE_KERNEL_MEMORY == 0
    for (i = 0; i < dmat->nrows; ++i)
	ker[i] = (mp_limb_t *) malloc(dmat->limbs_per_col * sizeof(mp_limb_t));
#endif

    dim = kernel(dmat->data, ker, dmat->nrows, dmat->ncols,
		 dmat->limbs_per_row, dmat->limbs_per_col);

    fprintf(stderr, "dim of kernel is %d\n", dim);

    fprintf(stderr, "printing 128 first vectors of kernel...\n");

    for (i = 0; i < (dim > 128 ? 128 : dim); ++i) {
	for (j = 0; j < dmat->limbs_per_col; j++) {
	    printf("%lx ", ker[i][j]);
	}
	printf("\n");
    }

    sparse_mat_clear(mat);
    fclose(file);

    return 0;
}
