#ifndef READMAT_H_
#define READMAT_H_

#include <gmp.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

// data is stored in a huge table:
//   [ nb_coeff c0 c1 ... ck nb_coeff c0 c1 ..... ]
typedef struct {
    unsigned int nrows;
    unsigned int ncols;
    unsigned int *data;

    unsigned int alloc;
    unsigned int wt;		// number of coeffs
} sparse_mat_t[1];

/* Reads the whole matrix */
extern void readmat(FILE *file, sparse_mat_t mat, int compact);

extern void read_matrix_header(FILE *file, sparse_mat_t mat);

/* Reads only one row of the matrix. Reallocates the sparse_mat_t
 * structure if needed. The data is put at the location pointed to by
 * dst. The position of the next matrix row is returned. */
extern unsigned int * read_matrix_row(FILE * file, sparse_mat_t mat, unsigned int * dst, int compact);

extern void sparse_mat_init(sparse_mat_t dst);
extern void sparse_mat_clear(sparse_mat_t dst);

#ifdef __cplusplus
}
#endif

#endif	/* READMAT_H_ */
