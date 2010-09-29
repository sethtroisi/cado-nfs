#ifndef CADO_LINALG_READMAT_H_
#define CADO_LINALG_READMAT_H_

#include <gmp.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

// data is stored in a huge table:
//   [ nb_coeff c0 c1 ... ck nb_coeff c0 c1 ..... ]
struct filter_matrix_s {
    unsigned int nrows;
    unsigned int ncols;
    unsigned int *data;

    unsigned int alloc;
    unsigned int wt;		// number of coeffs
};

typedef struct filter_matrix_s filter_matrix_t[1];

/* Reads the whole matrix */
extern void readmat(FILE *file, filter_matrix_t mat, int compact);

static inline void readmat_with_ab(FILE *file, filter_matrix_t mat)
{ readmat(file, mat, 0); }

extern void read_matrix_header(FILE *file, filter_matrix_t mat);

/* Reads only one row of the matrix. Reallocates the filter_matrix_t
 * structure if needed. The data is put at the location pointed to by
 * dst. A pointer to the end of the written area is returned.
 *
 * BEWARE: dst must be a pointer within mat->data. However, besides that,
 * the caller is free to place the data at an arbitrary place within that
 * buffer. An example use could be:
 *
 * dst = mat->data;     // beginning of buffer
 * dst0 = dst;          // save pointer.
 * dst = read_matrix_row(f,m,dst,flag);       // read first row
 * // note that the buffer is auto-extended to fit the matrix row.
 * // dst now points after the first row
 * // do something with [dst0..dst[
 * dst0 = dst;
 * dst = read_matrix_row(f,m,dst,flag);       // read second row
 * // now say the first two rows are completely unneeded.
 * // It is ok to trash the corresponding data.
 * dst = mat->data
 *
 */
extern unsigned int * read_matrix_row(FILE * file, filter_matrix_t mat, unsigned int * dst, int compact);

extern void filter_matrix_init(filter_matrix_t dst);
extern void filter_matrix_clear(filter_matrix_t dst);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_LINALG_READMAT_H_ */
