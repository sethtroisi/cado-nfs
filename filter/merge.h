#ifndef MERGE_H_
#define MERGE_H_

#define USE_TAB 1 // 1 for compact rows...

/* int32_t is defined in sparse.h */

#include "dclist.h"

/* rows correspond to relations, and columns to primes (or prime ideals) */
typedef struct {
  int nrows;
  int ncols;
  int rem_nrows;     /* number of remaining rows */
  int rem_ncols;     /* number of remaining columns */
#if USE_TAB == 0
  rel_t *data;
#else
  int32_t **rows;
#endif
  int *wt;           /* weight w of column j, if w <= cwmax,
                        else <= 1 for a deleted column
                        (trick: we store -w if w > cwmax) */
  unsigned long *ad;
  unsigned long weight;
  int cwmax;         /* bound on weight of j to enter the SWAR structure */
  int rwmax;         /* if a weight(row) > rwmax, kill that row */
  int keep;          /* target for nrows-ncols */
  int mergelevelmax; /* says it */
  dclist *S;         /* S[w] is a doubly-chained list with columns of weight w,
                        for w <= cwmax. For technical reasons (to avoid NULL
                        arguments that we can't modify), S[w] always contains
                        a first cell with value -1, thus the real list is
                        S[w]->next. */
  dclist *A;         /* A[j] points to the unique cell in S[w] containing j,
                        where w = weight(j), for w <= cwmax. A[j]=NULL for
                        w > cwmax. */
  int32_t **R;           /* R[j][k] contains the rows of the non-empty elements
                        of column j, 0 <= j < ncols, 1 <= k <= R[j][0], for
                        weight(j) <= cwmax.
                        R[j][k] = -1 if the corresponding row has been deleted.
                        R[j]=NULL for weight(j) > cwmax. */
} filter_matrix_t;

#if USE_TAB == 0
#define isRowNull(mat, i) ((mat)->data[(i)].val == NULL)
#define lengthRow(mat, i) (mat)->data[(i)].len
#define cell(mat, i, k) (mat)->data[(i)].val[(k)]
#define SPARSE_ITERATE(mat, i, k) for((k)=0; (k)<lengthRow((mat),(i)); (k)++)
#else
#define isRowNull(mat, i) ((mat)->rows[(i)] == NULL)
#define lengthRow(mat, i) (mat)->rows[(i)][0]
#define cell(mat, i, k) (mat)->rows[(i)][(k)]
#define SPARSE_ITERATE(mat, i, k) for((k)=1; (k)<=lengthRow((mat),(i)); (k)++)
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern void report1(FILE *outfile, int32_t i);
extern void removeCellSWAR(filter_matrix_t *mat, int i, int32_t j);

#ifdef USE_MPI
extern void mpi_send_inactive_cols(int i);
extern void mpi_add_rows(filter_matrix_t *mat, int m, int32_t j, int32_t *ind);
extern void mpi_load_rows_for_j(filter_matrix_t *mat, int m, int32_t j);
#endif

#ifdef __cplusplus
}
#endif

#endif	/* MERGE_H_ */
