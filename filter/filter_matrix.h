#ifndef SPARSE_MAT_H_
#define SPARSE_MAT_H_

#include "purgedfile.h"
#include "typedefs.h"

typedef struct {
  int32_t id;
  int32_t e;
} ideal_merge_ffs_t;

#define TRACE_COL -1 // 253224 // 231 // put to -1 if not...!
#define TRACE_ROW -1 // 59496 // put to -1 if not...!

#ifndef FOR_DL
#define typerow_t int32_t
#else
#define typerow_t ideal_merge_ffs_t
#endif


/* rows correspond to relations, and columns to primes (or prime ideals) */
typedef struct {
  int nrows;
  int ncols;
  int rem_nrows;     /* number of remaining rows */
  int rem_ncols;     /* number of remaining columns */
  typerow_t **rows;     /* rows[i][k] contains indices of an ideal of row[i] 
                         with 1 <= k <= rows[i][0] */
                        /* FOR_DL: struct containing also the exponent */
  int *wt;           /* weight w of column j, if w <= cwmax,
                        else <= 1 for a deleted column
                        (trick: we store -w if w > cwmax) */
  int nburied;     /* the number of buried columns */
  unsigned long weight;
  int cwmax;         /* bound on weight of j to enter the SWAR structure */
  int rwmax;         /* if a weight(row) > rwmax, kill that row */
  int keep;          /* target for nrows-ncols */
  int mergelevelmax; /* says it */
  int32_t **R;           /* R[j][k] contains the indices of the rows containing
                            the ideal of index j, 0 <= j < ncols,
                            1 <= k <= R[j][0], for weight(j) <= cwmax.
                        R[j][k] = -1 if the corresponding row has been deleted.
                        R[j]=NULL for weight(j) > cwmax. */
  int32_t *MKZQ;         /* priority queue for Markowitz stuff */    
  int32_t *MKZA;         /* MKZA[j] gives u s.t. MKZQ[2*u] = j and
                            MKZQ[2*u+1] is the Markowitz cost of column j,
                            otherwise it is MKZ_INF if the column is inactive
                            (either too heavy initially or deleted) */
  int wmstmax;
  int mkztype;       /* which type of count */
  int itermax;       /* used for performing some sampling */

} filter_matrix_t;

#ifdef __cplusplus
extern "C" {
#endif

extern void initMat(filter_matrix_t *mat);
extern void clearMat (filter_matrix_t *mat);
extern void filter_matrix_read_weights(filter_matrix_t *mat, purgedfile_stream_ptr);
extern void fillmat(filter_matrix_t *mat);
extern void filter_matrix_read (filter_matrix_t *, const char *);

extern void remove_j_from_row(filter_matrix_t *mat, int i, int j);
extern void print_row(filter_matrix_t *mat, int i);

#define isRowNull(mat, i) ((mat)->rows[(i)] == NULL)
#ifdef FOR_DL
#define matLengthRow(mat, i) (mat)->rows[(i)][0].id
#define matCell(mat, i, k) (mat)->rows[(i)][(k)].id
#define setCell(v, j, c) v = (typerow_t) {.id = j, .e = c}
#else
#define matLengthRow(mat, i) (mat)->rows[(i)][0]
#define matCell(mat, i, k) (mat)->rows[(i)][(k)]
#define setCell(v, j, e) v = j
#endif
#define SPARSE_ITERATE(mat, i, k) for((k)=1; (k)<=lengthRow((mat),(i)); (k)++)

extern void freeRj(filter_matrix_t *mat, int j);
extern void remove_i_from_Rj(filter_matrix_t *mat, int i, int j);
extern void add_i_to_Rj(filter_matrix_t *mat, int i, int j);
extern int decrS(int w);
extern int incrS(int w);
extern int weightSum(filter_matrix_t *mat, int i1, int i2, int32_t j);
extern int fillTabWithRowsForGivenj(int32_t *ind, filter_matrix_t *mat, int32_t j);
extern void checkData(filter_matrix_t *mat);
extern void destroyRow(filter_matrix_t *mat, int i);

#ifdef __cplusplus
}
#endif

#endif	/* SPARSE_MAT_H_ */
