#ifndef SPARSE_MAT_H_
#define SPARSE_MAT_H_

#include "purgedfile.h"
#include "typedefs.h"

#define TRACE_COL -1 // 253224 // 231 // put to -1 if not...!
#define TRACE_ROW -1 // 59496 // put to -1 if not...!

#ifndef FOR_DL
#define typerow_t index_t
#define cmp_typerow_t cmp_index
#else
#define typerow_t ideal_merge_t
#define cmp_typerow_t cmp_ideal_merge
#endif


/* rows correspond to relations, and columns to primes (or prime ideals) */
typedef struct {
  uint64_t nrows;
  uint64_t ncols;
  uint64_t rem_nrows;  /* number of remaining rows */
  uint64_t rem_ncols;  /* number of remaining columns */
  typerow_t **rows;    /* rows[i][k] contains indices of an ideal of row[i] 
                          with 1 <= k <= rows[i][0] */
                       /* FOR_DL: struct containing also the exponent */
  int32_t *wt;         /* weight w of column j, if w <= cwmax,
                          else <= 0 for a deleted column
                          (trick: we store -w if w > cwmax) */
                       /* 32 bits is sufficient as we only want precise weight
                          for column of low weight. If the weight exceed 2^32-1,
                          we saturate */
  uint64_t nburied;    /* the number of buried columns */
  uint64_t weight;     /* number of non-zero coefficients in the active part */
  uint64_t tot_weight; /* Initial total number of non-zero coefficients */
  int cwmax;           /* bound on weight of j to enter the SWAR structure */
  int64_t keep;        /* target for nrows-ncols */
  int mergelevelmax;   /* says it */
  index_t **R;         /* R[j][k] contains the indices of the rows containing
                          the ideal of index j, 0 <= j < ncols,
                          1 <= k <= R[j][0], for weight(j) <= cwmax.
                          R[j][k] = UMAX(index_t) if the corresponding row has
                                                                  been deleted.
                        R[j]=NULL for weight(j) > cwmax. */
  int32_t *MKZQ;         /* priority queue for Markowitz stuff */    
  index_t *MKZA;         /* MKZA[j] gives u s.t. MKZQ[2*u] = j and
                            MKZQ[2*u+1] is the Markowitz cost of column j,
                            otherwise it is MKZ_INF if the column is inactive
                            (either too heavy initially or deleted) */
  int wmstmax;
  int mkztype;       /* which type of count */
} filter_matrix_t;

#ifdef __cplusplus
extern "C" {
#endif

#define compute_WN(mat) ((mat)->rem_nrows * (mat)->weight)
#define compute_WoverN(mat) (((double)(mat)->weight)/((double)(mat)->rem_nrows))

extern void initMat(filter_matrix_t *, int, uint32_t, uint32_t);
extern void clearMat (filter_matrix_t *mat);
extern void fillmat(filter_matrix_t *mat);
extern void filter_matrix_read (filter_matrix_t *, const char *);

extern void remove_j_from_row(filter_matrix_t *mat, int i, int j);
extern void print_row(filter_matrix_t *mat, int i);

#define isRowNull(mat, i) ((mat)->rows[(i)] == NULL)
#ifdef FOR_DL
#define matLengthRow(mat, i) (mat)->rows[(i)][0].id
#define matCell(mat, i, k) (mat)->rows[(i)][(k)].id
#define setCell(cell, v, c) cell = (ideal_merge_t) {.id = v, .e = c}
#else
#define matLengthRow(mat, i) (mat)->rows[(i)][0]
#define matCell(mat, i, k) (mat)->rows[(i)][(k)]
#define setCell(cell, v, c) cell = v
#endif
#define SPARSE_ITERATE(mat, i, k) for((k)=1; (k)<=lengthRow((mat),(i)); (k)++)

extern void freeRj(filter_matrix_t *mat, int j);
extern void remove_i_from_Rj(filter_matrix_t *mat, int i, int j);
extern void add_i_to_Rj(filter_matrix_t *mat, int i, int j);
extern int decrS(int w);
extern int incrS(int w);
extern int weightSum(filter_matrix_t *mat, int i1, int i2, int32_t j);
extern int fillTabWithRowsForGivenj(int32_t *ind, filter_matrix_t *mat, int32_t j);
extern void destroyRow(filter_matrix_t *mat, int i);

#ifdef __cplusplus
}
#endif

#endif	/* SPARSE_MAT_H_ */
