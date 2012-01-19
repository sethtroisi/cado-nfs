#ifndef SPARSE_MAT_H_
#define SPARSE_MAT_H_

#include "purgedfile.h"
#include "dclist.h"

#define USE_TAB 1 // 1 for compact rows...

#define TRACE_COL -1 // 253224 // 231 // put to -1 if not...!
#define TRACE_ROW -1 // 59496 // put to -1 if not...!

/* we bury all ideals of density > 1/BURIED_MAX_DENSITY */
#define BURIED_MAX_DENSITY 2000.0

/* rows correspond to relations, and columns to primes (or prime ideals) */
typedef struct {
  int nrows;
  int ncols;
  int rem_nrows;     /* number of remaining rows */
  int rem_ncols;     /* number of remaining columns */
  int32_t jmin, jmax;    /* we are interested in columns [jmin..jmax[ */
  int32_t **rows;
  int *wt;           /* weight w of column j, if w <= cwmax,
                        else <= 1 for a deleted column
                        (trick: we store -w if w > cwmax) */
  int nburied;     /* the number of buried columns, hence an upper
			bound for wburied[i] */
  int *wburied;     /* wburied[i] counts the estimated weight of buried
			columns */
  unsigned long *ad;
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
#ifndef USE_MARKOWITZ
  dclist *S;         /* S[w] is a doubly-chained list with columns of weight w,
                        for w <= cwmax. For technical reasons (to avoid NULL
                        arguments that we can't modify), S[w] always contains
                        a first cell with value -1, thus the real list is
                        S[w]->next. */
  dclist *A;         /* A[j] points to the unique cell in S[w] containing j,
                        where w = weight(j), for w <= cwmax. A[j]=NULL for
                        w > cwmax. */
#else
  int32_t *MKZQ;         /* priority queue for Markowitz stuff */    
  int32_t *MKZA;         /* MKZA[j] gives u s.t. MKZQ[2*u] = j and
                            MKZQ[2*u+1] is the Markowitz cost of column j,
                            otherwise it is MKZ_INF if the colum is inactive
                            (either too heavy initially or deleted) */
  int wmstmax;
  int mkztype;       /* which type of count */
#endif
  int itermax;       /* used for performing some sampling */
} filter_matrix_t;

#ifdef __cplusplus
extern "C" {
#endif

#ifdef USE_MPI
#define GETJ(mat, j) ((j)-(mat)->jmin)
#else
#define GETJ(mat, j) (j)
#endif

extern void initMat(filter_matrix_t *mat, int32_t jmin, int32_t jmax);
extern void clearMat (filter_matrix_t *mat);
extern void filter_matrix_read_weights(filter_matrix_t *mat, purgedfile_stream_ptr);
extern void fillmat(filter_matrix_t *mat);
extern int filter_matrix_read (filter_matrix_t *mat, purgedfile_stream_ptr, int verbose);

extern void remove_j_from_row(filter_matrix_t *mat, int i, int j);
extern void print_row(filter_matrix_t *mat, int i);

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
extern void freeRj(filter_matrix_t *mat, int j);
extern void remove_i_from_Rj(filter_matrix_t *mat, int i, int j);
extern void add_i_to_Rj(filter_matrix_t *mat, int i, int j);
extern int decrS(int w);
extern int incrS(int w);
extern int weightSum(filter_matrix_t *mat, int i1, int i2);
extern void fillTabWithRowsForGivenj(int32_t *ind, filter_matrix_t *mat, int32_t j);
extern void checkData(filter_matrix_t *mat);
extern void destroyRow(filter_matrix_t *mat, int i);
  extern int buriedRowsWeight (filter_matrix_t *mat, int i1, int i2);

#ifdef __cplusplus
}
#endif

#endif	/* SPARSE_MAT_H_ */
