#ifndef SPARSE_MAT_H_
#define SPARSE_MAT_H_

#include "dclist.h"

#define USE_TAB 1 // 1 for compact rows...

#define TRACE_COL -1 // 253224 // 231 // put to -1 if not...!
#define TRACE_ROW -1 // 59496 // put to -1 if not...!

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
  int *wburried;     /* wburried[i] counts the estimated weight of burried
			columns */
  int nburried;     /* the number of burried columns, hence an upper
			bound for wburried[i] */
  unsigned long *ad;
  unsigned long weight;
  int cwmax;         /* bound on weight of j to enter the SWAR structure */
  int rwmax;         /* if a weight(row) > rwmax, kill that row */
  int keep;          /* target for nrows-ncols */
  int mergelevelmax; /* says it */
  int32_t **R;           /* R[j][k] contains the rows of the non-empty elements
                        of column j, 0 <= j < ncols, 1 <= k <= R[j][0], for
                        weight(j) <= cwmax.
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
  int32_t *MKZA;         /* MKZA[j] gives u s.t. MKZQ[u] = j */ 
  int wmstmax;
  int mkzrnd;        /* to randomize things */
  int mkztype;       /* which type of count */
#endif
  int itermax;       /* used for performing some sampling */
} sparse_mat_t;

#ifdef __cplusplus
extern "C" {
#endif

#define weightRow(mat, i) (mat->rows[(i)][0] + mat->wburried[(i)])

#ifdef USE_MPI
#define GETJ(mat, j) ((j)-(mat)->jmin)
#else
#define GETJ(mat, j) (j)
#endif

extern void initMat(sparse_mat_t *mat, int32_t jmin, int32_t jmax);
extern void clearMat (sparse_mat_t *mat);
extern void initWeightFromFile(sparse_mat_t *mat, FILE *purgedfile, int skipfirst);
extern void fillmat(sparse_mat_t *mat);
extern int readmat (sparse_mat_t *mat, FILE *file, int skipfirst,
                    int skipheavycols, int verbose);

extern void remove_j_from_row(sparse_mat_t *mat, int i, int j);
extern void print_row(sparse_mat_t *mat, int i);

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
extern void freeRj(sparse_mat_t *mat, int j);
extern void remove_i_from_Rj(sparse_mat_t *mat, int i, int j);
extern void add_i_to_Rj(sparse_mat_t *mat, int i, int j);
extern int decrS(int w);
extern int incrS(int w);
extern int weightSum(sparse_mat_t *mat, int i1, int i2);
extern void fillTabWithRowsForGivenj(int32_t *ind, sparse_mat_t *mat, int32_t j);
extern void checkData(sparse_mat_t *mat);
extern void destroyRow(sparse_mat_t *mat, int i);

#ifdef __cplusplus
}
#endif

#endif	/* SPARSE_MAT_H_ */
