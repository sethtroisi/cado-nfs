#define USE_TAB 1 // 1 for compact rows...

#define TRACE_COL -1 // 253224 // 231 // put to -1 if not...!
#define TRACE_ROW -1 // 59496 // put to -1 if not...!

/* rows correspond to relations, and columns to primes (or prime ideals) */
typedef struct {
  int nrows;
  int ncols;
  int rem_nrows;     /* number of remaining rows */
  int rem_ncols;     /* number of remaining columns */
  INT jmin, jmax;    /* we are interested in columns [jmin..jmax[ */
  INT **rows;
  int *wt;           /* weight w of column j, if w <= cwmax,
                        else <= 1 for a deleted column
                        (trick: we store -w if w > cwmax) */
  unsigned long *ad;
  unsigned long weight;
  int cwmax;         /* bound on weight of j to enter the SWAR structure */
  int rwmax;         /* if a weight(row) > rwmax, kill that row */
  int delta;         /* bound for nrows-ncols */
  int mergelevelmax; /* says it */
  INT **R;           /* R[j][k] contains the rows of the non-empty elements
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
  INT *MKZQ;         /* priority queue for Markowitz stuff */    
  INT *MKZA;         /* MKZA[j] gives u s.t. MKZQ[u] = j */ 
  int wmstmax;
#endif
} sparse_mat_t;

#ifdef USE_MPI
#define GETJ(mat, j) ((j)-(mat)->jmin)
#else
#define GETJ(mat, j) (j)
#endif

extern void addRowsWithWeight(sparse_mat_t *mat, int i1, int i2);
extern void removeWeightFromRow(sparse_mat_t *mat, int i);
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
extern int findAllRowsWithGivenj(INT *ind, sparse_mat_t *mat, INT j, int nb);
extern void fillTabWithRowsForGivenj(INT *ind, sparse_mat_t *mat, INT j);


