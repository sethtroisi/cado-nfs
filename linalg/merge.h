#define USE_TAB 1 // 1 for compact rows...

// doubly chained lists
typedef struct dclist{
    int j;
    struct dclist *prev, *next;
} *dclist;

typedef struct {
  int nrows;
  int ncols;
  int rem_nrows;     /* number of remaining rows */
  int rem_ncols;     /* number of remaining columns */
#if USE_TAB == 0
  rel_t *data;
#else
  int **rows;
#endif
  int *wt; /* weight of prime j, <= 1 for a deleted prime */
  unsigned long *ad;
  int weight;
  int cwmax;         /* bound on weight of j to enter the SWAR structure */
  int rwmax;         /* if a weight(row) > rwmax, kill that row */
  int delta;         /* bound for nrows-ncols */
  int mergelevelmax; /* says it */
  dclist *S, *A;
  int **R;
} sparse_mat_t;

#if USE_TAB == 0
#define isRowNull(mat, i) ((mat)->data[(i)].val == NULL)
#define lengthRow(mat, i) (mat)->data[(i)].len
#define cell(mat, i, k) (mat)->data[(i)].val[(k)]
#else
#define isRowNull(mat, i) ((mat)->rows[(i)] == NULL)
#define lengthRow(mat, i) (mat)->rows[(i)][0]
#define cell(mat, i, k) (mat)->rows[(i)][(k)]
#endif

extern void report1(int i);
extern void removeCellSWAR(sparse_mat_t *mat, int i, int j);
extern void destroyRow(sparse_mat_t *mat, int i);
extern int removeSingletons(sparse_mat_t *mat);
extern int deleteEmptyColumns(sparse_mat_t *mat);
