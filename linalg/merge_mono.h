#define USE_TAB 1 // 1 for compact rows...

#define MERGE_LEVEL_MAX 256 // maximum level for a merge; such a large value
                            // is only useful when not using BW

#define M_STRATEGY 3 // 0: finish that mergelevel
                     // 1: change if min weight < mergelevel
                     // 2: jump to minimal possible mergelevel
                     // 3: perform one merge, then check for next min weight

//#define USE_MARKOWITZ // says it...!

// 0 means dummy filtering using slow methods
// 1 means:
//  * fast with optimized-though-memory-consuming data structure;
//  * heavy columns are killed.
// 2 means:
//  * fast as above;
//  * heavy columns are present, but not in S[]; this yields more accurate
//    weights.
// 3 means:
//  * reintroduce columns whose weight is <= mergelevelmax (not suggested
//    for large numbers).
#define USE_MERGE_FAST 2

/* INT is defined in sparse.h */

// doubly chained lists
typedef struct dclist{
    INT j;
    struct dclist *prev, *next;
} *dclist;

/* rows correspond to relations, and columns to primes (or prime ideals) */
typedef struct {
  int nrows;
  int ncols;
  int rem_nrows;     /* number of remaining rows */
  int rem_ncols;     /* number of remaining columns */
  INT jmin, jmax;    /* we are interested in columns [jmin..jmax[ */
#if USE_TAB == 0
  rel_t *data;
#else
  INT **rows;
#endif
  int *wt;           /* weight w of column j, if w <= cwmax,
                        else <= 1 for a deleted column
                        (trick: we store -w if w > cwmax) */
  unsigned long *ad;
  unsigned long weight;
  int cwmax;         /* bound on weight of j to enter the SWAR structure */
  int rwmax;         /* if a weight(row) > rwmax, kill that row */
  int delta;         /* bound for nrows-ncols */
  int mergelevelmax; /* says it */
  dclist *S;         /* S[w] is a doubly-chained list with columns of weight w,
                        for w <= cwmax. For technical reasons (to avoid NULL
                        arguments that we can't modify), S[w] always contains
                        a first cell with value -1, thus the real list is
                        S[w]->next. */
  dclist *A;         /* A[j] points to the unique cell in S[w] containing j,
                        where w = weight(j), for w <= cwmax. A[j]=NULL for
                        w > cwmax. */
  INT **R;           /* R[j][k] contains the rows of the non-empty elements
                        of column j, 0 <= j < ncols, 1 <= k <= R[j][0], for
                        weight(j) <= cwmax.
                        R[j][k] = -1 if the corresponding row has been deleted.
                        R[j]=NULL for weight(j) > cwmax. */
#ifdef USE_MARKOWITZ
  INT *MKZQ;         /* priority queue for Markowitz stuff */    
#endif  
} sparse_mat_t;

// data structure for reporting actions during the merge; in standard mode
// (say mono proc), this is just a wrap around for a FILE; otherwise,
// it can be used to register things in an array that will be examined and
// flushed from time to time. See the MPI version for more.
typedef struct{
    char type;
    // '0' for the standard stuff
    FILE *outfile;
    // '1' for MPI
    INT **history;
    int mark;
} report_t;

#ifdef USE_MPI
#define GETJ(mat, j) ((j)-(mat)->jmin)
#else
#define GETJ(mat, j) (j)
#endif

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

extern void init_rep(report_t *rep, char *outname, sparse_mat_t *mat, int type);
extern void initMat(sparse_mat_t *mat, INT jmin, INT jmax);
extern void initWeightFromFile(sparse_mat_t *mat, FILE *purgedfile);
extern void fillSWAR(sparse_mat_t *mat);
extern void closeSWAR(/*sparse_mat_t *mat*/);

extern int readmat(sparse_mat_t *mat, FILE *file);
extern void report1(report_t *rep, INT i);
extern void report2(report_t *rep, INT i1, INT i2);
extern void removeCellSWAR(sparse_mat_t *mat, int i, INT j);
extern void addRowSWAR(sparse_mat_t *mat, int i);
extern void addRowsSWAR(sparse_mat_t *mat, int i1, int i2, int len);
extern void removeRowSWAR(sparse_mat_t *mat, int i);
extern void remove_j_from_SWAR(sparse_mat_t *mat, int j);
extern void destroyRow(sparse_mat_t *mat, int i);
extern int removeSingletons(report_t *rep, sparse_mat_t *mat);
extern int deleteEmptyColumns(sparse_mat_t *mat);
extern void removeRowDefinitely(report_t *rep, sparse_mat_t *mat, INT i);
extern int minColWeight(sparse_mat_t *mat);
extern void fillRowAddMatrix(int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], sparse_mat_t *mat, int m, INT *ind);
extern void MSTWithA(report_t *rep, sparse_mat_t *mat, int m, INT *ind, double *tMST, int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX]);
extern int deleteHeavyColumns(report_t *rep, sparse_mat_t *mat);

extern int cmp(const void *p, const void *q);
extern int number_of_superfluous_rows(sparse_mat_t *mat);
extern void merge(report_t *rep, sparse_mat_t *mat, int maxlevel, int verbose, int forbw);
extern void mergeOneByOne(report_t *rep, sparse_mat_t *mat, int maxlevel, int verbose, int forbw, double ratio, int coverNmax);
extern void doOneMerge(report_t *rep, sparse_mat_t *mat, int *njrem, double *totopt, double *totfill, double *totMST, double *totdel, int m, int maxdo, int useMST, int verbose);

extern void resume(report_t *rep, sparse_mat_t *mat, char *resumename);
