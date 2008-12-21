#define USE_TAB 1 // 1 for compact rows...

#define MERGE_LEVEL_MAX 256 // maximum level for a merge; such a large value
                            // is only useful when not using BW

#define M_STRATEGY 3 // 0: finish that mergelevel
                     // 1: change if min weight < mergelevel
                     // 2: jump to minimal possible mergelevel
                     // 3: perform one merge, then check for next min weight

#define TRACE_COL -1 // 253224 // 231 // put to -1 if not...!
#define TRACE_ROW -1 // put to -1 if not...!

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
// TODO: remove these one day, since the dependency is strange in swar.c
extern void destroyRj(sparse_mat_t *mat, int j);
extern void remove_i_from_Rj(sparse_mat_t *mat, int i, int j);
extern void add_i_to_Rj(sparse_mat_t *mat, int i, int j);
extern int decrS(int w);
extern int incrS(int w);
// TODO_END


extern void initMat(sparse_mat_t *mat, INT jmin, INT jmax);
extern void initWeightFromFile(sparse_mat_t *mat, FILE *purgedfile);

extern int readmat(sparse_mat_t *mat, FILE *file);
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
