extern void initMat(sparse_mat_t *mat, INT jmin, INT jmax);
extern void initWeightFromFile(sparse_mat_t *mat, FILE *purgedfile);

extern void fillmat(sparse_mat_t *mat);
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
