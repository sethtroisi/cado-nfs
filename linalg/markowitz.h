#define MkzIsQueueEmpty(Q) ((Q)[0] == 0)
#define MkzQueueCardinality(Q) ((Q)[0])

extern void MkzInit(sparse_mat_t *mat);
extern void MkzClose(sparse_mat_t *mat);

extern int MkzGetCount(INT *Q, INT *A, INT dj);
extern int MkzIsAlive(INT *A, INT dj);

extern void MkzPopQueue(INT *dj, INT *mkz, INT *Q, INT *A);
extern int MkzIncrCol(sparse_mat_t *mat, INT j);
extern void MkzUpdate(sparse_mat_t *mat, INT i, INT j);
extern void MkzDecreaseColWeight(sparse_mat_t *mat, INT j);
extern void MkzRemoveJ(sparse_mat_t *mat, INT j);
extern int MkzDeleteHeavyColumns(report_t *rep, sparse_mat_t *mat);
extern int MkzRemoveCols(report_t *rep, sparse_mat_t *mat, int wmin, int wmax);
