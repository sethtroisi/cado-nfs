#define MkzIsQueueEmpty(Q) ((Q)[0] == 0)
#define MkzQueueCardinality(Q) ((Q)[0])

extern void MkzInit(sparse_mat_t *mat);
extern void MkzClose(sparse_mat_t *mat);

extern int MkzGetCount(int32_t *Q, int32_t *A, int32_t dj);
extern int MkzIsAlive(int32_t *A, int32_t dj);

extern void MkzPopQueue(int32_t *dj, int32_t *mkz, int32_t *Q, int32_t *A);
extern void MkzRemove(int32_t *dj, int32_t *mkz, int32_t *Q, int32_t *A, int32_t k);
extern int MkzIncrCol(sparse_mat_t *mat, int32_t j);
extern void MkzUpdate(sparse_mat_t *mat, int32_t i, int32_t j);
extern void MkzDecreaseColWeight(sparse_mat_t *mat, int32_t j);
extern void MkzRemoveJ(sparse_mat_t *mat, int32_t j);
extern int MkzDeleteHeavyColumns(report_t *rep, sparse_mat_t *mat);
extern int MkzRemoveCols(report_t *rep, sparse_mat_t *mat, int wmin, int wmax);
