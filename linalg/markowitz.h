#define MkzIsQueueEmpty(Q) ((Q)[0] == 0)
#define MkzQueueCardinality(Q) ((Q)[0])

extern void MkzInit(sparse_mat_t *mat);
extern void MkzPopQueue(INT *dj, INT *mkz, INT *Q, INT *A);
extern int MkzIncrCol(sparse_mat_t *mat, INT j);
extern void MkzUpdate(sparse_mat_t *mat, INT j);
extern void MkzDecreaseColWeight(sparse_mat_t *mat, INT j);
extern void MkzRemoveJ(sparse_mat_t *mat, INT j);
