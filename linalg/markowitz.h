extern void MkzInit(sparse_mat_t *mat);
extern void MkzPopQueue(INT *dj, INT *mkz, INT *Q, INT *A);
extern int MkzAddCol(sparse_mat_t *mat, INT j);
extern void MkzUpdate(sparse_mat_t *mat, INT j);
extern void MkzRemoveCol(sparse_mat_t *mat, INT j);
extern void MkzRemoveJ(sparse_mat_t *mat, INT j);
