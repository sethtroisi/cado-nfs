/* INT is defined in sparse.h */

extern void initSWAR(sparse_mat_t *mat);
extern void closeSWAR(/*sparse_mat_t *mat*/);
extern void texSWAR(sparse_mat_t *mat);
extern int addColSWAR(sparse_mat_t *mat, INT j);
extern void printStatsSWAR(sparse_mat_t *mat);
extern int deleteEmptyColumnsSWAR(sparse_mat_t *mat);

