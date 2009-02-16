extern void initSWAR(sparse_mat_t *mat);
extern void closeSWAR(/*sparse_mat_t *mat*/);
extern void texSWAR(sparse_mat_t *mat);
extern int addColSWAR(sparse_mat_t *mat, int32_t j);
extern void printStatsSWAR(sparse_mat_t *mat);
extern int deleteEmptyColumnsSWAR(sparse_mat_t *mat);
extern void decreaseColWeightSWAR(sparse_mat_t *mat, int32_t j);
extern void remove_j_from_SWAR(sparse_mat_t *mat, int j);

