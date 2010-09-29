extern void initSWAR(filter_matrix_t *mat);
extern void closeSWAR(/*filter_matrix_t *mat*/);
extern void texSWAR(filter_matrix_t *mat);
extern int addColSWAR(filter_matrix_t *mat, int32_t j);
extern void printStatsSWAR(filter_matrix_t *mat);
extern int deleteEmptyColumnsSWAR(filter_matrix_t *mat);
extern void decreaseColWeightSWAR(filter_matrix_t *mat, int32_t j);
extern void remove_j_from_SWAR(filter_matrix_t *mat, int j);

