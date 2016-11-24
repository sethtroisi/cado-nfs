#ifndef CADO_MERGE_MARKOWITZ_H_
#define CADO_MERGE_MARKOWITZ_H_

#define MKZTYPE_CAVALLAR 0
#define MKZTYPE_PURE 1
#define MKZTYPE_LIGHT 2

#ifdef __cplusplus
extern "C" {
#endif

#define MkzIsQueueEmpty(Q) ((Q)[0] == 0)
#define MkzQueueCardinality(mat) ((mat->MKZQ)[0])

extern void MkzInit(filter_matrix_t *mat, int verbose);
extern void MkzClear(filter_matrix_t *mat, int verbose);

extern int MkzIsAlive(index_t *A, int32_t dj);

extern int  MkzPopQueue(int32_t *dj, int32_t *mkz, filter_matrix_t *mat);
extern void MkzRemove(int32_t *dj, int32_t *mkz, int32_t *Q, index_t *A, int32_t k);
extern int MkzIncrCol(filter_matrix_t *mat, int32_t j);
extern void MkzUpdate(filter_matrix_t *mat, int32_t i, int32_t j);
extern void MkzUpdateDown (filter_matrix_t *mat, int32_t j);
extern void MkzDecreaseColWeight(filter_matrix_t *mat, int32_t j);
extern void MkzRemoveJ(filter_matrix_t *mat, int32_t j);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_MERGE_MARKOWITZ_H_ */
