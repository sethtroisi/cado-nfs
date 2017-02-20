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

extern int MkzIsAlive(uint32_t *A, index_t dj);

extern int  MkzPopQueue(index_t *dj, index_signed_t *mkz, filter_matrix_t *mat);
extern void MkzRemove(index_t *dj, index_t *mkz, index_t *Q, uint32_t *A, index_t k);
extern int MkzIncrCol(filter_matrix_t *mat, index_t j);
extern void MkzUpdate(filter_matrix_t *mat, index_t i, index_t j);
extern void MkzUpdateDown (filter_matrix_t *mat, index_t j);
extern void MkzDecreaseColWeight(filter_matrix_t *mat, index_t j);
extern void MkzRemoveJ(filter_matrix_t *mat, index_t j);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_MERGE_MARKOWITZ_H_ */
