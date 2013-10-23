#ifndef CADO_MERGE_MARKOWITZ_H_
#define CADO_MERGE_MARKOWITZ_H_

#ifdef __cplusplus
extern "C" {
#endif

#define MkzIsQueueEmpty(Q) ((Q)[0] == 0)
#define MkzQueueCardinality(Q) ((Q)[0])

extern void MkzInit(filter_matrix_t *mat);
extern void MkzClose(filter_matrix_t *mat);

extern int MkzIsAlive(index_t *A, int32_t dj);

extern int  MkzPopQueue(int32_t *dj, int32_t *mkz, filter_matrix_t *mat);
extern void MkzRemove(int32_t *dj, int32_t *mkz, int32_t *Q, index_t *A, int32_t k);
extern int MkzIncrCol(filter_matrix_t *mat, int32_t j);
extern void MkzUpdate(filter_matrix_t *mat, int32_t i, int32_t j);
extern void MkzUpdateDown (filter_matrix_t *mat, int32_t i MAYBE_UNUSED, int32_t j);
extern void MkzDecreaseColWeight(filter_matrix_t *mat, int32_t j);
extern void MkzRemoveJ(filter_matrix_t *mat, int32_t j);
extern int MkzDeleteHeavyColumns(report_t *rep, filter_matrix_t *mat);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_MERGE_MARKOWITZ_H_ */
