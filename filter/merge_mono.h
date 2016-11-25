#ifndef MERGE_MONO_H_
#define MERGE_MONO_H_

#ifdef __cplusplus
extern "C" {
#endif

void removeRowDefinitely(report_t *rep, filter_matrix_t *mat, int32_t i);
int deleteHeavyColumns(report_t *rep, filter_matrix_t *mat);
int addFatherToSons(int history[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1], filter_matrix_t *mat, int m, int *ind, int32_t j,	int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], int *father, int *height MAYBE_UNUSED, int hmax MAYBE_UNUSED, int sons[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1]);

void mergeOneByOne(report_t *rep, filter_matrix_t *mat, int maxlevel, double target_density);

void resume(report_t *rep, filter_matrix_t *mat, const char *resumename);

#ifdef __cplusplus
}
#endif

#endif	/* MERGE_MONO_H_ */
