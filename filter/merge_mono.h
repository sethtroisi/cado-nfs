#ifndef MERGE_MONO_H_
#define MERGE_MONO_H_

#ifdef __cplusplus
extern "C" {
#endif

void mergeOneByOne(report_t *rep, filter_matrix_t *mat, int maxlevel, double target_density);

void resume(report_t *rep, filter_matrix_t *mat, const char *resumename);

#ifdef __cplusplus
}
#endif

#endif	/* MERGE_MONO_H_ */
