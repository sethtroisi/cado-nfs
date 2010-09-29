#ifndef MERGE_MONO_H_
#define MERGE_MONO_H_

#ifdef __cplusplus
extern "C" {
#endif

extern void addRowsAndUpdate(filter_matrix_t *mat, int i1, int i2, int len);
extern int deleteEmptyColumns(filter_matrix_t *mat);
extern void removeRowDefinitely(report_t *rep, filter_matrix_t *mat, int32_t i);
extern void MSTWithA(report_t *rep, filter_matrix_t *mat, int m, int32_t *ind, double *tMST, int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX]);
extern int deleteHeavyColumns(report_t *rep, filter_matrix_t *mat);
extern int addFatherToSons(int history[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1], filter_matrix_t *mat, int m, int *ind,	int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], int *father, int *height MAYBE_UNUSED, int hmax MAYBE_UNUSED, int sons[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1]);

extern int number_of_superfluous_rows(filter_matrix_t *mat);
extern void merge(report_t *rep, filter_matrix_t *mat, int maxlevel, int verbose, int forbw);
extern void mergeOneByOne(report_t *rep, filter_matrix_t *mat, int maxlevel, int verbose, int forbw, double ratio, int coverNmax);
extern void doOneMerge(report_t *rep, filter_matrix_t *mat, int *njrem, double *totopt, double *totfill, double *totMST, double *totdel, int m, int maxdo, int useMST, int verbose);

extern void resume(report_t *rep, filter_matrix_t *mat, char *resumename);

#ifdef __cplusplus
}
#endif

#endif	/* MERGE_MONO_H_ */
