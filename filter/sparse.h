#ifndef CADO_SPARSE_H_
#define CADO_SPARSE_H_

#include <stdint.h> /* for int32_t */

#ifndef FOR_FFS
#define typerow_t int32_t
#define rowLength(rows, i) rows[(i)][0]
#define rowCell(rows, i, k) rows[(i)][(k)]
#else
#define typerow_t ideal_merge_ffs_t
#define rowLength(rows, i) rows[(i)][0].id
#define rowCell(rows, i, k) rows[(i)][(k)].id
#endif

extern void fprintRow(FILE *file, typerow_t *row);
extern int32_t * copyRow(int32_t *row);
extern void removeWeight(int32_t **rows, int *wt, int i);
extern void addRows(typerow_t **rows, int i1, int i2, int32_t j);
extern int hasCol(int32_t **rows, int i, int32_t j);
extern int cmp(const void *p, const void *q);

extern int parse_hisfile_line (int32_t *ind, char *t, int32_t *j);
#endif  /* CADO_SPARSE_H_ */
