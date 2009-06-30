#ifndef CADO_SPARSE_H_
#define CADO_SPARSE_H_

#include <stdint.h> /* for int32_t */

extern void fprintRow(FILE *file, int32_t *row);
extern int32_t * copyRow(int32_t *row);
extern void removeWeight(int32_t **rows, int *wt, int i);
extern void addRows(int32_t **rows, int i1, int i2, int len0);
extern int hasCol(int32_t **rows, int i, int32_t j);
extern int cmp(const void *p, const void *q);

#endif  /* CADO_SPARSE_H_ */
