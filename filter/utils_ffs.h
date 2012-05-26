/* Contains functions used when doing filter step for FFS instead of NFS */
#ifndef UTILS_FFS_H_
#define UTILS_FFS_H_
#include "utils.h"

typedef struct {
  int32_t id; 
  int32_t e;   
} ideal_merge_ffs_t;

extern unsigned int weight_ffs (relation_t rel);
extern unsigned long findroot_ffs (long a, unsigned long b, unsigned long p);
extern void computeroots_ffs (relation_t * rel);
extern int ffs_poly_read(cado_poly poly, const char *filename);
//extern void sort_ffs (filter_matrix_t *mat, int i, int32_t size);

#endif
