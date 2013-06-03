/* Contains functions used when doing filter step for FFS instead of NFS */
#ifndef UTILS_FFS_H_
#define UTILS_FFS_H_
#include "utils.h"

extern index_t findroot_ffs (int64_t a, uint64_t b, index_t p);
extern int ffs_poly_read(cado_poly poly, const char *filename);

#endif
