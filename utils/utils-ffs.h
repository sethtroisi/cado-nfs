#ifndef	CADO_UTILS_FFS_H_
#define	CADO_UTILS_FFS_H_
/* Contains functions used when doing filter step for FFS instead of NFS */

#include "utils_with_io.h"

p_r_values_t ffs_relation_compute_r (int64_t a, uint64_t b, p_r_values_t p);
int ffs_poly_read(cado_poly poly, const char *filename);


#endif	/* CADO_UTILS_FFS_H_ */
