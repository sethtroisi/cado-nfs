/* Contains functions used when doing filter step for FFS instead of NFS */
#ifndef UTILS_FFS_H_
#define UTILS_FFS_H_

index_t findroot_ffs (int64_t a, uint64_t b, index_t p);
int ffs_poly_read(cado_poly poly, const char *filename);
int sq_is_irreducible(sq_srcptr p);

#endif
