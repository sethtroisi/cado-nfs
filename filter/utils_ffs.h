/* Contains functions used when doing filter step for FFS instead of NFS */
#ifndef UTILS_FFS_H_
#define UTILS_FFS_H_

extern unsigned long findroot_ffs (long a, unsigned long b, unsigned long p);
extern void computeroots_ffs (relation_t * rel);
extern int ffs_poly_read(cado_poly poly, const char *filename);


#endif
