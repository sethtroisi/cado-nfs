#ifndef LAS_FORWARDTYPES_H_
#define LAS_FORWARDTYPES_H_

/* These must be forward-declared, because tvarious header files use them */

struct sieve_info_s;
typedef struct sieve_info_s * sieve_info_ptr;
typedef const struct sieve_info_s * sieve_info_srcptr;
struct las_info_s;
typedef struct las_info_s * las_info_ptr;
typedef const struct las_info_s * las_info_srcptr;
struct where_am_I_s;
typedef struct where_am_I_s * where_am_I_ptr;
typedef const struct where_am_I_s * where_am_I_srcptr;

#endif
