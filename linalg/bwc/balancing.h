#ifndef BALANCING_H_
#define BALANCING_H_

#include "utils.h"

// #include "mf.h"

/* balancing structure */

#define FLAG_COLPERM    1
#define FLAG_ROWPERM    2
#define FLAG_PADDING    4       /* pad to largest dimension */
#define FLAG_REPLICATE  8       /* only balancing in one dimension */

struct balancing_header_s {
    // FIXME: add a magic number here ? This header is read directly in
    // binary format, so it might be a good idea.
    uint32_t nh;
    uint32_t nv;
    uint32_t nrows;
    uint32_t ncols;
    uint64_t ncoeffs;
    uint32_t checksum;
    uint32_t flags;
};
typedef struct balancing_header_s balancing_header[1];

struct balancing_s {
    balancing_header h;
    uint32_t trows;     // target number of rows. ==nrows if no padding
    uint32_t tcols;     // target number of cols. ==ncols if no padding
    uint32_t * rowperm; // row index for new mat. --> row index for old mat.
    uint32_t * colperm; // might be equal to colperm.
};
typedef struct balancing_s balancing[1];
typedef struct balancing_s * balancing_ptr;

#ifdef __cplusplus
extern "C" {
#endif

/* Once the flags and perm[] fields have been provided, the caller must
 * call _finalize() in order to 1) update the trows and tcols fields 2)
 * compute the checksum of the balancing.
 */
extern void balancing_finalize(balancing_ptr bal);
extern void balancing_write_inner(balancing_ptr bal, const char *);
extern void balancing_write(balancing_ptr bal, const char * , const char *);
extern void balancing_read(balancing_ptr bal, const char *);
extern void balancing_read_header(balancing_ptr bal, const char * filename);
extern void balancing_clear(balancing_ptr bal);
extern void balancing_init(balancing_ptr bal);

#ifdef __cplusplus
}
#endif

#endif	/* BALANCING_H_ */
