#ifndef MATMUL_H_
#define MATMUL_H_

#include "abase.h"

#ifdef __cplusplus
extern "C" {
#endif

struct matmul_public_s {
    /* The three fields here must be exposed by all implementations */
    unsigned int nrows;
    unsigned int ncols;
    unsigned long ncoeffs;
    /* The rest of the implementation-dependent storage comes right after
     * that, in memory. Therefore, only the pointer may be manipulated
     * freely.
     */
};

typedef struct matmul_public_s  * matmul_t;
typedef struct matmul_public_s  * matmul_ptr;

extern matmul_ptr matmul_build(abobj_ptr, const char * filename);
extern matmul_ptr matmul_reload_cache(abobj_ptr, const char * filename);
extern void matmul_save_cache(matmul_ptr, const char * filename);
extern void matmul_report(matmul_ptr);
extern void matmul(matmul_ptr, abt *, abt const *);
extern void matmul_clear(matmul_ptr mm);

#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_H_ */
