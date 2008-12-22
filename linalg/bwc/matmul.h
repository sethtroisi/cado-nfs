#ifndef MATMUL_H_
#define MATMUL_H_

/* This header is common to the different matrix product implementations
 * (note that it depends on abase, though). So the matrix product codes
 * are meant to be used as interchangeable .o files.
 */
#include "abase.h"

struct matmul_public_s {
    /* The fields here must be exposed by all implementations */
    unsigned int nrows;
    unsigned int ncols;
    unsigned long ncoeffs;
    size_t datasize;
    /* The rest of the implementation-dependent storage comes right after
     * that, in memory. Therefore, only the pointer may be manipulated
     * freely. The datasize field merely indicates the __on-disk__ data
     * size. The in-memory data size is allowed to involve pointers, and
     * is not required to match the on-disk structure exactly.
     */
};

typedef struct matmul_public_s  * matmul_t;
typedef struct matmul_public_s  * matmul_ptr;

#ifdef __cplusplus
extern "C" {
#endif


/* here are the functions which exists irrespective of the implementation */
extern matmul_ptr matmul_build(abobj_ptr, const char * filename);
extern matmul_ptr matmul_reload_cache(abobj_ptr, const char * filename);
extern void matmul_save_cache(matmul_ptr, const char * filename);

/* matmul() is the main entry point for the cpu-intensive critical loop.
 * besides the obvious, there's a direction argument which tells whether
 * we're concerned with the linear operator represented by the matrix
 * (d==1) or its transpose (d==0)
 *
 * Note that not all implementations properly handle both values of d. On
 * paper, both are dual to each other, but doing so equally efficiently
 * with the same in-memory structure is not easy.
 */
extern void matmul(matmul_ptr, abt *, abt const *, int);


extern void matmul_report(matmul_ptr);
extern void matmul_clear(matmul_ptr mm);
#ifdef __cplusplus
}
#endif

/* Now some cpp glue which sets up the different options */
#define MATMUL_NAME(kind,func) matmul_ ## kind ## _ ## func
#define MATMUL_NAME_(kind,func) MATMUL_NAME(kind,func)
#define MATMUL(func) MATMUL_NAME_(MATMUL_PREFERRED,func)

#define MATMUL_PREFERRED sliced
#include "matmul-sliced.h"

// #include "matmul-basic.h"

#endif	/* MATMUL_H_ */
