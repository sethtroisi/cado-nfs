#ifndef MATMUL_H_
#define MATMUL_H_

#include <stdarg.h>

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

/* These two functions are the means of initializing the mm layer. No
 * bare _init() function is publicly accessible */
extern matmul_ptr matmul_build(abobj_ptr, const char * filename);
extern matmul_ptr matmul_reload_cache(abobj_ptr, const char * filename);

/* Exit point */
extern void matmul_clear(matmul_ptr mm);

/* Convenience. Once the data has been built, it may be saved */
extern void matmul_save_cache(matmul_ptr, const char * filename);

/* matmul_mul() is the main entry point for the cpu-intensive critical loop.
 * besides the obvious, there's a direction argument which tells whether
 * we're concerned with the linear operator represented by the matrix
 * (d==1) or its transpose (d==0)
 *
 * Note that not all implementations properly handle both values of d. On
 * paper, both are dual to each other, but doing so equally efficiently
 * with the same in-memory structure is not easy.
 */
extern void matmul_mul(matmul_ptr, abt *, abt const *, int);

/* matmul_aux is used on some occasions for obtaining private
 * informations on the matmul structure. It's really a virtual method
 * hiding its name.
 *
 * _aux() is only the top-level. auxv is the one that does something.
 */
#define MATMUL_AUX_GET_READAHEAD        10
extern void matmul_auxv(matmul_ptr, int, va_list ap);
extern void matmul_aux(matmul_ptr, int, ...);


/* Used for statistics only */
extern void matmul_report(matmul_ptr);


#ifdef __cplusplus
}
#endif

/* Now some cpp glue which sets up the different options */
#define MATMUL_NAME(kind,func) matmul_ ## kind ## _ ## func
#define MATMUL_NAME_(kind,func) MATMUL_NAME(kind,func)
#define MATMUL(func) MATMUL_NAME_(MATMUL_PREFERRED,func)

#define MATMUL_PREFERRED threaded
#include "matmul-threaded.h"

// #include "matmul-basic.h"

#endif	/* MATMUL_H_ */
