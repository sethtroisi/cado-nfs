#ifndef MATMUL_H_
#define MATMUL_H_

#include <stdarg.h>

/* This header is common to the different matrix product implementations
 * (note that it depends on abase, though). So the matrix product codes
 * are meant to be used as interchangeable .o files.
 */
#include "abase.h"
#include "params.h"

struct matmul_public_s;

typedef struct matmul_public_s  * matmul_t;
typedef struct matmul_public_s  * matmul_ptr;

typedef matmul_ptr (*matmul_build_t)(abobj_ptr, const char * filename, param_list pl);
typedef matmul_ptr (*matmul_reload_cache_t)(abobj_ptr, const char * filename, param_list pl);
typedef void (*matmul_save_cache_t)(matmul_ptr, const char * filename);
typedef void (*matmul_mul_t)(matmul_ptr, abt *, abt const *, int);
typedef void (*matmul_report_t)(matmul_ptr);
typedef void (*matmul_clear_t)(matmul_ptr mm);
typedef void (*matmul_auxv_t)(matmul_ptr mm, int op, ...);
typedef void (*matmul_aux_t)(matmul_ptr mm, int op, va_list ap);

struct matmul_bindings_s {
	matmul_build_t		build;
	matmul_reload_cache_t	reload_cache;
	matmul_save_cache_t	save_cache;
	matmul_mul_t		mul;
	matmul_report_t		report;
	matmul_clear_t		clear;
	matmul_auxv_t		auxv;
	matmul_aux_t		aux;
};

struct matmul_public_s {
    /* The fields here must be exposed by all implementations */
    unsigned int dim[2];        /* dim[0] is nrows, dim[1] is ncols. The
                                   organization of the matrix in memory
                                   also obeys the ``transposed'' flag */
    unsigned long ncoeffs;

    /* Now the virtual method table */
    struct matmul_bindings_s bind[1];

    /* The rest of the implementation-dependent storage comes right after
     * that, in memory. Therefore, only the pointer may be manipulated
     * freely. The datasize field merely indicates the __on-disk__ data
     * size. The in-memory data size is allowed to involve pointers, and
     * is not required to match the on-disk structure exactly.
     */
};


#ifdef __cplusplus
extern "C" {
#endif

/* These two functions are the means of initializing the mm layer. No
 * bare _init() function is publicly accessible */
extern matmul_ptr matmul_build(abobj_ptr, const char * filename, const char * impl, param_list pl);
extern matmul_ptr matmul_reload_cache(abobj_ptr, const char * filename, const char * impl, param_list pl);

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

#endif	/* MATMUL_H_ */
