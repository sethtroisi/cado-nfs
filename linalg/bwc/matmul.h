#ifndef MATMUL_H_
#define MATMUL_H_

#include <stdarg.h>
#include <stdint.h>

/* This header is common to the different matrix product implementations
 * (note that it depends on abase, though). So the matrix product codes
 * are meant to be used as interchangeable .o files.
 */
#include "abase.h"
#include "params.h"

struct matmul_public_s;

typedef struct matmul_public_s  * matmul_t;
typedef struct matmul_public_s  * matmul_ptr;

typedef matmul_ptr (*matmul_init_t)(abobj_ptr, param_list pl, int);
typedef void (*matmul_build_cache_t)(matmul_ptr, uint32_t *);
typedef int (*matmul_reload_cache_t)(matmul_ptr);
typedef void (*matmul_save_cache_t)(matmul_ptr);
typedef void (*matmul_mul_t)(matmul_ptr, abt *, abt const *, int);
typedef void (*matmul_report_t)(matmul_ptr, double scale);
typedef void (*matmul_clear_t)(matmul_ptr mm);
typedef void (*matmul_auxv_t)(matmul_ptr mm, int op, ...);
typedef void (*matmul_aux_t)(matmul_ptr mm, int op, va_list ap);

struct matmul_bindings_s {
	matmul_build_cache_t	build_cache;
	matmul_reload_cache_t	reload_cache;
	matmul_save_cache_t	save_cache;
	matmul_mul_t		mul;
	matmul_report_t		report;
	matmul_clear_t		clear;
	matmul_init_t		init;
	matmul_auxv_t		auxv;
	matmul_aux_t		aux;
        const char *            impl;
};

struct matmul_public_s {
    /* The fields here must be exposed by all implementations */
    unsigned int dim[2];        /* dim[0] is nrows, dim[1] is ncols. The
                                   organization of the matrix in memory
                                   also obeys the ``transposed'' flag */
    uint64_t ncoeffs;

    unsigned int store_transposed;

    int iteration[2];

    const char * locfile;

    char * cachefile_name;
    char * local_cache_copy;

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

extern matmul_ptr matmul_init(abobj_ptr, unsigned int, unsigned int, const char * locfile, const char * impl, param_list pl, int);

/* These two functions are the means of initializing the mm layer. No
 * bare _init() function is publicly accessible */
extern void matmul_build_cache(matmul_ptr mm, uint32_t *);
extern int matmul_reload_cache(matmul_ptr mm);

/* Exit point */
extern void matmul_clear(matmul_ptr mm);

/* Convenience. Once the data has been built, it may be saved */
extern void matmul_save_cache(matmul_ptr);

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
#define MATMUL_AUX_ZERO_STATS        11
extern void matmul_auxv(matmul_ptr, int, va_list ap);
extern void matmul_aux(matmul_ptr, int, ...);


/* Used for statistics only */
extern void matmul_report(matmul_ptr, double scale);


#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_H_ */
