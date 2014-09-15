#ifndef VERBOSE_H_
#define VERBOSE_H_

#include <stdarg.h>
#include <stdio.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* This provides fine-grain control over which messages we want. In order
 * to add a new flag, you must define its name as a preprocessor integer
 * macro here, and define a corresponding text string in verbose.c
 */

#define CADO_VERBOSE_PRINT_CMDLINE                      0
#define CADO_VERBOSE_PRINT_MODIFIED_FILES               1
#define CADO_VERBOSE_PRINT_COMPILATION_INFO             2
#define CADO_VERBOSE_PRINT_BWC_DISPATCH_SLAVES          3
#define CADO_VERBOSE_PRINT_BWC_DISPATCH_MASTER          4
#define CADO_VERBOSE_PRINT_BWC_TIMING_GRIDS             5
#define CADO_VERBOSE_PRINT_BWC_ITERATION_TIMINGS        6
#define CADO_VERBOSE_PRINT_BWC_CACHE_BUILD              7

extern uint64_t verbose_flag_word;

/* This must be called in single-threaded context, preferably at program
 * start */
extern void verbose_set_enabled_flags(param_list pl);

/* returns true if the following verbose flag is enabled */
static int verbose_enabled(int flag) {
    return verbose_flag_word & (UINT64_C(1) << flag);
}

static inline int verbose_vfprintf(FILE * f, int flag, const char * fmt, va_list ap)
{
    if (verbose_enabled(flag)) {
        return vfprintf(f, fmt, ap);
    }
    return 1;
}

static inline int verbose_vprintf(int flag, const char * fmt, va_list ap)
{
    return verbose_vfprintf(stdout, flag, fmt, ap);
}
static inline int verbose_fprintf(FILE * f, int flag, const char * fmt, ...)
{
    va_list ap;
    int rc;
    va_start(ap, fmt);
    rc = verbose_vfprintf(f, flag, fmt, ap);
    va_end(ap);
    return rc;
}
static inline int verbose_printf(int flag, const char * fmt, ...)
{
    va_list ap;
    int rc;
    va_start(ap, fmt);
    rc = verbose_vprintf(flag, fmt, ap);
    va_end(ap);
    return rc;
}



#ifdef __cplusplus
}
#endif

#endif	/* VERBOSE_H_ */
