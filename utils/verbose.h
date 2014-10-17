#ifndef VERBOSE_H_
#define VERBOSE_H_

#include <stdarg.h>
#include <stdio.h>
#include "params.h"

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
#define CADO_VERBOSE_PRINT_BWC_DISPATCH_OUTER           8
#define CADO_VERBOSE_PRINT_BWC_CPUBINDING               9
#define CADO_VERBOSE_PRINT_BWC_CACHE_MAJOR_INFO         10
#define CADO_VERBOSE_PRINT_BWC_LOADING_MKSOL_FILES      11

/* This must be called in single-threaded context, preferably at program
 * start */
extern void verbose_set_enabled_flags(param_list pl);

extern int verbose_enabled(int flag);
extern int verbose_vfprintf(FILE * f, int flag, const char * fmt, va_list ap);
extern int verbose_vprintf(int flag, const char * fmt, va_list ap);
extern int verbose_fprintf(FILE * f, int flag, const char * fmt, ...);
extern int verbose_printf(int flag, const char * fmt, ...);
int verbose_output_init(size_t);
int verbose_output_clear();
int verbose_output_add(size_t, FILE *, int);
int verbose_output_print(size_t, int, const char *, ...);
int verbose_decl_usage(param_list pl);
FILE *verbose_output_get(size_t, int, size_t);

#ifdef __cplusplus
}
#endif

#endif	/* VERBOSE_H_ */
