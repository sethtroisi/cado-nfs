#ifndef CADO_UTILS_MISC_H_
#define CADO_UTILS_MISC_H_

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

extern char * cado_strndup(const char * a, size_t n);
extern int cado_posix_memalign(void** ptr, size_t alignment, size_t size);

void * malloc_check(const size_t x);
void * aligned_malloc(size_t size, size_t alignment);
void * malloc_pagealigned(size_t sz); 

#ifdef __cplusplus
}
#endif

#endif	/* CADO_UTILS_MISC_H_ */
