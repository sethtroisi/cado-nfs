#ifndef CADO_UTILS_MEMORY_H_
#define CADO_UTILS_MEMORY_H_

#include <stdlib.h>
#include "macros.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void * malloc_check(const size_t x);
extern void * physical_malloc(const size_t x, const int affect);
extern void physical_free(const void *, size_t);

void *malloc_hugepages(size_t);
extern long pagesize (void);
ATTRIBUTE((malloc)) extern void * malloc_aligned(size_t size, size_t alignment);
ATTRIBUTE((warn_unused_result)) void * realloc_aligned(void * p, 
        const size_t old_size, const size_t new_size, const size_t alignment);
extern void free_aligned(const void * ptr);

extern void * malloc_pagealigned(size_t sz);
extern void free_pagealigned(const void * ptr);

void *contiguous_malloc(size_t);
void contiguous_free(const void *);

#ifdef __cplusplus
}
#endif

#endif
