#ifndef CADO_UTILS_MEMORY_H_
#define CADO_UTILS_MEMORY_H_

#include <stdlib.h>

extern void * malloc_check(const size_t x);
extern void * physical_malloc(const size_t x, const int affect);

extern long pagesize (void);
extern void * malloc_aligned(size_t size, size_t alignment);
extern void free_aligned(const void * ptr, size_t alignment);

extern void * malloc_pagealigned(size_t sz);
extern void free_pagealigned(const void * ptr);

void *contiguous_malloc(size_t);
void contiguous_free(const void *);

#endif
