#ifndef CADO_UTILS_MEMORY_H_
#define CADO_UTILS_MEMORY_H_

#include <stdlib.h>

void *contiguous_malloc(size_t);
void contiguous_free(const void *, size_t);

#endif
