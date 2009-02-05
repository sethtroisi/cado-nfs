#define _POSIX_C_SOURCE 200112L

#include <stdlib.h>
#include "aligned_malloc.h"

void * aligned_malloc(size_t size, size_t alignment)
{
    void * res;
    int rc = posix_memalign(&res, alignment, size);
    return rc == 0 ? res : NULL;
}

