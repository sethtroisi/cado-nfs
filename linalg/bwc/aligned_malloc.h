#ifndef ALIGNED_MALLOC_H_
#define ALIGNED_MALLOC_H_

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

extern void * aligned_malloc(size_t size, size_t alignment);

#ifdef __cplusplus
}
#endif


#endif	/* ALIGNED_MALLOC_H_ */
