#ifndef DEBUG_H_
#define DEBUG_H_

#include "macros.h"

#ifdef __cplusplus
extern "C" {
#endif

void debug_write(const void * v, size_t n, const char * fmt, ...) ATTR_PRINTF(3,4);


#ifdef __cplusplus
}
#endif

#endif	/* DEBUG_H_ */
