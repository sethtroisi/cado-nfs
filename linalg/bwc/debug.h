#ifndef DEBUG_H_
#define DEBUG_H_

#include "abase.h"
#include "macros.h"

#ifdef __cplusplus
extern "C" {
#endif

void debug_write(abobj_ptr abase, const abt * v,
        unsigned int n, const char * fmt, ...) ATTR_PRINTF(4,5);


#ifdef __cplusplus
}
#endif

#endif	/* DEBUG_H_ */
