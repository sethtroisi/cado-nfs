#ifndef ABASE_COMMON_H_
#define ABASE_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

#define ABASE_CURRENT_INTERFACE abase_generic
typedef void * abase_generic_obj_ptr;
typedef void abase_generic_base_type;

#define ABASE_CODE_HEADER struct abase_function_pointers {
#define ABASE_CODE_TRAILER };
#define ABASE_MKCODE(t,n,a) t (*n) a;
#include "abase-api.h"

#ifdef __cplusplus
}
#endif

#endif	/* ABASE_COMMON_H_ */
