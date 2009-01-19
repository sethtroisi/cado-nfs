#ifndef ABASE_H_
#define ABASE_H_

#if defined(SELECT_ABASE_u64)
#include "abase-u64.h"
#elif defined(SELECT_ABASE_u64k)
#include "abase-u64k.h"
#elif defined(SELECT_ABASE_u64n)
#include "abase-u64n.h"
#elif defined(SELECT_ABASE_u128)
#include "abase-u128.h"
#else

#warning "Using default selection for abase"
#include "abase-u64.h"
#endif

#endif	/* ABASE_H_ */
