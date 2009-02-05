#ifndef ABASE_H_
#define ABASE_H_

/* This has to match the largest possible alignment requirement from
 * different abase layers.
 */
#define ALWAYS_ALIGN_LARGE_MALLOCS      16

#include "aligned_malloc.h"

#ifdef  ALWAYS_ALIGN_LARGE_MALLOCS
#define alignable_malloc(s)     aligned_malloc((s), ALWAYS_ALIGN_LARGE_MALLOCS)
#else
#define alignable_malloc(s)     malloc((s))
#endif
#define alignable_free(p,s)     free((p))


#include "abase-generic.h"

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
