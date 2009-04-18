#ifndef ABASE_COMMON_H_
#define ABASE_COMMON_H_

/* This has to match the largest possible alignment requirement from
 * different abase layers.
 */
#define ALWAYS_ALIGN_LARGE_MALLOCS      16

/* Make sure that large mallocs always have space for a multiple of this
 * value. Must be a power of 2 in any case.
 *
 * XXX Note that if ever there is a discrepancy here between the
 * ``global'' readahead value, and a ``local'' one defined in an abase
 * header, there could be some
 * problems.
 * --> This global readahead value must be a multiple of all local ones.
 * Otherwise, requirements of the extra abase may be unsatisfied.
 * --> If it is a strict multiple, one must ensure that alloc/free are
 * performed by symetric functions. Not the generic abase_generic_init
 * for one, and the specific abclear for the other...
 */
#define ABASE_UNIVERSAL_READAHEAD_ITEMS       2

#if 0   /* only for very desperate cases */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "electric_alloc.h"
#define alignable_malloc(s)     electric_alloc(next_multiple_of_powerof2(s,16))
#define alignable_free(p,s)     electric_free(p,next_multiple_of_powerof2(s,16))

#else

#ifdef  ALWAYS_ALIGN_LARGE_MALLOCS
#include "misc.h"       /* in aligned_malloc is in utils */
#define alignable_malloc(s)     aligned_malloc((s), ALWAYS_ALIGN_LARGE_MALLOCS)
#define alignable_free(p,s)     free((p))
#else
#define alignable_malloc(s)     malloc((s))
#define alignable_free(p,s)     free((p))
#endif

#endif

#endif	/* ABASE_COMMON_H_ */
