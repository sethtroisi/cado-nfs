#define _POSIX_C_SOURCE 200112L
#define _XOPEN_SOURCE   600
#define _DARWIN_C_SOURCE /* for getpagesize */

#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>

#include "cado_config.h"
#include "macros.h"
#include "misc.h"

/* Not every libc has this, and providing a workalike is very easy */

char *cado_strndup(const char *a, size_t n)
{
    char *r = malloc(n+1);
    memcpy(r, a, n+1);
    r[n] = '\0';
    return r;
}

/* Only in recent posix. Mac OS does not have it. Unfortunately,
 * providing something that works just the same isn't easy: even though
 * adjusting the pointer would be doable (albeit quirky), then free()
 * would choke. Therefore, if we know that the real posix_memalign, use
 * it, or else do something stupid.
 */

int cado_posix_memalign(void **ptr, size_t alignment, size_t size)
{
#ifdef HAVE_POSIX_MEMALIGN
    return posix_memalign(ptr, alignment, size);
#else
    /* Otherwise do something quite stupid. It's actually fairly
     * problematic, since we have to way to ensure that we get in return
     * something which is simultaneously suitably aligned _and_ freeable
     * with free(). In some cases it's a possible performance hit, or
     * even possibly a segmentation fault (sse-2 movdqa on a pentium4).
     * We prefer to abort early, because it's preferrable to be directed
     * here rather than being misled by a SEGV.
     */
    size_t r = size % alignment;
    if (r && alignment < size) {
	size += alignment - r;
    }
    *ptr = malloc(size);
    ASSERT_ALWAYS((((unsigned long) *ptr) % alignment) == 0);
    if (*ptr == NULL) {
	return ENOMEM;
    } else {
	return 0;
    }
#endif
}

void *malloc_check(const size_t x)
{
    void *p;
    p = malloc(x);
    ASSERT_ALWAYS(p != NULL);
    return p;
}


void *aligned_malloc(size_t size, size_t alignment)
{
    void *res;
    int rc = cado_posix_memalign(&res, alignment, size);
    return rc == 0 ? res : NULL;
}

void *malloc_pagealigned(size_t sz)
{
    void *p = aligned_malloc(sz, getpagesize());
    ASSERT_ALWAYS(p != NULL);
    return p;
}
