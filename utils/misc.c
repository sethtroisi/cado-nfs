#define _POSIX_C_SOURCE 200112L
#define _XOPEN_SOURCE   600     /* sometimes useful for posix_memalign */
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

void *malloc_check(const size_t x)
{
    void *p;
    p = malloc(x);
    ASSERT_ALWAYS(p != NULL);
    return p;
}

/* Not everybody has posix_memalign. In order to provide a viable
 * alternative, we need an ``aligned free'' matching the ``aligned
 * malloc''. We rely on posix_memalign if it is available, or else fall
 * back on ugly pointer arithmetic so as to guarantee alignment. Note
 * that not providing the requested alignment can have trouble some
 * consequences. At best, a performance hit, at worst a segv (sse-2
 * movdqa on a pentium4 causes a GPE if improperly aligned).
 */

void *malloc_aligned(size_t size, size_t alignment)
{
#ifdef HAVE_POSIX_MEMALIGN
    void *res = NULL;
    posix_memalign(&res, alignment, size);
    return res;
#else
    char * res;
    res = malloc(size + sizeof(size_t) + alignment);
    res += sizeof(size_t);
    size_t displ = alignment - ((unsigned long) res) % alignment;
    res += displ;
    memcpy(res - sizeof(size_t), &displ, sizeof(size_t));
    ASSERT_ALWAYS((((unsigned long) res) % alignment) == 0);
    return (void*) res;
#endif
}

void free_aligned(void * p, size_t size MAYBE_UNUSED, size_t alignment MAYBE_UNUSED)
{
#ifdef HAVE_POSIX_MEMALIGN
    free(p);
#else
    char * res = (char *) p;
    ASSERT_ALWAYS((((unsigned long) res) % alignment) == 0);
    size_t displ;
    memcpy(&displ, res - sizeof(size_t), sizeof(size_t));
    res -= displ;
    ASSERT_ALWAYS(displ == alignment - ((unsigned long) res) % alignment);
    res -= sizeof(size_t);
    free(res);
#endif
}

void *malloc_pagealigned(size_t sz)
{
    void *p = malloc_aligned(sz, getpagesize());
    ASSERT_ALWAYS(p != NULL);
    return p;
}

void free_pagealigned(void * p, size_t sz)
{
    free_aligned(p, sz, getpagesize());
}

