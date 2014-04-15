
#include "cado.h"
#include <stdio.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>
#include <limits.h>
#include <unistd.h>
#include <ctype.h>
#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif

#ifdef HAVE_SYS_MMAN_H
#include <sys/mman.h>
#endif

/* For MinGW Build */
#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#endif

#include "macros.h"
#include "portability.h"
#include "misc.h"
#include "memory.h"

#ifndef LARGE_PAGE_SIZE
#define LARGE_PAGE_SIZE (2UL*1024*1024)
#endif

void
*malloc_check (const size_t x)
{
    void *p;
    p = malloc (x);
    if (p == NULL)
      {
        fprintf (stderr, "Error, malloc of %zu bytes failed\n", x);
        fflush (stderr);
        abort ();
      }
    return p;
}


void *
malloc_hugepages(const size_t size)
{
#ifdef MADV_HUGEPAGE
  void *m = malloc_aligned(size, LARGE_PAGE_SIZE);
  int r;
  static int printed_error = 0;
  do {
    r = madvise(m, size, MADV_HUGEPAGE);
  } while (r == EAGAIN);
  if (r != 0 && !printed_error) {
    perror("madvise failed");
    printed_error = 1;
  }
  return m;
#else
  return malloc_pagealigned(size);
#endif
}

void
*physical_malloc (const size_t x, const int affect)
{
  void *p;
  p = malloc_hugepages(x);
  if (affect) {
    size_t i, m;
#ifdef HAVE_SSE2
    const __m128i a = (__m128i) {0, 0};
#endif    
    i = ((size_t) p + 15) & (~15ULL);
    m = ((size_t) p + x - 1) & (~15ULL);
    while (i < m) {
#ifdef HAVE_SSE2
      _mm_stream_si128((__m128i *)i, a);
#else
      *(unsigned char *) i = 0;
#endif
      i += pagesize ();
    }
  }
  return p;
}

/* Not everybody has posix_memalign. In order to provide a viable
 * alternative, we need an ``aligned free'' matching the ``aligned
 * malloc''. We rely on posix_memalign if it is available, or else fall
 * back on ugly pointer arithmetic so as to guarantee alignment. Note
 * that not providing the requested alignment can have some troublesome
 * consequences. At best, a performance hit, at worst a segv (sse-2
 * movdqa on a pentium4 causes a GPE if improperly aligned).
 */

void *malloc_aligned(size_t size, size_t alignment)
{
#ifdef HAVE_POSIX_MEMALIGN
    void *res = NULL;
    int rc = posix_memalign(&res, alignment, size);
    ASSERT_ALWAYS(rc == 0);
    return res;
#else
    char * res;
    res = malloc(size + sizeof(size_t) + alignment);
    res += sizeof(size_t);
    size_t displ = alignment - ((uintptr_t) res) % alignment;
    res += displ;
    memcpy(res - sizeof(size_t), &displ, sizeof(size_t));
    ASSERT_ALWAYS((((uintptr_t) res) % alignment) == 0);
    return (void*) res;
#endif
}

void free_aligned(const void * p, size_t alignment MAYBE_UNUSED)
{
#ifdef HAVE_POSIX_MEMALIGN
    free((void *) p);
#else
    const char * res = (const char *) p;
    ASSERT_ALWAYS((((uintptr_t) res) % alignment) == 0);
    size_t displ;
    memcpy(&displ, res - sizeof(size_t), sizeof(size_t));
    res -= displ;
    ASSERT_ALWAYS((displ + (uintptr_t) res) % alignment == 0);
    res -= sizeof(size_t);
    free((void *)res);
#endif
}

long pagesize (void)
{
#if defined(_WIN32) || defined(_WIN64)
  /* cf http://en.wikipedia.org/wiki/Page_%28computer_memory%29 */
  SYSTEM_INFO si;
  GetSystemInfo(&si);
  return si.dwPageSize;
#elif defined(HAVE_SYSCONF)
  return sysconf (_SC_PAGESIZE);
#else
  #error "Cannot determine page size"
#endif
}

void *malloc_pagealigned(size_t sz)
{
    void *p = malloc_aligned (sz, pagesize ());
    ASSERT_ALWAYS(p != NULL);
    return p;
}

void free_pagealigned(const void * p)
{
    free_aligned(p, pagesize ());
}

/* Functions for allocating contiguous physical memory, if large pages are available.
   If not, they just return malloc_pagealigned().
   
   We keep track of allocated memory with a linked list. No attempt is made to fill
   freed holes with later requests, new requests always get memory after the 
   furthest-back allocated chunk. We assume that all of the contiguous memory needed
   by a program is allocated at the start, and all freed at the end, so that this
   restriction has no impact. */

struct largepage_chunk
{
    size_t offset, size;
    struct largepage_chunk *next;
};

/* NULL means we never tried mmap() yet. MAP_FAILED means we tried and
   it failed. */
static void *largepages = NULL;
static size_t largepage_size = 0;
static struct largepage_chunk *chunks = NULL;

void *contiguous_malloc(const size_t size)
{
  if (size == 0) {
    return NULL;
  }

#if defined(HAVE_MMAP) && defined(MAP_HUGETLB)
  int dont_map = (getenv("LAS_NO_MAP_HUGETLB") != NULL);
  static int printed_dont_map = 0;
  if (dont_map && !printed_dont_map) {
    fprintf(stderr, "mmap()-ing huge page prevented by LAS_NO_MAP_HUGETLB\n");
    printed_dont_map = 1;
  }
  
  /* Try to mmap a large page, if we haven't already */
  if (largepages == NULL && !dont_map) {
#ifdef MAP_ANONYMOUS
    const int anon = MAP_ANONYMOUS;
#else
    const int anon = MAP_ANON;
#endif
    largepages = mmap (NULL, LARGE_PAGE_SIZE, PROT_READ|PROT_WRITE, anon | MAP_PRIVATE | MAP_HUGETLB, -1, 0);
    if (largepages != MAP_FAILED) {
      largepage_size = LARGE_PAGE_SIZE;
      // printf ("# mmap-ed large page at %p\n", largepages);
      /* Try to write to the mapping so we bomb out right away if it can't be accessed */
      ((char *)largepages)[0] = 0;
    } else {
      // Commented out because it's spammy
      // perror("mmap failed");
    }
  }
#endif

  /* Find end of linked list. new_chunk points to the address where the
     address of the new chunk should be written, i.e., either
     new_chunk = &chunks, of new_chunk = &(last_chunk->next) where
     last_chunk is the last entry in the linked list. */
  struct largepage_chunk **new_chunk = &chunks;
  size_t last_offset = 0, last_size = 0;
  while (*new_chunk != NULL) {
    last_offset = (*new_chunk)->offset;
    last_size = (*new_chunk)->size;
    new_chunk = &(*new_chunk)->next;
  }

  /* See if there is enough memory in the large page left to serve the
     current request */
  if (last_offset + last_size + size > largepage_size) {
    /* If not, return memory allocated by malloc_pagealigned */
    void *ptr = malloc_pagealigned(size);
    // printf ("# Not enough large page memory available, returning malloc_pagealigned() = %p\n", ptr);
    return ptr;
  }

  ASSERT_ALWAYS(largepages != NULL);
#ifdef MAP_FAILED
  ASSERT_ALWAYS(largepages != MAP_FAILED);
#endif

  size_t pgsize = pagesize();
  ASSERT_ALWAYS(largepage_size % pgsize == 0);
  size_t round_up_size = ((size - 1) / pgsize + 1) * pgsize;
  *new_chunk = malloc(sizeof(struct largepage_chunk));
  ASSERT_ALWAYS(*new_chunk != NULL);
  (*new_chunk)->offset = last_offset + last_size;
  (*new_chunk)->size = round_up_size;
  (*new_chunk)->next = NULL;

  void *ptr = largepages + (*new_chunk)->offset;  
  // printf ("# Returning large-page memory at %p\n", ptr);
  return ptr;
}

void contiguous_free(const void *ptr)
{
  struct largepage_chunk **next = &chunks;
  while (*next != NULL) {
    struct largepage_chunk *chunk = *next;
    if (largepages + chunk->offset == ptr) {
      *next = chunk->next;
      // printf ("# Freeing large-page memory at %p\n", ptr);
      free(chunk);
#if defined(HAVE_MMAP) && defined(MAP_HUGETLB)
      if (chunks == NULL && largepage_size != 0) {
        // printf ("# Last chunk freed, unmapping large page\n");
        if (munmap(largepages, LARGE_PAGE_SIZE) != 0) {
          perror("munmap() failed");
        } else {
          largepage_size = 0;
          largepages = NULL;
        }
      }
#endif
      return;
    }
    next = &(chunk->next);
  }
  // printf ("# large-page memory at %p not found, calling free_pagealigned()\n", ptr);
  free_pagealigned(ptr);
}
