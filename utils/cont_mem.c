
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
#ifdef HAVE_LINUX_BINFMTS_H
/* linux/binfmts.h defines MAX_ARG_STRLEN in terms of PAGE_SIZE, but does not
   include a header where PAGE_SIZE is defined, so we include sys/user.h
   as well. */
#include <sys/user.h>
#include <linux/binfmts.h>
#endif

#ifdef HAVE_SYS_MMAN_H
#include <sys/mman.h>
#endif

#include "macros.h"
#include "portability.h"
#include "misc.h"
#include "cont_mem.h"

#ifndef LARGE_PAGE_SIZE
#define LARGE_PAGE_SIZE (2UL*1024*1024)
#endif

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

#ifdef HAVE_MMAP
  /* Try to mmap a large page, if we haven't already */
  if (largepages == NULL) {
    largepages = mmap (NULL, LARGE_PAGE_SIZE, PROT_READ|PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB, -1, 0);
    if (largepages != MAP_FAILED) {
      largepage_size = LARGE_PAGE_SIZE;
      printf ("# mmap-ed large page at %p\n", largepages);
      /* Try to write to the mapping so we bomb out right away if it can't be accessed */
      ((char *)largepages)[0] = 0;
    } else {
      perror("mmap failed");
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
    printf ("# Not enough large page memory available, returning malloc_pagealigned() = %p\n", ptr);
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
  printf ("# Returning large-page memory at %p\n", ptr);
  return ptr;
}

void contiguous_free(const void *ptr, const size_t size)
{
  struct largepage_chunk **next = &chunks;
  while (*next != NULL) {
    struct largepage_chunk *chunk = *next;
    if (largepages + chunk->offset == ptr) {
      *next = chunk->next;
      printf ("# Freeing large-page memory at %p\n", ptr);
      free(chunk);
      if (chunks == NULL) {
        printf ("# Last chunk freed, unmapping large page\n");
        if (munmap(largepages, LARGE_PAGE_SIZE) != 0) {
          perror("munmap() failed");
        }
      }
      return;
    }
    next = &(chunk->next);
  }
  printf ("# large-page memory at %p not found, calling free_pagealigned()\n", ptr);
  free_pagealigned(ptr, size);
}
