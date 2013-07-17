#ifndef __VARLIST_H
#define __VARLIST_H

/*
   Macros for variable-length lists with iterator protocol.
*/

/* Need size_t */
#include <stdlib.h>


/*
    Type declaration for a variable-length list. This produces an anonymous
    struct. Use as, e.g.,

    VARLIST_DECLARE(int) buffer;
*/

#define VARLIST_DECLARE(__type) \
  struct { \
    __type *entries; \
    size_t alloc, len, typesize; \
  } \


/*
    Initialise a variable-length list to an empty list.
    Required, or VARLIST_APPEND() will write to a random memory location.
*/

#define VARLIST_INIT(__list) \
  do { \
    (__list).entries = NULL; \
    (__list).alloc = 0; \
    (__list).len = 0; \
    (__list).typesize = sizeof((__list).entries[0]); \
  } while (0) \


/*
    Append a value to a variable-length list. If the end of allocated memory is
    reached, doubles the allocation size. The first alloc is for 16 entries.
*/

#define VARLIST_APPEND(__list, __value) \
  do { \
    if ((__list).alloc == (__list).len) { \
      if ((__list).alloc == 0) \
        (__list).alloc = 16; \
      else \
        (__list).alloc *= 2; \
      (__list).entries = realloc((__list).entries, \
                                 (__list).alloc * (__list).typesize); \
      ASSERT_ALWAYS((__list).entries != NULL); \
    } \
    (__list).entries[(__list).len++] = (__value); \
  } while (0) \


/*
    Set one entry of a variable-length list. The index must be within existing
    elements of the list, or right at the end in which case the entry is
    appended.
*/

#define VARLIST_SET(__list, __index, __value) \
  do { \
    ASSERT((__index) <= (__list).len); \
    if ((__index) == (__list).len) { \
      VARLIST_APPEND(__list, __value); \
    } else { \
      ((__list).entries[__index]) = (__value); \
    } \
  } while (0) \


/* Get the entry at a given index from the list */

#define VARLIST_GET(__list, __index) \
  ((__list).entries[__index]) \


/* Free the memory allocated by this list */

#define VARLIST_CLEAR(__list) \
  do { \
    free((__list).entries); \
    (__list).entries = NULL; \
    (__list).alloc = 0; \
    (__list).len = 0; \
    (__list).typesize = 0; \
  } while (0)\


/*
    Type declaration for an iterator over a variable-length list of
    type __type. The type specified in __type here must agree with the
    type specified in VARLIST_DECLARE(), e.g.,

    VARLIST_DECLARE(int) buffer;
    VARLIST_ITERATOR_DECLARE(int) iterator;
*/

#define VARLIST_ITERATOR_DECLARE(__type) \
  struct { \
    __type *ptr; \
    size_t next, len; \
  } \


/* 
    Initialise an iterator to iterate over the entries of a variable-length
    list. Note that the length of the list is assumed to be constant during
    the life-time of the iterator.
*/

#define VARLIST_ITERATOR_INIT(__iter, __list) \
  do { \
    (__iter).ptr = (__list).entries; \
    (__iter).next = 0; \
    (__iter).len = (__list).len; \
  } while (0) \


/*
   Get next element from the iterator. The expression can be use as an lvalue.
*/

#define VARLIST_ITERATOR_NEXT(__iter) \
    ((__iter).ptr[(__iter).next++])


/* Test whether the iterator has reached the end of the variable-length list */

#define VARLIST_ITERATOR_FINISHED(__iter) \
    ((__iter).next >= (__iter).len)


/* Clear the iterator. Not strictly required, but good practice to avoid stale
   pointers to memory. */

#define VARLIST_ITERATOR_CLEAR(__iter) \
  do { \
    (__iter).ptr = NULL; \
    (__iter).next = 0; \
    (__iter).len = 0; \
  } while (0) \


/* Call a function __f(__data, entry) for each entry in the variable-length 
   list. The __data parameter is simply passed through and can be used to 
   transport auxiliary information. */
#define VARLIST_MAP(__type, __list, __f, __data) \
  do { \
    VARLIST_ITERATOR_DECLARE(__type) __iter; \
    VARLIST_ITERATOR_INIT(__iter, __list); \
    while (!VARLIST_ITERATOR_FINISHED(__iter)) { \
      (__f)((__data), VARLIST_ITERATOR_NEXT(__iter)); \
    } \
    VARLIST_ITERATOR_CLEAR(__iter); \
  } while (0) \


#if 0
/* Example code: */

#include "cado.h"
#include <stdio.h>

#include "macros.h"
#include "varlist.h"


int main() {
  /* Declare a list of int elements */
  VARLIST_DECLARE(int) intbuf;
  /* Declare a pointer to a list of strings */
  VARLIST_DECLARE(char *) *strbuf;
  /* Declare an iterator over a list of strings */
  VARLIST_ITERATOR_DECLARE(char *) striter;

  /* Simple example, varlist storage is allocated in declaration */
  VARLIST_INIT(intbuf);
  /* Append some int values */
  VARLIST_APPEND(intbuf, 1);
  VARLIST_APPEND(intbuf, 2);
  VARLIST_SET(intbuf, 2, 42);
  VARLIST_SET(intbuf, 2, 3);
  ASSERT_ALWAYS(VARLIST_GET(intbuf, 2) == 3);
  /* Done */
  VARLIST_CLEAR(intbuf);

  /* Varlist storage is allocated with malloc and is passed around as pointer */
  strbuf = malloc(sizeof(*strbuf));
  VARLIST_INIT(*strbuf);
  /* Append some strings */
  VARLIST_APPEND(*strbuf, "a");
  VARLIST_APPEND(*strbuf, "b");

  /* Iterate over the list entries. We simply print them */
  VARLIST_ITERATOR_INIT(striter, *strbuf);
  while (!VARLIST_ITERATOR_FINISHED(striter))
    printf ("%s\n", VARLIST_ITERATOR_NEXT(striter));
  VARLIST_ITERATOR_CLEAR(striter);

  /* The same functionality as the iterator above, but using the map 
     construct */
  VARLIST_MAP(char *, *strbuf, &printf, "%s\n");

  VARLIST_CLEAR(*strbuf);
  free(strbuf);

  return 0;
}

#endif

#endif /* __VARLIST_H */
