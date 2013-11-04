#include "cado.h"
#include <stdio.h>
#include "portability.h"
#include "utils.h"
#include "filter_memalloc.h"

/* relcompact_list is a list of blocks, each one of BLOCK_SIZE index_ts */
static index_t **relcompact_list = NULL;
static index_t *myrelcompact;
static unsigned int relcompact_used = BLOCK_SIZE; /* usage of current block */
static unsigned int relcompact_size = 0;  /* minimal size of relcompact_list */
static int relcompact_current = -1; /* index of current block */
static size_t my_malloc_bytes = 0;

/* return a pointer to an array of n (index_t) */
index_t *
my_malloc (unsigned int n)
{
  index_t *ptr;

  if (relcompact_used + n > BLOCK_SIZE) 
  {
    relcompact_used = 0;
    if (((unsigned int) (++relcompact_current)) == relcompact_size) 
    {
      my_malloc_bytes -= relcompact_size * sizeof (index_t *);
      relcompact_size = relcompact_size ? (relcompact_size << 1) : (1<<16);
      size_t tmp_size = relcompact_size * sizeof(index_t *);
      if (!(relcompact_list = (index_t **) realloc(relcompact_list, tmp_size))) 
      {
        fprintf (stderr, "my_malloc_int: realloc error: %s\n",strerror(errno));
        exit (1);
      }
      my_malloc_bytes += tmp_size;
    }
    SMALLOC(relcompact_list[relcompact_current], BLOCK_SIZE, "my_malloc_int 1");
    myrelcompact = relcompact_list[relcompact_current];
    my_malloc_bytes += BLOCK_SIZE * sizeof (index_t);
  }
  ptr = &(myrelcompact[relcompact_used]);
  relcompact_used += n;
  return ptr;
}

void
my_malloc_free_all (void)
{
  for ( ; relcompact_current >= 0; relcompact_current--)
  {
    SFREE(relcompact_list[relcompact_current]);
  }
  SFREE(relcompact_list);
  relcompact_list = NULL;
  relcompact_used = BLOCK_SIZE;
  relcompact_size = 0;
}

inline size_t 
get_my_malloc_bytes ()
{
  return my_malloc_bytes;
}
