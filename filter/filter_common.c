#include "cado.h"
#include <stdio.h>
#include <stdlib.h>

#include "portability.h"
#include "utils.h"
#include "filter_common.h"

int
cmp_index (const void *p, const void *q)
{
  index_t x = *((index_t *)p);
  index_t y = *((index_t *)q);
  return (x <= y ? -1 : 1);
}

int
cmp_ideal_merge (const void *p, const void *q)
{
  ideal_merge_t x = *((ideal_merge_t *)p);
  ideal_merge_t y = *((ideal_merge_t *)q);
  return (x.id <= y.id ? -1 : 1);
}

/* We also compare x[1] and y[1] to make the code deterministic
   since in case x[0] = y[0], qsort() may give different results on
   different machines */
int
cmp_int2 (const void *p, const void *q)
{
  int *x = (int*) p;
  int *y = (int*) q;

  if (x[0] < y[0])
    return -1;
  else if (x[0] > y[0])
    return 1;
  else
    return (x[1] < y[1]) ? 1 : -1;
}
