#include "cado.h"
#include <stdint.h>
#include <test_iter.h>
#include "macros.h"

void test_iter_init (test_iter_t iter, const unsigned long max)
{
  iter->next = 0;
  iter->max = max;
}

int test_iter_isdone (test_iter_t iter)
{
  ASSERT(iter->next <= iter->max);
  return iter->next == iter->max;
}

int64_t test_iter_int64_next(test_iter_t iter)
{
  unsigned long limit1 = (iter->max + 3) / 4;
  unsigned long limit2 = 2 * (iter->max - limit1) / 3;

  ASSERT(!test_iter_isdone(iter));
  
  unsigned long i = iter->next++;
  if (i < limit1)
    return INT64_MIN + (int64_t) i;
  i -= limit1;
  if (i < limit2)
    return -(int64_t) ((iter->max - limit1) / 3) + (int64_t) i;
  i -= limit2;
  return INT64_MAX - (int64_t) i;
}


uint64_t test_iter_uint64_next(test_iter_t iter)
{
  unsigned long limit1 = (iter->max + 1) / 2;

  ASSERT(!test_iter_isdone(iter));
  
  unsigned long i = iter->next++;
  if (i < limit1)
    return (uint64_t) i;
  i -= limit1;
  return UINT64_MAX - (int64_t) i;
}


#if 0
int main(int argc, char **argv)
{
  test_iter_t iter;
  test_iter_init(iter, atoi(argv[1]));
  while (!test_iter_isdone(iter)) {
    printf ("int64_t %lld\n", (long long) test_iter_int64_next(iter));
  }

  test_iter_init(iter, atoi(argv[1]));
  while (!test_iter_isdone(iter)) {
    printf ("uint64_t %llu\n", (unsigned long long) test_iter_uint64_next(iter));
  }

  return 0;
}
#endif
