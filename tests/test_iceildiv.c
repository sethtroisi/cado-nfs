#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "macros.h"

static int
do_test(int i, int j, double ref)
{
  int icd = iceildiv(i,j);
  if (ceil(ref) != icd) {
    printf("iceildiv(%d, %d) = %d != ceil(%f)\n", i, j, icd, ref);
    return 1;
  }
  return 0;
}

int main()
{
  int i, j, had_error = 0;
  for (i=-3; i <= 3; i++)
    for (j=-3; j <= 3; j++)
      if (j != 0) {
        had_error |= do_test(i, j, (double) i / (double) j);
      }
  had_error |= do_test(INT_MAX, INT_MAX, 1.0);
  had_error |= do_test(INT_MIN, INT_MIN, 1.0);
  if (had_error) {
    exit(EXIT_FAILURE);
  } else {
    exit(EXIT_SUCCESS);
  }
}
