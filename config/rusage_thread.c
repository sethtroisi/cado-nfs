/* this file checks that getrusage (RUSAGE_THREAD, ...) works */

#define _GNU_SOURCE
#include <sys/time.h>
#include <sys/resource.h>

int
main ()
{
  struct rusage res[1];
  getrusage (RUSAGE_THREAD, res);
  return 0;
}
