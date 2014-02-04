#include "cado.h"
#include <stdint.h>
#include <time.h>
#include "getprime.h"
#include "macros.h"

void
test_getprime ()
{
  unsigned long p;
  int i;

  /* there is no really something random we can test here, since
     getprime can start only at 2 */
  for (p = 2, i = 0; i < 1000000; p = getprime (p), i++);
  assert (p == 15485867);
  getprime (0);
  for (p = 2, i = 0; i < 1000000; p = getprime (p), i++);
  assert (p == 15485867);
  getprime (0);
}

int
main ()
{
  test_getprime ();
  exit (EXIT_SUCCESS);
}
