#include <sys/types.h>    /* for cputime */
#include <sys/resource.h> /* for cputime */

/* Common timing function: return user time since beginning of session
   in milliseconds. Example:

   int st = cputime ();
   foo (...);
   printf ("foo took %dms\", cputime () - st);
*/
int
cputime (void)
{
  struct rusage rus;

  getrusage (0, &rus);
  return rus.ru_utime.tv_sec * 1000 + rus.ru_utime.tv_usec / 1000;
}

