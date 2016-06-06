#include <stdio.h>
#include <stdlib.h>

int
main ()
{
  int i, v, count[2];

  srand (0);
  count[0] = count[1] = 0;
  for (i = 0; i < 1000; i++)
    {
      v = (rand () + i) & 1;
      count[v] ++;
    }
  if (count[0] == 0 || count[1] == 0)
    {
      printf ("bad random generator\n");
      return 1;
    }
  return 0;
}
