#include <stdlib.h>
#include <stdio.h>
#include "macros.h"
#include "sieving_interval.h"

void sieving_interval_init(sieving_interval_ptr H, unsigned int t)
{
  ASSERT(t != 0);

  H->t = t;
  H->h = (unsigned int * ) malloc(sizeof(unsigned int) * t);
}

void sieving_interval_set_hi(sieving_interval_ptr H, unsigned int i, unsigned
                             int value)

{
  ASSERT(i < H->t);

  H->h[i] = value;
}

void sieving_interval_number_element(uint64_t * nb, sieving_interval_srcptr H)
{
  * nb = 1;
  for (unsigned int i = 0; i < H->t - 1; i++) {
    * nb = * nb * (2 * (uint64_t)H->h[i] + 1);
  }
  * nb = *nb * ((uint64_t)H->h[H->t - 1] + 1);

  ASSERT(* nb != 1);
}

void sieving_interval_clear(sieving_interval_ptr H)
{
  free(H->h);
  H->t = 0;
}

void sieving_interval_fprintf(FILE * filew, sieving_interval_srcptr H)
{
  for (unsigned int i = 0 ; i < H->t - 1; i++) {
    fprintf(filew, "-H%u: -%u -- H%u: %u\n", i, H->h[i], i, H->h[i]);
  }
  fprintf(filew, "H%u: 0 -- H%u: %u\n", H->t - 1, H->t - 1, H->h[H->t - 1]);
}
