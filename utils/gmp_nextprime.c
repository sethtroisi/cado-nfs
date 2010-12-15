#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"

int
main (int argc, char *argv[])
{
  unsigned long p;
  if (argc != 2) {
      fprintf(stderr, "Usage: nextprime p\n");
      exit(1);
  }

  p = strtoul (argv[1], NULL, 10);
  mpz_t P;

  mpz_init_set_ui (P, p);
  mpz_nextprime (P, P);
  gmp_printf ("%Zd\n", P);
  mpz_clear (P);
  return 0;
}
