/* random_integer d --- generate a random integer of d digits
   random_integer d1 d2 --- generate a random integer of d1 to d2 digits */

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include "gmp.h"
#include "portability.h"

void usage_die(const char *argv0)
{
  fprintf (stderr, "Usage: %s [-p] d1 [d2]\n", argv0);
  fprintf (stderr, "       Generate a random integer of d1 digits.\n");
  fprintf (stderr, "       (resp. of d1 to d2 digits if d2 is given).\n");
  fprintf (stderr, "       With the -p option, generate prime integers.\n");
  exit (EXIT_FAILURE);
}

int
main (int argc, char *argv[])
{
  unsigned long d, d1, d2;
  gmp_randstate_t state;
  mpz_t n, b;
  int prime = 0;
  char * argv0 = argv[0];
 
  if (argc == 1)
    usage_die(argv0);
  if (argv[1][0] == '-') {
    if (argv[1][1] != 'p')
      usage_die(argv0);
    prime = 1;
    argc--;
    argv++;
  }
  if (argc != 2 && argc != 3)
    usage_die(argv0);

  d1 = strtoul (argv[1], NULL, 0);
  d2 = (argc == 2) ? d1 : strtoul (argv[2], NULL, 0);

  /* first choose d1 <= d <= d2 */
  d = d1 + (getpid () % (d2 - d1 + 1));

  gmp_randinit_default (state);
  gmp_randseed_ui (state, getpid ());
  mpz_init (b);
  mpz_init (n);
  mpz_ui_pow_ui (b, 10, d);
  do {
      do
          mpz_urandomm (n, state, b);
      while (mpz_sizeinbase (n, 10) != d);
  } while (prime && !mpz_probab_prime_p(n, 25));
  gmp_printf ("%Zd\n", n);
  mpz_clear (n);
  mpz_clear (b);
  gmp_randclear (state);

  return 0;
}
