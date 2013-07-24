#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

// Note: we use GMP for primality test and random functions; not really
// for multiprecision arithmetic.

void usage(const char *argv0) {
  fprintf(stderr, "usage: ./%s <bitsize> <mod12_class>\n", argv0);
  fprintf(stderr, "  bitsize must be in [15..60]\n");
  fprintf(stderr, "  mod12class must be in {1,5,7,11}\n");
  exit(EXIT_FAILURE);
}

int main(int argc, char **argv) {
  int bitsize, mod12;
  if (argc != 3)
    usage(argv[0]);
  bitsize = atoi(argv[1]);
  mod12 = atoi(argv[2]);
  if (bitsize < 15 || bitsize > 60 || 
      !(mod12 == 1 || mod12 == 5 || mod12 == 7 || mod12 == 11))
    usage(argv[0]);

  gmp_randstate_t state;
  gmp_randinit_default(state);
  mpz_t p, r;
  mpz_init(p);
  mpz_init(r);
  
  mpz_set_ui(p, 1);
  mpz_mul_2exp(p, p, bitsize-1);

  for (int i = 0; i < 5000; ++i) {
    if (bitsize <= 19) {
      // exhaustive list
      do {
        mpz_nextprime(p, p);
      } while (mpz_mod_ui(r, p, 12) != mod12);
      if (mpz_sizeinbase(p, 2) > bitsize)
        break;
      gmp_printf("%Zd\n", p);
    } else {
      // randomize
      do {
        do {
          mpz_urandomb(p, state, bitsize);
        } while (mpz_sizeinbase(p, 2) != bitsize);
        mpz_nextprime(p, p);
      } while (mpz_sizeinbase(p, 2) != bitsize || mpz_mod_ui(r, p, 12) != mod12);
      gmp_printf("%Zd\n", p);
    }
  }
  mpz_clear(p);
  mpz_clear(r);
  gmp_randclear(state);
  return EXIT_SUCCESS;
}


  


