#include "cado.h"
#include <gmp.h>
#include <stdio.h>

int main (int argc, char **argv)
{	
	if (argc != 2) {
		fprintf(stderr, "usage: %s integer\n", argv[0]);
		return 1;
	}

    mpz_t x;
    mpz_init_set_str (x,argv[1],10);
	if (mpz_probab_prime_p (x, 5))
		printf ("1\n");
	else
		printf ("0\n");

	mpz_clear (x);
	return 0;
}

