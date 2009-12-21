#include "gmp.h"
#include "stdio.h"

int main(int argc, char **argv)
{
	if (argc < 3) {
		fprintf(stderr, "usage: %s integer [factor factor ...]\n", argv[0]);
		return 1;
	}

	mpz_t N, product_primes;

	mpz_init_set_str (N,argv[1],10);
	mpz_init_set_ui (product_primes, 1);
	int i;
	for (i=2; i<argc; i++) {
		mpz_t prime_factor;
		mpz_init_set_str (prime_factor,argv[i],10);
		mpz_mul (product_primes, prime_factor, product_primes);
		mpz_clear (prime_factor);
	}

	printf ("%d\n", mpz_cmp (N, product_primes));
	mpz_clear (N);
	mpz_clear (product_primes);
	return 0;
}

