#include "gmp.h"
#include "stdio.h"

int main(int argc, char **argv)
{
	if (argc < 3) {
		fprintf(stderr, "usage: %s integer [factor factor ...]\n", argv[0]);
		return 1;
	}

	mpz_t N, product;

	mpz_init_set_str (N,argv[1],10);
	mpz_init_set_ui (product, 1);
	int i;
	for (i=2; i<argc; i++) {
		mpz_t factor;
		mpz_init_set_str (factor,argv[i],10);
		mpz_mul (product, factor, product);
		mpz_clear (factor);
	}

	printf ("%d\n", mpz_cmp (N, product));
	mpz_clear (N);
	mpz_clear (product);
	return 0;
}

