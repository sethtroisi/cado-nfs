#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <errno.h>

void usage()
{
	fprintf(stderr, "Usage: ./dec2bin <field width> [ files... ]\n");
	exit(1);
}

void dec2bin(unsigned int w, FILE ** io, unsigned int n)
{
	mpz_t z, pad, r;


	mpz_init(z);
	mpz_init(r);
	mpz_init(pad);

	mpz_set_ui(pad, 1);
	mpz_mul_2exp(pad, pad, w);
	
	for(;;) {
		unsigned int f;
		int stop = 0;
		for(f = 0 ; f < n ; f++) {
			unsigned int k;
			if (feof(io[f])) { stop++; continue; }
			gmp_fscanf(io[f], "%Zd", z);
			if (mpz_sgn(z) < 0) mpz_add(z, z, pad);
			for(k = 0 ; k < w ; k+= 8) {
				mp_limb_t foo;
				mpz_fdiv_r_2exp(r, z, 8);
				foo = mpz_get_ui(r);
				fwrite(&foo, 1, 1, stdout);
				mpz_fdiv_q_2exp(z, z, 8);
			}
		}
		if (stop != 0 && stop != n) {
			fprintf(stderr, "bad stride\n");
			abort();
		}
		if (stop == n)
			break;
	}

	mpz_clear(r);
	mpz_clear(z);
	mpz_clear(pad);
}

void bin2dec(unsigned int w, FILE ** io, unsigned int n)
{
	mpz_t z,a;

	mpz_init(z);
	mpz_init(a);

	for(;;) {
		unsigned int f;
		int stop = 0;
		for(f = 0 ; f < n ; f++) {
			unsigned int k;
			if (feof(io[f])) { stop++; continue; }
			mpz_set_ui(z,0);
			for(k = 0 ; k < w ; k+= 8) {
				mp_limb_t foo = 0;
				fread(&foo, 1, 1, io[f]);
				mpz_set_ui(a,foo);
				mpz_mul_2exp(a, a, k);
				mpz_add(z, z, a);
			}
			gmp_printf("%Zd%c", z, f == (n-1) ? '\n' : ' ');
		}
		if (stop != 0 && stop != n) {
			fprintf(stderr, "bad stride\n");
			abort();
		}
		if (stop == n)
			break;
	}

	mpz_clear(a);
	mpz_clear(z);
}

int main(int argc, char * argv[])
{
	unsigned int w;
	FILE ** ios;
	unsigned int n;

	if (argc > 2) {
		if (strcmp(argv[1], "--dec2bin") == 0) {
			argv++, argc--;
		} else if (strcmp(argv[1], "--bin2dec") == 0) {
			argv++, argc--;
		}
	}

	if (argc < 2) {
		usage();
	}
	w = atoi(argv[1]);
	if (w == 0 || (w % 32 != 0)) {
		usage();
	}

	if (w % mp_bits_per_limb != 0) {
		usage();
	}

	if (argc == 2) {
		n = 1;
		ios = malloc(n * sizeof(FILE*));
		ios[0] = stdin;
	} else {
		unsigned int i;
		n = argc - 2;
		ios = malloc(n * sizeof(FILE*));
		for(i = 0 ; i < n ; i++) {
			ios[i] = fopen(argv[i+2], "r");
			if (ios[i] == NULL) {
				fprintf(stderr, "%s: %s\n", 
						argv[i+2], strerror(errno));
			}
		}
	}

	if (strstr(argv[0], "bin2dec") != NULL) {
		bin2dec(w, ios, n);
	} else {
		dec2bin(w, ios, n);
	}
	return 0;
}
