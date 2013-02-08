#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "rootfinder.h"
#include "portability.h"

int cmp(mpz_t * a, mpz_t * b)
{
    return mpz_cmp(*a, *b);
}

int main(int argc, char *argv[])
{
    mpz_t p;
    int i;
    int d;
    mpz_t *f;
    mpz_t *r;

    d = argc - 3;

    if (d < 1) {
	fprintf(stderr, "Usage: test-rootfinder p a_d [...] a_1 a_0\n");
	exit(1);
    }

    mpz_init_set_str(p, argv[1], 0);

    f = (mpz_t *) malloc((d + 1) * sizeof(mpz_t));
    r = (mpz_t *) malloc((d + 1) * sizeof(mpz_t));
    for (i = 0; i <= d; i++) {
	mpz_init_set_str(f[i], argv[2 + d - i], 0);
	mpz_init(r[i]);
    }

    int n = poly_roots(r, f, d, p);

    qsort(r, n, sizeof(mpz_t), (int(*)(const void *, const void*)) cmp);

    if (n > 0) {
	gmp_printf("%Zu", r[0]);
	for (i = 1; i < n; i++)
	    gmp_printf(" %Zu", r[i]);
	printf("\n");
    }
    for (i = 0; i <= d; i++) {
	mpz_clear(f[i]);
	mpz_clear(r[i]);
    }
    free(f);
    free(r);
    mpz_clear(p);
    return 0;
}

#if 0
// magma code for producing test cases.
s:=1.2;            
p:=10;
while p lt 2^200 do
    for i in [1..100] do
        p:=NextPrime(p);
        d:=Random([2..7]);
        coeffs:=[Random(GF(p)):i in [0..d]];
        F:=PolynomialRing(GF(p))!coeffs;
        printf "in %o", p;
        for c in Reverse(coeffs) do printf " %o", c; end for;
        printf "\n";
        r:=Sort([x[1]: x in Roots(F) ]);
        printf "out";
        for c in r do printf " %o", c; end for;
        printf "\n";
    end for;
    p := Ceiling(p*s);
end while;


#endif

