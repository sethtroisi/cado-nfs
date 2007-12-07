#define _BSD_SOURCE
#include <stdlib.h>
#include <string.h>
#include "gf2x.h"

void random_wordstring(unsigned long *a, long n)
{
    for (long i = 0; i < n; i++)
	a[i] = random();
}

int main(int argc, char *argv[])
{
    if (argc != 2 && argc != 3) {
	fprintf(stderr, "Usage: ./test <n1> [<n2>]\n");
	exit(1);
    }

    long n1 = atol(argv[1]);
    long n2 = n1;
    if (argc == 3)
	n2 = atol(argv[2]);

    unsigned long *a = malloc(n1 * sizeof(unsigned long));
    unsigned long *b = malloc(n2 * sizeof(unsigned long));
    unsigned long *c = malloc((n1 + n2) * sizeof(unsigned long));
    unsigned long *d = malloc((n1 + n2) * sizeof(unsigned long));
    for (int i = 0; i < 100; i++) {
	srandom(1);
	random_wordstring(a, n1);
	random_wordstring(b, n2);
	mul_basecase(c, a, n1, b, n2);
	mul_gf2x(d, a, n1, b, n2);
	if (memcmp(c, d, (n1 + n2) * sizeof(unsigned long)) == 0) {
	    continue;
	}
	printf("test%d:=[", i);
	printf("[");
	int j;
	for (j = 0; j < n1; j++) {
	    if (j)
		printf(", ");
	    printf("%lu", a[j]);
	}
	printf("],");
	printf("[");
	for (j = 0; j < n2; j++) {
	    if (j)
		printf(", ");
	    printf("%lu", b[j]);
	}
	printf("],");
	printf("[");
	for (j = 0; j < (n1 + n2); j++) {
	    if (j)
		printf(", ");
	    printf("%lu", c[j]);
	}
	printf("],");
	printf("[");
	for (j = 0; j < (n1 + n2); j++) {
	    if (j)
		printf(", ");
	    printf("%lu", d[j]);
	}
	printf("]");
	printf("];\n");

	fprintf(stderr, "// Error on test %d\n", i);
	abort();
    }
    free(a);
    free(b);
    free(c);
    free(d);
    return 0;
}
