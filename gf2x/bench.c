#include <sys/resource.h>
#include <sys/types.h>
#include <stdlib.h>
#include "gf2x.h"

/* Usage: ./bench <limit size> */

double runtime(void)
{
   struct rusage used; 

   getrusage(RUSAGE_SELF, &used);
   return (used.ru_utime.tv_sec + used.ru_stime.tv_sec +
      (used.ru_utime.tv_usec + used.ru_stime.tv_usec) / 1e6);
}

int main(int argc, char * argv[])
{
    unsigned int i;
    unsigned long n0 = 1;
    unsigned long n = 200;
    if (argc > 1) {
        n = atoi(argv[1]);
    }
    if (argc > 2) {
        n0 = n;
        n = atoi(argv[2]);
    }

    setbuf(stdout,NULL);

    unsigned long * a, * b, * c;

    a = malloc(n * sizeof(unsigned long));
    b = malloc(n * sizeof(unsigned long));
    c = malloc(2 * n * sizeof(unsigned long));

    for(i = 0 ; i < n ; i++) {
        a[i] = rand();
        b[i] = rand();
    }

    unsigned int jmax = 100000;

    for(i = n0 ; i < n ; i++) {
        /* warm up */
        mul_gf2x(c,a,i,b,i);
        mul_gf2x(c,a,i,b,i);
        mul_gf2x(c,a,i,b,i);
        double t = runtime();
        unsigned int j;
        double d;
#if 0
        for(j = 0 ; j < jmax ; j++) {
            mul_gf2x(c,a,i,b,i);
        }
        d = runtime() - t;

        if (d > 1 && jmax > 8) {
            jmax -= jmax / 3;
        }
#else
        for(j = 0 ; (d=runtime()-t) < 1 && (j < 2000 || d < 0.2) ; j++) {
            mul_gf2x(c,a,i,b,i);
        }
#endif
        printf("%u %e\t%u %.2f\n", i, d/j,j,d);

    }

    free(a);
    free(b);
    free(c);
}

