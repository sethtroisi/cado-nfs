#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

inline void addmul(uint64_t a, uint64_t b, uint64_t * w)
{
    unsigned int i;
#if 0
    /* Dans un sens */
    for (i = 0; i < 64; i++) {
	*w++ ^= b & -(a & 1);
	a >>= 1;
    }
#endif
#if 0
    /* Dans l'autre -- va plus vite. */
    for (i = 0; i < 64; i++) {
	*w++ ^= b & (((int64_t) a) >> 63);
	a <<= 1;
    }
#endif
#if 0
    /* Ã€ peu prÃ¨s comme la mÃ©thode 1, mais pas mieux */
    typedef uint64_t mvec_t[2];
    mvec_t mb[4] = {
	{0, 0}, {b, 0}, {0, b}, {b, b},
    };
    for (i = 0; i < 64; i += 2) {
	const uint64_t *y = mb[a & 3];
	*w++ ^= y[0];
	*w++ ^= y[1];
	a >>= 2;
    }
#endif
#if 1
    /* Avec des sse-2 */
    typedef uint64_t sse2_t __attribute__ ((vector_size(16)));
    sse2_t mb[4] = {
	(sse2_t) {0, 0},
	(sse2_t) {b, 0},
	(sse2_t) {0, b},
	(sse2_t) {b, b},
    };
    sse2_t *sw = (sse2_t *) w;
    for (i = 0; i < 64; i += 2) {
	*sw++ ^= mb[a & 3];
	a >>= 2;
    }
#endif
}


int main()
{
    clock_t t0, t1;

    uint64_t *a;
    uint64_t *b;
    uint64_t *w;

    unsigned int n = 2 * 1000 * 1000;

    int j;

    w = (uint64_t *) malloc(64 * sizeof(uint64_t));
    a = (uint64_t *) malloc(n * sizeof(uint64_t));
    b = (uint64_t *) malloc(n * sizeof(uint64_t));

    t0 = clock();

    for (j = 0; j < 16; j++) {
	unsigned int i;
	memset(w, 0, 64 * sizeof(uint64_t));
	for (i = 0; i < n; i++) {
	    addmul(a[i], b[i], w);
	}
    }

    t1 = clock() - t0;

    printf("%d times in %.4f s each\n", j, t1 * 1.0 / CLOCKS_PER_SEC / j);

    free(w);
    free(a);
    free(b);

    return 0;
}
