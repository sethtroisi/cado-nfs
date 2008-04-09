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

void bit_transpose_mat(unsigned long *dst, const unsigned long *src)
{
    int i, j;
    for (i = 0; i < 64; i++) {
	dst[i] = 0;
	for (j = 0; j < 64; j++) {
	    dst[i] ^= ((src[j] >> i) & 1UL) << j;
	}
    }
}

void mul_vec_mat(unsigned long *C,
		 const unsigned long *A,
		 const unsigned long *B, unsigned long m)
{
    int i;
    unsigned long j;
    memset(C, 0, m * sizeof(unsigned long));

    j = 0;

#if 1				/* a la main */

#if 1				/* sse */
    typedef uint64_t sse2_t __attribute__ ((vector_size(16)));

    sse2_t *Cw = (sse2_t *) C;
    sse2_t *Aw = (sse2_t *) A;

    for (j = 0; j < m; j += 2) {
	sse2_t c = { 0, 0 };
	sse2_t a = *Aw++;

	sse2_t one = { 1, 1, };
#if 1
#define SHR(x,r) (sse2_t)__builtin_ia32_psrlqi128   ((x),(r))
	for (i = 0; i < 64; i++) {
	    sse2_t bw = { B[i], B[i], };

	    c ^= (bw & -(a & one));
	    a = SHR(a, 1);
	}
#else
#endif
	*Cw++ = c;
    }
    C += j;
    A += j;
#endif
    for (; j < m; j++) {
	unsigned long c = 0UL;
	unsigned long a = *A++;
#if 1
	for (i = 0; i < 64; i++) {
	    c ^= (B[i] & -(a & 1UL));
	    a >>= 1UL;
	}
#else
	for (i = 64 - 1; i >= 0; i--) {
	    c ^= (B[i] & (((long) a) >> (64 - 1)));
	    a <<= 1UL;
	}
#endif
	*C++ = c;
    }


#else				/* parity */

    unsigned long *tb = malloc(64 * sizeof(unsigned long));

    bit_transpose_mat(tb, B);

    for (j = 0; j < m; j++) {
	unsigned long a = *A++;
	unsigned long c = 0UL;
	for (i = 0; i < 64; i++) {
	    c ^= (((unsigned long) __builtin_parityl(a & tb[i])) << i);
	}
	*C++ = c;
    }

    free(tb);

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

#if 0
    t0 = clock();

    for (j = 0; j < 4; j++) {
	unsigned int i;
	memset(w, 0, 64 * sizeof(uint64_t));
	for (i = 0; i < n; i++) {
	    addmul(a[i], b[i], w);
	}
    }

    t1 = clock() - t0;

    printf("%d times in %.4f s each\n", j, t1 * 1.0 / CLOCKS_PER_SEC / j);
#endif



    t0 = clock();

    for (j = 0; j < 16; j++) {
	mul_vec_mat(a, b, w, n);
    }

    t1 = clock() - t0;

    printf("mul_vec_mat : %d times in %.4f s each\n", j,
	   t1 * 1.0 / CLOCKS_PER_SEC / j);

    free(w);
    free(a);
    free(b);

    return 0;
}
