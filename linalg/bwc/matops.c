#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <emmintrin.h>

/* This is **only** for 64 * 64 matrices */

#define WBITS   64
typedef uint64_t mat64[64];
typedef uint64_t * mat64_ptr;
typedef const uint64_t * mat64_srcptr;

#define mul_6464_6464 mul_6464_6464_sse
#define add_6464_6464 add_6464_6464_C
#define mul_N64_6464 mul_N64_6464_sse
#define mul_N64_T6464 mul_N64_6464_vec
#define addmul_o64_64o addmul_o64_64o_lsb_sse_v1
#define mul_o64_6464 mul_o64_6464_C_lsb
#define mul_o64_T6464 mul_o64_T6464_C_parity

/* level 1 */

static inline void mul_6464_6464_sse(mat64_ptr C, mat64_srcptr A, mat64_srcptr B)
{
    int i;
    memset(C, 0, sizeof(mat64));
 
    /* As per emmintrin.h, only the __m128i type is declared with
     * attribute __may_alias__, meaning that accesses through this
     * pointer may alias other types (e.g. uint64_t's, in our cases). For
     * this reason, we _must_ use this pointer type, and not __v2di, for
     * our pointer types. Reading from a __m128 * into a __v2di, or the
     * converse, are legal operations.
     */
    __m128i *Cw = (__m128i *) C;
    __m128i *Aw = (__m128i *) A;

    for (int j = 0; j < 64; j += 2) {
	__v2di c = { 0, 0 };
	__v2di a = *Aw++;

	__v2di one = { 1, 1, };
#define SHR(x,r) _mm_srli_epi64((x),(r))
	for (i = 0; i < 64; i++) {
	    __v2di bw = { B[i], B[i], };

	    c ^= (bw & -(a & one));
	    a = SHR(a, 1);
	}
#undef  SHR
	*Cw++ = c;
    }
}

static inline void add_6464_6464_C(mat64_ptr C, mat64_srcptr A, mat64_srcptr B)
{
    for (int j = 0; j < 64; j++) {
        C[j] = A[j] ^ B[j];
    }
}


static inline void addmul_o64_64o_lsb(uint64_t a, uint64_t b, uint64_t * w)
{
    /* Dans un sens */
    for (unsigned int i = 0; i < 64; i++) {
	*w++ ^= b & -(a & 1);
	a >>= 1;
    }
}

static inline void addmul_o64_64o_msb(uint64_t a, uint64_t b, uint64_t * w)
{
    /* Dans l'autre -- va un poil plus vite. */
    for (unsigned int i = 0; i < 64; i++) {
	*w++ ^= b & (((int64_t) a) >> 63);
	a <<= 1;
    }
}

static inline void addmul_o64_64o_lsb_sse_v0(uint64_t a, uint64_t b, uint64_t * w)
{
    /* À peu près comme la méthode 1, mais pas mieux */
    typedef uint64_t mvec_t[2];
    mvec_t mb[4] = {
	{0, 0}, {b, 0}, {0, b}, {b, b},
    };
    for (int i = 0; i < 64; i += 2) {
	const uint64_t *y = mb[a & 3];
	*w++ ^= y[0];
	*w++ ^= y[1];
	a >>= 2;
    }
}

static inline void addmul_o64_64o_lsb_sse_v1(uint64_t a, uint64_t b, uint64_t * w)
{
    /* Avec des sse-2 */
    __v2di mb[4] = {
	(__v2di) {0, 0},
	(__v2di) {b, 0},
	(__v2di) {0, b},
	(__v2di) {b, b},
    };
    __m128i *sw = (__m128i *) w;
    for (int i = 0; i < 64; i += 2) {
	*sw++ ^= mb[a & 3];
	a >>= 2;
    }
}

static inline void mul_o64_6464_C_lsb(uint64_t * w, uint64_t a, mat64_srcptr b)
{
    unsigned long c = 0;
    for (unsigned int i = 0; i < WBITS; i++) {
	c ^= (b[i] & -(a & 1UL));
	a >>= 1UL;
    }
    *w = c;
}

static inline void mul_o64_6464_C_msb(uint64_t *w,
                   uint64_t a,
                   mat64_srcptr B)
{
    uint64_t c = 0UL;
    for (int i = 64 - 1; i >= 0; i--) {
        c ^= (B[i] & (((long) a) >> (64 - 1)));
        a <<= 1UL;
    }
    *w = c;
}


static inline void mul_o64_T6464_C_parity(uint64_t * w, uint64_t a, mat64_srcptr b)
{
    // Uses unoptimized __builtin function -- maybe with gcc 4.3 is beter
    unsigned long c = 0;
    for (unsigned int i = 0; i < WBITS; i++) {
	c ^= ((__builtin_parityl(a & b[i])/* & 1UL*/) << i);
    }
    *w = c;
}

static inline void transp_6464(mat64_ptr dst, mat64_srcptr src)
{
    int i, j;
    for (i = 0; i < 64; i++) {
	dst[i] = 0;
	for (j = 0; j < 64; j++) {
	    dst[i] ^= ((src[j] >> i) & 1UL) << j;
	}
    }
}

static inline void copy_6464(mat64_ptr dst, mat64_srcptr src)
{
    memcpy(dst, src, sizeof(mat64));
}


/* level 2 */

static inline void copy_N64(uint64_t * dst, const uint64_t * src, unsigned long m)
{
    memcpy(dst, src, m * sizeof(uint64_t));
}

static inline int cmp_N64(const uint64_t * dst, const uint64_t * src, unsigned long m)
{
    return memcmp(dst, src, m * sizeof(uint64_t));
}


static inline void mul_N64_6464_vec(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, unsigned long m)
{

    memset(C, 0, m * sizeof(unsigned long));
    for (unsigned long i = 0; i < m; i++) {
        mul_o64_6464(C++, *A++, B);
    }
}

static inline void mul_N64_T6464_vec(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, unsigned long m)
{

    memset(C, 0, m * sizeof(unsigned long));
    for (unsigned long i = 0; i < m; i++) {
        mul_o64_T6464(C++, *A++, B);
    }
}

static inline void mul_N64_6464_transB(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, unsigned long m)
{
    uint64_t *tb = malloc(64 * sizeof(uint64_t));
    transp_6464(tb, B);
    mul_N64_T6464(C, A, B, m);
    free(tb);
}


static inline void mul_N64_6464_sse(uint64_t *C,
		 const uint64_t *A,
		 const uint64_t *B, unsigned long m)
{
    unsigned long j;
    memset(C, 0, m * sizeof(uint64_t));

    __m128i *Cw = (__m128i *) C;
    __m128i *Aw = (__m128i *) A;

    for (j = 0; j < m; j += 2) {
	__v2di c = { 0, 0 };
	__v2di a = *Aw++;

	__v2di one = { 1, 1, };
#define SHR(x,r) _mm_srli_epi64((x),(r))
	for (int i = 0; i < 64; i++) {
	    __v2di bw = { B[i], B[i], };

	    c ^= (bw & -(a & one));
	    a = SHR(a, 1);
	}
#undef  SHR
	*Cw++ = c;
    }
    C += j;
    A += j;
    for (; j < m; j++) {
	uint64_t c = 0UL;
	uint64_t a = *A++;
	for (int i = 0; i < 64; i++) {
	    c ^= (B[i] & -(a & 1UL));
	    a >>= 1UL;
	}
	*C++ = c;
    }
}

static inline void mul_64N_N64_addmul(uint64_t *w, uint64_t *a, uint64_t *b, unsigned long n)
{
    memset(w, 0, 64 * sizeof(uint64_t));
    for (unsigned long i = 0; i < n; i++) {
        addmul_o64_64o(a[i], b[i], w);
    }
}


static inline void mul_TN32_N64_C(uint64_t * b, uint32_t * A, uint64_t * x, unsigned int ncol)
{
    uint32_t idx, i, rA;
    uint64_t rx;

    memset(b, 0, 32 * sizeof(uint64_t));
    for(idx = 0; idx < ncol; idx++) {
        rA = A[idx];
        rx = x[idx];
        for(i = 0; i < 32; i++) {
            b[i] ^= rx & -(rA & 1);
            rA >>= 1;
        }
    }
}

static inline void mul_TN64_N64_C(uint64_t * b, uint64_t * A, uint64_t * x, unsigned int ncol)
{
    uint64_t idx, i, rA;
    uint64_t rx;

    memset(b, 0, 64 * sizeof(uint64_t));
    for(idx = 0; idx < ncol; idx++) {
        rA = A[idx];
        rx = x[idx];
        for(i = 0; i < 64; i++) {
            b[i] ^= rx & -(rA & 1);
            rA >>= 1;
        }
    }
}




#if 0   /* haven't checked yet what the funny-named functions actually do... *//*{{{*/


static inline void TVUBit_v2(unsigned long m,
	       unsigned long n,
	       const unsigned long *A, const unsigned long *B,
	       unsigned long *C)
{
    unsigned long i, P, k;

    memset(C, 0, WBITS * sizeof(unsigned long));

    for (i = 0; i < m; i++) {

	P = *A++;
	for (k = 0; k < WBITS; k++) {

	    //if (P & 1UL) C[k] ^= B[i];

	    C[k] ^= (B[i] & -(P & 1UL));

	    P >>= 1UL;
	}

    }


}


#if 1
static inline void VUBit_v2(unsigned long m,
	      unsigned long n,
	      const unsigned long *A, const unsigned long *B,
	      unsigned long *C)
{
    unsigned long i, P, k;
    memset(C, 0, m * sizeof(unsigned long));
    for (i = 0; i < m; i++) {
	P = *A++;
	for (k = 0; k < n; k++) {
	    C[i] ^= (B[k] & -(P & 1UL));
	    P >>= 1UL;
	}

    }


}
#endif


#if 1

#endif




#if 0


static inline void VUBit(unsigned long m,
	   unsigned long n,
	   const unsigned long *A, const unsigned long *B, unsigned long *C)
{

#if 1



#endif


}




#endif











#if 0
static inline void VUBit(unsigned long m,
	      unsigned long n,
	      const unsigned long *A, const unsigned long *B,
	      unsigned long *C)
{
    unsigned long i, P, k, j;
    memset(C, 0, m * sizeof(unsigned long));
    for (i = 0; i < m; i++) {
	P = *A++;
	j = C[i];
	for (k = 0; k < n; k++) {
	    j ^= (B[k] & -(P & 1UL));
	    P >>= 1UL;
	}
	C[i] = j;
    }
}
#endif

#endif/*}}}*/


/* polynomials */

/* lengths of a1 and a2 are n1 and n2 */
typedef uint64_t (*m64pol_ptr)[64];
typedef uint64_t (*const m64pol_srcptr)[64];

void m64pol_addmul(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2)
{
    assert(r != a1 && r != a2);
    for(unsigned int i = 0 ; i < n1 ; i++) {
        for(unsigned int j = 0 ; j < n2 ; j++) {
            mat64 x;
            mul_6464_6464(x, a1[i], a2[j]);
            add_6464_6464(r[i+j], r[i+j], x);
        }
    }
}

void m64pol_add(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n)
{
    for(unsigned int i = 0 ; i < n ; i++) {
        add_6464_6464(r[i], a1[i], a2[i]);
    }
}

void m64pol_mul(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2)
{
    memset(r, 0, (n1 + n2 - 1) * sizeof(mat64));
    m64pol_addmul(r, a1, a2, n1, n2);
}

void m64pol_mul_kara(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2)
{
    assert(r != a1 && r != a2);
    assert(n1 == n2);
    assert((n1 & (n1 - 1)) == 0);
    if (n1 == 1) {
        m64pol_mul(r, a1, a2, n1, n2);
        return;
    }
    unsigned int h = n1 >> 1;
    memset(r, 0, (n1 + n2 - 1) * sizeof(mat64));

    m64pol_add(r, a1, a1 + h, h);
    m64pol_add(r + 2 * h, a2, a2 + h, h);

    mat64 * t = (mat64 *) malloc((2*h-1) * sizeof(mat64));
    m64pol_mul_kara(t, r, r + 2 * h, h, h);

    m64pol_mul_kara(r, a1, a2, h, h);
    m64pol_mul_kara(r + 2 * h, a1 + h, a2 + h, h, h);
    m64pol_add(t, t, r, 2 * h - 1);
    m64pol_add(t, t, r + 2 * h, 2 * h - 1);
    m64pol_add(r + h, r + h, t, 2 * h - 1);
    free(t);
}

void m64pol_addmul_kara(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2)
{
    mat64 * t = (mat64 *) malloc((n1 + n2 -1) * sizeof(mat64));
    m64pol_mul_kara(t, a1, a2, n1, n2);
    m64pol_add(r, r, t, n1 + n2 - 1);
    free(t);
}


// 2.11 -- 512 512 based on basecase @ 512
// 1.49 -- 512 512 based on basecase @ 256
// 1.07 -- 512 512 based on basecase @ 128
// 0.76 -- 512 512 based on basecase @ 64
// 0.56 -- 512 512 based on basecase @ 32
// 0.41 -- 512 512 based on basecase @ 16
// 0.31 -- 512 512 based on basecase @ 8
// 0.24 -- 512 512 based on basecase @ 4
// 0.19 -- 512 512 based on basecase @ 2
// 0.14 -- 512 512 based on basecase @ 1

#define TIME1(maxtime, what, args) do {			        	\
    clock_t measuring_time = maxtime * CLOCKS_PER_SEC;			\
    clock_t t0, t1;							\
    int j;								\
    t0 = clock();							\
    for (j = 0; ; j++) {						\
        what args;							\
        t1 = clock() - t0;						\
        if (j && t1 > measuring_time)					\
            break;							\
    }									\
    double t = t1;							\
    t /= CLOCKS_PER_SEC;						\
    t /= j;								\
    char * unit = "s";							\
    if (t < 1.0e-7) { unit = "ns"; t *= 1.0e9;			\
    } else if (t < 1.0e-4) { unit = "micros"; t *= 1.0e6;		\
    } else if (t < 1.0e-1) { unit = "micros"; t *= 1.0e3; }		\
    printf(#what " \t%d times in %.4f %s each\n",       		\
            j, t, unit);		                        	\
} while (0)

/* Same spirit, but treat multiplication of 64K by 64K matrices (of
 * polynomials).
 * 
 * We assume that there is no pointer aliasing, and that matrices are
 * stored row-major, with all polynomials contiguous (and of fixed
 * lengths n1, resp n2).
 */
void m64polblock_mul(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2, unsigned int K)
{
    assert(r != a1 && r != a2);
    memset(r, 0, (n1 + n2 - 1) * K * K * sizeof(mat64));
    for(unsigned int i = 0 ; i < K ; i++) {
        m64pol_srcptr ra1 = a1 + i * n1 * K;
        m64pol_srcptr rr = r + i * (n1 + n2 - 1) * K;
    for(unsigned int j = 0 ; j < K ; j++) {
        m64pol_srcptr ca2 = a2 + j * n2;
        m64pol_srcptr pr = rr + j * (n1 + n2 - 1);
    for(unsigned int k = 0 ; k < K ; k++) {
        m64pol_srcptr pa1 = ra1 + k * n1;
        m64pol_srcptr pa2 = ca2 + k * n2 * K;
        m64pol_addmul(pr, pa1, pa2, n1, n2);
    }
    }
    }
}

void m64polblock_mul_kara(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2, unsigned int K)
{
    assert(r != a1 && r != a2);
    memset(r, 0, (n1 + n2 - 1) * K * K * sizeof(mat64));
    for(unsigned int i = 0 ; i < K ; i++) {
        m64pol_srcptr ra1 = a1 + i * n1 * K;
        m64pol_srcptr rr = r + i * (n1 + n2 - 1) * K;
    for(unsigned int j = 0 ; j < K ; j++) {
        m64pol_srcptr ca2 = a2 + j * n2;
        m64pol_srcptr pr = rr + j * (n1 + n2 - 1);
    for(unsigned int k = 0 ; k < K ; k++) {
        m64pol_srcptr pa1 = ra1 + k * n1;
        m64pol_srcptr pa2 = ca2 + k * n2 * K;
        m64pol_addmul_kara(pr, pa1, pa2, n1, n2);
    }
    }
    }
}

int main()
{
    uint64_t *a, *b, *c, *w;

    unsigned int n = 2 * 1000 * 1000;

    w = (uint64_t *) malloc(64 * sizeof(uint64_t));
    a = (uint64_t *) malloc(n * sizeof(uint64_t));
    b = (uint64_t *) malloc(n * sizeof(uint64_t));
    c = (uint64_t *) malloc(n * sizeof(uint64_t));

    for(unsigned int i = 0 ; i < n ; i++) {
        a[i] = rand();
        b[i] = rand();
        c[i] = rand();
    }
    for(unsigned int i = 0 ; i < 64 ; i++) {
        w[i] = rand();
    }

    if (0) {
        printf("-- level-1 benches --\n");
        TIME1(1, transp_6464, (a,b));
        TIME1(1, copy_6464, (a,b));
        TIME1(1, mul_6464_6464_sse, (w, a, b));
        TIME1(1, add_6464_6464_C, (w, a, b));
        TIME1(1, mul_o64_T6464_C_parity, (w, *a, b));
        TIME1(1, mul_o64_6464_C_lsb, (w, *a, b));
        TIME1(1, mul_o64_6464_C_msb, (w, *a, b));
        TIME1(1, addmul_o64_64o_lsb, (*a,*b,w));
        TIME1(1, addmul_o64_64o_msb, (*a,*b,w));
        TIME1(1, addmul_o64_64o_lsb_sse_v0, (*a,*b,w));
        TIME1(1, addmul_o64_64o_lsb_sse_v1, (*a,*b,w));
        TIME1(1, addmul_o64_64o, (*a,*b,w));
    }

    if (0) {
        printf("-- level-2 benches (N=%u) --\n", n);
        TIME1(2, mul_N64_6464_vec, (a,b,w,n));

        uint64_t * mul_N64_6464_ref = (uint64_t *) malloc(n * sizeof(uint64_t));
        copy_N64(mul_N64_6464_ref, a, n);

        TIME1(2, mul_N64_6464_transB, (a,b,w,n));
        TIME1(2, mul_N64_6464_sse, (a,b,w,n));

        TIME1(1, mul_64N_N64_addmul, (w,a,b,n));
        TIME1(5, mul_TN32_N64_C, (w,(uint32_t*)a,b,n));
        TIME1(5, mul_TN64_N64_C, (w,a,b,n));
    }

    if (0) {
        size_t n = 512;
        mat64 * A = (mat64 *) malloc(n * sizeof(mat64));
        mat64 * B = (mat64 *) malloc(n * sizeof(mat64));
        mat64 * C = (mat64 *) malloc(2 *n * sizeof(mat64));
        printf("-- polynomials (N=%zu) --\n", n);
        TIME1(5, m64pol_mul, (C,A,B,n,n));
        TIME1(5, m64pol_mul_kara, (C,A,B,n,n));
        free(A);
        free(B);
        free(C);
    }

    {
        size_t n = 512;
        mat64 * A = (mat64 *) malloc(2 * 2 * n * sizeof(mat64));
        mat64 * B = (mat64 *) malloc(2 * 2 * n * sizeof(mat64));
        mat64 * C = (mat64 *) malloc(2 * 2 * 2 *n * sizeof(mat64));
        printf("-- polynomials, larger matrices (N=%zu) --\n", n);
        TIME1(5, m64polblock_mul, (C,A,B,n,n,2));
        TIME1(5, m64polblock_mul_kara, (C,A,B,n,n,2));
        free(A);
        free(B);
        free(C);
    }

    free(w);
    free(a);
    free(b);

    return 0;
}
