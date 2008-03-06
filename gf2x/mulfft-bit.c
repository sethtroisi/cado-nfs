/* Multiplication over GF(2)[x] using Fast Fourier Transform.
   Bit-aligned version.

   Copyright Paul Zimmermann 2007.
   
  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or (at your
  option) any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
  more details.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
  02111-1307, USA.

   Reference: Schnelle Multiplikation von Polynomen u"ber Ko"rpern der
   Charakteristik 2, Arnold Scho"nhage, Acta Informatica, vol. 7, 1977.
*/

#include "gf2x.h"

// #define DEBUG
// #define DEBUG_LSHIFT
// #define DEBUG_MULMOD
// #define VERBOSE
// #define TIMING

/* The check in tunefft finds that FFTMul does not work correctly for
   balanced operands of size <= 13 words.  Hence this threshhold. */


#ifndef MUL_FFT_THRESHOLD
#define MUL_FFT_THRESHOLD 28
#endif

#if (MUL_FFT_THRESHOLD < 28)
#error "MUL_FFT_THRESHOLD too small, should be at least 28"
#endif

/* Assume wordlength  WLEN is 32 or 64 */

# define WLEN WORDSIZE

/* CEIL(a,b) = ceiling(a/b) */
#define CEIL(a,b) (1 + ((a)-1)/(b))

/* W(b) is the number of words needed to store b bits */
#define W(b) CEIL(b,WLEN)

/* I(b) is the index word of bit b, assuming bits 0..WLEN-1
   have index 0 */
#define I(b) ((b) / WLEN)

#define R(b) ((b) % WLEN)
#define R2(b) ((WLEN - R(b)) % WLEN)	/* remaining bits */

#define MASK(x) ((1UL << (x)) - 1UL)

/* GETBIT(a,i)   gets the i-th bit of the bit-array starting at a[0],
   XORBIT(a,i,x) xors this bit with the bit x, where x = 0 or 1.   */

#define GETBIT(a,i) ((a[(i)/WLEN] >> ((i)%WLEN))&1)

#define XORBIT(a,i,x) (a[(i)/WLEN] ^= ((x)) << ((i)%WLEN))

static void Copy(unsigned long *a, const unsigned long *b, long n)
{
    ASSERT(n >= 0);
    for (long i = 0; i < n; i++)
	a[i] = b[i];
}

static void Zero(unsigned long *a, long n)
{
    ASSERT(n >= 0);
    for (long i = 0; i < n; i++)
	a[i] = 0;
}

static void Clear(unsigned long *a, long low, long high)
{
    for (long i = low; i < high; i++)
	a[i] = 0;
}

/* a <- b + c */

static void AddMod(unsigned long *a, unsigned long *b, unsigned long *c,
		   long n)
{
    ASSERT(n >= 0);
    for (long i = 0; i < n; i++)
	a[i] = b[i] ^ c[i];
}

/* c <- a * x^k, return carry out, 0 <= k < WLEN */

static unsigned long Lsh(unsigned long *c, unsigned long *a, long n, long k)
{
    unsigned long t, cy = 0;

    ASSERT(n >= 0);
    if (k == 0) {
	if (c != a)
	    Copy(c, a, n);
	return 0;
    }

    /* {c, n} and {a, n} should not overlap */
    ASSERT(c <= a || a + n <= c);

    for (long i = 0; i < n; i++) {
	t = (a[i] << k) | cy;
	cy = a[i] >> (WLEN - k);
	c[i] = t;
    }
    return cy;
}

/* c <- c + a * x^k, return carry out, 0 <= k < WLEN */

static unsigned long AddLsh(unsigned long *c, unsigned long *a, long n,
			    long k)
{
    unsigned long t, cy = 0;

    ASSERT(n >= 0);
    if (k == 0) {
	AddMod(c, c, a, n);
	return 0;
    }

    /* {c, n} and {a, n} should not overlap */
    ASSERT(c <= a || a + n <= c);

    for (long i = 0; i < n; i++) {
	t = (a[i] << k) | cy;
	cy = a[i] >> (WLEN - k);
	c[i] ^= t;
    }
    return cy;
}

/* c <- a / x^k, return carry out, 0 <= k < WLEN */

static unsigned long Rsh(unsigned long *c, const unsigned long *a, long n,
			 long k)
{
    unsigned long t, cy = 0;

    ASSERT(n >= 0);
    if (k == 0) {
	if (c != a)
	    Copy(c, a, n);
	return 0;
    }

    for (long i = n - 1; i >= 0; i--) {
	t = (a[i] >> k) | cy;
	cy = a[i] << (WLEN - k);
	c[i] = t;
    }
    return cy;
}

/* c <- c + a / x^k, return carry out, 0 <= k < WLEN */

static unsigned long AddRsh(unsigned long *c, unsigned long *a, long n,
			    long k)
{
    unsigned long t, cy = 0;

    ASSERT(n >= 0);
    if (k == 0) {
	AddMod(c, c, a, n);
	return 0;
    }

    for (long i = n - 1; i >= 0; i--) {
	t = (a[i] >> k) | cy;
	cy = a[i] << (WLEN - k);
	c[i] ^= t;
    }
    return cy;
}

#if 0
// unused
/* c <- x^j * a / x^k */

static unsigned long LshRsh(unsigned long *c, unsigned long *a, long n,
			    long j, long k)
{
    ASSERT(n >= 0);
    if (j > k)
	Lsh(c, a, n, j - k);
    else
	Rsh(c, a, n, k - j);
}

#endif

#if (defined(DEBUG) || defined(DEBUG_LSHIFT) || defined(DEBUG_MULMOD))

static void dump(const unsigned long *a, long n)
{
    for (long i = 0; i < n; i++) {
	printf("+%lu*X^%lu", a[i], i);
	if ((i + 1) % 3 == 0)
	    printf("\n");
    }
    printf(":\n");
}
#endif

/* a <- b * x^k mod x^(2*N)+x^N+1.
   Assume a and b do not overlap.
*/
static void Lshift(unsigned long *a, unsigned long *b, long k, long N)
{
    long r, h, l, n, ih, il;
    unsigned long s1, s2=0; // hush gcc

    ASSERT(k >= 0);

    n = W(2 * N);

#ifdef DEBUG_LSHIFT
    printf("R:=x^%u+x^%u+1: k:=%u:\nb:=", 2 * N, N, k);
    dump(b, n);
#endif

    k = k % (3 * N);

    ASSERT(0 <= k && k < 3 * N);
    ASSERT(n >= 0);

    if (k == 0) {
	if (a != b)
	    Copy(a, b, n);
    } else if (k <= N) {
	/*  ------------------------------------------
	   |  L0  |      L1      |  L2  |      H      |
	   ------------------------------------------ 
	   ------------------------------------------
	   |      H      |  L0  |    L1+H      |  L2  |
	   ------------------------------------------ 
	   L0 has l bits, L1 has h bits, L2 has l bits, H has h bits
	 */
	h = k;			/* 1 <= h <= N */
	l = N - h;		/* 0 <= l < N */
	/* A <- 0:L0:L1:L2 */
	ASSERT(W(N + l) >= 0);
	s1 = Lsh(a + I(h), b, W(N + l), R(h));
	/* b[W(N + l)-1] has R2(N+l) bits from H, thus s1 has no bit from L
	   if R(h) <= R2(N+l), and R(h)-R2(N+l) bits from L otherwise */
	if (R(h) > R2(N + l))
	    a[I(h) + W(N + l)] = s1 & MASK(R(h) - R2(N + l));
	if (n - I(N + l) > I(h))
	    s1 = a[I(h)];
	if (n - I(N + l) > I(h) + 1)
	    s2 = a[I(h) + 1];
	ASSERT(n - I(N + l) <= I(h) + 2);
	/* copy H in low words: will clobber a[0]..a[n-I(N+l)-1] */
	Rsh(a, b + I(N + l), n - I(N + l), R(N + l));
	if (n - I(N + l) > I(h))
	    a[I(h)] ^= s1;
	if (n - I(N + l) > I(h) + 1)
	    a[I(h) + 1] ^= s2;
	/* we now have a = H:L0:L1:L2 */
	/* add H in high words */
	if (R(N + l) <= R(N)) {
	    s1 = AddLsh(a + I(N), b + I(N + l), n - I(N + l),
			R(N) - R(N + l));
	    /* we have to remove the high R(N+l) bits from L2 */
	    a[I(N)] ^= (b[I(N + l)] & MASK(R(N + l))) << (R(N) - R(N + l));
	    if (s1)
		a[I(N) + n - I(N + l)] ^= s1;
	} else {
	    /* the low word from H in b contains R(N + l) bits from L2, and the
	       low word from L1+H in a contains R(N) bits from L0. After the
	       shift, it remains R(N) bits from L2 */
	    AddRsh(a + I(N), b + I(N + l), n - I(N + l), R(N + l) - R(N));
	    if (R(N))
		a[I(N)] ^=
		    (b[I(N + l)] & MASK(R(N + l))) >> (R(N + l) - R(N));
	}
    } else if (k <= 2 * N) {
	/*  ------------------------------------------
	   |  L   |      H0      |  H1  |     H2      |
	   ------------------------------------------ 
	   ------------------------------------------
	   |     H0+H2   |  H1  |     H0      |  L+H1 |
	   ------------------------------------------ 
	   L has l bits, H0 has h bits, H1 has l bits, H2 has h bits
	 */
	h = k - N;		/* 1 <= h <= N */
	l = N - h;		/* 0 <= l < N */
	r = R(l);
	Rsh(a, b + I(l), W(N + l) - I(l), r);
	/* mask high bits */
	if (R(N))
	    a[I(N)] &= MASK(R(N));
	/* now we have a = H0:H1 */
	if (R(N) > 0)
	    s1 = a[I(N)];
	ASSERT(I(N) >= 0);
	s2 = Lsh(a + I(N), a, I(N), R(N));
	if (R(N) > 0) {
	    a[I(N)] ^= s1;	/* restore high R(N) bits */
	    a[2 * I(N)] = s2 ^ (s1 << R(N));
	    if (2 * R(N) > WLEN)
		a[2 * I(N) + 1] = s1 >> (WLEN - R(N));
	}
	/* now we have a = H0:H1:H0:H1 */
	AddRsh(a, b + I(N + l), W(2 * N) - I(N + l), R(N + l));
	/* now we have a = H0+H2:H1:H0:H1 */
	s1 = AddLsh(a + I(N + h), b, W(l), R(N + h));
	if (r > 0)		/* mask shifted low bits from H0 */
	    a[I(N + h) + I(l)] ^= (b[I(l)] & ~MASK(R(l))) << R(N + h);
	/* b[W(l)-1] contains R2(l) bits from H0, thus if R(N+h) > R2(l),
	   then s1 contains R(N+h)-R2(l) bits from L */
	if (R(N + h) > R2(l))
	    a[I(N + h) + W(l)] ^= s1 & MASK(R(N + h) - R2(l));
    } else {			/* 2*N < k < 3*N */

	/*  ------------------------------------------
	   |  L   |      H0      |  H1  |     H2      |
	   ------------------------------------------ 
	   ------------------------------------------
	   |     H0+H2   |  H1  |     H0      |  L+H1 | after k-N
	   ------------------------------------------ 
	   ------------------------------------------
	   |    H0       | L+H1 |     H2      |   L   |
	   ------------------------------------------ 
	   L has l bits, H0 has h bits, H1 has l bits, H2 has h bits
	 */
	h = k - 2 * N;		/* 1 <= h < N */
	l = N - h;		/* 0 < l < N */
	il = W(2 * N) - I(l);
	Rsh(a, b + I(l), il, R(l));
	/* we now have a = H0:H1:H2:0 */
	s1 = AddLsh(a + I(h), b, W(l), R(h));	/* H0:H1+L:H2:0 */
	if (R(h) > R2(l))
	    /* add least R(h) - R2(l) bits from s1 */
	    a[I(h) + W(l)] ^= s1 & MASK(R(h) - R2(l));
	else if (R2(l) > R(h))
	    /* mask low bits from H0 */
	    a[I(h) + W(l) - 1] ^= (b[W(l) - 1] & ~MASK(R(l))) << R(h);
	/* we now have a = H0:H1+L:H2:0 */
	ih = I(N + h);		/* index of lowest word from L in a */
	if (R(N + h) > 0)	/* H2 and L share a common word a[ih] */
	    s2 = a[ih];
	ASSERT(W(l) >= 0);
	s1 = Lsh(a + ih, b, W(l), R(N + h));	/* H0:H1+L:H2:L */
	/* b[W(l)-1] contains R2(l) bits of H0, thus if R(N+h) > R2(l),
	   s1 contains R(N+h) - R2(l) bits of L */
	if (R(N + h) > R2(l))
	    a[ih + W(l)] = s1 & MASK(R(N + h) - R2(l));
	if (R(N + h) > 0)
	    a[ih] ^= s2;
    }
    /* mask high bits of result */
    r = R(2 * N);
    if (r > 0)
	a[n - 1] &= MASK(r);
#ifdef DEBUG_LSHIFT
    printf("a:=");
    dump(a, n);
#endif
}

/* a <- b * c mod x^(2*N)+x^(N)+1.
   Assumes t has space for 2n words, and u for toomspace(n) words,
   where n = ceil(2N/WLEN).
   a and b may be equal.
   a must have space for n words.
*/

static void MulMod(unsigned long *a, unsigned long *b, unsigned long *c,
		   long N, unsigned long *t, unsigned long *u)
{
    long n = W(2 * N), sh, l;

    /* FIXME: in practice N is divisible by a multiple of 3, thus if WLEN is
       a power of two, R(N) > 0 and W(N) = I(N) + 1, thus we can avoid a few
       tests below. */

#ifdef DEBUG_MULMOD
    printf("b:=");
    dump(b, n);
    printf("c:=");
    dump(c, n);
#endif
    mul_toom(t, b, c, n, u);	/* t = | L0:N | L1:N | H0:N | H1:N | */
    /* t mod x^(2*N)+x^(N)+1 = | L0+H0+H1:N | L1+H0:N | */
    l = W(N);
    Rsh(a, t + I(2 * N), W(3 * N) - I(2 * N), R(2 * N));	/* H0 */
    if (R(N) > 0)		/* mask low bits from H1 if necessary */
	a[I(N)] &= MASK(R(N));
    if (I(N) < l)		/* necessarily l = I(N)+1 */
	u[0] = a[l - 1];
    ASSERT(I(N) >= 0);
    u[1] = Lsh(a + I(N), a, I(N), R(N));	/* H0:H0 */
    /* u[0] and a[1] have R(N) bits */
    if (I(N) < l) {
	ASSERT(2 * l - 2 < n);
	a[2 * l - 2] = u[1] ^ (u[0] << R(N));	/* 2R(N) bits */
	if (2 * R(N) > WLEN) {
	    ASSERT(2 * l - 1 < n);
	    a[2 * l - 1] = u[0] >> (WLEN - R(N));
	}
	a[l - 1] ^= u[0];	/* restore low R(N) bits */
    }
    AddMod(a, a, t, W(2 * N));	/* H0+L0:H0+L1 */
    AddRsh(a, t + I(3 * N), W(4 * N) - I(3 * N), R(3 * N));	/* L0+H0+H1:H0+L1 */
    /* mask high bits */
    sh = R(2 * N);
    if (sh > 0)
	a[n - 1] &= MASK(sh);
#ifdef DEBUG_MULMOD
    printf("a:=");
    dump(a, W(2 * N));
#endif
}

/* starting from i=0, j=0, K and Z=1, enumerates all pairs (i,j)
   for 0 <= i < K such that j = bitrev(i).
*/

static void bitrev(long i, long j, long K, long Z, long *perm)
{
    if (K == 1)
	perm[i] = j;
    else {
	bitrev(i, j, K / 3, 3 * Z, perm);
	bitrev(i + K / 3, j + Z, K / 3, 3 * Z, perm);
	bitrev(i + 2 * K / 3, j + 2 * Z, K / 3, 3 * Z, perm);
    }
}

/* performs an FFT of length K on A[0], A[stride], A[(K-1)*stride] with
   omega=x^j as root of unity, where all computations are
   done modulo x^(2Np) + x^Np + 1.
   Each A[i] has space for 2np words where np=ceil(Np/WLEN).
   Assume omega^K = 1 mod x^(2Np) + x^Np + 1, i.e., mod x^(3Np)+1.
   t1, t2, t3 are buffers of 2np words.
   Assumes 0 <= j < 3*Np.
   For 0 <= i < K, p[i] = tritrev(i), i.e., the integer j such that the
   base-3 representation of i and j are reverse one from each other.
*/

static void fft(unsigned long **A, long K, long j, long Np, long stride,
		unsigned long *t1, unsigned long *t2, unsigned long *t3,
		long *p)
{
    ASSERT(0 <= j && j < 3 * Np);

    if (K == 1)
	return;

    long i, k = K / 3, twonp = W(2 * Np), ii;

    fft(A, k, (3 * j) % (3 * Np), Np, 3 * stride, t1, t2, t3, p);
    fft(A + stride, k, (3 * j) % (3 * Np), Np, 3 * stride, t1, t2, t3, p);
    fft(A + 2 * stride, k, (3 * j) % (3 * Np), Np, 3 * stride, t1, t2, t3, p);

#define a A[3*i*stride]
#define b A[(3*i+1)*stride]
#define c A[(3*i+2)*stride]
    for (i = 0; i < k; i++) {
	ii = p[3 * stride * i];	/* bitrev(i,K/3) = bitrev(3*stride*i,K) */
	Lshift(t1, b, ii * j, Np);	/* t1 = w^ii*b */
	Lshift(t2, b, (ii + 2 * K / 3) * j, Np);	/* t2 = w^(ii+2K/3)*b */
	Lshift(t3, b, (ii + K / 3) * j, Np);
	AddMod(b, a, t3, twonp);
	Lshift(t3, c, (2 * ii + 2 * K / 3) * j, Np);
	AddMod(b, b, t3, twonp);
	Lshift(t3, c, (2 * ii + 4 * K / 3) * j, Np);
	AddMod(t2, t2, t3, twonp);
	Lshift(t3, c, 2 * ii * j, Np);	/* t3 = w^(2ii)*c */
	AddMod(c, a, t2, twonp);
	AddMod(a, a, t1, twonp);
	AddMod(a, a, t3, twonp);
    }
#undef a
#undef b
#undef c
}

/* allocate A[0]...A[K-1], and put there {a, an} cut into K chunks of M bits;
   return pointer to block of memory containing A[0]...A[K-1] (to be freed
   by the calling routine) */

static unsigned long *decompose(unsigned long **A, const unsigned long *a,
				long an, long M, long K, long np)
{
    long i, j, k, l, sh = R(M);
    unsigned long *A0;

    // Allocate space for K buffers each of size 2*np words
    // Since this can not be freed in decompose, A0 is returned.
    // The calling routine should free(A0) when A[0]..A[K-1] are
    // no longer needed.

    A0 = (unsigned long *) malloc(2 * np * K * sizeof(unsigned long));

    for (i = 0; i < K; i++)
	A[i] = A0 + 2 * i * np;

    for (i = 0, j = 0, k = 0; i < K; i++) {
	/* invariant: we have already used k bits of a[j] */
	l = W(k + M);		/* number of words for k+M bits */
	if (j + l > an)
	    l = (an > j) ? an - j : 0;	/* remains l words a[j]..a[j+l-1] */
	if (l > 0) {
	    ASSERT(0 <= j && j + l <= an);
	    Rsh(A[i], a + j, l, k);
	    /* mask last bits */
	    if (sh) {
		ASSERT(I(M) < 2 * np);
		A[i][I(M)] &= MASK(sh);
	    }
	    if (l > W(M))
		l = W(M);
	}
	ASSERT(0 <= l && l < 2 * np);
	Zero(A[i] + l, 2 * np - l);
	k += M;
	j += k / WLEN;
	k %= WLEN;
    }
    return (A0);
}

/* multiplies {a, an} by {b, bn} using an FFT of length K,
   and stores the result into {c, an+bn}. 
   The result is computed mod (x^N+1) where N = K*M. 
   Thus for a full product K*M should be >= (an+bn)*WLEN, 
   the size of product in bits. For the result mod (x^N+1) 
   it is only required that 2*K*M >= (an+bn)*WLEN */

void FFTMul0(unsigned long *c, const unsigned long *a, long an,
	     const unsigned long *b, long bn, long K, long M)
{
    long cn = an + bn;		// size in words of a*b
    // long N = K * M;		// N >= bits(a*b) is a multiple of K
    long Mp = CEIL(M, K / 3);	// ceil(M/(K/3))
    long Np = Mp * (K / 3);	// Np >= M, Np multiple of K/3
    long np = W(Np);		// Words to store Np bits
    long i, j, k, j1, k1, l, ltmp, z;
    long *perm;
    unsigned long **A, **B, *A0, *B0, *tmp1, *tmp2, *tmp3;

    ASSERT((K % 3) == 0);
    ASSERT(2 * K * M >= cn * WLEN);

    A = (unsigned long **) malloc(K * sizeof(unsigned long *));
    B = (unsigned long **) malloc(K * sizeof(unsigned long *));
    A0 = decompose(A, a, an, M, K, np);
    B0 = decompose(B, b, bn, M, K, np);
    i = toomspace(2 * np);
    if (i < 2 * np)
	i = 2 * np;
    ltmp = 4 * np + i;
    tmp1 = (unsigned long *) malloc(ltmp * sizeof(unsigned long));
    tmp2 = tmp1 + 2 * np;
    tmp3 = tmp2 + 2 * np;	/* space for max(2np,toomspace(2np)) words */
    perm = (long *) malloc(K * sizeof(long));
    bitrev(0, 0, K, 1, perm);

    fft(A, K, Mp, Np, 1, tmp1, tmp2, tmp3, perm);
    fft(B, K, Mp, Np, 1, tmp1, tmp2, tmp3, perm);

#ifdef TIMING
    double st = GetTime();
#endif
    for (i = 0; i < K; i++)
	MulMod(B[i], B[i], A[i], Np, tmp1, tmp3);
#ifdef TIMING
    printf("   FFT(%ld): pointwise products on Np=%ld took %f\n", K, Np,
	   GetTime() - st);
#endif

    /* bit reverse: A[i] <- B[bitrev(i)] */
    for (i = 0; i < K; i++)
	A[i] = B[perm[i]];

    fft(A, K, 3 * Np - Mp, Np, 1, tmp1, tmp2, tmp3, perm);

    /* bit reverse: B[i] <- A[bitrev(i)] */
    for (i = 0; i < K; i++)
	B[i] = A[perm[i]];

    /* reconstruct C = sum(B[i]*X^(M*i) mod x^N+1.
       We first compute sum(B[i]*X^(M*i), then reduce it using the wrap function.
       Since we know the result has at most cn words, any value exceeding cn
       words is necessarily zero.
       Each B[i] has 2*Np bits, thus the full C has (K-1)*M+2*Np
       = N - M + 2*Np >= N + Np >= 2*n*WLEN + Np bits.
       Thus exactly 2*Np-M bits wrap around mod x^N+1.
     */

    l = 2 * Np - M;		/* number of overlapping bits with previous B[i] */
    for (i = 0, j = 0, k = 0, j1 = I(l), k1 = R(l); i < K; i++) {
	// unsigned long cy;
	/* invariants:
	   - first bit of B[i] is bit k of c[j]
	   - first bit of B[i] non overlapping with B[i-1] is bit k1 of c[j1] */
	/* add B[i] shifted by i*M to c, where K*M = N, and Np >= M */
	if (i == 0)		/* FIXME: we could set B[0] to c to avoid this copy */
	    Copy(c, B[i], W(2 * Np) < cn ? W(2 * Np) : cn);
	else {
	    /* we have already set bit k of c[j] up to bit k1 of c[j1]
	       (excluded), i.e., words c[j] up to c[j1 - (k1 == 0)] */
	    z = j1 + (k1 != 0) - j;	/* number of overlapping words */
	    /* first treat the high (non overlapping) words of B[i], i.e.,
	       {B[i] + z, W(2*Np) - z} */
	    if (j + W(2 * Np) < cn) {
		if (z < W(2 * Np)) {
		    ASSERT(W(2 * Np) - z >= 0);
		    ASSERT((0 <= j + z) && (j + W(2 * Np) < cn));
		    c[j + W(2 * Np)] =
			Lsh(c + j + z, B[i] + z, W(2 * Np) - z, k);
		} else if (z == W(2 * Np)) {	/* all words overlap with B[i-1] */
		    ASSERT((0 <= j + W(2 * Np)) && (j + W(2 * Np) < cn));
		    c[j + W(2 * Np)] = 0UL;
		}
	    } else if (j + z < cn) {
		ASSERT(0 <= j + z);
		Lsh(c + j + z, B[i] + z, cn - j - z, k);
	    }
	    /* then deal with the low bits of B[i], overlapping with B[i-1] */
	    if (j + z < cn)
		c[j + z] ^= AddLsh(c + j, B[i], z, k);
	    else if (j < cn)
		AddLsh(c + j, B[i], cn - j, k);
	}

	k += M;
	j += k / WLEN;
	k %= WLEN;
	k1 += M;
	j1 += k1 / WLEN;
	k1 %= WLEN;
    }
    free(perm);			// With some implementations of malloc it is
    free(tmp1);			// most efficient to free in the reverse order
    free(B);			// to the mallocs (i.e. using a stack discipline)
    free(A);
    free(B0);
    free(A0);
}

// Wraps the polynomial represented by c mod x^N + 1
// Assumes wraps at most once, i.e. deg(c) < 2N.
// The high part of c (bits N to WLEN*cn) are cleared.
// RPB 20070429

void wrap(unsigned long *c, long cn, long N)
{
    long i;
    long Nw = I(N);
    long Nb = R(N);
    long Nbc = WLEN - Nb;

    // Perhaps most of this could be done by a call to AddLsh ?

    if (N < WLEN * cn) {	// xor bits N .. WLEN*cn of c to c[0...]
	if (Nb == 0) {
	    for (i = 0; i < cn - Nw - 1; i++)
		c[i] ^= c[i + Nw];
	} else {
	    for (i = 0; i < cn - Nw - 1; i++)
		c[i] ^= (c[i + Nw] >> Nb) | (c[i + Nw + 1] << Nbc);
	}
	ASSERT(cn > Nw);
	c[cn - Nw - 1] ^= (c[cn - 1] >> Nb);

	// Now clear remaining bits of c

	c[Nw] &= MASK(Nb);	// Clear high Nbc bits of c[Nw]
	Clear(c, Nw + 1, cn);
    }
}

/* multiplies {a, an} by {b, bn} using an FFT of length K,
   and stores the result into {c, an+bn}. If an+bn is too small
   then Toom-Cook is used.  */

void FFTMul(unsigned long *c, const unsigned long *a, long an,
	    const unsigned long *b, long bn, long K)
{
    if (an + bn < MUL_FFT_THRESHOLD) {
	printf("FFTMul: arguments (%ld, %ld) too small\n", an, bn);
	exit(1);
    }

    long M = CEIL((an + bn) * WLEN, K);	// ceil(bits(product)/K)

    FFTMul0(c, a, an, b, bn, K, M);	// FFTMul0 does the work
}

// Multiplies {a, an} by {b, bn} mod (x^N + 1) and stores the
// result in {c, cn}. Here N = K*M and the FFT is performed on K points
// using O(M)-bit multiplications at each point.
// cn should be at least W(2*N) even though the result is of size W(N).
// The high part of c (bits N to cn*WLEN) are cleared.
// 2*K*M should be >= (an+bn)*WLEN, the size of product in bits. */

void FFTMul1(unsigned long *c, long cn,
	     const unsigned long *aa, long an,
	     const unsigned long *bb, long bn, long K, long M)
{
    long N = K * M;

    ASSERT(cn >= W(2 * N));
    ASSERT(2 * K * M >= (an + bn) * WLEN);

    unsigned long *a, *b;
    // long i;
    long sa = W(N);
    if (sa < an)
	sa = an;		// max (W(N), an, bn) is the
    if (sa < bn)
	sa = bn;		// maximum space needed for temporaries

// Allocate temporaries with enough space and copy inputs.
// In principle we could avoid some of this overhead.

    a = (unsigned long *) malloc(sa * sizeof(unsigned long));
    b = (unsigned long *) malloc(sa * sizeof(unsigned long));

    Copy(a, aa, an);		// Copy aa to a
    Clear(a, an, sa);		// Clear upper part of a

    Copy(b, bb, bn);		// Similarly for bb
    Clear(b, bn, sa);		// Clear upper part of b

    Clear(c, I(N), cn);		// Clear upper part of c

    if (WLEN * an > N)		// Wrap if necessary
    {
	wrap(a, an, N);		// Wrap a mod x^N + 1
	an = W(N);		// New size of a
    }

    if (WLEN * bn > N) {
	wrap(b, bn, N);		// Wrap b mod x^N + 1
	bn = W(N);		// New size of b
    }

    ASSERT(an + bn <= cn);	// Check space for result

    FFTMul0(c, a, an, b, bn, K, M);	// Do the multiplication

    wrap(c, an + bn, N);	// Wrap c mod x^N + 1

    free(b);			// Free temporaries
    free(a);
}

/* Multiplies {a, an} by {b, bn} using (one or) two FFTs of length K,
   and stores the result into {cc, an+bn}. RPB 20070429 */

void FFTMul2(unsigned long *c, const unsigned long *a, long an,
	     const unsigned long *b, long bn, long K)
{


// Avoid splitting for small cases (this should never 
// occur as TC3 is faster than FFTMul for K < 81) 

    if (K < WLEN)		// FFTMul1 needs K >= WLEN
    {				// so in this case do just one
	FFTMul(c, a, an, b, bn, K);	// multiplication using FFTMul
	return;
    }

    long cn, j, m1, m2, n1, n2;
    long cn2 = CEIL(an + bn, 2);	// Space for half product

    unsigned long *c1, *c2;

    m2 = CEIL(cn2 * WLEN, K);	// m2 = ceil(cn2*WLEN/K)
    n2 = K * m2;		// n2 smallest possible multiple of K
    m1 = m2 + 1;		// next possible M
    n1 = K * m1;		// next possible multiple of K 
    cn = W(2 * n1);		// n1 > n2 so cn words is enough 
    // space for temporaries 

// Sometimes (cn > an+bn) so need temporary c1 (as well as c2)

    c1 = (unsigned long *) malloc(cn * sizeof(unsigned long));
    c2 = (unsigned long *) malloc(cn * sizeof(unsigned long));

    FFTMul1(c1, cn, a, an, b, bn, K, m1);	// multiplication mod x^n1 + 1
    FFTMul1(c2, cn, a, an, b, bn, K, m2);	// multiplication mod x^n2 + 1

    long n = WLEN * (an + bn);	// Max bit-size of full product
    long delta = K;		// delta = n1 - n2;
    long jw, jn1w, jn1b, jn1bc, jdw, jdb, jdbc;
    unsigned long t, next;

// Now extract the result. First do a partial word bit-by-bit.

    for (j = n - n1 - 1; (j % WLEN) != (WLEN - 1); j--) {
	t = GETBIT(c1, j + delta) ^ GETBIT(c2, j + delta);
	XORBIT(c1, j + n1, t);	// XOR assumes high part of c1 was zero
	XORBIT(c1, j, t);
    }

// Now do the rest using full-word operations.

    j -= WLEN - 1;
    jdb = R(j + delta);
    jdbc = WLEN - 1 - jdb;
    jn1b = R(j + n1);
    jw = I(j);
    jdw = I(j + delta);
    jn1w = I(j + n1);
    next = c1[jdw + 1] ^ c2[jdw + 1];

    if (jn1b == 0) {		// Unusual case
	for (; jw >= 0; jw--, jdw--, jn1w--) {
	    t = (next << 1) << jdbc;
	    next = c1[jdw] ^ c2[jdw];
	    t ^= next >> jdb;
	    c1[jw] ^= t;
	    c1[jn1w] = t;
	}
    } else			// Usual case
    {
	for (jn1bc = WLEN - jn1b; jw >= 0; jw--, jdw--, jn1w--) {
	    t = (next << 1) << jdbc;
	    next = c1[jdw] ^ c2[jdw];
	    t ^= next >> jdb;
	    c1[jw] ^= t;
	    c1[jn1w] ^= t << jn1b;
	    c1[jn1w + 1] ^= t >> jn1bc;
	}
    }

// Do a consistency check. This is cheap and detects most errors.
// If DEBUG defined we check the first delta bits, otherwise we only
// check the first WLEN bits.

#ifdef DEBUG
    for (j = 0; j < delta; j++) {
	if ((GETBIT(c2, j) ^ GETBIT(c1, j) ^ GETBIT(c1, j + n2)) != 0) {
	    printf("Consistency check failed in FFTMul2, bit %ld\n", j);
	    exit(1);
	}
    }
#endif

    t = c2[0] ^ c1[0] ^ (c1[n2 / WLEN] >> n2 % WLEN) ^
	((c1[n2 / WLEN + 1] << 1) << (WLEN - 1 - n2 % WLEN));
    if (t != 0) {
	printf("Consistency check failed in FFTMul2, low word %lx\n", t);
	exit(1);
    }

    Copy(c, c1, an + bn);	// Copy result

    free(c2);			// Free temporaries
    free(c1);

}

#undef ASSERT
#undef CHECK_ASSERT
