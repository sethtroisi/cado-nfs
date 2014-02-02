#include <string.h>
#include <unistd.h>

#include "flint-fft/fft.h"
/*{{{ macros */
#ifndef PTR
#define PTR(x) ((x)->_mp_d)
#endif
#ifndef SIZ
#define SIZ(x) ((x)->_mp_size)
#endif
#ifndef ALLOC
#define ALLOC(x) ((x)->_mp_alloc)
#endif

#ifndef MPN_NORMALIZE
#define MPN_NORMALIZE(DST, NLIMBS) \
  do {                                                                  \
    while ((NLIMBS) > 0)                                                \
      {                                                                 \
        if ((DST)[(NLIMBS) - 1] != 0)                                   \
          break;                                                        \
        (NLIMBS)--;                                                     \
      }                                                                 \
  } while (0)
#endif

#define BEGIN_BLOCK	do {
#define END_BLOCK	} while (0)

#define MPZ_GROW_ALLOC(DST, NLIMBS)					\
	BEGIN_BLOCK							\
	if (ALLOC(DST) < (int) (NLIMBS)) {				\
		ALLOC(DST) = (NLIMBS);					\
		PTR(DST)=(mp_limb_t *) realloc(PTR(DST),		\
				(NLIMBS) * sizeof(mp_limb_t));		\
	}								\
	END_BLOCK

#define MPZ_INIT_SET_MPN(DST, SRC, NLIMBS)				\
	BEGIN_BLOCK							\
	ALLOC(DST) = (NLIMBS);						\
	SIZ(DST) = (NLIMBS);						\
	PTR(DST) = (mp_limb_t *)malloc((NLIMBS) * sizeof(mp_limb_t));	\
	memcpy(PTR(DST),(SRC),(NLIMBS) * sizeof(mp_limb_t));		\
	MPN_NORMALIZE(PTR(DST),SIZ(DST));				\
	END_BLOCK
	
#define MPZ_SET_MPN(DST, SRC, NLIMBS)					\
	BEGIN_BLOCK							\
	MPZ_GROW_ALLOC(DST, NLIMBS);					\
	SIZ(DST) = (NLIMBS);						\
	memcpy(PTR(DST),(SRC),(NLIMBS) * sizeof(mp_limb_t));		\
	MPN_NORMALIZE(PTR(DST),SIZ(DST));				\
	END_BLOCK
	   
#define MPN_SET_MPZ(DST,NLIMBS,SRC)					\
	BEGIN_BLOCK							\
	memcpy((DST),PTR(SRC),SIZ(SRC) * sizeof(mp_limb_t));		\
	memset((DST)+SIZ(SRC),0,((NLIMBS)-SIZ(SRC)) * sizeof(mp_limb_t));\
	END_BLOCK

#ifndef ABS
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#endif

#ifndef MIN
#define MIN(l,o) ((l) < (o) ? (l) : (o))
#endif

#ifndef MAX
#define MAX(h,i) ((h) > (i) ? (h) : (i))
#endif

#ifndef iceildiv
/* unfortunately this fails miserably if x+y-1 overflows */
#define iceildiv(x,y)	(((x)+(y)-1)/(y))
#endif
/*}}}*/

void get_ft_hash(mpz_t h, int bits_per_coeff, void * data, struct fft_transform_info * fti);

#define xxPARI
/*{{{ display macros */
#ifdef PARI
#define ppol(_name, _x, _nx) do {					\
    mpz_t tmp;								\
    mpz_init(tmp);							\
    printf("v"_name"=[");						\
    for(int i = 0 ; i < _nx ; i++) {					\
        MPZ_SET_MPN(tmp, _x + (_nx-1-i) * np, np);			\
        if (i) gmp_printf(", ");					\
        gmp_printf("%Zd", tmp);	        		        	\
    }									\
    printf("];\n");							\
    printf(_name"=Pol(vector(#v"_name",i,Mod(v"_name"[i],p)));\n");	\
    mpz_clear(tmp);							\
} while (0)
#define pint(_name, _x, _nx) gmp_printf(_name "=%Nd;\n", _x, _nx)
#else
#define ppol(_name, _x, _nx) do {					\
    mpz_t tmp;								\
    mpz_init(tmp);							\
    printf(_name ":=Polynomial(GF(p),[");				\
    for(int i = 0 ; i < _nx ; i++) {					\
        MPZ_SET_MPN(tmp, _x + i * np, np);				\
        if (i) gmp_printf(", ");					\
        gmp_printf("%Zd", tmp);	        		                \
    }									\
    printf("]);\n");							\
    mpz_clear(tmp);							\
} while (0)
#define pint(_name, _x, _nx) gmp_printf(_name ":=%Nd;\n", _x, _nx)
#endif
/*}}}*/

/*{{{ setup operand sizes */
int operand_sizes(int * xbits, int * ybits, int s, gmp_randstate_t rstate)
{
    int base;
    int rs;

    /* We can't handle too small examples */
    do {
        rs = 10 + gmp_urandomm_ui(rstate, 200);
        if (s == 0) s = rs;
        base = 120 * s + gmp_urandomm_ui(rstate, 20 * s);
        *xbits = base + gmp_urandomm_ui(rstate, 5 * s) - 2*s;
        *ybits = base + gmp_urandomm_ui(rstate, 5 * s) - 2*s;
        if (*xbits < 10) *xbits = 10;
        if (*ybits < 10) *ybits = 10;
        /* Do this so that we don't loop forever on easy cases */
        s++;
    } while (*xbits + *ybits < 4000);
    if (*xbits > *ybits) { int a; a = *xbits; *xbits = *ybits; *ybits = a; }
    return 1;
}

int operand_sizes_fppol(int * nx, int * ny, mpz_t p, int s, gmp_randstate_t rstate)
{
    size_t bits_of_p;
    int n;
    bits_of_p = 32 + gmp_urandomm_ui(rstate, 512);

    mp_size_t np;

    do {
        mpz_ui_pow_ui(p, 2, bits_of_p);
        mpz_sub_ui(p, p, 1);
        for( ; !mpz_probab_prime_p(p, 2) ; mpz_sub_ui(p, p, 2));
        np = mpz_size(p);
        n = 120 * s + gmp_urandomm_ui(rstate, 20 * s);
        *nx = n + gmp_urandomm_ui(rstate, 5 * s) - 2*s;
        *ny = n + gmp_urandomm_ui(rstate, 5 * s) - 2*s;
        if (*nx < 10) *nx = 10;
        if (*ny < 10) *ny = 10;
    } while ((*nx+*ny-1) * np * FLINT_BITS < 4000);

    // finer-grain control can go here. But it does not change the
    // picture much.
    if (*nx > *ny) { int a; a = *nx; *nx = *ny; *ny = a; }
    return 1;
}
/*}}}*/

/*{{{ setting entries to random */
void bitrandom(mp_limb_t * x, int xbits, int longstrings, gmp_randstate_t rstate)
{
    int nx = iceildiv(xbits, FLINT_BITS);
    if (longstrings) {
        mpn_rrandom(x, rstate, nx);
    } else {
        mpn_randomb(x, rstate, nx);
    }
    if (xbits % FLINT_BITS) { x[nx-1] &= (1UL<<(xbits%FLINT_BITS))-1; }
}

void bitrandom_fppol(mp_limb_t * x, int nx, mpz_srcptr p, int longstrings, gmp_randstate_t rstate)
{
    mpz_t tmp;
    mpz_init(tmp);
    int np = mpz_size(p);
    int xbits = nx * mpz_size(p) * FLINT_BITS;
    if (longstrings) {
        mpn_rrandom(x, rstate, nx);
    } else {
        mpn_randomb(x, rstate, nx);
    }
    if (xbits % FLINT_BITS) { x[nx-1] &= (1UL<<(xbits%FLINT_BITS))-1; }
    for(int i = 0 ; i < nx ; i++) {
        MPZ_SET_MPN(tmp, x + i * np, np);
        mpz_mod(tmp, tmp, p);
        MPN_SET_MPZ(x + i * np, np, tmp);
    }
    mpz_clear(tmp);
}
/*}}}*/

/* test multiplication of integers */
int test_mul(gmp_randstate_t rstate) /*{{{*/
{
    int seed = getpid();
    int s = 0;
    int longstrings = 0;
    int xbits;
    int ybits;

    // s=84; seed=12682; longstrings=0; 
    // s=93; seed=13156; longstrings=0;
    // s=10; seed=22442; longstrings=0;

    gmp_randseed_ui(rstate, seed);
    operand_sizes(&xbits, &ybits, s, rstate);

    fprintf(stderr, "/* s=%d; seed=%d; longstrings=%d; */\n", s, seed, longstrings);
    fprintf(stderr, "xbits:=%d; ybits:=%d;\n", xbits, ybits);
    int zbits = xbits + ybits;

    mp_size_t nx = iceildiv(xbits, FLINT_BITS);
    mp_size_t ny = iceildiv(ybits, FLINT_BITS);
    mp_size_t nz = iceildiv(zbits, FLINT_BITS);
    mp_limb_t * x = malloc(nx * sizeof(mp_limb_t));
    mp_limb_t * y = malloc(ny * sizeof(mp_limb_t));
    mp_limb_t * z = malloc(nz * sizeof(mp_limb_t));

    bitrandom(x, xbits, longstrings, rstate);
    bitrandom(y, ybits, longstrings, rstate);

    struct fft_transform_info fti[1];
    size_t fft_alloc_sizes[3];

    /* 3 is the maximum number of products we intend to accumulate */
    fft_get_transform_info(fti, xbits, ybits, 3);

#ifndef PARI
    printf("fti_bits:=%lu; fti_ks_coeff_bits:=%lu; fti_depth:=%zu;\n",
            fti->bits, fti->ks_coeff_bits, fti->depth);
    printf("fti_trunc0:=%lu;\n", fti->trunc0);
    printf("fti_w:=%lu;\n", fti->w);
#endif

    fft_get_transform_allocs(fft_alloc_sizes, fti);

    void * tx = malloc(fft_alloc_sizes[0]);
    void * ty = malloc(fft_alloc_sizes[0]);
    void * tz = malloc(fft_alloc_sizes[0]);
    void * tt = malloc(MAX(fft_alloc_sizes[1], fft_alloc_sizes[2]));
    fft_transform_prepare(tx, fti);
    fft_transform_prepare(ty, fti);
    fft_transform_prepare(tz, fti);

#ifdef PARI
    printf("allocatemem(800000000)\n");
#endif

    pint("A0", x, nx);
    pint("A1", y, ny);

    fft_do_dft(tx, x, nx, tt, fti);
    rename("/tmp/before_dft.m", "/tmp/P0_before_dft.m");
    rename("/tmp/after_dft.m", "/tmp/P0_after_dft.m");
    fft_do_dft(ty, y, ny, tt, fti);
    rename("/tmp/before_dft.m", "/tmp/P1_before_dft.m");
    rename("/tmp/after_dft.m", "/tmp/P1_after_dft.m");

    fft_mul(tz, tx, ty, tt, fti);
    fft_add(tz, tz, tz, fti);
    fft_add(tz, tz, tx, fti);
    fft_add(tz, tz, ty, fti);
    fft_do_ift(z, nz, tz, tt, fti);
    rename("/tmp/before_ift.m", "/tmp/P2_before_ift.m");
    rename("/tmp/after_ift.m", "/tmp/P2_after_ift.m");
    pint("A2", z, nz);

#ifdef  PARI
    printf("print(A2==2*A0*A1+A0+A1)\n");
    // printf("assert P2 eq P0*P1;\n");
    printf("quit\n");
#else
    printf("A2 eq 2*A0*A1+A0+A1;\n");
#endif

    /* magma check code
Z:=Integers();
ZP<T>:=PolynomialRing(Z);
n:=2^fti_depth;
R:=Integers(2^(n*fti_w) + 1);
load "/tmp/P0_before_dft.m"; Q0:=Polynomial([R!Seqint(x,2^64):x in data]);
Q0 eq Polynomial(R,Intseq(A0,2^fti_bits));
load "/tmp/P1_before_dft.m"; Q1:=Polynomial([R!Seqint(x,2^64):x in data]);
Q1 eq Polynomial(R,Intseq(A1,2^fti_bits));
load "/tmp/P0_after_dft.m"; tQ0:=([R!Seqint(x,2^64):x in data]);
load "/tmp/P1_after_dft.m"; tQ1:=([R!Seqint(x,2^64):x in data]);
bitrev:=func<x,n|Seqint(Reverse(Intseq(x,2,n)),2)>;
bitrevseq:=func<n|[bitrev(i,n):i in [0..2^n-1]]>;
rho:=R!2^(fti_w div 2);
seqmatch:=func<S,T,n|S[1..n] eq T[1..n]>;
seqmatch([Evaluate(Q0,rho^i):i in bitrevseq(fti_depth+2)], tQ0, fti_trunc0);
seqmatch([Evaluate(Q1,rho^i):i in bitrevseq(fti_depth+2)], tQ1, fti_trunc0);
load "/tmp/P2_before_ift.m"; tQ2:=([R!Seqint(x,2^64):x in data]);
Q2:=2*Q0*Q1+Q0+Q1 mod (T^(4*n)-1);
seqmatch([Evaluate(Q2,rho^i):i in bitrevseq(fti_depth+2)], tQ2, fti_trunc0);
A2 eq Evaluate(ChangeRing(Q2,Z),2^fti_bits);
    */

    free(tx);
    free(ty);
    free(tz);
    free(tt);
    free(x);
    free(y);
    free(z);
    return 0;
}/*}}}*/

/* test wrapped product of integers. This computea A*B mod base^n\pm1 */
int test_mulmod(gmp_randstate_t rstate) /*{{{*/
{
    int seed = getpid();
    int s = 0;
    int longstrings = 0;
    int xbits;
    int ybits;

    // s=84; seed=12682; longstrings=0; 
    // s=93; seed=13156; longstrings=0;
    // s=10; seed=22442; longstrings=0;
    // s=0; seed=8412; longstrings=0;

    gmp_randseed_ui(rstate, seed);
    operand_sizes(&xbits, &ybits, s, rstate);
    int wrap = gmp_urandomm_ui(rstate, 128);

    

    fprintf(stderr, "/* s=%d; seed=%d; longstrings=%d; */\n", s, seed, longstrings);
    fprintf(stderr, "xbits:=%d; ybits:=%d;\n", xbits, ybits);

    // int zbits = ybits - xbits + 1;


    struct fft_transform_info fti[1];
    /* 3 is the maximum number of products we intend to accumulate */
    fft_get_transform_info_mulmod(fti, xbits, ybits, 3, ybits + wrap);

    int a;
    int zbits = (int) fft_get_mulmod(fti, &a);

    mp_size_t nx = iceildiv(xbits, FLINT_BITS);
    mp_size_t ny = iceildiv(ybits, FLINT_BITS);
    mp_size_t nz = fft_get_mulmod_output_minlimbs(fti);
    mp_limb_t * x = malloc(nx * sizeof(mp_limb_t));
    mp_limb_t * y = malloc(ny * sizeof(mp_limb_t));
    mp_limb_t * z = malloc(nz * sizeof(mp_limb_t));

    bitrandom(x, xbits, longstrings, rstate);
    bitrandom(y, ybits, longstrings, rstate);


    size_t fft_alloc_sizes[3];

#ifndef PARI
    printf("fti_bits:=%lu; fti_ks_coeff_bits:=%lu; fti_depth:=%zu;\n",
            fti->bits, fti->ks_coeff_bits, fti->depth);
    printf("fti_trunc0:=%lu;\n", fti->trunc0);
    printf("fti_w:=%lu;\n", fti->w);
#endif

    fft_get_transform_allocs(fft_alloc_sizes, fti);

    void * tx = malloc(fft_alloc_sizes[0]);
    void * ty = malloc(fft_alloc_sizes[0]);
    void * tz = malloc(fft_alloc_sizes[0]);
    void * tt = malloc(MAX(fft_alloc_sizes[1], fft_alloc_sizes[2]));
    fft_transform_prepare(tx, fti);
    fft_transform_prepare(ty, fti);
    fft_transform_prepare(tz, fti);

#ifdef PARI
    printf("allocatemem(800000000)\n");
#endif

    pint("A0", x, nx);
    pint("A1", y, ny);

    fft_do_dft(tx, x, nx, tt, fti);
    rename("/tmp/before_dft.m", "/tmp/P0_before_dft.m");
    rename("/tmp/after_dft.m", "/tmp/P0_after_dft.m");
    fft_do_dft(ty, y, ny, tt, fti);
    rename("/tmp/before_dft.m", "/tmp/P1_before_dft.m");
    rename("/tmp/after_dft.m", "/tmp/P1_after_dft.m");

    fft_mul(tz, tx, ty, tt, fti);
    fft_add(tz, tz, tz, fti);
    fft_add(tz, tz, tx, fti);
    fft_add(tz, tz, ty, fti);
    /*
    */
    fft_do_ift(z, nz, tz, tt, fti);
    rename("/tmp/before_ift.m", "/tmp/P2_before_ift.m");
    rename("/tmp/after_ift.m", "/tmp/P2_after_ift.m");
    pint("A2", z, nz);

#ifdef  PARI
    printf("print(A2==2*A0*A1+A0+A1)\n");
    // printf("assert P2 eq P0*P1;\n");
    printf("quit\n");
#else
    printf("A2 eq (2*A0*A1+A0+A1) mod (2^%d-%d);\n", zbits, a);
#endif

    /* TODO: we must take into account the fact that the wraparound adds
     * extra summands to each polynomial coefficient */
    /* magma check code

Z:=Integers();
ZP<T>:=PolynomialRing(Z);
n:=2^fti_depth;
R:=Integers(2^(n*fti_w) + 1);
load "/tmp/P0_before_dft.m"; Q0:=Polynomial([R!Seqint(x,2^64):x in data]);
Q0 eq Polynomial(R,Intseq(A0,2^fti_bits));
load "/tmp/P1_before_dft.m"; Q1:=Polynomial([R!Seqint(x,2^64):x in data]);
Q1 eq Polynomial(R,Intseq(A1,2^fti_bits));
load "/tmp/P0_after_dft.m"; tQ0:=([R!Seqint(x,2^64):x in data]);
load "/tmp/P1_after_dft.m"; tQ1:=([R!Seqint(x,2^64):x in data]);
bitrev:=func<x,n|Seqint(Reverse(Intseq(x,2,n)),2)>;
bitrevseq:=func<n|[bitrev(i,n):i in [0..2^n-1]]>;
rho:=R!2^(fti_w div 2);
seqmatch:=func<S,T,n|S[1..n] eq T[1..n]>;
seqmatch([Evaluate(Q0,rho^i):i in bitrevseq(fti_depth+2)], tQ0, fti_trunc0);
seqmatch([Evaluate(Q1,rho^i):i in bitrevseq(fti_depth+2)], tQ1, fti_trunc0);
load "/tmp/P2_after_ift.m"; cQ2:=Polynomial([R!Seqint(x,2^64):x in data]);
load "/tmp/P2_before_ift.m"; tQ2:=([R!Seqint(x,2^64):x in data]);
Q2:=(2*Q0*Q1+Q0+Q1) mod (T^(4*n)-1);
// Q2:=Q0*Q1 mod (T^(4*n)-1);
Q2 eq cQ2;
seqmatch([Evaluate(Q2,rho^i):i in bitrevseq(fti_depth+2)], tQ2, fti_trunc0);
A2 eq Evaluate(ChangeRing(Q2,Z),2^fti_bits) mod (2^(4*n*fti_bits)-1);

    */

    free(tx);
    free(ty);
    free(tz);
    free(tt);
    free(x);
    free(y);
    free(z);
    return 0;
}/*}}}*/

/* multiplication of polynomials */
int test_mul_fppol(gmp_randstate_t rstate) /*{{{*/
{
    mpz_t p;
    int seed = getpid();
    int s = 0;
    int longstrings = 0;
    mpz_init(p);
    int nx;
    int ny;

    // seed=6286; longstrings=0;

    gmp_randseed_ui(rstate, seed);
    operand_sizes_fppol(&nx, &ny, p, s, rstate);
    mp_size_t np = mpz_size(p);


    fprintf(stderr, "/* s=%d; seed=%d; longstrings=%d; */\n", s, seed, longstrings);
    fprintf(stderr, "nx:=%d; ny:=%d; bits_of_p:=%zu;\n", nx, ny, mpz_sizeinbase(p, 2));
    int nz = nx + ny - 1;


    mp_limb_t * x = malloc(nx * np * sizeof(mp_limb_t));
    mp_limb_t * y = malloc(ny * np * sizeof(mp_limb_t));
    mp_limb_t * z = malloc(nz * np * sizeof(mp_limb_t));
    bitrandom_fppol(y, ny, p, longstrings, rstate);
    bitrandom_fppol(x, nx, p, longstrings, rstate);

    // mul_mfa_truncate_sqrt2(z, x, nx*np, y, ny*np, 10, 6);

    struct fft_transform_info fti[1];
    size_t fft_alloc_sizes[3];

    /* 3 is the maximum number of products we intend to accumulate */
    fft_get_transform_info_fppol_mp(fti, p, nx, ny, 3);

#ifndef PARI
    printf("fti_bits:=%lu; fti_ks_coeff_bits:=%lu; fti_depth:=%zu;\n",
            fti->bits, fti->ks_coeff_bits, fti->depth);
    printf("fti_trunc0:=%lu;\n", fti->trunc0);
    printf("fti_w:=%lu;\n", fti->w);
#endif

    fft_get_transform_allocs(fft_alloc_sizes, fti);

    void * tx = malloc(fft_alloc_sizes[0]);
    void * ty = malloc(fft_alloc_sizes[0]);
    void * tz = malloc(fft_alloc_sizes[0]);
    void * tt = malloc(MAX(fft_alloc_sizes[1], fft_alloc_sizes[2]));
    fft_transform_prepare(tx, fti);
    fft_transform_prepare(ty, fti);
    fft_transform_prepare(tz, fti);

#ifdef PARI
    printf("allocatemem(800000000)\n");
    gmp_printf("p=%Zd;\n", p);
#else
    gmp_printf("p:=%Zd;\n", p);
    printf("KP<x>:=PolynomialRing(GF(p));\n");
#endif

    ppol("P0", x, nx);
    ppol("P1", y, ny);

    fft_do_dft_fppol(tx, x, nx * np, tt, fti, p);
    // get_ft_hash(tmp, 1, tx, fti);
    // gmp_fprintf(stderr, "%Zx\n", tmp);
    rename("/tmp/before_dft.m", "/tmp/P0_before_dft.m");
    rename("/tmp/after_dft.m", "/tmp/P0_after_dft.m");
    fft_do_dft_fppol(ty, y, ny * np, tt, fti, p);
    rename("/tmp/before_dft.m", "/tmp/P1_before_dft.m");
    rename("/tmp/after_dft.m", "/tmp/P1_after_dft.m");

    /* /tmp/P0_before_dft.m must be such that the following holds (this
     * is with PARI=0:
     
     n:=2^fti_depth;
     R:=Integers(2^(n*fti_w) + 1);
     load "/tmp/P0_before_dft.m"; Q0:=Polynomial([R!Seqint(x,2^64):x in data]);
     Q0 eq Polynomial(R,Intseq(Evaluate(Polynomial(ChangeUniverse(Eltseq(P0),Integers())),2^fti_ks_coeff_bits),2^fti_bits));
     load "/tmp/P1_before_dft.m"; Q1:=Polynomial([R!Seqint(x,2^64):x in data]);
     Q1 eq Polynomial(R,Intseq(Evaluate(Polynomial(ChangeUniverse(Eltseq(P1),Integers())),2^fti_ks_coeff_bits),2^fti_bits));

     load "/tmp/P0_after_dft.m"; tQ0:=([R!Seqint(x,2^64):x in data]);
     load "/tmp/P1_after_dft.m"; tQ1:=([R!Seqint(x,2^64):x in data]);
    bitrev:=func<x,n|Seqint(Reverse(Intseq(x,2,n)),2)>;
bitrevseq:=func<n|[bitrev(i,n):i in [0..2^n-1]]>;
     rho:=R!2^(fti_w div 2);
     [Evaluate(Q0,rho^i):i in bitrevseq(fti_depth+2)] eq tQ0;
     [Evaluate(Q1,rho^i):i in bitrevseq(fti_depth+2)] eq tQ1;

    */
    fft_mul(tz, tx, ty, tt, fti);
    fft_add(tz, tz, tz, fti);
    fft_add(tz, tz, tx, fti);
    fft_add(tz, tz, ty, fti);
    fft_do_ift_fppol(z, nz * np, tz, tt, fti, p);
    rename("/tmp/before_ift.m", "/tmp/P2_before_ift.m");
    rename("/tmp/after_ift.m", "/tmp/P2_after_ift.m");

    ppol("P2", z, nz);

#ifdef  PARI
    printf("print(P2==2*P0*P1+P0+P1)\n");
    // printf("assert P2 eq P0*P1;\n");
    printf("quit\n");
#endif

    free(tx);
    free(ty);
    free(tz);
    free(tt);
    free(x);
    free(y);
    free(z);
    mpz_clear(p);
    return 0;
}/*}}}*/

#if 0
/* middle product of polynomials */
int test_mp_fppol(gmp_randstate_t rstate)/*{{{*/
{
    mpz_t p;
    int seed = getpid();
    int s = 0;
    int longstrings = 0;
    mpz_init(p);
    int nx;
    int ny;

    // Here is the list of setup bugs I had to cover, successively.
    // s=3; seed = 17769; longstrings=1;
    // s=12; seed=1010; longstrings=0;
    s=24; seed=6931; longstrings=0;

    gmp_randseed_ui(rstate, seed);
    operand_sizes_fppol(&nx, &ny, p, s, rstate);
    mp_size_t np = mpz_size(p);

    /* We'll do the transpose of
     * MUL(nx, ny) == nz ; which is MP(nx, nz) == ny.
     * But we'll rewrite this as MP(nx, ny) == nz by swapping ny and nz.
     */
    int nz = ny;
    ny = nx + nz - 1;

    fprintf(stderr, "/* s=%d; seed=%d; longstrings=%d; */\n", s, seed, longstrings);
    fprintf(stderr,
            "/* MP(degree %d, degree %d) -> degree %d */\n",
            nx - 1, ny - 1, nz - 1);


    mp_limb_t * x = malloc(nx * np * sizeof(mp_limb_t));
    mp_limb_t * y = malloc(ny * np * sizeof(mp_limb_t));
    mp_limb_t * z = malloc(nz * np * sizeof(mp_limb_t));
    bitrandom_fppol(y, ny, p, longstrings, rstate);
    bitrandom_fppol(x, nx, p, longstrings, rstate);

    mpz_init(p);
    mpz_ui_pow_ui(p, 2, bits_of_p);
    mpz_sub_ui(p, p, 1);
    for( ; !mpz_probab_prime_p(p, 2) ; mpz_sub_ui(p, p, 2));
    mp_size_t np = mpz_size(p);

    mp_limb_t * x = malloc(nx * np * sizeof(mp_limb_t));
    mp_limb_t * z = malloc(nz * np * sizeof(mp_limb_t));
    mp_limb_t * y = malloc(ny * np * sizeof(mp_limb_t));

    // mul_mfa_truncate_sqrt2(y, x, nx*np, z, nz*np, 10, 6);

    struct fft_transform_info fti[1];
    size_t fft_alloc_sizes[3];

    /* 3 is the maximum number of products we intend to accumulate */
    fft_get_transform_info_fppol_mp(fti, p, nx, nz, 3);

#ifndef PARI
    printf("fti_bits:=%lu; fti_ks_coeff_bits:=%lu; fti_depth:=%zu;\n",
            fti->bits, fti->ks_coeff_bits, fti->depth);
    printf("fti_trunc0:=%lu;\n", fti->trunc0);
    printf("fti_w:=%lu;\n", fti->w);
#endif

    fft_get_transform_allocs(fft_alloc_sizes, fti);

    void * tx = malloc(fft_alloc_sizes[0]);
    void * tz = malloc(fft_alloc_sizes[0]);
    void * ty = malloc(fft_alloc_sizes[0]);
    void * tt = malloc(MAX(fft_alloc_sizes[1], fft_alloc_sizes[2]));
    fft_transform_prepare(tx, fti);
    fft_transform_prepare(tz, fti);
    fft_transform_prepare(ty, fti);

#ifdef PARI
    printf("allocatemem(800000000)\n");
    gmp_printf("p=%Zd;\n", p);
#else
    gmp_printf("p:=%Zd;\n", p);
    printf("KP<x>:=PolynomialRing(GF(p));\n");
#endif

    ppol("P0", x, nx);
    ppol("P1", z, nz);

    fft_do_dft_fppol(tx, x, nx * np, tt, fti, p);
    // get_ft_hash(tmp, 1, tx, fti);
    // gmp_fprintf(stderr, "%yx\n", tmp);
    rename("/tmp/before_dft.m", "/tmp/P0_before_dft.m");
    rename("/tmp/after_dft.m", "/tmp/P0_after_dft.m");

    fft_do_dft_fppol(tz, z, nz * np, tt, fti, p);
    rename("/tmp/before_dft.m", "/tmp/P1_before_dft.m");
    rename("/tmp/after_dft.m", "/tmp/P1_after_dft.m");

#if 0

Here is some magma code for checking.
     
n:=2^fti_depth;
R:=Integers(2^(n*fti_w) + 1);
RP<T>:=PolynomialRing(R);
Z:=Integers();
// integers
zP0:=Evaluate(ChangeRing(P0,Z),2^fti_ks_coeff_bits);
zP1:=Evaluate(ChangeRing(P1,Z),2^fti_ks_coeff_bits);
load "/tmp/P0_before_dft.m"; Q0:=Polynomial([R!Seqint(x,2^64):x in data]);
load "/tmp/P1_before_dft.m"; Q1:=Polynomial([R!Seqint(x,2^64):x in data]);
// now polynomials over R. We pick blocks of fti_bits, which is less
// than half of n*fti_w
Q0 eq Polynomial(R,Intseq(zP0,2^fti_bits));
Q1 eq Polynomial(R,Intseq(zP1,2^fti_bits));
load "/tmp/P0_after_dft.m"; tQ0:=([R!Seqint(x,2^64):x in data]);
load "/tmp/P1_after_dft.m"; tQ1:=([R!Seqint(x,2^64):x in data]);
bitrev:=func<x,n|Seqint(Reverse(Intseq(x,2,n)),2)>;
bitrevseq:=func<n|[bitrev(i,n):i in [0..2^n-1]]>;
rho:=R!2^(fti_w div 2);
[Evaluate(Q0,rho^i):i in bitrevseq(fti_depth+2)][1..fti_trunc0] eq tQ0[1..fti_trunc0];
[Evaluate(Q1,rho^i):i in bitrevseq(fti_depth+2)][1..fti_trunc0] eq tQ1[1..fti_trunc0];
load "/tmp/P2_before_ift.m"; tQ2:=([R!Seqint(x,2^64):x in data]);
Q2:=(Q0*Q1) mod (T^(4*n)-1);
tQ2[1..fti_trunc0] eq [Evaluate(Q2,rho^i):i in bitrevseq(fti_depth+2)][1..fti_trunc0];
zP2:=Evaluate(ChangeRing(Q2,Z),2^fti_bits);
// check that product in RP gives the info to recover MP(P0,P1) as
// desired. We're checking all this between the integers on one side, and
// GF(p) on the other side. There are other ways to write the same
// (explicitly writing "mod p", e.g.).
MP:=func<P0,P1|[Coefficient(P,i) :i in [Min(Degree(P0),Degree(P1))..Max(Degree(P0),Degree(P1))]] where P is P0*P1>;
P2:=Intseq(zP2,2^fti_ks_coeff_bits)[Degree(P0)+1..Degree(P1)+1];
MP(P0,P1) eq P2;
#endif
    /* TODO: we're doing a middle product, which means that we must pick
     * our result coefficients precisely from the point where they sit in
     * the in-memory data. It's not the same as ift_fppol. Or if we
     * insist on using fft_do_ift_fppol, then we first have to rotate the
     * data appropriately (this will be clumsy). */
    fft_mul(ty, tx, tz, tt, fti);

    fft_do_ift_fppol(y, ny * np, ty, tt, fti, p);
    rename("/tmp/before_ift.m", "/tmp/P2_before_ift.m");
    rename("/tmp/after_ift.m", "/tmp/P2_after_ift.m");

    ppol("P2", y, ny);

#ifdef  PARI
    /* P2 should contain the middle product of P0 and P1, but we do not
     * have complete code at this point */
    printf("print(-1)\n");
    printf("quit\n");
#endif

    free(tx);
    free(tz);
    free(ty);
    free(tt);
    free(x);
    free(z);
    free(y);
    return 0;
}/*}}}*/
#endif


int main()
{
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    // test_mul(rstate);
    test_mulmod(rstate);
    // test_mul_fppol(rstate);
    // test_mp_fppol(rstate);
    gmp_randclear(rstate);
    return 0;
}

