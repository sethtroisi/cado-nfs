#include <string.h>
#include <unistd.h>

#include "flint-fft/fft.h"

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

void get_ft_hash(mpz_t h, int bits_per_coeff, void * data, struct fft_transform_info * fti);

#define xxPARI

#ifdef PARI
#define ppol(_name, _x, _nx) do {					\
    printf("v"_name"=[");						\
    for(int i = 0 ; i < _nx ; i++) {					\
        MPZ_SET_MPN(tmp, _x + (_nx-1-i) * np, np);			\
        if (i) gmp_printf(", ");					\
        gmp_printf("%Zd", tmp);	        		        	\
    }									\
    printf("];\n");							\
    printf(_name"=Pol(vector(#v"_name",i,Mod(v"_name"[i],p)));\n");			\
} while (0)
#else
#define ppol(_name, _x, _nx) do {					\
    printf(_name ":=Polynomial(GF(p),[");				\
    for(int i = 0 ; i < _nx ; i++) {					\
        MPZ_SET_MPN(tmp, _x + i * np, np);				\
        if (i) gmp_printf(", ");					\
        gmp_printf("%Zd", tmp);	        		                \
    }									\
    printf("]);\n");							\
} while (0)
#endif

int test_mul_fppol(gmp_randstate_t rstate) /*{{{*/
{
    mpz_t p;
    mpz_t tmp;

    int seed = getpid();
    int s = 0;
    int longstrings = 0;

    // seed=6286; longstrings=0;

    gmp_randseed_ui(rstate, seed);
    int rs = 10 + gmp_urandomm_ui(rstate, 200);
    if (s == 0) s = rs;

    fprintf(stderr, "/* s=%d; seed=%d; longstrings=%d; */\n", s, seed, longstrings);

    size_t bits_of_p = 32 + gmp_urandomm_ui(rstate, 512);


    int n = 120 * s + gmp_urandomm_ui(rstate, 20 * s);
    int nx = n + gmp_urandomm_ui(rstate, 5 * s) - 2*s;
    int ny = n + gmp_urandomm_ui(rstate, 5 * s) - 2*s;
    if (nx < 10) nx = 10;
    if (ny < 10) ny = 10;

    // finer-grain control can go here. But it does not change the
    // picture much.

#ifndef PARI
    fprintf(stderr, "nx:=%d; ny:=%d; bits_of_p:=%zu;\n", nx, ny, bits_of_p);
#else
    fprintf(stderr, "nx=%d; ny=%d; bits_of_p=%zu;\n", nx, ny, bits_of_p);
#endif

    int nz = nx + ny - 1;

    mpz_init(p);
    mpz_ui_pow_ui(p, 2, bits_of_p);
    mpz_sub_ui(p, p, 1);
    for( ; !mpz_probab_prime_p(p, 2) ; mpz_sub_ui(p, p, 2));
    mp_size_t np = mpz_size(p);

    mpz_init(tmp);
    mp_limb_t * x = malloc(nx * np * sizeof(mp_limb_t));
    mp_limb_t * y = malloc(ny * np * sizeof(mp_limb_t));
    mp_limb_t * z = malloc(nz * np * sizeof(mp_limb_t));

    if (longstrings) {
        mpn_rrandom(x, rstate, nx*np);
        mpn_rrandom(y, rstate, ny*np);
    } else {
        mpn_randomb(x, rstate, nx*np);
        mpn_randomb(y, rstate, ny*np);
    }

    // mul_mfa_truncate_sqrt2(z, x, nx*np, y, ny*np, 10, 6);

    for(int i = 0 ; i < nx ; i++) {
        MPZ_SET_MPN(tmp, x + i * np, np);
        mpz_mod(tmp, tmp, p);
        MPN_SET_MPZ(x + i * np, np, tmp);
    }
    for(int i = 0 ; i < ny ; i++) {
        MPZ_SET_MPN(tmp, y + i * np, np);
        mpz_mod(tmp, tmp, p);
        MPN_SET_MPZ(y + i * np, np, tmp);
    }

    struct fft_transform_info fti[1];
    size_t fft_alloc_sizes[3];

    /* 3 is the maximum number of products we intend to accumulate */
    fft_get_transform_info_fppol(fti, p, nx, ny, 3);

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
     [Evaluate(Q0,rho^i):i in bitrevseq(n)] eq tQ0;
     [Evaluate(Q1,rho^i):i in bitrevseq(n)] eq tQ1;

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
    return 0;
}/*}}}*/

void fft_do_ift_fppol_hack_debug_mp(mp_limb_t * x, mp_size_t nx, void * y, void * temp, struct fft_transform_info * fti, mpz_srcptr p);

int test_mp_fppol(gmp_randstate_t rstate)/*{{{*/
{
    mpz_t p;
    mpz_t tmp;

    int seed = getpid();
    int s = 0;
    int longstrings = 0;

    // seed=6286; longstrings=0;
    // s=3; seed = 17769; longstrings=1;
    // s=12; seed=1010; longstrings=0;


    gmp_randseed_ui(rstate, seed);
    int rs = 10 + gmp_urandomm_ui(rstate, 20);
    if (s == 0) s = rs;

    fprintf(stderr, "/* s=%d; seed=%d; longstrings=%d; */\n", s, seed, longstrings);


    size_t bits_of_p = 32 + gmp_urandomm_ui(rstate, 64);


    int n = 20 * s + gmp_urandomm_ui(rstate, 20 * s);
    int nx = n + gmp_urandomm_ui(rstate, 5 * s) - 2*s;
    int nz = n + gmp_urandomm_ui(rstate, 5 * s) + 2*s;
    if (nx < 10) nx = 10;
    if (nz < 10) nz = 10;
    if (nx > nz) { int a; a = nx; nx = nz; nz = a; }

    // finer-grain control can go here. But it does not change the
    // picture much.

#ifndef PARI
    fprintf(stderr, "nx:=%d; nz:=%d; bits_of_p:=%zu;\n", nx, nz, bits_of_p);
#else
    fprintf(stderr, "nx=%d; nz=%d; bits_of_p=%zu;\n", nx, nz, bits_of_p);
#endif

    int ny = MAX(nx, nz) - MIN(nx, nz) + 1;

    fprintf(stderr,
            "/* MP(degree %d, degree %d) -> degree %d */\n",
            nx - 1, nz - 1, ny - 1);


    mpz_init(p);
    mpz_ui_pow_ui(p, 2, bits_of_p);
    mpz_sub_ui(p, p, 1);
    for( ; !mpz_probab_prime_p(p, 2) ; mpz_sub_ui(p, p, 2));
    mp_size_t np = mpz_size(p);

    mpz_init(tmp);
    mp_limb_t * x = malloc(nx * np * sizeof(mp_limb_t));
    mp_limb_t * z = malloc(nz * np * sizeof(mp_limb_t));
    mp_limb_t * y = malloc(ny * np * sizeof(mp_limb_t));

    if (longstrings) {
        mpn_rrandom(x, rstate, nx*np);
        mpn_rrandom(z, rstate, nz*np);
    } else {
        mpn_randomb(x, rstate, nx*np);
        mpn_randomb(z, rstate, nz*np);
    }

    // mul_mfa_truncate_sqrt2(y, x, nx*np, z, nz*np, 10, 6);

    for(int i = 0 ; i < nx ; i++) {
        MPZ_SET_MPN(tmp, x + i * np, np);
        mpz_mod(tmp, tmp, p);
        MPN_SET_MPZ(x + i * np, np, tmp);
    }
    for(int i = 0 ; i < nz ; i++) {
        MPZ_SET_MPN(tmp, z + i * np, np);
        mpz_mod(tmp, tmp, p);
        MPN_SET_MPZ(z + i * np, np, tmp);
    }

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

    fft_do_ift_fppol_hack_debug_mp(y, ny * np, ty, tt, fti, p);
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


int test_mul()
{
    mpz_t x,y,z;
    mpz_init(x);
    mpz_init(y);
    mpz_init(z);
    mpz_ui_pow_ui(x, 10, 250000);

    mpz_realloc(z, 2*mpz_size(x));
    SIZ(z) = 2 * mpz_size(x);

    struct fft_transform_info fti[1];
    size_t fft_alloc_sizes[3];

    fft_get_transform_info(fti, mpz_sizeinbase(x, 2), mpz_sizeinbase(x, 2), 1);
    fft_get_transform_allocs(fft_alloc_sizes, fti);

    void * tt = malloc(MAX(fft_alloc_sizes[1], fft_alloc_sizes[2]));

    // For just computing a square, we may use only one space for the
    // transform.
    void * tx = malloc(fft_alloc_sizes[0]);
    // void * ty = malloc(fft_alloc_sizes[0]);
    // void * tz = malloc(fft_alloc_sizes[0]);

    fft_transform_prepare(tx, fti);
    fft_do_dft(tx, PTR(x), SIZ(x), tt, fti);
    fft_mul(tx, tx, tx, tt, fti);
    // fft_add(tx, tx, tx, fti);
    // fft_add(tx, tx, tx, fti);
    // fft_add(tx, tx, tx, fti);
    fft_do_ift(PTR(z), SIZ(z), tx, tt, fti);

    MPN_NORMALIZE(PTR(z), SIZ(z));

    gmp_printf("%Zd\n", z);
    mpz_clear(x);
    mpz_clear(z);
    return 0;
}



int main()
{
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    // test_mul_fppol(rstate);
    test_mp_fppol(rstate);
    gmp_randclear(rstate);
    return 0;
}

