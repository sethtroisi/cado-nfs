#define __STDC_LIMIT_MACROS
#include "cado.h"
#include <string.h>
#include <limits.h>
#include <math.h>               /* ceil */

/* #undef HAVE_SSE2 */ /* Only for tests */ 

#ifdef HAVE_SSE41
#include <smmintrin.h>
#else
#ifdef HAVE_SSSE3
#include <tmmintrin.h>
#else
#ifdef HAVE_SSE3
#include <pmmintrin.h>
#else
#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif
#endif
#endif
#endif

#include "las-config.h"
#include "las-debug.h"
#include "las-norms.h"
#include "utils.h"
#include "mpz_poly.h"
#include "portability.h"

/* Input: i, double. i >= 0 needed!
   Output: o , trunc(o) == trunc(log2(i)) && o <= log2(i) < o + 0.0861.
   Careful: o ~= log2(i) iif add = 0x3FF00000 & scale = 1/0x100000.
   Add & scale are need to compute o'=f(log2(i)) where f is an affine function.
*/
static inline uint8_t inttruncfastlog2(double i, double add, double scale) {
#ifdef HAVE_SSE2
  __asm__ ( "psrlq $0x20,  %0    \n"
	    "cvtdq2pd      %0, %0\n" /* Mandatory in packed double even it's non packed! */
	    : "+x" (i));             /* Really need + here! i is not modified in C code! */
  return (uint8_t) ((i-add)*scale);
#else
  /* Same function, but in x86 gcc needs to transfer the input i from a
     xmm register to a classical register. No other way than use memory.
     So this function needs at least 6 to 8 cycles more than the previous,
     which uses ~3 cycles.
     NOTE: tg declaration is mandatory: it's the place where gcc use memory
     to do the transfert. Without it, a warning appears but the code is false!
  */
  void *tg = &i;
  return (uint8_t) ((((double) (*((uint64_t *)tg) >> 0x20)) - add) * scale);
#endif
}

/* Same than previous, but the result is duplicated 8 times in a "false" double
   in SSE version, i.e. in a xmm register, or in a 64 bits general register in
   non SSE version, and written at addr[decal].
   The SSE version has 4 interests: no memory use, no General Register use, 
   no xmm -> GR conversion, no * 0x0101010101010101 (3 to 4 cycles) but only
   2 P-instructions (1 cycle).
*/
static inline void uint64truncfastlog2(double i, double add, double scale, uint8_t *addr, ssize_t decal) {
#ifdef HAVE_SSE2
  __asm__ ( "psrlq $0x20,  %0                 \n"
	    "cvtdq2pd      %0,              %0\n"
	    : "+x" (i));
  i = (i - add) * scale;
  __asm__ ( 
	    "cvttpd2dq     %0,       %0       \n" /* 0000 0000 0000 000Y */
	    "punpcklbw     %0,       %0       \n" /* 0000 0000 0000 00YY */
	    "pshuflw    $0x00,       %0,    %0\n" /* 0000 0000 YYYY YYYY */
	    : "+x" (i));
  *(double *)&addr[decal] = i;
#else
  void *tg = &i;
  *(uint64_t *)&addr[decal] = 0x0101010101010101 * (uint64_t) (((double)(*(uint64_t *)tg >> 0x20) - add) * scale);
#endif
}

#ifdef HAVE_SSE2
/* These functions are for the SSE2 init algebraics.
   I prefer use X86 ASM directly and avoid intrinsics because the trick of
   cvtdq2pd (I insert 2 doubles values and do a cvtdq2pd on them in order to
   compute their log2).
*/
/* 2 values are computed + fabs at the beginning. These values are duplicated and
   written.
*/
static inline void w128itruncfastlog2fabs(__m128d i, __m128d add, __m128d scale, uint8_t *addr, ssize_t decal) {
  __asm__ (
	   "shufps    $0xED,    %0,    %0\n"
	   "pslld     $0x01,    %0       \n" /* Dont use pabsd! */
	   "psrld     $0x01,    %0       \n"
	   "cvtdq2pd  %0,       %0       \n"
	   : "+x"(i));
  i = _mm_mul_pd(_mm_sub_pd(i, add), scale);
  __asm__ (
	   "cvttpd2dq %0,       %0       \n" /* 0000 0000 000X 000Y */
	   "packssdw  %0,       %0       \n" /* 0000 0000 0000 0X0Y */
	   "punpcklbw %0,       %0       \n" /* 0000 0000 00XX 00YY */
	   "pshuflw   $0xA0,    %0,    %0\n" /* 0000 0000 XXXX YYYY */
	   "shufps    $0x50,    %0,    %0\n" /* XXXX XXXX YYYY YYYY */
	   : "+x"(i));
  *(__m128d *)&addr[decal] = i; /* addr and decal are 16 bytes aligned: MOVAPD */
}

/* 2 values are computed + fabs at the beginning. These values are written at scale+addr
   (2 bytes).
   NB: the double fabs is done by pslld 1 + psrld 1. */
static inline void w16itruncfastlog2fabs(__m128d i, __m128d add, __m128d scale, uint8_t *addr, ssize_t decal) {
  __asm__ (
	   "shufps    $0xED,    %0,    %0\n"
	   "pslld     $0x01,    %0       \n"
	   "psrld     $0x01,    %0       \n"
	   "cvtdq2pd  %0,       %0       \n"
	   : "+x" (i));
  i = _mm_mul_pd(_mm_sub_pd(i, add), scale);
  addr[decal] = (uint8_t) _mm_cvttsd_si32(i);
  i = _mm_unpackhi_pd(i,i);
  (&addr[decal])[1] = (uint8_t) _mm_cvttsd_si32(i);
}

/* Same than previous, but the 2 values are written at different places */
static inline void w8ix2truncfastlog2fabs(__m128d i, __m128d add, __m128d scale, uint8_t *addr, ssize_t decal1, ssize_t decal2) {
  __asm__ (
	   "shufps    $0xED,    %0,    %0\n"
	   "pslld     $0x01,    %0       \n"
	   "psrld     $0x01,    %0       \n"
	   "cvtdq2pd  %0,       %0       \n"
	   : "+x" (i));
  i = _mm_mul_pd(_mm_sub_pd(i, add), scale);
  addr[decal1] = (uint8_t) _mm_cvttsd_si32(i);
  i = _mm_unpackhi_pd(i,i);
  addr[decal2] = (uint8_t) _mm_cvttsd_si32(i);
}
#endif

/* Initialize lognorms on the rational side for the bucket_region
 * number N.
 * For the moment, nothing clever, wrt discarding (a,b) pairs that are
 * not coprime, except for the line j=0.
 */
void init_rat_norms_bucket_region(unsigned char *S,
                                 /* no condition array here ! */
                                 unsigned int j,
                                 sieve_info_ptr si)
{
#define COMPUTE_Y(G) ((LIKELY((G) > 1.0)) ? inttruncfastlog2 ((G), add, scale) : GUARD)
/* #define COMPUTE_Y(G) (inttruncfastlog2 ((G), add, scale)) */

  /* #define DEBUG_INIT_RAT 1 */ /* For internal debug: trace all */
  /* #define CHECK_INIT_RAT 1 */ /* For internal debug: control all */

  sieve_side_info_ptr rat = si->sides[RATIONAL_SIDE];
  int halfI = (si->I)>>1,
    int_i;
  uint8_t oy, y;
  double						\
    u0 = si->sides[RATIONAL_SIDE]->fijd[0], // gj
    u1 = si->sides[RATIONAL_SIDE]->fijd[1], // gi
    invu1 = 1.0/u1,
    u0j,
    d0_init,
    scale = rat->scale * (1.0/0x100000),
    add = 0x3FF00000 - GUARD / scale,
    g, rac, d0, d1, i;
  size_t ts;
  unsigned int						\
    j1 = LOG_BUCKET_REGION - si->conf->logI,
    inc;

  j = j << j1;
  j1 = (1U << j1) + j;
  u0j = u0 * j;
  d0_init = rat->cexp2[((unsigned int)GUARD) - 1U];
  if (!j) {
    // compute only the norm for i = 1. Everybody else is 255.
    memset(S, 255, halfI<<1);
    S[halfI + 1] = COMPUTE_Y(fabs(u1));
    S+= halfI<<1;
    j++;
    u0j += u0;
  }
  for( ; j < j1 ; j++, u0j += u0) {
#ifdef CHECK_INIT_RAT
    unsigned char *cS = S + halfI;
    memset (S, 0, halfI<<1);
#endif
    int_i = -halfI;
    g = u0j + u1 * int_i;
    rac = u0j * (-invu1);
    d0 = d0_init;
    d1 = rac - d0 * rac;
    /* g sign is mandatory at the beginning of the line intialization.
       If g ~= 0, the sign of g is not significant; so the sign of g
       is the sign of u1 (g+u1 is the next g value). CAREFUL:
       g ~== 0 is coded g + u0j == u0j. It's possible it's not suffisant;
       in this case, the right test will be g >= fabs(u0j)*(1.0/(1ULL<<51)) */
    if (LIKELY(g + u0j != u0j))
      if (signbit(g)) {
	g = -g;
	y = COMPUTE_Y(g);
	if (rac >= -halfI) goto cas3; else goto cas2;
      }
      else {
	y = COMPUTE_Y(g);
	if (rac >= -halfI) goto cas1; else goto cas4;
      }
    else {
      y = GUARD;
      if (signbit(u1)) goto cas2; else goto cas4;
    }
  cas1:
    /* In this case, we exit from the loop when ts == 0, at the exception
       of the first iteration. In this special case, old_i = -halfI and
       int_i = trunc (i), where i=[inverse of the function g](trunc(y)) and
       y=g(old_i).
       So, it's possible if y is very trunc trunc(y), old_i == int_i, so ts == 0.
       We have to iterate at least one time to avoid this case => this is the
       use of inc here. */
    for (i = rac + rat->cexp2[y] * invu1, inc = 1;; y--) {
      ts = -int_i;
      int_i = (int) trunc(i); 
      if (UNLIKELY(int_i >= halfI)) {
	ts += halfI;
#ifdef DEBUG_INIT_RAT
	fprintf (stderr, "A1.END : i1=%ld i2=%d, ts=%ld, y=%u, rac=%e\n", halfI - ts, halfI, ts, y, rac);
#endif
	memset(S, y, ts);
	S += ts;
	goto nextj;
      }
      ts += int_i;
#ifdef DEBUG_INIT_RAT
      fprintf (stderr, "A1 : i1=%ld, i2=%d, ts=%ld, y=%u, rac=%e\n", int_i - ts, int_i, ts, y, rac);
#endif
      if (LIKELY(ts <= MEMSET_MIN)) {
	if (!(ts + inc)) goto np1;
	memset(S, y, MEMSET_MIN);
      }
      else memset(S, y, ts);
      S += ts;
      i = i * d0 + d1;
      inc = 0;
    }
  np1:
    g = u0j + u1 * int_i;
    if (UNLIKELY(trunc(rac) >= halfI - 1)) {
      for ( ; int_i < halfI; int_i++) {
	y = COMPUTE_Y(g);
#ifdef DEBUG_INIT_RAT
	fprintf (stderr, "A2.1 : i=%d, y=%u, rac=%e\n", int_i, y, rac);
#endif
	*S++ = y;
	g += u1;
      }
      goto nextj;
    }
    for (inc = 0; g > 0; g += u1) {
      y = COMPUTE_Y(g);
#ifdef DEBUG_INIT_RAT
      fprintf (stderr, "A2.2 : i=%d, y=%u, rac=%e\n", int_i + inc, y, rac);
#endif
      S[inc++] = y;
    }
    int_i += inc;
    S += inc;
    g = -g;
    y = COMPUTE_Y(g);
  cas2:
    do {
#ifdef DEBUG_INIT_RAT
      fprintf (stderr, "A3 : i=%d, y=%u, rac=%e\n", int_i, y, rac);
#endif
      *S++ = y;
      if (++int_i >= halfI) goto nextj;
      oy = y;
      g -= u1;
      y = COMPUTE_Y(g);
    } while (oy != y);
    d0 = 1.0/d0;
    d1 = rac - d0 * rac;
    y++;
    i = rac - rat->cexp2[(unsigned int)y + 1] * invu1;
    for (;; y++) {
      ts = -int_i;
      int_i = (int) trunc(i);
      if (UNLIKELY(int_i >= halfI)) {
	ts += halfI;
#ifdef DEBUG_INIT_RAT
	fprintf (stderr, "A4.END : i1=%ld, i2=%d, ts=%ld, y=%u, rac=%e\n", halfI - ts, halfI, ts, y, rac);
#endif
	memset(S, y, ts);
	S += ts;
	goto nextj;
      }
      ts += int_i;
#ifdef DEBUG_INIT_RAT
      fprintf (stderr, "A4 : i1=%ld, i2=%d, ts=%ld, y=%u, rac=%e\n", int_i - ts, int_i, ts, y, rac);
#endif
      if (LIKELY(ts <= MEMSET_MIN))
	memset (S, y, MEMSET_MIN);
      else
	memset(S, y, ts);
      S += ts;
      i = i * d0 + d1;
    }

    /* Now, the same from cas1 but log2(-g): CAREFUL, not the same formula */
  cas3:
    for (i = rac - rat->cexp2[y] * invu1, inc = 1;; y--) {
      ts = -int_i;
      int_i = (int) trunc(i); 
      if (UNLIKELY(int_i >= halfI)) {
	ts += halfI;
#ifdef DEBUG_INIT_RAT
	fprintf (stderr, "B1.END : i1=%ld, i2=%d, ts=%ld, y=%u, rac=%e\n", halfI - ts, halfI, ts, y, rac);
#endif
	memset(S, y, ts);
	S += ts;
	goto nextj;
      }
      ts += int_i;
#ifdef DEBUG_INIT_RAT
      fprintf (stderr, "B1 : i1=%ld, i2=%d, ts=%ld, y=%u, rac=%e\n", int_i - ts, int_i, ts, y, rac);
#endif
      if (LIKELY(ts <= MEMSET_MIN)) {
	if (!(ts + inc)) goto np2;
	memset(S, y, MEMSET_MIN);
      }
      else memset(S, y, ts);
      S += ts;
      i = i * d0 + d1;
      inc = 0;
    }
  np2:
    g = -(u0j + u1 * int_i);
    if (UNLIKELY(trunc(rac) >= halfI - 1)) {
      for ( ; int_i < halfI; int_i++) {
	y = COMPUTE_Y(g);
#ifdef DEBUG_INIT_RAT
	fprintf (stderr, "B2.1 : i=%d, y=%u, rac=%e\n", int_i, y, rac);
#endif
	*S++ = y;
	g -= u1;
      }
      goto nextj;
    }
    for (inc = 0; g > 0; g -= u1) {
      y = COMPUTE_Y(g);
#ifdef DEBUG_INIT_RAT
      fprintf (stderr, "B2.2 : i=%d, y=%u, rac=%e\n", int_i + inc, y, rac);
#endif
      S[inc++] = y;
    }
    int_i += inc;
    S += inc;
    g = -g;
    y = COMPUTE_Y(g);
  cas4:
    do {
#ifdef DEBUG_INIT_RAT
      fprintf (stderr, "B3 : i=%d, y=%u, rac=%e\n", int_i, y, rac);
#endif
      *S++ = y;
      if (++int_i == halfI) goto nextj;
      oy = y;
      g += u1;
      y = COMPUTE_Y(g);
    } while (oy != y);
    d0 = 1.0/d0;
    d1 = rac - d0 * rac;
    y++;
    i = rac + rat->cexp2[(unsigned int)y + 1] * invu1;
    for (;; y++) {
      ts = -int_i;
      int_i = (int) trunc(i);
      if (UNLIKELY(int_i >= halfI)) {
	ts += halfI;
#ifdef DEBUG_INIT_RAT
	fprintf (stderr, "B4.END : i1=%ld i2=%d, ts=%ld, y=%u, rac=%e\n", halfI - ts, halfI, ts, y, rac);
#endif
	memset(S, y, ts);
	S += ts;
	goto nextj;
      }
      ts += int_i;
#ifdef DEBUG_INIT_RAT
      fprintf (stderr, "B4 : i1=%ld i2=%d, ts=%ld, y=%u, rac=%e\n", int_i - ts, int_i, ts, y, rac);
#endif
      if (LIKELY(ts <= MEMSET_MIN)) 
	memset(S, y, MEMSET_MIN);
      else
	memset(S, y, ts);
      S += ts;
      i = i * d0 + d1;
    }
  nextj:
    for (;0;); /* gcc needs something after a label */
#ifdef CHECK_INIT_RAT 
    /* First MANDATORY condition. The exact line must be initialised. */
    if (UNLIKELY(cS + halfI != S)) {
      fprintf (stderr, "init_rat_norms_bucket_region: S control Error: OldS(%p) + I(%d) != S(%p)", cS - halfI, si->I, S);
      exit (1);
    }
    /* Not really mandatory: when g ~= 0, the two formula (log2 & fastlog2) could really differ */
    int_i = -halfI;
    g = u0j + u1 * int_i;
    while (int_i < halfI) {
      y = LIKELY(fabs(g) > 1.0) ? log2(fabs(g))*rat->scale+GUARD : GUARD;
      if (UNLIKELY(fabs(cS[int_i] - y) > 1.)) {
	fprintf (stderr, "\
init_rat_norms_bucket_region: possible problem in S, offset %d:\n	\
   real value=%d, S init value=%d. If g + uj0 ~= uj0,\n			\
   it could be OK: g=%e, u0j=%e\n", int_i, y, cS[int_i], g, u0j);
      }
      int_i++;
      g += u1;
    }
#endif
  }
}

static inline void fpoly_scale(double * u, const double * t, unsigned int d, double h)
{
  double hpow;
  u[d] = t[d];
  for (hpow = h; --d != UINT_MAX; hpow *= h) u[d] = t[d] * hpow;
}


/**************************************************************************
       8 algorithms for the initialization of the algebraics:
                 init_alg_norms_bucket_region.
       All the combinations with HAVE_SSE2, ALG_LAZY, ALG_RAT.
       All these 8 algorithms use a big switch for the degree of
       the algebraics; all the cases are independant and end by
                      a return for lisibility.
       All have a default which could replace all others cases.
     NB: Somes SSE algorithms needs __x86_64; all __x86_64 are SSE2.
**************************************************************************/

/**************** used #define from SSE2 to -41 **************************/
#ifdef HAVE_SSE2

#ifdef HAVE_SSE41
#define TSTZXMM(A) _mm_testz_si128(A,A)
#else
#define TSTZXMM(A) ((uint16_t)_mm_movemask_epi8(_mm_cmpeq_epi8(_mm_add_epi8(A,tt),tt))==0xFFFF)
#endif

/* This macro avoids a stupid & boring C types control */
#define _MM_SHUFFLE_EPI32(A,B,C) __asm__ ("pshufd $" #C ", %1, %0\n":"=x"(A):"x"(B))

/* CAREFUL! __m128 a=x|y => x is the strongest quadword, higher in memory.
   INITALGN: uN= (j^N)*(alg->fijd[N]) | (j^N)*(alg->fijd[N])  
   INITALG0 to 2 are simple. 3,5,7,9 & 4,6,8 have similar algorithms. */
#define INITALG0 							\
  uu0=_mm_set1_epi64((__m64)(0x0101010101010101*			\
			     inttruncfastlog2(fabs(alg->fijd[0]),*(double *)&add,*(double *)&scale)))

#define INITALG1(A)				\
  u1=_mm_load1_pd(alg->fijd+1);			\
  u0=_mm_set1_pd((double)(A)*alg->fijd[0])

#define INITALG2(A)				\
  h=_mm_set1_pd((double)(A));			\
  u2=_mm_load1_pd(alg->fijd+2);			\
  u1=_mm_load1_pd(alg->fijd+1);			\
  u0=_mm_load1_pd(alg->fijd);			\
  u1=_mm_mul_pd(u1,h);				\
  u0=_mm_mul_pd(u0,h);				\
  u0=_mm_mul_pd(u0,h)

#define INITALG3(A)						\
  h=_mm_set_sd((double)(A));     /* h=0|j */			\
  _MM_SHUFFLE_EPI32(g, h, 0x44);				\
  g=_mm_mul_pd(g,g);             /* g=j2|j2 */			\
  u2=_mm_load_pd(alg->fijd+2);   /* u2=d3|d2 */			\
  u0=_mm_load_pd(alg->fijd);     /* u0=d1|d0 */			\
  u2=_mm_mul_sd(u2,h);           /* u2=d3|d2*j */		\
  u0=_mm_mul_sd(u0,h);						\
  u0=_mm_mul_pd(u0,g);						\
  _MM_SHUFFLE_EPI32(u1, u0, 0xEE);				\
  u0=_mm_unpacklo_pd(u0,u0);					\
  _MM_SHUFFLE_EPI32(u3, u2, 0xEE);				\
  u0=_mm_unpacklo_pd(u2,u2)

#define INITALG4(A)						\
  h=_mm_set1_pd((double)(A));					\
  h=_mm_mul_sd(h,h);             /* h=j|j2 */			\
  _MM_SHUFFLE_EPI32(g, h, 0x44); /* g=j2|j2 */			\
  u4=_mm_load1_pd(alg->fijd+4);					\
  u2=_mm_load_pd(alg->fijd+2);					\
  u0=_mm_load_pd(alg->fijd);					\
  u2=_mm_mul_pd(u2,h);						\
  u0=_mm_mul_pd(u0,h);						\
  u0=_mm_mul_pd(u0,g);						\
  _MM_SHUFFLE_EPI32(u1, u0, 0xEE);				\
  u0=_mm_unpacklo_pd(u0,u0);					\
  _MM_SHUFFLE_EPI32(u3, u2, 0xEE);				\
  u2=_mm_unpacklo_pd(u2,u2)

#define INITALG5(A)						\
  h=_mm_set_sd((double)(A));					\
  _MM_SHUFFLE_EPI32(g, h, 0x44);				\
  g=_mm_mul_pd(g,g);						\
  u4=_mm_load_pd(alg->fijd+4);					\
  u2=_mm_load_pd(alg->fijd+2);					\
  u0=_mm_load_pd(alg->fijd);					\
  u4=_mm_mul_sd(u4,h);						\
  u2=_mm_mul_sd(u2,h);						\
  u0=_mm_mul_sd(u0,h);						\
  u2=_mm_mul_pd(u2,g);						\
  u0=_mm_mul_pd(u0,g);						\
  u0=_mm_mul_pd(u0,g);						\
  _MM_SHUFFLE_EPI32(u1, u0, 0xEE);				\
  u0=_mm_unpacklo_pd(u0,u0);					\
  _MM_SHUFFLE_EPI32(u3, u2, 0xEE);				\
  u2=_mm_unpacklo_pd(u2,u2);					\
  _MM_SHUFFLE_EPI32(u5, u4, 0xEE);				\
  u4=_mm_unpacklo_pd(u4,u4)

#define INITALG6(A)						\
  h=_mm_set1_pd((double)(A));					\
  h=_mm_mul_sd(h,h);						\
  _MM_SHUFFLE_EPI32(g, h, 0x44);				\
  u6=_mm_load1_pd(alg->fijd+6);					\
  u4=_mm_load_pd(alg->fijd+4);					\
  u2=_mm_load_pd(alg->fijd+2);					\
  u0=_mm_load_pd(alg->fijd);					\
  u4=_mm_mul_pd(u4,h);						\
  u2=_mm_mul_pd(u2,h);						\
  u0=_mm_mul_pd(u0,h);						\
  u2=_mm_mul_pd(u2,g);						\
  u0=_mm_mul_pd(u0,g);						\
  u0=_mm_mul_pd(u0,g);						\
  _MM_SHUFFLE_EPI32(u1, u0, 0xEE);				\
  u0=_mm_unpacklo_pd(u0,u0);					\
  _MM_SHUFFLE_EPI32(u3, u2, 0xEE);				\
  u2=_mm_unpacklo_pd(u2,u2);					\
  _MM_SHUFFLE_EPI32(u5, u4, 0xEE);				\
  u4=_mm_unpacklo_pd(u4,u4)

#define INITALG7(A)						\
  h=_mm_set_sd((double)(A));					\
  _MM_SHUFFLE_EPI32(g, h, 0x44);				\
  g=_mm_mul_pd(g,g);						\
  u6=_mm_load_pd(alg->fijd+6);					\
  u4=_mm_load_pd(alg->fijd+4);					\
  u2=_mm_load_pd(alg->fijd+2);					\
  u0=_mm_load_pd(alg->fijd);					\
  u6=_mm_mul_sd(u6,h);						\
  u4=_mm_mul_sd(u4,h);						\
  u2=_mm_mul_sd(u2,h);						\
  u0=_mm_mul_sd(u0,h);						\
  u4=_mm_mul_pd(u4,g);						\
  u0=_mm_mul_pd(u0,g);						\
  g=_mm_mul_pd(g,g);						\
  u2=_mm_mul_pd(u2,g);						\
  u0=_mm_mul_pd(u0,g);						\
  _MM_SHUFFLE_EPI32(u1, u0, 0xEE);				\
  u0=_mm_unpacklo_pd(u0,u0);					\
  _MM_SHUFFLE_EPI32(u3, u2, 0xEE);				\
  u2=_mm_unpacklo_pd(u2,u2);					\
  _MM_SHUFFLE_EPI32(u5, u4, 0xEE);				\
  u4=_mm_unpacklo_pd(u4,u4);					\
  _MM_SHUFFLE_EPI32(u7, u6, 0xEE);				\
  u6=_mm_unpacklo_pd(u6,u6)

#define INITALG8(A)						\
  h=_mm_set1_pd((double)(A));					\
  h=_mm_mul_sd(h,h);						\
  _MM_SHUFFLE_EPI32(g, h, 0x44);				\
  u8=_mm_load1_pd(alg->fijd+8);					\
  u6=_mm_load_pd(alg->fijd+6);					\
  u4=_mm_load_pd(alg->fijd+4);					\
  u2=_mm_load_pd(alg->fijd+2);					\
  u0=_mm_load_pd(alg->fijd);					\
  u6=_mm_mul_pd(u6,h);						\
  u4=_mm_mul_pd(u4,h);						\
  u2=_mm_mul_pd(u2,h);						\
  u0=_mm_mul_pd(u0,h);						\
  u4=_mm_mul_pd(u4,g);						\
  u0=_mm_mul_pd(u0,g);						\
  g=_mm_mul_pd(g,g);						\
  u2=_mm_mul_pd(u2,g);						\
  u0=_mm_mul_pd(u0,g);						\
  _MM_SHUFFLE_EPI32(u1, u0, 0xEE);				\
  u0=_mm_unpacklo_pd(u0,u0);					\
  _MM_SHUFFLE_EPI32(u3, u2, 0xEE);				\
  u2=_mm_unpacklo_pd(u2,u2);					\
  _MM_SHUFFLE_EPI32(u5, u4, 0xEE);				\
  u4=_mm_unpacklo_pd(u4,u4);					\
  _MM_SHUFFLE_EPI32(u7, u6, 0xEE);				\
  u6=_mm_unpacklo_pd(u6,u6)

#define INITALG9(A)						\
  h=_mm_set_sd((double)(A));					\
  _MM_SHUFFLE_EPI32(g, h, 0x44);				\
  g=_mm_mul_pd(g,g);						\
  u8=_mm_load_pd(alg->fijd+8);					\
  u6=_mm_load_pd(alg->fijd+6);					\
  u4=_mm_load_pd(alg->fijd+4);					\
  u2=_mm_load_pd(alg->fijd+2);					\
  u0=_mm_load_pd(alg->fijd);					\
  u8=_mm_mul_sd(u8,h);						\
  u6=_mm_mul_sd(u6,h);						\
  u4=_mm_mul_sd(u4,h);						\
  u2=_mm_mul_sd(u2,h);						\
  u0=_mm_mul_sd(u0,h);						\
  u6=_mm_mul_pd(u6,g);						\
  u2=_mm_mul_pd(u2,g);						\
  g=_mm_mul_pd(g,g);						\
  u4=_mm_mul_pd(u4,g);						\
  u2=_mm_mul_pd(u2,g);						\
  g=_mm_mul_pd(g,g);						\
  u0=_mm_mul_pd(u0,g);						\
  _MM_SHUFFLE_EPI32(u1, u0, 0xEE);				\
  u0=_mm_unpacklo_pd(u0,u0);					\
  _MM_SHUFFLE_EPI32(u3, u2, 0xEE);				\
  u2=_mm_unpacklo_pd(u2,u2);					\
  _MM_SHUFFLE_EPI32(u5, u4, 0xEE);				\
  u4=_mm_unpacklo_pd(u4,u4);					\
  _MM_SHUFFLE_EPI32(u7, u6, 0xEE);				\
  u6=_mm_unpacklo_pd(u6,u6);					\
  _MM_SHUFFLE_EPI32(u9, u8, 0xEE);				\
  u8=_mm_unpacklo_pd(u8,u8)

#define INITALGD(A) do {					\
    double du[d+1];						\
    fpoly_scale(du, alg->fijd, d, (double) (A));		\
    for (uint32_t k=0; k<=d; k++) u[k] = _mm_set1_pd(du[k]);	\
  } while (0)

#define G1 g=_mm_add_pd(_mm_mul_pd(g,u1),u0)

#define G2 g=_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(g,u2),u1),h),u0)

#define G3 g=_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(g,u3),u2),h),u1),h),u0)

#define G4 g=_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(g,u4),u3),h),u2),h),u1),h),u0)

#define G5 g=_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(g,u5),u4),h),u3),h),u2),h),u1),h),u0)

#define G6 g=_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(g,u6),u5),h),u4),h),u3),h),u2),h),u1),h),u0)

#define G7 g=_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(g,u7),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0)

#define G8 g=_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(g,u8),u7),h),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0)

#define G9 g=_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(g,u9),u8),h),u7),h),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0)

#define GD g=u[d-1];for (unsigned int k=d;--k!=UINT_MAX;g=_mm_add_pd(_mm_mul_pd(g,h),u[k]))

#define ALG_INIT_SCALE_ADD_A						\
  __m128d scale = _mm_set1_pd(alg->scale * (1.0/0x100000)),             \
    add = _mm_set1_pd(0x3FF00000 - GUARD / (alg->scale * (1.0/0x100000))); \
  const __m128d a MAYBE_UNUSED = _mm_set1_pd(16.)

#define ALG_INIT_SCALE_ADD_A_TT						\
  ALG_INIT_SCALE_ADD_A;							\
  __m128i tt = _mm_set1_epi64((__m64)(0x0101010101010101* (uint64_t)    \
				      ((si->sides[RATIONAL_SIDE]->bound + \
					(si->sides[RATIONAL_SIDE]->bound != 255 ? 1 : 0)))))

#endif /* HAVE_SSE2 */
/**************** End of used #define SSE ********************************/

/**************** 1: HAVE_SSE2, ALG_LAZY, ALG_RAT ************************/
#ifdef HAVE_SSE2
#ifdef ALG_LAZY 
#ifdef ALG_RAT
/* Smart initialization of the algebraics. In one time, 16 bytes of
   the rationals are tested. If one byte is > rat->bound, computes two
   initialization (3th and 11th bytes) and initializes (8+8) algrebraics.
   At the beginning, all are initialised to 255.
   2nd fastest initialization, but maybe the smartest. SSE version.
   32 bits compatible.
   NOTA : the SSE and non SSE versions differs! */
void init_alg_norms_bucket_region(unsigned char *S,
				  unsigned char *xS, unsigned int j,
				  sieve_info_ptr si)
{
  sieve_side_info_ptr alg = si->sides[ALGEBRAIC_SIDE];
  ssize_t ih, Idiv2 = (si->I >> 1);
  unsigned int ej = 1U << (LOG_BUCKET_REGION - si->conf->logI),
    d = si->cpoly->alg->degree;
  __m128i m;
  ALG_INIT_SCALE_ADD_A_TT;

  memset (S, 255, (Idiv2<<1) * ej);
  j *= ej;
  ej += j;
  S += Idiv2;
  xS += Idiv2;
  switch (d) {
  case 0 : {
    __m128i uu0;
    INITALG0;
    for (; j < ej; j++, S += (Idiv2<<1), xS += (Idiv2<<1)) {
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
        __asm__("movdqa %1,%0\npsubusb (%2,%3),%0\n":"=&x"(m):"x"(tt),"r"(xS),"r"(ih));
        if (UNLIKELY(!TSTZXMM(m))) {
          __asm__( "movdqa %0,(%1,%2)\n"::"x"(uu0),"r"(S),"r"(ih));
        }}}}
    return;
  case 1 : {
    __m128d g, h, u0, u1;
    for (; j < ej; j++, S += (Idiv2<<1), xS += (Idiv2<<1)) {
      INITALG1(j);
      h =_mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
        __asm__("movdqa %1,%0\npsubusb (%2,%3),%0\n":"=&x"(m):"x"(tt),"r"(xS),"r"(ih));
        if (UNLIKELY(!TSTZXMM(m))) {
	  G1;
          w128itruncfastlog2fabs(g, add, scale, S, ih);
        }
	h = _mm_add_pd(h, a);
      }}}
    return;
  case 2 : {
    __m128d g, h, u0, u1, u2;
    for (; j < ej; j++, S += (Idiv2<<1), xS += (Idiv2<<1)) {
      INITALG2(j);
      h =_mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
        __asm__("movdqa %1,%0\npsubusb (%2,%3),%0\n":"=&x"(m):"x"(tt),"r"(xS),"r"(ih));
        if (UNLIKELY(!TSTZXMM(m))) {
          g = h; G2;
	  w128itruncfastlog2fabs(g, add, scale, S, ih);
	}
	h = _mm_add_pd(h, a);
      }}}
    return;
  case 3 : {
    __m128d g, h, u0, u1, u2, u3;
    for (; j < ej; j++, S += (Idiv2<<1), xS += (Idiv2<<1)) {
      INITALG3(j);
      h =_mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	__asm__("movdqa %1,%0\npsubusb (%2,%3),%0\n":"=&x"(m):"x"(tt),"r"(xS),"r"(ih));
	if (UNLIKELY(!TSTZXMM(m))) {
	  g = h; G3;
	  w128itruncfastlog2fabs(g, add, scale, S, ih);
	}
	h = _mm_add_pd(h, a);
      }}}
    return;
  case 4 : {
    __m128d g, h, u0, u1, u2, u3, u4;
    for (; j < ej; j++, S += (Idiv2<<1), xS += (Idiv2<<1)) {
      INITALG4(j);
      h =_mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	__asm__("movdqa %1,%0\npsubusb (%2,%3),%0\n":"=&x"(m):"x"(tt),"r"(xS),"r"(ih));
	if (UNLIKELY(!TSTZXMM(m))) {
	  g = h; G4;
	  w128itruncfastlog2fabs(g, add, scale, S, ih);
	}
	h = _mm_add_pd(h, a);
      }}}
    return;
  case 5 : {
    __m128d g, h, u0, u1, u2, u3, u4, u5;
    for (; j < ej; j++, S += (Idiv2<<1), xS += (Idiv2<<1)) {
      INITALG5(j);
      h =_mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	__asm__("movdqa %1,%0\npsubusb (%2,%3),%0\n":"=&x"(m):"x"(tt),"r"(xS),"r"(ih));
	if (UNLIKELY(!TSTZXMM(m))) {
	  g = h; G5;
	  w128itruncfastlog2fabs(g, add, scale, S, ih);
	}
	h = _mm_add_pd(h, a);
      }}}
    return;
  case 6 : {
    __m128d g, h, u0, u1, u2, u3, u4, u5, u6;
    for (; j < ej; j++, S += (Idiv2<<1), xS += (Idiv2<<1)) {
      INITALG6(j);
      h =_mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	__asm__("movdqa %1,%0\npsubusb (%2,%3),%0\n":"=&x"(m):"x"(tt),"r"(xS),"r"(ih));
	if (UNLIKELY(!TSTZXMM(m))) {
 	  g = h; G6;
	  w128itruncfastlog2fabs(g, add, scale, S, ih);
	}
	h = _mm_add_pd(h, a);
      }}}
    return;
  case 7 : {
    __m128d g, h, u0, u1, u2, u3, u4, u5, u6, u7;
    for (; j < ej; j++, S += (Idiv2<<1), xS += (Idiv2<<1)) {
      INITALG7(j);
      h =_mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	__asm__("movdqa %1,%0\npsubusb (%2,%3),%0\n":"=&x"(m):"x"(tt),"r"(xS),"r"(ih));
	if (UNLIKELY(!TSTZXMM(m))) {
	  g = h; G7;
	  w128itruncfastlog2fabs(g, add, scale, S, ih);
	}}}}
    return;
  case 8 : {
    __m128d g, h, u0, u1, u2, u3, u4, u5, u6, u7, u8;
    for (; j < ej; j++, S += (Idiv2<<1), xS += (Idiv2<<1)) {
      INITALG8(j);
      h =_mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	__asm__("movdqa %1,%0\npsubusb (%2,%3),%0\n":"=&x"(m):"x"(tt),"r"(xS),"r"(ih));
	if (UNLIKELY(!TSTZXMM(m))) {
	  g = h; G8;
	  w128itruncfastlog2fabs(g, add, scale, S, ih);
	}
	h = _mm_add_pd(h, a);
      }}}
    return;
  case 9 : {
    __m128d g, h, u0, u1, u2, u3, u4, u5, u6, u7, u8, u9;
    for (; j < ej; j++, S += (Idiv2<<1), xS += (Idiv2<<1)) {
      INITALG9(j);
      h =_mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	__asm__("movdqa %1,%0\npsubusb (%2,%3),%0\n":"=&x"(m):"x"(tt),"r"(xS),"r"(ih));
	if (UNLIKELY(!TSTZXMM(m))) {
	  g = h; G9;
	  w128itruncfastlog2fabs(g, add, scale, S, ih);
	}
	h = _mm_add_pd(h, a);
      }}}
    return;
  default: {
    __m128d g, h, u[d+1];
    for (; j < ej; j++, S += (Idiv2<<1), xS += (Idiv2<<1)) {
      INITALGD(j);
      h =_mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	__asm__("movdqa %1,%0\npsubusb (%2,%3),%0\n":"=&x"(m):"x"(tt),"r"(xS),"r"(ih));
	if (UNLIKELY(!TSTZXMM(m))) {
	  g = h; GD;
	  w128itruncfastlog2fabs(g, add, scale, S, ih);
	}
	h = _mm_add_pd(h, a);
      }}}
    return;
  }
}
#endif
#endif
#endif

/**************** 2: HAVE_SSE2, ALG_LAZY, !ALG_RAT ***********************/
#ifdef HAVE_SSE2
#ifdef ALG_LAZY
#ifndef ALG_RAT
/* Smart initialization of the algebraics. Computes the central
   initialization value of a box of (horizontail=8,vertical=p) and
   propagates it in the box. 8 is locked by the code and p is the
   minimum of VERT_NORM_STRIDE (#define, =4) and ej, with ej =
   2^(LOG_BUCKET_REGION-si->conf->logI).
   So, the fastest code is for logI <= 14.
   It's the fastest initialization of the algebraics.
   SSE version, 32bits compatible. */
void init_alg_norms_bucket_region(unsigned char *S,
				  unsigned char *xS MAYBE_UNUSED, unsigned int j,
				  sieve_info_ptr si)
{
  sieve_side_info_ptr alg = si->sides[ALGEBRAIC_SIDE];
  ssize_t ih, Idiv2 = (si->I >> 1);
  unsigned int ej = 1U << (LOG_BUCKET_REGION - si->conf->logI),
    d = si->cpoly->alg->degree,
    p = ej < VERT_NORM_STRIDE ? ej : VERT_NORM_STRIDE;
  ALG_INIT_SCALE_ADD_A;

#define FILL_OTHER_LINES						\
  if (LIKELY(p > 1)) {							\
    unsigned char *mS = S + Idiv2;					\
    S -= Idiv2;								\
    for (unsigned int k = 1; k < p; k++, mS += (Idiv2<<1))		\
      memcpy (mS, S, (Idiv2<<1));					\
    S = mS + Idiv2;							\
    j += p;								\
  }									\
  else do {								\
      S += (Idiv2<<1);							\
      j++;								\
    } while(0)
  
  j *= ej;
  ej += j;
  S += Idiv2;
  switch (d) {
  case 0 :
    memset(S-Idiv2,inttruncfastlog2(fabs(alg->fijd[0]),*(double *)&add,*(double *)&scale),(Idiv2<<1)*(ej-j));
    return;
  case 1 : {
    __m128d h, g, u0, u1;
    while (j < ej) {
      INITALG1(j + (p >> 1));
      h = _mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	g = h; G1;
	w128itruncfastlog2fabs(g, add, scale, S, ih);
	h = _mm_add_pd(h, a);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 2 : {
    __m128d h, g, u0, u1, u2;
    while (j < ej) {
      INITALG2(j + (p >> 1));
      h = _mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	g = h; G2;
	w128itruncfastlog2fabs(g, add, scale, S, ih);
	h = _mm_add_pd(h, a);
      }
      FILL_OTHER_LINES;
    }}	
    return;
  case 3 : {
    __m128d h, g, u0, u1, u2, u3;
    while (j < ej) {
      INITALG3(j + (p >> 1));
      h = _mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	g = h; G3;
	w128itruncfastlog2fabs(g, add, scale, S, ih);
	h = _mm_add_pd(h, a);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 4 : {
    __m128d h, g, u0, u1, u2, u3, u4;
    while (j < ej) {
      INITALG4(j + (p >> 1));
      h = _mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	g = h; G4;
	w128itruncfastlog2fabs(g, add, scale, S, ih);
	h = _mm_add_pd(h, a);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 5 : {
    __m128d h, g, u0, u1, u2, u3, u4, u5;
    while (j < ej) {
      INITALG5(j + (p >> 1));
      h = _mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	g = h; G5;
	w128itruncfastlog2fabs(g, add, scale, S, ih);
	h = _mm_add_pd(h, a);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 6 : {
    __m128d h, g, u0, u1, u2, u3, u4, u5, u6;
    while (j < ej) {
      INITALG6(j + (p >> 1));
      h =_mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	g = h; G6;
	w128itruncfastlog2fabs(g, add, scale, S, ih);
	h = _mm_add_pd(h, a);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 7 : {
    __m128d h, g, u0, u1, u2, u3, u4, u5, u6, u7;
    while (j < ej) {
      INITALG7(j + (p >> 1));
      h =_mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	g = h; G7;
	w128itruncfastlog2fabs(g, add, scale, S, ih);
	h = _mm_add_pd(h, a);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 8 : {
    __m128d h, g, u0, u1, u2, u3, u4, u5, u6, u7, u8;
    while (j < ej) {
      INITALG8(j + (p >> 1));
      h =_mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	g = h; G8;
	w128itruncfastlog2fabs(g, add, scale, S, ih);
	h = _mm_add_pd(h, a);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 9 : {
    __m128d h, g, u0, u1, u2, u3, u4, u5, u6, u7, u8, u9;
    while (j < ej) {
      INITALG9(j + (p >> 1));
      h =_mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	g = h; G9;
	w128itruncfastlog2fabs(g, add, scale, S, ih);
	h = _mm_add_pd(h, a);
      }
      FILL_OTHER_LINES;
    }}
    return;
  default: {
    __m128d h, g, u[d+1];
    while (j < ej) {
      INITALGD(j + (p >> 1));
      h =_mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	g = h; GD;
	w128itruncfastlog2fabs(g, add, scale, S, ih);
	h = _mm_add_pd(h, a);
      }
      FILL_OTHER_LINES;
    }}
    return;
  }
}
#undef FILL_OTHER_LINES
#endif 
#endif
#endif

/**************** 3: HAVE_SSE2, !ALG_LAZY, ALG_RAT **************************/
#ifdef NOTDEFINED /* __x86_64 */ /* Not active: code slower than non SSE version */
#ifndef ALG_LAZY
#ifdef ALG_RAT
/* Smart initialization of the algebraics : only if the corresponding
   rationnal is <= rat->bound, other 255 (by first memset), so ~ x 22 less
   computation than systematic computation; but requires of course many
   tests. SSE version.

   CAREFUL: only x86-64 version! It's not possible to code this algorithm
   with the intrinsics because it uses the carry flag with the instruction
   rcr (rotates on 8/16/32/64 + carry flag) & immediatly a jump if not carry.
   GCC cannot manipulate carry flag - it's a well known limitation.

   NB: this algorithm is LESS EFFICIENT on my machine than the non SSE
   version. Why ? Because it needs to find 2 interesting rationals
   before the vectorised computation of the both corresponding algebraics.
   So, I need 5 jumps to find these 2 algebraics in this code (1 when found,
   1 to code its emplacement)x2, and one to return after computation.
   In non SSE version, only 2 jumps by algebraic : one when found and one
   to return after computation.
   
   I try to code it without any jump except loop with conditionnal move.
   The idea is to fill an temp array with all the i value (from [-Idiv2
   to Idiv2[) but the pointer is << conditionaly >> incremented. The
   main asm code :
   
#define ALG_1BYTE                                               \
  "xor %%ebx, %%ebx\n"                                          \
  "mov %%eax, (%0)\n"                                           \
  "rcr $1, %%ecx\n"                                             \
  "rcl $3, %%ebx\n"                                             \
  "inc %%rax\n"                                                 \
  "add %%rbx, %0\n"
           
#define ALG_IL_P                                                \
           __asm__ (                                            \
           "pxor %%xmm1, %%xmm1\n"                              \
           "movq %2, %0\n"                                      \
           "xor %%rbx, %%rbx\n"                                 \
           "xor %%rax, %%rax\n"                                 \
           "sub %3, %%rax\n"                                    \
           ".balign 8\n"                                        \
           "0: movdqa %4, %%xmm0\n"                             \
           "psubusb (%1, %%rax), %%xmm0\n"                      \
           "pcmpeqb %%xmm1, %%xmm0\n"                           \
           "pmovmskb %%xmm0, %%ecx\n"                           \
           "xor $0xFFFFFFFF, %%ecx\n"                           \
           ALG_1BYTE                                            \
           ALG_1BYTE                                            \
           ALG_1BYTE                                            \
           ALG_1BYTE                                            \
           ALG_1BYTE                                            \
           ALG_1BYTE                                            \
           ALG_1BYTE                                            \
           ALG_1BYTE                                            \
           ALG_1BYTE                                            \
           ALG_1BYTE                                            \
           ALG_1BYTE                                            \
           ALG_1BYTE                                            \
           ALG_1BYTE                                            \
           ALG_1BYTE                                            \
           ALG_1BYTE                                            \
           ALG_1BYTE                                            \
           "cmp %3, %%rax\n"                                    \
           "jl 0b\n"                                            \
           : "=&r"(cP)                                          \
           : "r"(xS), "r"(P), "r"(Idiv2), "x"(tt)               \
           : "%rax", "%rbx", "%rcx", "%xmm0", "%xmm1");

   But the result is x1.5 slower than this version, because it needs
   6 instructions to look at only one rational (ALG_1BYTE).
*/
   
void init_alg_norms_bucket_region (unsigned char *S,
				   unsigned char *xS, unsigned int j,
				   sieve_info_ptr si)
{
  sieve_side_info_ptr alg = si->sides[ALGEBRAIC_SIDE];
  ssize_t Idiv2 = (si->I >> 1);
  unsigned int ej = 1U << (LOG_BUCKET_REGION - si->conf->logI),
    d = si->cpoly->alg->degree;
  ALG_INIT_SCALE_ADD_A_TT;

  /******************* Beginning of local define ****************/
#define ALG_IL0A "movapd %6, %%xmm1\n"
#define ALG_IL1A "movapd %7, %%xmm1\n"
#define ALG_IL2A "movapd %8, %%xmm1\n"
#define ALG_IL3A "movapd %9, %%xmm1\n"
#define ALG_IL4A "movapd %10, %%xmm1\n"
#define ALG_IL5A "movapd %11, %%xmm1\n"
#define ALG_IL6A "movapd %12, %%xmm1\n"
#define ALG_IL7A "movapd %13, %%xmm1\n"
#define ALG_IL8A "movapd %14, %%xmm1\n"
#define ALG_IL9A "movapd %15, %%xmm1\n"

#define ALG_IL1B "mulpd %%xmm0, %%xmm1\n addpd %6, %%xmm1\n"
#define ALG_IL2B "mulpd %%xmm0, %%xmm1\n addpd %7, %%xmm1\n" ALG_IL1B 
#define ALG_IL3B "mulpd %%xmm0, %%xmm1\n addpd %8, %%xmm1\n" ALG_IL2B 
#define ALG_IL4B "mulpd %%xmm0, %%xmm1\n addpd %9, %%xmm1\n" ALG_IL3B 
#define ALG_IL5B "mulpd %%xmm0, %%xmm1\n addpd %10, %%xmm1\n" ALG_IL4B 
#define ALG_IL6B "mulpd %%xmm0, %%xmm1\n addpd %11, %%xmm1\n" ALG_IL5B 
#define ALG_IL7B "mulpd %%xmm0, %%xmm1\n addpd %12, %%xmm1\n" ALG_IL6B 
#define ALG_IL8B "mulpd %%xmm0, %%xmm1\n addpd %13, %%xmm1\n" ALG_IL7B 
#define ALG_IL9B "mulpd %%xmm0, %%xmm1\n addpd %14, %%xmm1\n" ALG_IL8B 

#define ALG_IL1C "mulsd %%xmm0, %%xmm1\n addsd %6, %%xmm1\n"
#define ALG_IL2C "mulsd %%xmm0, %%xmm1\n addsd %7, %%xmm1\n" ALG_IL1B 
#define ALG_IL3C "mulsd %%xmm0, %%xmm1\n addsd %8, %%xmm1\n" ALG_IL2B 
#define ALG_IL4C "mulsd %%xmm0, %%xmm1\n addsd %9, %%xmm1\n" ALG_IL3B 
#define ALG_IL5C "mulsd %%xmm0, %%xmm1\n addsd %10, %%xmm1\n" ALG_IL4B 
#define ALG_IL6C "mulsd %%xmm0, %%xmm1\n addsd %11, %%xmm1\n" ALG_IL5B 
#define ALG_IL7C "mulsd %%xmm0, %%xmm1\n addsd %12, %%xmm1\n" ALG_IL6B 
#define ALG_IL8C "mulsd %%xmm0, %%xmm1\n addsd %13, %%xmm1\n" ALG_IL7B 
#define ALG_IL9C "mulsd %%xmm0, %%xmm1\n addsd %14, %%xmm1\n" ALG_IL8B 

#define ALG_ILD							\
  "mov %7, %%r13\n"						\
  "movapd (%6,%%r13), %%xmm1\n"					\
  ".balign 8\n"							\
  "7: test %%r13d, %%r13d\n"					\
  "jz 8f\n"							\
  "sub $0x10, %%r13d\n"						\
  "mulpd %%xmm0, %%xmm1\n"					\
  "addpd (%6,%%r13), %%xmm1\n"					\
  "jmp 7b\n"							\
  "8:\n"

#define ALG_ILP1							\
  __asm__ ( "pxor %%xmm2, %%xmm2\n"					\
  "xor %%r12, %%r12\n"							\
  "sub %2, %%r12\n"							\
  "jmp 1f\n"								\
									\
  ".balign 8\n"								\
  "0: add $0x10, %%r12\n"						\
  "cmp %2, %%r12\n"							\
  "jnl 6f\n"                               /* Exit */			\
  "1: movdqa %3, %%xmm0\n"						\
  "psubusb (%0, %%r12), %%xmm0\n"					\
  "pcmpeqb %%xmm2, %%xmm0\n"						\
  "pmovmskb %%xmm0, %%eax\n"						\
  ".balign 4\n"								\
									\
  "rcr $0x1, %%eax\n jnc 100f\n"					\
  "2: rcr $0x1, %%eax\n jnc 101f\n"					\
  "rcr $0x1, %%eax\n jnc 102f\n"					\
  "rcr $0x1, %%eax\n jnc 103f\n"					\
  "rcr $0x1, %%eax\n jnc 104f\n"					\
  "rcr $0x1, %%eax\n jnc 105f\n"					\
  "rcr $0x1, %%eax\n jnc 106f\n"					\
  "rcr $0x1, %%eax\n jnc 107f\n"					\
  "rcr $0x1, %%eax\n jnc 108f\n"					\
  "rcr $0x1, %%eax\n jnc 109f\n"					\
  "rcr $0x1, %%eax\n jnc 110f\n"					\
  "rcr $0x1, %%eax\n jnc 111f\n"					\
  "rcr $0x1, %%eax\n jnc 112f\n"					\
  "rcr $0x1, %%eax\n jnc 113f\n"					\
  "rcr $0x1, %%eax\n jnc 114f\n"					\
  "rcr $0x1, %%eax\n jnc 115f\n"					\
									\
  "jmp 0b\n"								\
  									\
  "100: movq %%r12, %%r14\n     jmp 200f\n"				\
  "101: lea 0x1(%%r12), %%r14\n jmp 201f\n"				\
  "102: lea 0x2(%%r12), %%r14\n jmp 202f\n"				\
  "103: lea 0x3(%%r12), %%r14\n jmp 203f\n"				\
  "104: lea 0x4(%%r12), %%r14\n jmp 204f\n"				\
  "105: lea 0x5(%%r12), %%r14\n jmp 205f\n"				\
  "106: lea 0x6(%%r12), %%r14\n jmp 206f\n"				\
  "107: lea 0x7(%%r12), %%r14\n jmp 207f\n"				\
  "108: lea 0x8(%%r12), %%r14\n jmp 208f\n"				\
  "109: lea 0x9(%%r12), %%r14\n jmp 209f\n"				\
  "110: lea 0xA(%%r12), %%r14\n jmp 210f\n"				\
  "111: lea 0xB(%%r12), %%r14\n jmp 211f\n"				\
  "112: lea 0xC(%%r12), %%r14\n jmp 212f\n"				\
  "113: lea 0xD(%%r12), %%r14\n jmp 213f\n"				\
  "114: lea 0xE(%%r12), %%r14\n jmp 214f\n"				\
  "115: lea 0xF(%%r12), %%r14\n jmp 215f\n"				\
  									\
  ".balign 8\n"								\
  "3: psubusb (%0, %%r12), %%xmm0\n"					\
  "pcmpeqb %%xmm2, %%xmm0\n"						\
  "pmovmskb %%xmm0, %%eax\n"						\
  "     rcr $0x1, %%eax\n jnc 300f\n"					\
  "200: rcr $0x1, %%eax\n jnc 301f\n"					\
  "201: rcr $0x1, %%eax\n jnc 302f\n"					\
  "202: rcr $0x1, %%eax\n jnc 303f\n"					\
  "203: rcr $0x1, %%eax\n jnc 304f\n"					\
  "204: rcr $0x1, %%eax\n jnc 305f\n"					\
  "205: rcr $0x1, %%eax\n jnc 306f\n"					\
  "206: rcr $0x1, %%eax\n jnc 307f\n"					\
  "207: rcr $0x1, %%eax\n jnc 308f\n"					\
  "208: rcr $0x1, %%eax\n jnc 309f\n"					\
  "209: rcr $0x1, %%eax\n jnc 310f\n"					\
  "210: rcr $0x1, %%eax\n jnc 311f\n"					\
  "211: rcr $0x1, %%eax\n jnc 312f\n"					\
  "212: rcr $0x1, %%eax\n jnc 313f\n"					\
  "213: rcr $0x1, %%eax\n jnc 314f\n"					\
  "214: rcr $0x1, %%eax\n jnc 315f\n"					\
									\
  "215: add $0x10, %%r12\n"						\
  "cmp %2, %%r12\n"							\
  "jnl 5f\n"                   /* Only 1 alg. init to compute & end */	\
  "movdqa %3, %%xmm0\n"							\
  "jmp 3b\n"								\
  									\
  "300: movq %%r12, %%r15\n     jmp 4f\n"				\
  "301: lea 0x1(%%r12), %%r15\n jmp 4f\n"				\
  "302: lea 0x2(%%r12), %%r15\n jmp 4f\n"				\
  "303: lea 0x3(%%r12), %%r15\n jmp 4f\n"				\
  "304: lea 0x4(%%r12), %%r15\n jmp 4f\n"				\
  "305: lea 0x5(%%r12), %%r15\n jmp 4f\n"				\
  "306: lea 0x6(%%r12), %%r15\n jmp 4f\n"				\
  "307: lea 0x7(%%r12), %%r15\n jmp 4f\n"				\
  "308: lea 0x8(%%r12), %%r15\n jmp 4f\n"				\
  "309: lea 0x9(%%r12), %%r15\n jmp 4f\n"				\
  "310: lea 0xA(%%r12), %%r15\n jmp 4f\n"				\
  "311: lea 0xB(%%r12), %%r15\n jmp 4f\n"				\
  "312: lea 0xC(%%r12), %%r15\n jmp 4f\n"				\
  "313: lea 0xD(%%r12), %%r15\n jmp 4f\n"				\
  "314: lea 0xE(%%r12), %%r15\n jmp 4f\n"				\
  "315: lea 0xF(%%r12), %%r15\n"					\
									\
  ".balign 8\n"								\
  "4: cvtsi2sd %%r15d, %%xmm0\n"					\
  "unpcklpd %%xmm0, %%xmm0\n"						\
  "cvtsi2sd %%r14d, %%xmm0\n"

#define ALG_ILP0						\
  "shufps $0xED, %%xmm1, %%xmm1\n"				\
  "pslld $0x01, %%xmm1\n"					\
  "psrld $0x01, %%xmm1\n"					\
  "cvtdq2pd %%xmm1, %%xmm1\n"					\
  "subpd %4, %%xmm1\n"						\
  "mulpd %5, %%xmm1\n"						\
  "cvtsd2si %%xmm1, %%r13d\n"					\
  "mov %%r13b, (%1, %%r14)\n"

#define ALG_ILP2						\
  ALG_ILP0							\
  "unpckhpd %%xmm1, %%xmm1\n"					\
  "cvtsd2si %%xmm1, %%r13d\n"					\
  "mov %%r13b, (%1, %%r15)\n"					\
  								\
  "sub %%r12, %%r15\n"						\
  "lea 2b(,%%r15,4), %%r15\n"					\
  "jmp * %%r15\n"						\
  								\
  ".balign 8\n"							\
  "5: cvtsi2sd %%r14d, %%xmm0\n"
  
#define ALG_ILP3							\
  ALG_ILP0								\
  									\
  "6:"									\
  :: "r"(xS), "r"(S), "r"(Idiv2), "x"(tt), "x"(add),			\
    "x"(scale), "x"(u0), "x"(u1), "x"(u2), "x"(u3),			\
    "x"(u4), "x"(u5), "x"(u6), "x"(u7), "x"(u8),			\
    "x"(u9)								\
    : "%rax", "%r12", "%r13", "%r14", "%r15", "%xmm0", "%xmm1", "%xmm2")

#define ALG_ILP3D							\
  ALG_ILP0								\
  									\
  "6:"									\
  :: "r"(xS), "r"(S), "r"(Idiv2), "x"(tt), "x"(add),			\
    "x"(scale), "r"(u), "r"((uint64_t) (d<<4))				\
    : "%rax", "%r12", "%r13", "%r14", "%r15", "%xmm0", "%xmm1", "%xmm2")
  /****************** End of local define ******************************/

  memset (S, 255, (Idiv2<<1) * ej);
  j *= ej;
  ej += j;
  xS += Idiv2;
  S += Idiv2;
  /* The next pragma is needed because u1...u9 are not "really" used 
     when the keyword MAYBE_UNUSED is used: there are only used as false
     inputs of the main X86 asm. But gcc warns there are not initialized!
     So I desactive the warning. */
#pragma GCC diagnostic ignored "-Wuninitialized"
  switch(d) {
  case 0: { /* Not optimal, but so useless... */
    __m128d u0=_mm_set1_pd(alg->fijd[0]), u1 MAYBE_UNUSED, u2 MAYBE_UNUSED,
      u3 MAYBE_UNUSED, u4 MAYBE_UNUSED, u5 MAYBE_UNUSED, u6 MAYBE_UNUSED,
      u7 MAYBE_UNUSED, u8 MAYBE_UNUSED, u9 MAYBE_UNUSED;
    for (; j < ej; S += (Idiv2<<1), xS += (Idiv2<<1), j++) {
      ALG_ILP1 ALG_IL0A ALG_ILP2 ALG_IL0A ALG_ILP3;
    }}
    return;
  case 1: {
    __m128d u0, u1,
      u2 MAYBE_UNUSED, u3 MAYBE_UNUSED, u4 MAYBE_UNUSED, u5 MAYBE_UNUSED,
      u6 MAYBE_UNUSED, u7 MAYBE_UNUSED, u8 MAYBE_UNUSED, u9 MAYBE_UNUSED;
    for (; j < ej; S += (Idiv2<<1), xS += (Idiv2<<1), j++) {
      INITALG1(j);
      ALG_ILP1 ALG_IL1A ALG_IL1B ALG_ILP2 ALG_IL1A ALG_IL1C ALG_ILP3;
    }}
    return;
  case 2: {
    __m128d h, u0, u1, u2,
      u3 MAYBE_UNUSED, u4 MAYBE_UNUSED, u5 MAYBE_UNUSED,
      u6 MAYBE_UNUSED, u7 MAYBE_UNUSED, u8 MAYBE_UNUSED, u9 MAYBE_UNUSED;
    for (; j < ej; S += (Idiv2<<1), xS += (Idiv2<<1), j++) {
      INITALG2(j);
      ALG_ILP1 ALG_IL2A ALG_IL2B ALG_ILP2 ALG_IL2A ALG_IL2C ALG_ILP3;
    }}
    return;
  case 3: {
    __m128d g, h, u0, u1, u2, u3,
      u4 MAYBE_UNUSED, u5 MAYBE_UNUSED, u6 MAYBE_UNUSED,
      u7 MAYBE_UNUSED, u8 MAYBE_UNUSED, u9 MAYBE_UNUSED;
    for (; j < ej; S += (Idiv2<<1), xS += (Idiv2<<1), j++) {
      INITALG3(j);
      ALG_ILP1 ALG_IL3A ALG_IL3B ALG_ILP2 ALG_IL3A ALG_IL3C ALG_ILP3;
    }}
    return;
  case 4: {
    __m128d g, h, u0, u1, u2, u3, u4,
      u5 MAYBE_UNUSED, u6 MAYBE_UNUSED,
      u7 MAYBE_UNUSED, u8 MAYBE_UNUSED, u9 MAYBE_UNUSED;
    for (; j < ej; S += (Idiv2<<1), xS += (Idiv2<<1), j++) {
      INITALG4(j);
      ALG_ILP1 ALG_IL4A ALG_IL4B ALG_ILP2 ALG_IL4A ALG_IL4C ALG_ILP3;
    }}
    return;
  case 5: {
    __m128d g, h, u0, u1, u2, u3, u4, u5,
      u6 MAYBE_UNUSED, u7 MAYBE_UNUSED, u8 MAYBE_UNUSED, u9 MAYBE_UNUSED;
    for (; j < ej; S += (Idiv2<<1), xS += (Idiv2<<1), j++) {
      INITALG5(j);
      ALG_ILP1 ALG_IL5A ALG_IL5B ALG_ILP2 ALG_IL5A ALG_IL5C ALG_ILP3;
    }}
    return;
  case 6: {
    __m128d g, h, u0, u1, u2, u3, u4, u5, u6,
      u7 MAYBE_UNUSED, u8 MAYBE_UNUSED, u9 MAYBE_UNUSED;
    for (; j < ej; S += (Idiv2<<1), xS += (Idiv2<<1), j++) {
      INITALG6(j);
      ALG_ILP1 ALG_IL6A ALG_IL6B ALG_ILP2 ALG_IL6A ALG_IL6C ALG_ILP3;
    }}
    return;
  case 7: {
    __m128d g, h, u0, u1, u2, u3, u4, u5, u6, u7,
      u8 MAYBE_UNUSED, u9 MAYBE_UNUSED;
    for (; j < ej; S += (Idiv2<<1), xS += (Idiv2<<1), j++) {
      INITALG7(j);
      ALG_ILP1 ALG_IL7A ALG_IL7B ALG_ILP2 ALG_IL7A ALG_IL7C ALG_ILP3;
    }}
    return;
  case 8: {
    __m128d g, h, u0, u1, u2, u3, u4, u5, u6, u7, u8,
      u9 MAYBE_UNUSED;
    for (; j < ej; S += (Idiv2<<1), xS += (Idiv2<<1), j++) {
      INITALG8(j);
      ALG_ILP1 ALG_IL8A ALG_IL8B ALG_ILP2 ALG_IL8A ALG_IL8C ALG_ILP3;
    }}
    return;
#pragma GCC diagnostic warning "-Wuninitialized"			
  case 9: {
    __m128d g, h, u0, u1, u2, u3, u4, u5, u6, u7, u8, u9;
    for (; j < ej; S += (Idiv2<<1), xS += (Idiv2<<1), j++) {
      INITALG9(j);
      ALG_ILP1 ALG_IL9A ALG_IL9B ALG_ILP2 ALG_IL9A ALG_IL9C ALG_ILP3;
    }}
    return;
  default: {
    __m128d u[d+1];
    for (; j < ej; S += (Idiv2<<1), xS += (Idiv2<<1), j++) {
      INITALGD(j);
      ALG_ILP1 ALG_ILD ALG_ILP2 ALG_ILD ALG_ILP3D;
    }}
    return;
  }
}
#undef ALG_IL0A
#undef ALG_IL1A
#undef ALG_IL2A
#undef ALG_IL3A
#undef ALG_IL4A
#undef ALG_IL5A
#undef ALG_IL6A
#undef ALG_IL7A
#undef ALG_IL8A
#undef ALG_IL9A
#undef ALG_IL1B
#undef ALG_IL2B
#undef ALG_IL3B
#undef ALG_IL4B
#undef ALG_IL5B
#undef ALG_IL6B
#undef ALG_IL7B
#undef ALG_IL8B
#undef ALG_IL9B
#undef ALG_IL1C
#undef ALG_IL2C
#undef ALG_IL3C
#undef ALG_IL4C
#undef ALG_IL5C
#undef ALG_IL6C
#undef ALG_IL7C
#undef ALG_IL8C
#undef ALG_IL9C
#undef ALG_ILD
#undef ALG_ILP1
#undef ALG_ILP0
#undef ALG_ILP2
#undef ALG_ILP3
#undef ALG_ILP3D
#endif
#endif
#endif

/**************** 4: HAVE_SSE2, !ALG_LAZY, !ALG_RAT *************************/
#ifdef HAVE_SSE2
#ifndef ALG_LAZY
#ifndef ALG_RAT

/* Full initialization of the algebraics with SSE. 32bits compatible.
   No so slow! About x6 slower than the best. */
void init_alg_norms_bucket_region (unsigned char *S,
				   unsigned char *xS MAYBE_UNUSED, unsigned int j,
				   sieve_info_ptr si)
{
  sieve_side_info_ptr alg = si->sides[ALGEBRAIC_SIDE];
  ssize_t ih, Idiv2 = (si->I >> 1);
  unsigned int ej = 1U << (LOG_BUCKET_REGION - si->conf->logI),
    d = si->cpoly->alg->degree;
  __m128d two = _mm_set1_pd(2.0);
  ALG_INIT_SCALE_ADD_A;
  
  j *= ej;
  ej += j;
  S += Idiv2;
  switch(d) {
  case 0:
    memset(S-Idiv2,inttruncfastlog2(fabs(alg->fijd[0]),*(double *)&add,*(double *)&scale),(Idiv2<<1)*(ej-j));
    return;
  case 1: {
    __m128d g, h, u0, u1;
    for (; j < ej; S += (Idiv2<<1), j++) {
      INITALG1(j);
      g = h = _mm_set_pd(1.0 - Idiv2, (double) -Idiv2);
      for (ih = -Idiv2; ih < Idiv2;) {
	G1;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G1;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G1;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G1;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
      }}}
    return;
  case 2: {
    __m128d g, h, u0, u1, u2;
    for (; j < ej; S += (Idiv2<<1), j++) {
      INITALG2(j);
      g = h = _mm_set_pd(1.0 - Idiv2, (double) -Idiv2);
      for (ih = -Idiv2; ih < Idiv2;) {
	G2;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G2;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G2;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G2;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
      }}}
    return;
  case 3: {
    __m128d g, h, u0, u1, u2, u3;
    for (; j < ej; S += (Idiv2<<1), j++) {
      INITALG3(j);
      g = h = _mm_set_pd(1.0 - Idiv2, (double) -Idiv2);
      for (ih = -Idiv2; ih < Idiv2;) {
	G3;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G3;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G3;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G3;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
      }}}
    return;
  case 4: {
    __m128d g, h, u0, u1, u2, u3, u4;
    for (; j < ej; S += (Idiv2<<1), j++) {
      INITALG4(j);
      g = h = _mm_set_pd(1.0 - Idiv2, (double) -Idiv2);
      for (ih = -Idiv2; ih < Idiv2;) {
	G4;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G4;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G4;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G4;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
      }}}
    return;
  case 5: {
    __m128d g, h, u0, u1, u2, u3, u4, u5;
    for (; j < ej; S += (Idiv2<<1), j++) {
      INITALG5(j);
      g = h = _mm_set_pd(1.0 - Idiv2, (double) -Idiv2);
      for (ih = -Idiv2; ih < Idiv2;) {
	G5;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G5;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G5;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G5;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
      }}}
    return;
  case 6: {
    __m128d g, h, u0, u1, u2, u3, u4, u5, u6;
    for (; j < ej; S += (Idiv2<<1), j++) {
      INITALG6(j);
      g = h = _mm_set_pd(1.0 - Idiv2, (double) -Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 2) {
	G6;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G6;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G6;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G6;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
      }}}
    return;
  case 7: {
    __m128d g, h, u0, u1, u2, u3, u4, u5, u6, u7;
    for (; j < ej; S += (Idiv2<<1), j++) {
      INITALG7(j);
      g = h = _mm_set_pd(1.0 - Idiv2, (double) -Idiv2);
      for (ih = -Idiv2; ih < Idiv2;) {
	G7;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G7;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G7;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G7;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
      }}}
    return;
  case 8: {
    __m128d g, h, u0, u1, u2, u3, u4, u5, u6, u7, u8;
    for (; j < ej; S += (Idiv2<<1), j++) {
      INITALG8(j);
      g = h = _mm_set_pd(1.0 - Idiv2, (double) -Idiv2);
      for (ih = -Idiv2; ih < Idiv2;) {
	G8;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G8;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G8;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G8;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
      }}}
    return;
  case 9: {
    __m128d g, h, u0, u1, u2, u3, u4, u5, u6, u7, u8, u9;
    for (; j < ej; S += (Idiv2<<1), j++) {
      INITALG9(j);
      g = h = _mm_set_pd(1.0 - Idiv2, (double) -Idiv2);
      for (ih = -Idiv2; ih < Idiv2;) {
	G9;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G9;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G9;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	G9;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
      }}}
    return;
  default: {
    __m128d g, h, u[d+1];
    for (; j < ej; S += (Idiv2<<1), j++) {
      INITALGD(j);
      h = _mm_set_pd(1.0 - Idiv2, (double) -Idiv2);
      for (ih = -Idiv2; ih < Idiv2;) {
	GD;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	GD;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	GD;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
	GD;
	w16itruncfastlog2fabs (g, add, scale, S, ih); g = h =_mm_add_pd (h, two); ih += 2;
      }}}
    return;
  }
}
#endif
#endif
#endif

/************************* undef of all SSE define *****************************/
#ifdef HAVE_SSE2
#undef TSTZXMM
#undef _MM_SHUFFLE_EPI32
#undef INITALG0
#undef INITALG1
#undef INITALG2
#undef INITALG3
#undef INITALG4
#undef INITALG5
#undef INITALG6
#undef INITALG7
#undef INITALG8
#undef INITALG9
#undef INITALGD
#undef G1
#undef G2
#undef G3
#undef G4
#undef G5
#undef G6
#undef G7
#undef G8
#undef G9
#undef GD
#undef ALG_INIT_SCALE_ADD_A
#undef ALG_INIT_SCALE_ADD_A_TT
#endif
/******************** end of undef of all SSE #define ************************/

/***************** Beginning of all non SSE #define **************************/
/* I cannot do a #ifndef HAVE_SSE here because the 7th algorithm is faster   */
/* than its SSE version, 3th algorithm. So 7th is the default => all these   */
/* defines may be used with a SSE machine.                                   */
/*****************************************************************************/
/* #ifdef HAVE_SSE2 */ 

#define INITALG0 uu0=0x0101010101010101*inttruncfastlog2(fabs(alg->fijd[0]),add,scale)

#define INITALG1(A)				\
  u1 = alg->fijd[1]    ;			\
  u0 = alg->fijd[0] * ((double) (A))

#define INITALG2(A)				\
  u2 = alg->fijd[2]    ; h  = ((double) (A));	\
  u1 = alg->fijd[1] * h; h *= h;		\
  u0 = alg->fijd[0] * h

#define INITALG3(A)				 \
  u3 = alg->fijd[3]    ; h  = ((double) (A));	 \
  u2 = alg->fijd[2] * h; g  = h * h;		 \
  u1 = alg->fijd[1] * g; g *= h;		 \
  u0 = alg->fijd[0] * g

#define INITALG4(A)				 \
  u4 = alg->fijd[4]    ; h  = ((double) (A));	 \
  u3 = alg->fijd[3] * h; g  = h * h;		 \
  u2 = alg->fijd[2] * g; g *= h;		 \
  u1 = alg->fijd[1] * g; g *= h;		 \
  u0 = alg->fijd[0] * g

#define INITALG5(A)				 \
  u5 = alg->fijd[5]    ; h  = ((double) (A));	 \
  u4 = alg->fijd[4] * h; g  = h * h;		 \
  u3 = alg->fijd[3] * g; g *= h;		 \
  u2 = alg->fijd[2] * g; g *= h;		 \
  u1 = alg->fijd[1] * g; g *= h;		 \
  u0 = alg->fijd[0] * g

#define INITALG6(A)					\
  u6 = alg->fijd[6]    ; h  = ((double) (A));		\
  u5 = alg->fijd[5] * h; g  = h * h;			\
  u4 = alg->fijd[4] * g; g *= h;			\
  u3 = alg->fijd[3] * g; g *= h;			\
  u2 = alg->fijd[2] * g; g *= h;			\
  u1 = alg->fijd[1] * g; g *= h;			\
  u0 = alg->fijd[0] * g

#define INITALG7(A)					\
  u7 = alg->fijd[7]    ; h  = ((double) (A));		\
  u6 = alg->fijd[6] * h; g  = h * h;			\
  u5 = alg->fijd[5] * g; g *= h;			\
  u4 = alg->fijd[4] * g; g *= h;			\
  u3 = alg->fijd[3] * g; g *= h;			\
  u2 = alg->fijd[2] * g; g *= h;			\
  u1 = alg->fijd[1] * g; g *= h;			\
  u0 = alg->fijd[0] * g

#define INITALG8(A)					\
  u8 = alg->fijd[8]    ; h  = ((double) (A));		\
  u7 = alg->fijd[7] * h; g  = h * h;			\
  u6 = alg->fijd[6] * g; g *= h;			\
  u5 = alg->fijd[5] * g; g *= h;			\
  u4 = alg->fijd[4] * g; g *= h;			\
  u3 = alg->fijd[3] * g; g *= h;			\
  u2 = alg->fijd[2] * g; g *= h;			\
  u1 = alg->fijd[1] * g; g *= h;			\
  u0 = alg->fijd[0] * g

#define INITALG9(A)					\
  u9 = alg->fijd[9]    ; h  = ((double) (A));		\
  u8 = alg->fijd[8] * h; g  = h * h;			\
  u7 = alg->fijd[7] * g; g *= h;			\
  u6 = alg->fijd[6] * g; g *= h;			\
  u5 = alg->fijd[5] * g; g *= h;			\
  u4 = alg->fijd[4] * g; g *= h;			\
  u3 = alg->fijd[3] * g; g *= h;			\
  u2 = alg->fijd[2] * g; g *= h;			\
  u1 = alg->fijd[1] * g; g *= h;			\
  u0 = alg->fijd[0] * g

#define INITALGD(A)				\
  fpoly_scale(u, alg->fijd, d, (double) (A))

#define ABSG1 g=fabs(u0+h*u1)

#define ABSG2 g=fabs(u0+h*(u1+h*u2))

#define ABSG3 g=fabs(u0+h*(u1+h*(u2+h*u3)))

#define ABSG4 g=fabs(u0+h*(u1+h*(u2+h*(u3+h*u4))))

#define ABSG5 g=fabs(u0+h*(u1+h*(u2+h*(u3+h*(u4+h*u5)))))

#define ABSG6 g=fabs(u0+h*(u1+h*(u2+h*(u3+h*(u4+h*(u5+h*u6))))))

#define ABSG7 g=fabs(u0+h*(u1+h*(u2+h*(u3+h*(u4+h*(u5+h*(u6+h*u7)))))))

#define ABSG8 g=fabs(u0+h*(u1+h*(u2+h*(u3+h*(u4+h*(u5+h*(u6+h*(u7+h*u8))))))))

#define ABSG9 g=fabs(u0+h*(u1+h*(u2+h*(u3+h*(u4+h*(u5+h*(u6+h*(u7+h*(u8+h*u9)))))))))

#define ABSGD g=u[d]; for(unsigned int k=d;--k!=UINT_MAX;g=g*h+u[k]); g=fabs(g)

#define ALG_INIT_SCALE_ADD			\
  double scale = alg->scale * (1.0/0x100000),	\
    add = 0x3FF00000 - GUARD / scale

/* Careful! non SSE t and SSE tt are not the same formula, but the
   same use. */
#define ALG_INIT_SCALE_ADD_T						\
  ALG_INIT_SCALE_ADD;							\
  uint64_t t = (si->sides[RATIONAL_SIDE]->bound < 255 ?			\
		((si->sides[RATIONAL_SIDE]->bound+1)>>1):127)*0x0101010101010101

/* #endif !HAVE_SSE2 */
/**************** End of used non SSE #define *****************************/

/**************** 5: !HAVE_SSE2, ALG_LAZY, ALG_RAT ************************/
#ifndef HAVE_SSE2
#ifdef ALG_LAZY
#ifdef ALG_RAT
/* Smart initialization of the algebraics. In one test, 8 bytes of
   the rationals are tested. If one is > rat->bound, computes the
   central initialization value and duplicates it 8 times, and writes it.
   At the beginning, all are initialised to 255.
   2nd fastest initialization, but maybe the smartest.
   Non SSE version, 32bits compatible.
   
   NOTA : the SSE and non SSE versions differs!
   
   NB: the void *tg is mandatory, because x86 cannot transfer without
   memory a xmm register in a general register. And for this, gcc
   needs dedicated memory: it doesn't use the stack!

   NB2: 2 possibilities:
   1. The code sets h before the for (ih = -Idiv2...) and increase h
   with a in the loop.
   2. The code sets h only when EIGHT_BYTES_TEST is true -> a int to
   double conversion is needed.
   On my machine, same time for the both codes. So, I prefer don't use
   a, because the limited floating register number.
*/
void init_alg_norms_bucket_region(unsigned char *S,
				  unsigned char *xS, unsigned int j,
				  sieve_info_ptr si)
{
  sieve_side_info_ptr alg = si->sides[ALGEBRAIC_SIDE];
  ssize_t ih, Idiv2 = (si->I >> 1);
  unsigned int ej = 1U << (LOG_BUCKET_REGION - si->conf->logI),
    d = si->cpoly->alg->degree;
  ALG_INIT_SCALE_ADD_T;

  /* True if one byte of the uint64_t *(xS + ih) are ~ >
     si->sides[RATIONAL_SIDE]->bound */
#define EIGHT_BYTES_TEST				\
  ((((*(uint64_t *)&xS[ih]) >> 1)			\
    & 0x7F7F7F7F7F7F7F7F) - t) & 0x8080808080808080

  memset (S, 255, (Idiv2<<1) * ej);
  j *= ej;
  ej += j;
  S += Idiv2;
  xS += Idiv2;
  switch (d) {
  case 0: {
    uint64_t uu0;
    INITALG0;
    for (; j < ej; j++, S += (Idiv2<<1), xS += (Idiv2<<1)) {
      for (ih = -Idiv2; ih < Idiv2; ih += 8)
	if (UNLIKELY(EIGHT_BYTES_TEST)) {
	  *(uint64_t *)&S[ih] = uu0;
	}}}
    return;
  case 1 : {
    double g, h, u0, u1;
    for (; j < ej; j++, S += (Idiv2<<1), xS += (Idiv2<<1)) {
      INITALG1(j);
      for (ih = -Idiv2; ih < Idiv2; ih += 8)
	if (UNLIKELY(EIGHT_BYTES_TEST)) {
	  h = (double) (ih + 3);
	  ABSG1;
	  uint64truncfastlog2(g,add,scale,S,ih);
	}}}
    return;
  case 2 : {
    double g, h, u0, u1, u2;
    for (; j < ej; j++, S += (Idiv2<<1), xS += (Idiv2<<1)) {
      INITALG2(j);
      for (ih = -Idiv2; ih < Idiv2; ih += 8)
	if (UNLIKELY(EIGHT_BYTES_TEST)) {
	  h = (double) (ih + 3);
	  ABSG2;
	  uint64truncfastlog2(g,add,scale,S,ih);
	}}}
    return;
  case 3 : {
    double g, h, u0, u1, u2, u3;
    for (; j < ej; j++, S += (Idiv2<<1), xS += (Idiv2<<1)) {
      INITALG3(j);
      for (ih = -Idiv2; ih < Idiv2; ih += 8)
	if (UNLIKELY(EIGHT_BYTES_TEST)) {
	  h = (double) (ih + 3);
	  ABSG3;
	  uint64truncfastlog2(g,add,scale,S,ih);
	}}}
    return;
  case 4 : {
    double g, h, u0, u1, u2, u3, u4;
    for (; j < ej; j++, S += (Idiv2<<1), xS += (Idiv2<<1)) {
      INITALG4(j);
      for (ih = -Idiv2; ih < Idiv2; ih += 8)
	if (UNLIKELY(EIGHT_BYTES_TEST)) {
	  h = (double) (ih + 3);
	  ABSG4;
	  uint64truncfastlog2(g,add,scale,S,ih);
	}}}
    return;
  case 5 : {
    double g, h, u0, u1, u2, u3, u4, u5;
    for (; j < ej; j++, S += (Idiv2<<1), xS += (Idiv2<<1)) {
      INITALG5(j);
      for (ih = -Idiv2; ih < Idiv2; ih += 8)
	if (UNLIKELY(EIGHT_BYTES_TEST)) {
	  h = (double) (ih + 3);
	  ABSG5;
	  uint64truncfastlog2(g,add,scale,S,ih);
	}}}
    return;
  case 6 : {
    double g, h, u0, u1, u2, u3, u4, u5, u6;
    for (; j < ej; j++, S += (Idiv2<<1), xS += (Idiv2<<1)) {
      INITALG6(j);
      for (ih = -Idiv2; ih < Idiv2; ih += 8)
	if (UNLIKELY(EIGHT_BYTES_TEST)) {
	  h = (double) (ih + 3);
	  ABSG6;
	  uint64truncfastlog2(g,add,scale,S,ih);
	}}}
    return;
  case 7 : {
    double g, h, u0, u1, u2, u3, u4, u5, u6, u7;
    for (; j < ej; j++, S += (Idiv2<<1), xS += (Idiv2<<1)) {
      INITALG7(j);
      for (ih = -Idiv2; ih < Idiv2; ih += 8)
	if (UNLIKELY(EIGHT_BYTES_TEST)) {
	  h = (double) (ih + 3);
	  ABSG7;
	  uint64truncfastlog2(g,add,scale,S,ih);
	}}}
    return;
  case 8 : {
    double g, h, u0, u1, u2, u3, u4, u5, u6, u7, u8;
    for (; j < ej; j++, S += (Idiv2<<1), xS += (Idiv2<<1)) {
      INITALG8(j);
      for (ih = -Idiv2; ih < Idiv2; ih += 8)
	if (UNLIKELY(EIGHT_BYTES_TEST)) {
	  h = (double) (ih + 3);
	  ABSG8;
	  uint64truncfastlog2(g,add,scale,S,ih);
	}}}
    return;
  case 9 : {
    double g, h, u0, u1, u2, u3, u4, u5, u6, u7, u8, u9;
    for (; j < ej; j++, S += (Idiv2<<1), xS += (Idiv2<<1)) {
      INITALG9(j);
      for (ih = -Idiv2; ih < Idiv2; ih += 8)
	if (UNLIKELY(EIGHT_BYTES_TEST)) {
	  h = (double) (ih + 3);
	  ABSG9;
	  uint64truncfastlog2(g,add,scale,S,ih);
	}}}
    return;
  default : {
    double g, h, u[d+1];
    for (; j < ej; j++, S += (Idiv2<<1), xS += (Idiv2<<1)) {
      INITALGD(j);
      for (ih = -Idiv2; ih < Idiv2; ih += 8)
	if (UNLIKELY(EIGHT_BYTES_TEST)) {
	  h = (double) (ih + 3);
	  ABSGD;
	  uint64truncfastlog2(g,add,scale,S,ih);
	}}}
    return;
  }
}
#endif
#endif
#endif

/**************** 6: !HAVE_SSE2, ALG_LAZY, !ALG_RAT ***********************/
#ifndef HAVE_SSE2
#ifdef ALG_LAZY
#ifndef ALG_RAT
/* Smart initialization of the algebraics. Computes the central initialization
   value of a box of (horizontail=8,vertical=p) and propagates it in the box.
   8 is locked by the code and p is the minimum of VERT_NORM_STRIDE (#define, = 4)
   and ej, with ej = 2^(LOG_BUCKET_REGION-si->conf->logI).
   So, the fastest code is for logI <= 14.
   It's the fastest initialization of the algebraics. Non SSE version. */
void init_alg_norms_bucket_region(unsigned char *S,
				  unsigned char *xS MAYBE_UNUSED, unsigned int j,
				  sieve_info_ptr si)
{
  sieve_side_info_ptr alg = si->sides[ALGEBRAIC_SIDE];
  ssize_t ih, Idiv2 = (si->I >> 1);
  unsigned int ej = 1U << (LOG_BUCKET_REGION - si->conf->logI),
    d = si->cpoly->alg->degree,
    p = ej < VERT_NORM_STRIDE ? ej : VERT_NORM_STRIDE;
  ALG_INIT_SCALE_ADD;

#define FILL_OTHER_LINES						\
  if (LIKELY(p > 1)) {							\
    unsigned char *mS = S + Idiv2;					\
    S -= Idiv2;								\
    for (unsigned int k = 1; k < p; k++, mS += (Idiv2<<1))		\
      memcpy (mS, S, (Idiv2<<1));					\
    S = mS + Idiv2;							\
    j += p;								\
  }									\
  else do {								\
      S += (Idiv2<<1);							\
      j++;								\
    } while(0)
  
  j *= ej;
  ej += j;
  S += Idiv2;
  switch (d) {
  case 0: {
    unsigned int uu0;
    INITALG0;
    memset(S - Idiv2, uu0, (Idiv2 << 1) * (ej - j));
  }
    return;
  case 1 : {
    double g, h, u0, u1;
    while (j < ej) {
      INITALG1(j + (p >> 1));
      h = (double) (3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 8, h += 8.) {
	ABSG1;
	uint64truncfastlog2(g,add,scale,S,ih);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 2 : {
    double g, h, u0, u1, u2;
    while (j < ej) {
      INITALG2(j + (p >> 1));
      h = (double) (3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 8, h += 8.) {
	ABSG2;
	uint64truncfastlog2(g,add,scale,S,ih);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 3 : {
    double g, h, u0, u1, u2, u3;
    while (j < ej) {
      INITALG3(j + (p >> 1));
      h = (double) (3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 8, h += 8.) {
	ABSG3;
	uint64truncfastlog2(g,add,scale,S,ih);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 4 : {
    double g, h, u0, u1, u2, u3, u4;
    while (j < ej) {
      INITALG4(j + (p >> 1));
       h = (double) (3 - Idiv2);
       for (ih = -Idiv2; ih < Idiv2; ih += 8, h += 8.) {
	ABSG4;
	uint64truncfastlog2(g,add,scale,S,ih);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 5 : {
    double g, h, u0, u1, u2, u3, u4, u5;
    while (j < ej) {
      INITALG5(j + (p >> 1));
      h = (double) (3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 8, h += 8.) {
	ABSG5;
	uint64truncfastlog2(g,add,scale,S,ih);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 6 : {
    double g, h, u0, u1, u2, u3, u4, u5, u6;
    while (j < ej) {
      INITALG6(j + (p >> 1));
      h = (double) (3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 8, h += 8.) {
	ABSG6;
	uint64truncfastlog2(g,add,scale,S,ih);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 7 : {
    double g, h, u0, u1, u2, u3, u4, u5, u6, u7;
    while (j < ej) {
      INITALG7(j + (p >> 1));
      h = (double) (3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 8, h += 8.) {
	ABSG7;
	uint64truncfastlog2(g,add,scale,S,ih);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 8 : {
    double g, h, u0, u1, u2, u3, u4, u5, u6, u7, u8;
    while (j < ej) {
      INITALG8(j + (p >> 1));
      h = (double) (3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 8, h += 8.) {
	ABSG8;
	uint64truncfastlog2(g,add,scale,S,ih);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 9 : {
    double g, h, u0, u1, u2, u3, u4, u5, u6, u7, u8, u9;
    while (j < ej) {
      INITALG9(j + (p >> 1));
      h = (double) (3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 8, h += 8.) {
	ABSG9;
	uint64truncfastlog2(g,add,scale,S,ih);
      }
      FILL_OTHER_LINES;
    }}
    return;
  default : {
    double g, h, u[d];
    while (j < ej) {
      INITALGD(j + (p >> 1));
      h = (double) (3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 8, h += 8.) {
	ABSGD;
	uint64truncfastlog2(g,add,scale,S,ih);
      }
      FILL_OTHER_LINES;
    }}
    return;
  }
}
#undef FILL_OTHER_LINES
#endif
#endif
#endif

/********** 7: !(HAVE_SSE2 & x86-64), !ALG_LAZY, ALG_RAT **************/
#ifndef NOTDEFINED /* __x86_64 */ /* Active by default: SSE code slower */
#ifndef ALG_LAZY
#ifdef ALG_RAT
/* Smart initialization of the algebraics : only if the corresponding
   rationnal is < rat->bound, other 255 (by first memset), so ~ x 22 less
   computation; but requires of course many tests. Without SSE.
   NB: this procedure is used if !ALG_LAZY, ALG_RAT & !x86-64. */
void init_alg_norms_bucket_region(unsigned char *S,
				  unsigned char *xS, unsigned int j,
				  sieve_info_ptr si)
{
  sieve_side_info_ptr alg = si->sides[ALGEBRAIC_SIDE];
  ssize_t ih, Idiv2 = (si->I >> 1);
  unsigned int ej = 1U << (LOG_BUCKET_REGION - si->conf->logI),
    d = si->cpoly->alg->degree;
  uint64_t m;
  uint8_t ratbound = si->sides[RATIONAL_SIDE]->bound;
  ALG_INIT_SCALE_ADD;

#define ALG0(A) if (UNLIKELY((uint8_t) m <= ratbound))	\
    S[ih+A]=(uint8_t)uu0
#define ALG1(A) if (UNLIKELY((uint8_t) m <= ratbound)) do {	\
      h = ((double) (ih + A));					\
      ABSG1;							\
      S[ih+A]=inttruncfastlog2(g,add,scale); } while(0)
#define ALG2(A) if (UNLIKELY((uint8_t) m <= ratbound)) do {	\
      h = ((double) (ih + A));					\
      ABSG2;							\
      S[ih+A]=inttruncfastlog2(g,add,scale); } while(0)
#define ALG3(A) if (UNLIKELY((uint8_t) m <= ratbound)) do {	\
      h = ((double) (ih + A));					\
      ABSG3;							\
      S[ih+A]=inttruncfastlog2(g,add,scale); } while(0)
#define ALG4(A) if (UNLIKELY((uint8_t) m <= ratbound)) do {	\
      h = ((double) (ih + A));					\
      ABSG4;							\
      S[ih+A]=inttruncfastlog2(g,add,scale); } while(0)
#define ALG5(A) if (UNLIKELY((uint8_t) m <= ratbound)) do {	\
      h = ((double) (ih + A));					\
      ABSG5;							\
      S[ih+A]=inttruncfastlog2(g,add,scale); } while(0)
#define ALG6(A) if (UNLIKELY((uint8_t) m <= ratbound)) do {	\
      h = ((double) (ih + A));					\
      ABSG6;							\
      S[ih+A]=inttruncfastlog2(g,add,scale); } while(0)
#define ALG7(A) if (UNLIKELY((uint8_t) m <= ratbound)) do {	\
      h = ((double) (ih + A));					\
      ABSG7;							\
      S[ih+A]=inttruncfastlog2(g,add,scale); } while(0)
#define ALG8(A) if (UNLIKELY((uint8_t) m <= ratbound)) do {	\
      h = ((double) (ih + A));					\
      ABSG8;							\
      S[ih+A]=inttruncfastlog2(g,add,scale); } while(0)
#define ALG9(A) if (UNLIKELY((uint8_t) m <= ratbound)) do {	\
      h = ((double) (ih + A));					\
      ABSG9;							\
      S[ih+A]=inttruncfastlog2(g,add,scale); } while(0)
#define ALGD(A) if (UNLIKELY((uint8_t) m <= ratbound)) do {	\
      h = ((double) (ih + A));					\
      ABSGD;							\
      S[ih+A]=inttruncfastlog2(g,add,scale); } while(0)

  memset (S, 255, (Idiv2<<1) * ej);
  j *= ej;
  ej += j;
  S += Idiv2;
  xS += Idiv2;
  switch (d) {
  case 0 : {
    uint8_t uu0;
    INITALG0;
    for (; j < ej; S += (Idiv2<<1), xS += (Idiv2<<1), j++) {
      for (ih = -Idiv2; ih < Idiv2; ih += 8) {
	m = *(uint64_t *)&xS[ih];
	ALG0(0); m >>= 8; ALG0(1); m >>= 8; ALG0(2); m >>= 8; ALG0(3); m >>= 8; 
	ALG0(4); m >>= 8; ALG0(5); m >>= 8; ALG0(6); m >>= 8; ALG0(7); 
      }}}
    return;
  case 1 : {
    double g, h, u0, u1;
    for (; j < ej; S += (Idiv2<<1), xS += (Idiv2<<1), j++) {
      INITALG1(j);
      for (ih = -Idiv2; ih < Idiv2; ih += 8) {
	m = *(uint64_t *)&xS[ih];
	ALG1(0); m >>= 8; ALG1(1); m >>= 8; ALG1(2); m >>= 8; ALG1(3); m >>= 8; 
	ALG1(4); m >>= 8; ALG1(5); m >>= 8; ALG1(6); m >>= 8; ALG1(7); 
      }}}
    return;
  case 2 : {
    double g, h, u0, u1, u2;
    for (; j < ej; S += (Idiv2<<1), xS += (Idiv2<<1), j++) {
      INITALG2(j);
      for (ih = -Idiv2; ih < Idiv2; ih += 8) {
	m = *(uint64_t *)&xS[ih];
	ALG2(0); m >>= 8; ALG2(1); m >>= 8; ALG2(2); m >>= 8; ALG2(3); m >>= 8; 
	ALG2(4); m >>= 8; ALG2(5); m >>= 8; ALG2(6); m >>= 8; ALG2(7); 
      }}}
    return;
  case 3 : {
    double g, h, u0, u1, u2, u3;
    for (; j < ej; S += (Idiv2<<1), xS += (Idiv2<<1), j++) {
      INITALG3(j);
      for (ih = -Idiv2; ih < Idiv2; ih += 8) {
	m = *(uint64_t *)&xS[ih];
	ALG3(0); m >>= 8; ALG3(1); m >>= 8; ALG3(2); m >>= 8; ALG3(3); m >>= 8; 
	ALG3(4); m >>= 8; ALG3(5); m >>= 8; ALG3(6); m >>= 8; ALG3(7); 
      }}}
    return;
  case 4 : {
    double g, h, u0, u1, u2, u3, u4;
    for (; j < ej; S += (Idiv2<<1), xS += (Idiv2<<1), j++) {
      INITALG4(j);
      for (ih = -Idiv2; ih < Idiv2; ih += 8) {
	m = *(uint64_t *)&xS[ih];
	ALG4(0); m >>= 8; ALG4(1); m >>= 8; ALG4(2); m >>= 8; ALG4(3); m >>= 8; 
	ALG4(4); m >>= 8; ALG4(5); m >>= 8; ALG4(6); m >>= 8; ALG4(7); 
      }}}
    return;
  case 5 : {
    double g, h, u0, u1, u2, u3, u4, u5;
    for (; j < ej; S += (Idiv2<<1), xS += (Idiv2<<1), j++) {
      INITALG5(j);
      for (ih = -Idiv2; ih < Idiv2; ih += 8) {
	m = *(uint64_t *)&xS[ih];
	ALG5(0); m >>= 8; ALG5(1); m >>= 8; ALG5(2); m >>= 8; ALG5(3); m >>= 8; 
	ALG5(4); m >>= 8; ALG5(5); m >>= 8; ALG5(6); m >>= 8; ALG5(7); 
      }}}
    return;
  case 6 : {
    double g, h, u0, u1, u2, u3, u4, u5, u6;
    for (; j < ej; S += (Idiv2<<1), xS += (Idiv2<<1), j++) {
      INITALG6(j);
      for (ih = -Idiv2; ih < Idiv2; ih += 8) {
	m = *(uint64_t *)&xS[ih];
	ALG6(0); m >>= 8; ALG6(1); m >>= 8; ALG6(2); m >>= 8; ALG6(3); m >>= 8; 
	ALG6(4); m >>= 8; ALG6(5); m >>= 8; ALG6(6); m >>= 8; ALG6(7); 
      }}}
    return;
  case 7 : {
    double g, h, u0, u1, u2, u3, u4, u5, u6, u7;
    for (; j < ej; S += (Idiv2<<1), xS += (Idiv2<<1), j++) {
      INITALG7(j);
      for (ih = -Idiv2; ih < Idiv2; ih += 8) {
	m = *(uint64_t *)&xS[ih];
	ALG7(0); m >>= 8; ALG7(1); m >>= 8; ALG7(2); m >>= 8; ALG7(3); m >>= 8; 
	ALG7(4); m >>= 8; ALG7(5); m >>= 8; ALG7(6); m >>= 8; ALG7(7); 
      }}}
    return;
  case 8 : {
    double g, h, u0, u1, u2, u3, u4, u5, u6, u7, u8;
    for (; j < ej; S += (Idiv2<<1), xS += (Idiv2<<1), j++) {
      INITALG8(j);
      for (ih = -Idiv2; ih < Idiv2; ih += 8) {
	m = *(uint64_t *)&xS[ih];
	ALG8(0); m >>= 8; ALG8(1); m >>= 8; ALG8(2); m >>= 8; ALG8(3); m >>= 8; 
	ALG8(4); m >>= 8; ALG8(5); m >>= 8; ALG8(6); m >>= 8; ALG8(7); 
      }}}
    return;
  case 9 : {
    double g, h, u0, u1, u2, u3, u4, u5, u6 ,u7, u8, u9;
    for (; j < ej; S += (Idiv2<<1), xS += (Idiv2<<1), j++) {
      INITALG9(j);
      for (ih = -Idiv2; ih < Idiv2; ih += 8) {
	m = *(uint64_t *)&xS[ih];
	ALG9(0); m >>= 8; ALG9(1); m >>= 8; ALG9(2); m >>= 8; ALG9(3); m >>= 8; 
	ALG9(4); m >>= 8; ALG9(5); m >>= 8; ALG9(6); m >>= 8; ALG9(7); 
      }}}
    return;
  default : {
    double g, h, u[d+1];
    for (; j < ej; S += (Idiv2<<1), xS += (Idiv2<<1), j++) {
      INITALGD(j);
      for (ih = -Idiv2; ih < Idiv2; ih += 8) {
	m = *(uint64_t *)&xS[ih];
	ALGD(0); m >>= 8; ALGD(1); m >>= 8; ALGD(2); m >>= 8; ALGD(3); m >>= 8; 
	ALGD(4); m >>= 8; ALGD(5); m >>= 8; ALGD(6); m >>= 8; ALGD(7); 
      }}}
  }
#undef ALG0
#undef ALG1
#undef ALG2
#undef ALG3
#undef ALG4
#undef ALG5
#undef ALG6
#undef ALG7
#undef ALG8
#undef ALG9
#undef ALGD
}
#endif
#endif
#endif

/**************** 8: !HAVE_SSE2, !ALG_LAZY, !ALG_RAT ************************/
#ifndef HAVE_SSE2
#ifndef ALG_LAZY
#ifndef ALG_RAT
/* Full initialization of the algebraics without SSE. Slowest */
void init_alg_norms_bucket_region(unsigned char *S,
				  unsigned char *xS MAYBE_UNUSED, unsigned int j,
				  sieve_info_ptr si)
{
  sieve_side_info_ptr alg = si->sides[ALGEBRAIC_SIDE];
  ssize_t ih, Idiv2 = (si->I >> 1);
  unsigned int ej = 1U << (LOG_BUCKET_REGION - si->conf->logI),
    d = si->cpoly->alg->degree;
  ALG_INIT_SCALE_ADD;
  
#define ALG1(A)					\
  ABSG1;					\
  (&S[ih])[A]=inttruncfastlog2(g,add,scale)
#define ALG2(A)					\
  ABSG2;					\
  (&S[ih])[A]=inttruncfastlog2(g,add,scale)
#define ALG3(A)					\
  ABSG3;					\
  (&S[ih])[A]=inttruncfastlog2(g,add,scale)
#define ALG4(A)					\
  ABSG4;					\
  (&S[ih])[A]=inttruncfastlog2(g,add,scale)
#define ALG5(A)					\
  ABSG5;					\
  (&S[ih])[A]=inttruncfastlog2(g,add,scale)
#define ALG6(A)					\
  ABSG6;					\
  (&S[ih])[A]=inttruncfastlog2(g,add,scale)
#define ALG7(A)					\
  ABSG7;					\
  (&S[ih])[A]=inttruncfastlog2(g,add,scale)
#define ALG8(A)					\
  ABSG8;					\
  (&S[ih])[A]=inttruncfastlog2(g,add,scale)
#define ALG9(A)					\
  ABSG9;					\
  (&S[ih])[A]=inttruncfastlog2(g,add,scale)
#define ALGD(A)					\
  ABSGD;					\
  (&S[ih])[A]=inttruncfastlog2(g,add,scale)

  j *= ej;
  ej += j;
  S += Idiv2;
  switch (d) {
  case 0 : {
    unsigned int uu0;
    INITALG0;
    memset(S - Idiv2, uu0, (Idiv2 << 1) * (ej - j));
  }
    return;
  case 1 : {
    double g, h, u0, u1;
    for (; j < ej; S += (Idiv2<<1), j++) {
      INITALG1(j);
      h = (double) -Idiv2;
      for (ih = -Idiv2; ih < Idiv2; ih += 4) {
	ALG1(0); h++; ALG1(1); h++; ALG1(2); h++; ALG1(3); h++;
      }}}
    return;
  case 2 : {
    double g, h, u0, u1, u2;
    for (; j < ej; S += (Idiv2<<1), j++) {
      INITALG2(j);
      h = (double) -Idiv2;
      for (ih = -Idiv2; ih < Idiv2; ih += 4) {
	ALG2(0); h++; ALG2(1); h++; ALG2(2); h++; ALG2(3); h++;
      }}}
    return;
  case 3 : {
    double g, h, u0, u1, u2, u3;
    for (; j < ej; S += (Idiv2<<1), j++) {
      INITALG3(j);
      h = (double) -Idiv2;
      for (ih = -Idiv2; ih < Idiv2; ih += 4) {
	ALG3(0); h++; ALG3(1); h++; ALG3(2); h++; ALG3(3); h++;
      }}}
    return;
  case 4 : {
    double g, h, u0, u1, u2, u3, u4;
    for (; j < ej; S += (Idiv2<<1), j++) {
      INITALG4(j);
      h = (double) -Idiv2;
      for (ih = -Idiv2; ih < Idiv2; ih += 4) {
	ALG4(0); h++; ALG4(1); h++; ALG4(2); h++; ALG4(3); h++;
      }}}
    return;
  case 5 : {
    double g, h, u0, u1, u2, u3, u4, u5;
    for (; j < ej; S += (Idiv2<<1), j++) {
      INITALG5(j);
      h = (double) -Idiv2;
      for (ih = -Idiv2; ih < Idiv2; ih += 4) {
	ALG5(0); h++; ALG5(1); h++; ALG5(2); h++; ALG5(3); h++;
      }}}
    return;
  case 6 : {
    double g, h, u0, u1, u2, u3, u4, u5, u6;
    for (; j < ej; S += (Idiv2<<1), j++) {
      INITALG6(j);
      h = (double) -Idiv2;
      for (ih = -Idiv2; ih < Idiv2; ih += 4) {
	ALG6(0); h++; ALG6(1); h++; ALG6(2); h++; ALG6(3); h++;
      }}}
    return;
  case 7 : {
    double g, h, u0, u1, u2, u3, u4, u5, u6, u7;
    for (; j < ej; S += (Idiv2<<1), j++) {
      INITALG7(j);
      h = (double) -Idiv2;
      for (ih = -Idiv2; ih < Idiv2; ih += 4) {
	ALG7(0); h++; ALG7(1); h++; ALG7(2); h++; ALG7(3); h++;
      }}}
    return;
  case 8 : {
    double g, h, u0, u1, u2, u3, u4, u5, u6, u7, u8;
    for (; j < ej; S += (Idiv2<<1), j++) {
      INITALG8(j);
      h = (double) -Idiv2;
      for (ih = -Idiv2; ih < Idiv2; ih += 4) {
	ALG8(0); h++; ALG8(1); h++; ALG8(2); h++; ALG8(3); h++;
      }}}
    return;
  case 9 : {
    double g, h, u0, u1, u2, u3, u4, u5, u6, u7, u8, u9;
    for (; j < ej; S += (Idiv2<<1), j++) {
      INITALG9(j);
      h = (double) -Idiv2;
      for (ih = -Idiv2; ih < Idiv2; ih += 4) {
	ALG9(0); h++; ALG9(1); h++; ALG9(2); h++; ALG9(3); h++;
      }}}
    return;
  default : {
    double g, h, u[d+1];
    for (; j < ej; S += (Idiv2<<1), j++) {
      INITALGD(j);
      h = (double) -Idiv2;
      for (ih = -Idiv2; ih < Idiv2; ih += 4) {
	ALGD(0); h++; ALGD(1); h++; ALGD(2); h++; ALGD(3); h++;
      }}}
    return;
  }
#undef ALG1
#undef ALG2
#undef ALG3
#undef ALG4
#undef ALG5
#undef ALG6
#undef ALG7
#undef ALG8
#undef ALG9
#undef ALGD
}
#endif
#endif
#endif

/********************* Undef of all non SSE defines ********************/
#undef INITALG0
#undef INITALG1
#undef INITALG2
#undef INITALG3
#undef INITALG4
#undef INITALG5
#undef INITALG6
#undef INITALG7
#undef INITALG8
#undef INITALG9
#undef INITALGD
#undef ABSG1
#undef ABSG2
#undef ABSG3
#undef ABSG4
#undef ABSG5
#undef ABSG6
#undef ABSG7
#undef ABSG8
#undef ABSG9
#undef ABSGD
#undef ALG_INIT_SCALE_ADD
#undef ALG_INIT_SCALE_ADD_T
/***************** end of undef of all non SSE defines *****************/
/********* End of the 8 procedures init_alg_norms_bucket_region ********/

/* return max |g(x)| for x in (0, s) where s can be negative,
   and g(x) = g[d]*x^d + ... + g[1]*x + g[0] */
static double
get_maxnorm_aux (double *g, const unsigned int d, double s)
{
  unsigned int k, l, sign_change, new_sign_change;
  double **dg;    /* derivatives of g */
  double a, va, b, vb;
  double *roots, gmax;

  dg = (double**) malloc (d * sizeof (double*));
  FATAL_ERROR_CHECK(dg == NULL, "malloc failed");
  dg[0] = g;
  for (k = 1; k < d; k++) { /* dg[k] is the k-th derivative, thus has
                             degree d-k, i.e., d-k+1 coefficients */
    dg[k] = (double*) malloc ((d - k + 1) * sizeof (double));
    FATAL_ERROR_CHECK(dg[k] == NULL, "malloc failed");
  }
  roots = (double*) malloc (d * sizeof (double));
  FATAL_ERROR_CHECK(roots == NULL, "malloc failed");
  for (k = 1; k < d; k++)
    for (l = 0; l <= d - k; l++)
      dg[k][l] = (l + 1) * dg[k - 1][l + 1];
  /* now dg[d-1][0]+x*dg[d-1][1] is the (d-1)-th derivative: it can have at
     most one sign change in (0, s), this happens iff dg[d-1][0] and
     dg[d-1][0]+s*dg[d-1][1] have different signs */
  if (dg[d-1][0] * (dg[d-1][0] + s * dg[d-1][1]) < 0)
    {
      sign_change = 1;
      roots[0] = - dg[d-1][0] / dg[d-1][1]; /* root of (d-1)-th derivative */
    }
  else
    sign_change = 0;
  roots[sign_change] = s; /* end of interval */
  for (k = d - 1; k-- > 1;)
    {
      /* invariant: sign_change is the number of sign changes of the
         (k+1)-th derivative, with corresponding roots in roots[0]...
         roots[sign_change-1], and roots[sign_change] = s. */
      a = 0.0;
      va = dg[k][0]; /* value of dg[k] at x=0 */
      new_sign_change = 0;
      for (l = 0; l <= sign_change; l++)
        {
          b = roots[l]; /* root of dg[k+1], or end of interval */
          vb = fpoly_eval (dg[k], d - k, b);
          if (va * vb < 0) /* root in interval */
            roots[new_sign_change++] = fpoly_dichotomy (dg[k], d - k,
                                                        a, b, va, 20);
          a = b;
          va = vb;
        }
      roots[new_sign_change] = s; /* end of interval */
      sign_change = new_sign_change;
    }
  /* now all extrema of g are 0, roots[0], ..., roots[sign_change] = s */
  gmax = fabs (g[0]);
  for (k = 0; k <= sign_change; k++)
    {
      va = fabs (fpoly_eval (g, d, roots[k]));
      if (va > gmax)
        gmax = va;
    }
  free (roots);
  for (k = 1; k < d; k++)
    free (dg[k]);
  free (dg);
  return gmax;
}

#if 1
/* returns the maximal value of log2|F(a,b)/q|, or log2|F(a,b)| if ratq=0,
   for a = a0 * i + a1 * j, b = b0 * i + b1 * j and q >= q0,
   -I/2 <= i <= I/2, 0 <= j <= I/2*min(s*B/|a1|,B/|b1|)
   where B >= sqrt(2*q/s/sqrt(3)) for all special-q in the current range
   (s is the skewness, and B = l_infty_snorm_u).

   Since |a0| <= s*B and |b0| <= B, then
   |a0 * i + a1 * j| <= s*B*I and |b0 * i + b1 * j| <= B*I,
   thus it suffices to compute M = max |F(x,y)| in the rectangle
   -s <= x <= s, 0 <= y <= 1, and to multiply M by (B*I)^deg(F).

   Since F is homogeneous, we know M = max |F(x,y)| is attained on the border
   of the rectangle, i.e.:
   (a) either on F(s,y) for 0 <= y <= 1
   (b) either on F(x,1) for -s <= x <= s
   (c) either on F(-s,y) for 0 <= y <= 1
   (d) or on F(x,0) for -s <= x <= s, but this maximum is f[d]*s^d,
       and is attained in (a) or (c).
*/
static double
get_maxnorm_alg (cado_poly cpoly, sieve_info_ptr si, double q0d, double l_infty_snorm_u, int qside)
{
  unsigned int d = cpoly->alg->degree, k;
  double *fd; /* double-precision coefficients of f */
  double norm, max_norm, pows, tmp;
  double B = l_infty_snorm_u;

  fd = (double*) malloc ((d + 1) * sizeof (double));
  FATAL_ERROR_CHECK(fd == NULL, "malloc failed");
  for (k = 0; k <= d; k++)
    fd[k] = mpz_get_d (cpoly->alg->f[k]);

  /* (b1) determine the maximum of |f(x)| for 0 <= x <= s */
  max_norm = get_maxnorm_aux (fd, d, cpoly->skew);

  /* (b2) determine the maximum of |f(-x)| for 0 <= x <= s */
  norm = get_maxnorm_aux (fd, d, -cpoly->skew);
  if (norm > max_norm)
    max_norm = norm;

  for (pows = 1.0, k = 0; k <= d; k++)
    {
      fd[k] *= pows;
      pows *= cpoly->skew;
    }
  /* swap coefficients; if d is odd, we need to go up to k = floor(d/2) */
  for (k = 0; k <= d / 2; k++)
    {
      tmp = fd[k];
      fd[k] = fd[d - k];
      fd[d - k] = tmp;
    }

  /* (a) determine the maximum of |g(y)| for 0 <= y <= 1, with g(y) = F(s,y) */
  norm = get_maxnorm_aux (fd, d, 1.0);
  if (norm > max_norm)
    max_norm = norm;

  /* (c) determine the maximum of |g(-y)| for 0 <= y <= 1 */
  norm = get_maxnorm_aux (fd, d, -1.0);
  if (norm > max_norm)
    max_norm = norm;

  free (fd);

  /* multiply by (B*I)^d: don't use the pow() function since it does not
     work properly when the rounding mode is not to nearest (at least under
     Linux with the glibc). Moreover for a small exponent d a direct loop
     as follows should not be much slower (if any), anyway the efficiency of
     that function is not critical */
  for (tmp = max_norm, k = 0; k < d; k++)
    tmp *= B * (double) si->I;
  /* divide by q0 if sieving on alg side */
  if (qside == ALGEBRAIC_SIDE)
      tmp /= q0d;
  return log2(tmp);
}
#else
/* simpler but less accurate version */
static double
get_maxnorm_alg (cado_poly cpoly, sieve_info_ptr si, uint64_t q0)
{
  unsigned int d = cpoly->alg->degree, k;
  double *fd, maxnorm;

  /* if F(a,b) = f[d]*a^d + f[d-1]*a^(d-1)*b + ... + f[0]*b^d,
     then |F(a,b)| <= |f[d]|*|a|^d + |f[d-1]|*|a|^(d-1)*b + ... + f[0]*b^d
     and the maximum is attained for a=s and b=1 (see above) */

  fd = (double*) malloc ((d + 1) * sizeof (double));
  FATAL_ERROR_CHECK(fd == NULL, "malloc failed");
  for (k = 0; k <= d; k++)
    {
      fd[k] = fabs (mpz_get_d (cpoly->alg->f[k]));
      fprintf (stderr, "f[%d]=%e\n", k, fd[k]);
    }
  maxnorm = fpoly_eval (fd, d, cpoly->skew);
  fprintf (stderr, "s=%f maxnorm=%e\n", cpoly->skew, maxnorm);
  free (fd);
  for (k = 0; k < d; k++)
    maxnorm *= si->B * (double) si->I;
  if (!si->ratq)
    maxnorm /= (double) q0;
  return log2 (maxnorm);
}
#endif

/* this function initializes the scaling factors and report bounds on the
   rational and algebraic sides */
void sieve_info_init_norm_data(FILE * output, sieve_info_ptr si, double q0d, int qside)
{
  double step, begin;
  for (int side = 0; side < 2; side++)
    {
      int d = si->cpoly->pols[side]->degree;
      si->sides[side]->fij = (mpz_t *) malloc((d + 1) * sizeof(mpz_t));
      FATAL_ERROR_CHECK(si->sides[side]->fij == NULL, "malloc failed");
      si->sides[side]->fijd = (double *) malloc_aligned((d + 1) * sizeof(double), 16);
      FATAL_ERROR_CHECK(si->sides[side]->fijd == NULL, "malloc failed");
      for (int k = 0; k <= d; k++)
        mpz_init (si->sides[side]->fij[k]);
    }

  double r, maxlog2;
  sieve_side_info_ptr rat = si->sides[RATIONAL_SIDE];
  sieve_side_info_ptr alg = si->sides[ALGEBRAIC_SIDE];

  /* Let u = (a0, b0) and v = (a1, b1) the two vectors obtained from skew
     Gauss reduction, and u' = (a0, s*b0), v' = (a1, s*b1). Assume |u'|<=|v'|.
     We know from Gauss reduction that the vectors u' and v' form an
     angle of at least pi/3, and their determinant is qs, thus:
     |u'|^2 <= |u'|*|v'| <= qs/sin(pi/3) = 2/sqrt(3)*q*s.
     Define B := sqrt(2/sqrt(3)*q/s), then |a0| <= s*B and |b0| <= B. */

  /* this was called si->B earlier. This bounds the skewed-L\infty norm
   * of u=(a0,b0), following |a0| <= s*B and |b0| <= B.
   * (IMO we have another notational convention elsewhere regarding what
   * the skewed norm is. But it's the one which is used here).
   */
  double B = sqrt (2.0 * q0d / (si->cpoly->skew * sqrt (3.0)));

  /************************** rational side **********************************/

  /* If J is chosen such that J<=I/2*s*B/max(|a1|,s*|b1|),
   * then |j*a1| <= I/2*s*B, and |j*b1| <= I/2*B
   * J is set to honour this requirement in sieve_info_adjust_IJ (but see
   * also bug #15617).
   * Now (a,b)=i*(a0,b0)+j*(a1,b1) with |i|<I/2 and |j|<=J gives
   * |a| <= s*I*B and |b| <= I*B, whence |G(a,b)| <= (|g[1]|*s+|g[0]|) * I*B */
  r = fabs (mpz_get_d (si->cpoly->rat->f[1])) * si->cpoly->skew
    + fabs (mpz_get_d (si->cpoly->rat->f[0]));
  r *= B * (double) si->I;

  /* if the special-q is on the rational side, divide by it */
  if (qside == RATIONAL_SIDE)
    r /= q0d;

  rat->logmax = log2 (r);

  /* we increase artificially 'logmax', to allow larger values of J */
  rat->logmax += 1.0;

  /* we know that |G(a,b)| < 2^(rat->logmax) when si->ratq = 0,
     and |G(a,b)/q| < 2^(rat->logmax) when si->ratq <> 0 */

  maxlog2 = rat->logmax;
  fprintf (output, "# Rat. side: log2(maxnorm)=%1.2f logbase=%1.6f",
           maxlog2, exp2 (maxlog2 / ((double) UCHAR_MAX - GUARD)));
  /* we want to map 0 <= x < maxlog2 to GUARD <= y < UCHAR_MAX,
     thus y = GUARD + x * (UCHAR_MAX-GUARD)/maxlog2 */
  rat->scale = ((double) UCHAR_MAX - GUARD) / maxlog2;
  step = 1 / rat->scale;
  begin = -step * GUARD;
  for (unsigned int inc = 0; inc < 257; begin += step) rat->cexp2[inc++] = exp2(begin);
  /* we want to select relations with a cofactor of less than r bits on the
     rational side */
  r = MIN(si->conf->sides[RATIONAL_SIDE]->lambda * (double) si->conf->sides[RATIONAL_SIDE]->lpb, maxlog2 - GUARD / rat->scale);
  rat->bound = (unsigned char) (r * rat->scale + GUARD);
  fprintf (output, " bound=%u\n", rat->bound);
  double max_rlambda = (maxlog2 - GUARD / rat->scale) / si->cpoly->rat->lpb;
  if (si->cpoly->rat->lambda > max_rlambda) {
      fprintf(output, "# Warning, rlambda>%.1f does not make sense (capped to limit)\n", max_rlambda);
  }
  /* Obsolete: rat->Bound is replaced by a single threshold alg->bound */
  /*
  sieve_info_init_lognorm (rat->Bound, rat->bound, si->conf->sides[RATIONAL_SIDE]->lim,
                           si->conf->sides[RATIONAL_SIDE]->lpb, rat->scale);
  */

  /************************** algebraic side *********************************/

  alg->logmax = get_maxnorm_alg (si->cpoly, si, q0d, B, qside); /* log2(max norm) */
  /* we know that |F(a,b)/q| < 2^(alg->logmax) when si->ratq = 0,
     and |F(a,b)| < 2^(alg->logmax) when si->ratq <> 0 */

  /* we increase artificially 'logmax', to allow larger values of J */
  alg->logmax += 2.0;

  /* on the algebraic side, we want that the non-reports on the rational
     side, which are set to 255, remain larger than the report bound 'r',
     even if the algebraic norm is totally smooth. For this, we artificially
     increase by 'r' the maximal range.
     If lambda * lpb < logmax, which is the usual case, then non-reports on
     the rational side will start from logmax + lambda * lpb or more, and
     decrease to lambda * lpb or more, thus above the threshold as we want.
     If logmax < lambda * lpb, non-reports on the rational side will start
     from 2*logmax or more, and decrease to logmax or more. */
  r = MIN(si->conf->sides[ALGEBRAIC_SIDE]->lambda * (double) si->conf->sides[ALGEBRAIC_SIDE]->lpb, alg->logmax);
  maxlog2 = alg->logmax + r;

  fprintf (output, "# Alg. side: log2(maxnorm)=%1.2f logbase=%1.6f",
           alg->logmax, exp2 (maxlog2 / ((double) UCHAR_MAX - GUARD)));
  /* we want to map 0 <= x < maxlog2 to GUARD <= y < UCHAR_MAX,
     thus y = GUARD + x * (UCHAR_MAX-GUARD)/maxlog2 */
  alg->scale = ((double) UCHAR_MAX - GUARD) / maxlog2;
  step = 1 / alg->scale;
  begin = -step * GUARD;
  for (unsigned int inc = 0; inc < 257; begin += step) alg->cexp2[inc++] = exp2(begin);
  /* we want to report relations with a remaining log2-norm after sieving of
     at most lambda * lpb, which corresponds in the y-range to
     y >= GUARD + lambda * lpb * scale */
  alg->bound = (unsigned char) (r * alg->scale + GUARD);
  fprintf (output, " bound=%u\n", alg->bound);
  double max_alambda = (alg->logmax) / si->cpoly->alg->lpb;
  if (si->cpoly->alg->lambda > max_alambda) {
      fprintf(output, "# Warning, alambda>%.1f does not make sense (capped to limit)\n", max_alambda);
  }
  /* Obsolete: alg->Bound is replaced by a single threshold alg->bound */
  /*
  sieve_info_init_lognorm (alg->Bound, alg->bound, si->conf->sides[ALGEBRAIC_SIDE]->lim,
                           si->conf->sides[ALGEBRAIC_SIDE]->lpb, alg->scale);
  */
}

void sieve_info_clear_norm_data(sieve_info_ptr si)
{
    for(int side = 0 ; side < 2 ; side++) {
        sieve_side_info_ptr s = si->sides[side];
        cado_poly_side_ptr ps = si->cpoly->pols[side];
        for (int k = 0; k <= ps->degree; k++)
            mpz_clear(s->fij[k]);
        free(s->fij);
        free(s->fijd);
    }
}

/* return largest possible J by simply bounding the Fij and Gij polynomials
   using the absolute value of their coefficients */
double
sieve_info_update_norm_data_Jmax (sieve_info_ptr si)
{
  double Iover2 = (double) (si->I >> 1);
  double Jmax = Iover2;
  double F[MAXDEGREE + 1];
  for (int side = 0; side < 2; side++)
    {
      sieve_side_info_ptr s = si->sides[side];
      cado_poly_side_ptr ps = si->cpoly->pols[side];
      double maxnorm = pow (2.0, s->logmax), v, powIover2 = 1.0;
      for (int k = 0; k <= ps->degree; k++)
        {
          /* reverse the coefficients since fij[k] goes with i^k but j^(d-k) */
          F[ps->degree - k] = fabs (s->fijd[k]) * powIover2;
          powIover2 *= Iover2;
        }
      v = fpoly_eval (F, ps->degree, Jmax);
      if (v > maxnorm)
        { /* use dichotomy to determine largest Jmax */
          double a, b, c;
          a = 0.0;
          b = Jmax;
          while (trunc (a) != trunc (b))
            {
              c = (a + b) * 0.5;
              v = fpoly_eval (F, ps->degree, c);
              if (v < maxnorm)
                a = c;
              else
                b = c;
            }
          Jmax = trunc (a) + 1; /* +1 since we don't sieve for j = Jmax */
        }
    }
  return (unsigned int) Jmax;
}

void
sieve_info_update_norm_data (sieve_info_ptr si, int nb_threads)
{
    int64_t H[4] = { si->a0, si->b0, si->a1, si->b1 };

    /* Update floating point version of algebraic poly (do both, while
     * we're at it...) */
    for (int side = 0; side < 2; side++) {
        sieve_side_info_ptr s = si->sides[side];
        cado_poly_side_ptr ps = si->cpoly->pols[side];
        mp_poly_homography (s->fij, ps->f, ps->degree, H);
        /* On the special-q side, divide all the coefficients of the 
           transformed polynomial by q */
        if (si->conf->side == side) {
            for (int i = 0; i <= ps->degree; i++) {
                ASSERT_ALWAYS(mpz_divisible_p(s->fij[i], si->doing->p));
                mpz_divexact(s->fij[i], s->fij[i], si->doing->p);
            }
        }
        for (int k = 0; k <= ps->degree; k++)
            s->fijd[k] = mpz_get_d (s->fij[k]);
    }

    /* improve bound on J if possible */
    unsigned int Jmax;
    Jmax = sieve_info_update_norm_data_Jmax (si);
    if (Jmax > si->J)
      {
        /* see sieve_info_adjust_IJ */
        ASSERT_ALWAYS(LOG_BUCKET_REGION >= si->conf->logI);
        uint32_t i = 1U << (LOG_BUCKET_REGION - si->conf->logI);
        i *= nb_threads;
        si->J = (Jmax / i) * i; /* cannot be zero since the previous value
                                   of si->J was already a multiple of i,
                                   and this new value is larger */
      }
}

