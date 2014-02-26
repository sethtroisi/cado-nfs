#define __STDC_LIMIT_MACROS
#include "cado.h"
#include <string.h>
#include <limits.h>
#include <math.h>               /* ceil */

#ifdef HAVE_SSE41
#include <smmintrin.h>
#elif defined(HAVE_SSSE3)
#include <tmmintrin.h>
#elif defined(HAVE_SSE3)
#include <pmmintrin.h>
#elif defined(HAVE_SSE2)
#include <emmintrin.h>
#endif

#include "las-config.h"
#include "las-debug.h"
#include "las-norms.h"
#include "utils.h"
#include "portability.h"

/* Input: i, double. i >= 0 needed!
   Output: o , trunc(o) == trunc(log2(i)) && o <= log2(i) < o + 0.0861.
   Careful: o ~= log2(i) iif add = 0x3FF00000 & scale = 1/0x100000.
   Add & scale are need to compute o'=f(log2(i)) where f is an affine function.
*/
static inline uint8_t inttruncfastlog2(double i, double add, double scale) {
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
  __asm__ __volatile__ (
            "psrlq $0x20,  %0    \n"
	    "cvtdq2pd      %0, %0\n" /* Mandatory in packed double even it's non packed! */
	    : "+&x" (i));             /* Really need + here! i is not modified in C code! */
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
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
  __asm__ __volatile__ (
            "psrlq $0x20,  %0                 \n"
	    "cvtdq2pd      %0,              %0\n"
	    : "+&x" (i));
  i = (i - add) * scale;
  __asm__ __volatile__ ( 
	    "cvttpd2dq     %0,       %0       \n" /* 0000 0000 0000 000Y */
	    "punpcklbw     %0,       %0       \n" /* 0000 0000 0000 00YY */
	    "pshuflw    $0x00,       %0,    %0\n" /* 0000 0000 YYYY YYYY */
	    : "+&x" (i));
  *(double *)&addr[decal] = i;
#else
  void *tg = &i;
  *(uint64_t *)&addr[decal] = 0x0101010101010101 * (uint64_t) (((double)(*(uint64_t *)tg >> 0x20) - add) * scale);
#endif
}

#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
/* This function is for the SSE2 init algebraics.
   I prefer use X86 ASM directly and avoid intrinsics because the trick of
   cvtdq2pd (I insert 2 doubles values and do a cvtdq2pd on them in order to
   compute their log2).
   CAREFUL! This function adds (1., 1.) to i to avoid negative value in the
   exposant (numbers < 1.) !
*/
static inline void w128itruncfastlog2fabs(__m128d i, __m128d add, __m128d scale, uint8_t *addr, ssize_t decal, __m128d un) {
  __asm__ __volatile__ (
	   "psllq     $0x01,    %0       \n" /* Dont use pabsd! */
	   "psrlq     $0x01,    %0       \n"
	   "addpd     %1,       %0       \n"
	   "shufps    $0xED,    %0,    %0\n"
	   "cvtdq2pd  %0,       %0       \n"
	   : "+&x"(i):"x"(un));
  i = _mm_mul_pd(_mm_sub_pd(i, add), scale);
  __asm__ __volatile__ (
	   "cvttpd2dq %0,       %0       \n" /* 0000 0000 000X 000Y */
	   "packssdw  %0,       %0       \n" /* 0000 0000 0000 0X0Y */
	   "punpcklbw %0,       %0       \n" /* 0000 0000 00XX 00YY */
	   "pshuflw   $0xA0,    %0,    %0\n" /* 0000 0000 XXXX YYYY */
	   "shufps    $0x50,    %0,    %0\n" /* XXXX XXXX YYYY YYYY */
	   : "+&x"(i));
  *(__m128d *)&addr[decal] = i; /* addr and decal are 16 bytes aligned: MOVAPD */
}
#endif

/* Initialize lognorms on the rational side for the bucket_region
 * number N.
 * For the moment, nothing clever, wrt discarding (a,b) pairs that are
 * not coprime, except for the line j=0.
 */
void init_rat_norms_bucket_region(unsigned char *S,
                                 unsigned int j,
                                 sieve_info_ptr si)
{
  /* #define COMPUTE_Y(G) ((LIKELY((G) > 1.)) ? inttruncfastlog2 ((G), add, scale) : GUARD) */
#define COMPUTE_Y(G) (inttruncfastlog2 ((G) + 1., add, scale))
  /* #define COMPUTE_Y(G) (inttruncfastlog2 ((G), add, scale)) */

  /* #define DEBUG_INIT_RAT 1 */ /* For internal debug: trace all */
  /* #define CHECK_INIT_RAT 1 */ /* For internal debug: control all */

  sieve_side_info_ptr rat = si->sides[RATIONAL_SIDE];
  const int halfI = (si->I)>>1;
  const double halfI_double = (double) halfI,
    halfI_double_minus_one = halfI_double - 1.0,
    u0 = si->sides[RATIONAL_SIDE]->fijd[0], // gj
    u1 = si->sides[RATIONAL_SIDE]->fijd[1], // gi
    scale = rat->scale * (1./0x100000),
    add = 0x3FF00000 - GUARD / scale;
  ASSERT_ALWAYS(u1 != 0.);
  const double invu1 = 1./u1;
  double u0j, d0_init, g, rac, d0, d1, i;
  size_t ts;
  unsigned int j1, inc;
  int int_i;
  uint8_t oy, y;

  ASSERT_ALWAYS(isfinite(invu1));
  j1 = LOG_BUCKET_REGION - si->conf->logI;
  j = j << j1;
  j1 = (1U << j1) + j;
  u0j = u0 * j;
  d0_init = rat->cexp2[((unsigned int)GUARD) - 1U];
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
       in this case, the right test will be g >= fabs(u0j)*(1./(1ULL<<51)) */
    /* Bug #16388 :
       u0 = 5.19229685853482763e+33, u1 = 4.27451782396958897e+33.
       So, u0j = 8.75421250348971938e+36, rac = -2.04800000000000000e+03.
       At the beginning, int_i = -2048, so the true value of g is 0.
       In fact, the computed value of g is 1.18059162071741130e+21.
       And u0j + g = 8.75421250348972056e+36; so u0j + g != u0j.
       => This test is false!
       => We have to use the slower but correct test,
       fabs(g) * (double) (((uint64_t) 1)<<51) >= fabs(u0j).
    */
    if (LIKELY(fabs(g) * (double) (((uint64_t) 1)<<51) >= fabs(u0j)))
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
       So, it's possible if y is very near trunc(y), old_i == int_i, so ts == 0.
       We have to iterate at least one time to avoid this case => this is the
       use of inc here. */
    for (i = rac + rat->cexp2[y] * invu1, inc = 1;; y--) {
      ts = -int_i;
      if (UNLIKELY(i >= halfI_double)) {
	ts += halfI;
#ifdef DEBUG_INIT_RAT
	fprintf (stderr, "A1.END : i1=%ld i2=%d, ts=%ld, y=%u, rac=%e\n", halfI - ts, halfI, ts, y, rac);
#endif
	memset(S, y, ts);
	S += ts;
	goto nextj;
      }
      int_i = (int) i; /* Overflow is not possible here */
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
    if (UNLIKELY(trunc(rac) >= halfI_double_minus_one)) {
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
    d0 = 1./d0;
    d1 = rac - d0 * rac;
    y++;
    i = rac - rat->cexp2[(unsigned int)y + 1] * invu1;
    for (;; y++) {
      ts = -int_i;
      if (UNLIKELY(i >= halfI_double)) {
	ts += halfI;
#ifdef DEBUG_INIT_RAT
	fprintf (stderr, "A4.END : i1=%ld, i2=%d, ts=%ld, y=%u, rac=%e\n", halfI - ts, halfI, ts, y, rac);
#endif
	memset(S, y, ts);
	S += ts;
	goto nextj;
      }
      int_i = (int) i; /* Overflow is not possible here */
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
      if (UNLIKELY(i >= halfI_double)) {
	ts += halfI;
#ifdef DEBUG_INIT_RAT
	fprintf (stderr, "B1.END : i1=%ld, i2=%d, ts=%ld, y=%u, rac=%e\n", halfI - ts, halfI, ts, y, rac);
#endif
	memset(S, y, ts);
	S += ts;
	goto nextj;
      }
      int_i = (int) i; /* Overflow is not possible here */
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
    if (UNLIKELY(trunc(rac) >= halfI_double_minus_one)) {
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
    d0 = 1./d0;
    d1 = rac - d0 * rac;
    y++;
    i = rac + rat->cexp2[(unsigned int)y + 1] * invu1;
    for (;; y++) {
      ts = -int_i;
      if (UNLIKELY(i >= halfI_double)) {
	ts += halfI;
#ifdef DEBUG_INIT_RAT
	fprintf (stderr, "B4.END : i1=%ld i2=%d, ts=%ld, y=%u, rac=%e\n", halfI - ts, halfI, ts, y, rac);
#endif
	memset(S, y, ts);
	S += ts;
	goto nextj;
      }
      int_i = (int) i; /* Overflow is not possible here */
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
      y = LIKELY(fabs(g) > 1.) ? log2(fabs(g))*rat->scale+GUARD : GUARD;
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

static inline void
poly_scale (double *u, const double *t, unsigned int d, double h)
{
  double hpow;
  u[d] = t[d];
  for (hpow = h; --d != UINT_MAX; hpow *= h) u[d] = t[d] * hpow;
}

/* Smart initialization of the algebraics. Computes the central
   initialization value of a box of (horizontail=8,vertical=p) and
   propagates it in the box. 8 is locked by the code and p is the
   minimum of VERT_NORM_STRIDE (#define, =4) and ej, with ej =
   2^(LOG_BUCKET_REGION-si->conf->logI).
   So, the fastest code is for logI <= 14.
   2 versions: SSE version (32bits compatible) & non SSE versions. */

/* Internal function, only with simple types, for unit/integration testing */
void init_alg_norms_bucket_region_internal (unsigned char *S, unsigned int j, size_t I, unsigned int d, double scale, double *fijd) { 
  
/* Macro to fill the others lines, used on SSE and non SSE versions.
   The first line is now OK; this macro copies it on the
   VERT_NORM_STRIDE-1 next lines. */
#define FILL_OTHER_LINES do {						\
    if (LIKELY(p > 1)) {						\
      unsigned char *mS = S + Idiv2;					\
      S -= Idiv2;							\
      for (unsigned int k = p; --k > 0; mS += I) aligned_medium_memcpy (mS, S, I); \
      S = mS + Idiv2; j += p;						\
    }									\
    else { S += I; j++; }						\
  } while (0)

  unsigned int ej = 1U << (LOG_BUCKET_REGION - ctzll((unsigned long long) I));
  ssize_t ih, Idiv2 = ((ssize_t) I) >> 1;
  unsigned int p = ej < VERT_NORM_STRIDE ? ej : VERT_NORM_STRIDE;
  j *= ej;
  ej += j;
  S += Idiv2;

#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM

  /*************** Some SSE macros for this function ***********************/
  /* This macro avoids a stupid & boring C types control */
#define _MM_SHUFFLE_EPI32(A,B,C) __asm__ __volatile__ ("pshufd $" #C ", %1, %0\n":"=x"(A):"x"(B))

  /* Initialisation of the data with intrinsics. */
#define ALG_INIT_SCALE_ADD_ONE_SIXTEEN const __m128d			\
    _one MAYBE_UNUSED = _mm_set1_pd(1.),				\
    _sixteen MAYBE_UNUSED = _mm_set1_pd(16.),				\
    _scale = _mm_set1_pd(scale * (1./0x100000)),			\
    _add = _mm_set1_pd(0x3FF00000 - GUARD / (scale * (1./0x100000)))	\
  /*************** End of the macros for this function *****************/
  
  /* This function is only a big switch in fact */
  switch (d) {
  case 0 : {
    ALG_INIT_SCALE_ADD_ONE_SIXTEEN;
    memset (S-Idiv2, inttruncfastlog2(fabs(fijd[0]), *(double *)&_add, *(double *)&_scale),
	    I * (ej - j)); 
  }
    return;
  case 1 : {
    __m128d h, g, u0, u1;
    ALG_INIT_SCALE_ADD_ONE_SIXTEEN;
    while (j < ej) {
      u1 = _mm_load1_pd(fijd+1);    
      u0 = _mm_set1_pd((double)(j + (p >> 1))*fijd[0]);
      h = _mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	g = h;
	g = _mm_add_pd(_mm_mul_pd(g,u1),u0);
	w128itruncfastlog2fabs(g, _add, _scale, S, ih, _one);
	h = _mm_add_pd(h, _sixteen);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 2 : {
    __m128d h, g, u0, u1, u2;
    ALG_INIT_SCALE_ADD_ONE_SIXTEEN;
    while (j < ej) {
      h = _mm_set1_pd((double)(j + (p >> 1)));
      u2 = _mm_load1_pd(fijd+2);
      u1 = _mm_load1_pd(fijd+1);
      u0 = _mm_load1_pd(fijd);
      u1 = _mm_mul_pd(u1,h);
      u0 = _mm_mul_pd(_mm_mul_pd(u0,h),h);
      h = _mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	g = h;
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(g,u2),u1),h),u0);
	w128itruncfastlog2fabs(g, _add, _scale, S, ih, _one);
	h = _mm_add_pd(h, _sixteen);
      }
      FILL_OTHER_LINES;
    }}	
    return;
  case 3 : {
    __m128d h, g, u0, u1, u2, u3;
    ALG_INIT_SCALE_ADD_ONE_SIXTEEN;
    while (j < ej) {
      h = _mm_set_sd((double)(j + (p >> 1)));
      _MM_SHUFFLE_EPI32(g, h, 0x44);
      g = _mm_mul_pd(g,g);
      u2 = _mm_load_pd(fijd+2);
      u0 = _mm_load_pd(fijd);
      u2 = _mm_mul_sd(u2,h);
      u0 = _mm_mul_pd(_mm_mul_sd(u0,h),g);
      _MM_SHUFFLE_EPI32(u1, u0, 0xEE);
      u0 = _mm_unpacklo_pd(u0,u0);
      _MM_SHUFFLE_EPI32(u3, u2, 0xEE);
      u0 = _mm_unpacklo_pd(u2,u2);
      h = _mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	g = h;
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(g,u3),u2),h),u1),h),u0);
	w128itruncfastlog2fabs(g, _add, _scale, S, ih, _one);
	h = _mm_add_pd(h, _sixteen);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 4 : {
    __m128d h, g, u0, u1, u2, u3, u4;
    ALG_INIT_SCALE_ADD_ONE_SIXTEEN;
    while (j < ej) {
      h = _mm_set1_pd((double)(j + (p >> 1)));
      h = _mm_mul_sd(h,h);
      _MM_SHUFFLE_EPI32(g, h, 0x44);
      u4 = _mm_load1_pd(fijd+4);
      u2 = _mm_load_pd(fijd+2);
      u0 = _mm_load_pd(fijd);
      u2 = _mm_mul_pd(u2,h);
      u0 = _mm_mul_pd(_mm_mul_pd(u0,h),g);
      _MM_SHUFFLE_EPI32(u1, u0, 0xEE);
      u0 = _mm_unpacklo_pd(u0,u0);
      _MM_SHUFFLE_EPI32(u3, u2, 0xEE);
      u2 = _mm_unpacklo_pd(u2,u2);
      h = _mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	g = h;
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(g,u4),u3),h),u2),h),u1),h),u0);
	w128itruncfastlog2fabs(g, _add, _scale, S, ih, _one);
	h = _mm_add_pd(h, _sixteen);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 5 : {
    __m128d h, g, u0, u1, u2, u3, u4, u5;
    ALG_INIT_SCALE_ADD_ONE_SIXTEEN;
    while (j < ej) {
      h = _mm_set_sd((double)(j + (p >> 1)));
      _MM_SHUFFLE_EPI32(g, h, 0x44);
      g = _mm_mul_pd(g,g);
      u4 = _mm_load_pd(fijd+4);
      u2 = _mm_load_pd(fijd+2);
      u0 = _mm_load_pd(fijd);
      u4 = _mm_mul_sd(u4,h);
      u2 = _mm_mul_pd(_mm_mul_sd(u2,h),g);
      u0 = _mm_mul_pd(_mm_mul_pd(_mm_mul_sd(u0,h),g),g);
      _MM_SHUFFLE_EPI32(u1, u0, 0xEE);
      u0 = _mm_unpacklo_pd(u0,u0);
      _MM_SHUFFLE_EPI32(u3, u2, 0xEE);
      u2 = _mm_unpacklo_pd(u2,u2);
      _MM_SHUFFLE_EPI32(u5, u4, 0xEE);
      u4 = _mm_unpacklo_pd(u4,u4);
      h = _mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	g = h; 
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(g,u5),u4),h),u3),h),u2),h),u1),h),u0);
	w128itruncfastlog2fabs(g, _add, _scale, S, ih, _one);
	h = _mm_add_pd(h, _sixteen);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 6 : {
    __m128d h, g, u0, u1, u2, u3, u4, u5, u6;
    ALG_INIT_SCALE_ADD_ONE_SIXTEEN;
    while (j < ej) {
      h = _mm_set1_pd((double)(j + (p >> 1)));
      h = _mm_mul_sd(h,h);
      _MM_SHUFFLE_EPI32(g, h, 0x44);
      u6 = _mm_load1_pd(fijd+6);
      u4 = _mm_load_pd(fijd+4);
      u2 = _mm_load_pd(fijd+2);
      u0 = _mm_load_pd(fijd);
      u4 = _mm_mul_pd(u4,h);
      u2 = _mm_mul_pd(_mm_mul_pd(u2,h),g);
      u0 = _mm_mul_pd(_mm_mul_pd(_mm_mul_pd(u0,h),g),g);
      _MM_SHUFFLE_EPI32(u1, u0, 0xEE);
      u0 = _mm_unpacklo_pd(u0,u0);
      _MM_SHUFFLE_EPI32(u3, u2, 0xEE);
      u2 = _mm_unpacklo_pd(u2,u2);
      _MM_SHUFFLE_EPI32(u5, u4, 0xEE);
      u4 = _mm_unpacklo_pd(u4,u4);
      h = _mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	g = h;
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(g,u6),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	w128itruncfastlog2fabs(g, _add, _scale, S, ih, _one);
	h = _mm_add_pd(h, _sixteen);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 7 : {
    __m128d h, g, u0, u1, u2, u3, u4, u5, u6, u7;
    ALG_INIT_SCALE_ADD_ONE_SIXTEEN;
    while (j < ej) {
      h = _mm_set_sd((double)(j + (p >> 1)));
      _MM_SHUFFLE_EPI32(g, h, 0x44);
      g = _mm_mul_pd(g,g);
      u6 = _mm_load_pd(fijd+6);
      u4 = _mm_load_pd(fijd+4);
      u2 = _mm_load_pd(fijd+2);
      u0 = _mm_load_pd(fijd);
      u6 = _mm_mul_sd(u6,h);
      u2 = _mm_mul_sd(u2,h);
      u4 = _mm_mul_pd(_mm_mul_sd(u4,h),g);
      u0 = _mm_mul_pd(_mm_mul_sd(u0,h),g);
      g = _mm_mul_pd(g,g);
      u2 = _mm_mul_pd(u2,g);
      u0 = _mm_mul_pd(u0,g);
      _MM_SHUFFLE_EPI32(u1, u0, 0xEE);
      u0 = _mm_unpacklo_pd(u0,u0);
      _MM_SHUFFLE_EPI32(u3, u2, 0xEE);
      u2 = _mm_unpacklo_pd(u2,u2);
      _MM_SHUFFLE_EPI32(u5, u4, 0xEE);
      u4 = _mm_unpacklo_pd(u4,u4);
      _MM_SHUFFLE_EPI32(u7, u6, 0xEE);
      u6 = _mm_unpacklo_pd(u6,u6);
      h = _mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	g = h; 
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(g,u7),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	w128itruncfastlog2fabs(g, _add, _scale, S, ih, _one);
	h = _mm_add_pd(h, _sixteen);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 8 : {
    __m128d h, g, u0, u1, u2, u3, u4, u5, u6, u7, u8;
    ALG_INIT_SCALE_ADD_ONE_SIXTEEN;
    while (j < ej) {
      h = _mm_set1_pd((double)(j + (p >> 1)));
      h = _mm_mul_sd(h,h);
      _MM_SHUFFLE_EPI32(g, h, 0x44);
      u8 = _mm_load1_pd(fijd+8);
      u6 = _mm_load_pd(fijd+6);
      u4 = _mm_load_pd(fijd+4);
      u2 = _mm_load_pd(fijd+2);
      u0 = _mm_load_pd(fijd);
      u6 = _mm_mul_pd(u6,h);
      u2 = _mm_mul_pd(u2,h);
      u4 = _mm_mul_pd(_mm_mul_pd(u4,h),g);
      u0 = _mm_mul_pd(_mm_mul_pd(u0,h),g);
      g = _mm_mul_pd(g,g);
      u2 = _mm_mul_pd(u2,g);
      u0 = _mm_mul_pd(u0,g);
      _MM_SHUFFLE_EPI32(u1, u0, 0xEE);
      u0 = _mm_unpacklo_pd(u0,u0);
      _MM_SHUFFLE_EPI32(u3, u2, 0xEE);
      u2 = _mm_unpacklo_pd(u2,u2);
      _MM_SHUFFLE_EPI32(u5, u4, 0xEE);
      u4 = _mm_unpacklo_pd(u4,u4);
      _MM_SHUFFLE_EPI32(u7, u6, 0xEE);
      u6 = _mm_unpacklo_pd(u6,u6);
      h = _mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	g = h; 
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(g,u8),u7),h),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	w128itruncfastlog2fabs(g, _add, _scale, S, ih, _one);
	h = _mm_add_pd(h, _sixteen);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 9 : {
    __m128d h, g, u0, u1, u2, u3, u4, u5, u6, u7, u8, u9;
    ALG_INIT_SCALE_ADD_ONE_SIXTEEN;
    while (j < ej) {
      h = _mm_set_sd((double)(j + (p >> 1)));
      _MM_SHUFFLE_EPI32(g, h, 0x44);
      g = _mm_mul_pd(g,g);
      u8 = _mm_load_pd(fijd+8);
      u6 = _mm_load_pd(fijd+6);
      u4 = _mm_load_pd(fijd+4);
      u2 = _mm_load_pd(fijd+2);
      u0 = _mm_load_pd(fijd);
      u8 = _mm_mul_sd(u8,h);
      u4 = _mm_mul_sd(u4,h);
      u0 = _mm_mul_sd(u0,h);
      u6 = _mm_mul_pd(_mm_mul_sd(u6,h),g);
      u2 = _mm_mul_pd(_mm_mul_sd(u2,h),g);
      g = _mm_mul_pd(g,g);
      u4 = _mm_mul_pd(u4,g);
      u2 = _mm_mul_pd(u2,g);
      g = _mm_mul_pd(g,g);
      u0 = _mm_mul_pd(u0,g);
      _MM_SHUFFLE_EPI32(u1, u0, 0xEE);
      u0 = _mm_unpacklo_pd(u0,u0);
      _MM_SHUFFLE_EPI32(u3, u2, 0xEE);
      u2 = _mm_unpacklo_pd(u2,u2);
      _MM_SHUFFLE_EPI32(u5, u4, 0xEE);
      u4 = _mm_unpacklo_pd(u4,u4);
      _MM_SHUFFLE_EPI32(u7, u6, 0xEE);
      u6 = _mm_unpacklo_pd(u6,u6);
      _MM_SHUFFLE_EPI32(u9, u8, 0xEE);
      u8 = _mm_unpacklo_pd(u8,u8);
      h =_mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	g = h;
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(g,u9),u8),h),u7),h),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	w128itruncfastlog2fabs(g, _add, _scale, S, ih, _one);
	h = _mm_add_pd(h, _sixteen);
      }
      FILL_OTHER_LINES;
    }}
    return;
  default: { /* This default could remplace all (optimized) previous cases. */
    __m128d h, g, u[d+1];
    ALG_INIT_SCALE_ADD_ONE_SIXTEEN;
    while (j < ej) {
      { /* Generic initialization of all the u?. */
	double du[d+1];
	poly_scale(du, fijd, d, (double) (j + (p >> 1)));
	for (unsigned int k = 0; k <= d; k++) u[k] = _mm_set1_pd(du[k]);
      }
      h =_mm_set_pd(11 - Idiv2, 3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 16) {
	g = u[d];
	for (unsigned int k = d; --k != UINT_MAX; g = _mm_add_pd(_mm_mul_pd(g,h),u[k]));
	w128itruncfastlog2fabs(g, _add, _scale, S, ih, _one);
	h = _mm_add_pd(h, _sixteen);
      }
      FILL_OTHER_LINES;
    }}
    return;
  } /* End of the big switch of the SSE version - and end of the function */
#undef _MM_SHUFFLE_EPI32
#undef ALG_INIT_SCALE_ADD_ONE_SIXTEEN

#else /* !HAVE_GCC_STYLE_AMD64_INLINE_ASM, end of SSE version of the function */

  /* Initialisation of the double data for the non SSE version. */
#define ALG_INIT_SCALE_ADD_ONE						\
  const double one MAYBE_UNUSED = 1., add;				\
  scale *= (1./0x100000); add = ((double) 0x3FF00000) - GUARD / scale
	
  /* This function is only a big switch in fact */
  switch (d) {
  case 0: {
    ALG_INIT_SCALE_ADD_ONE;
    memset(S - Idiv2, 0x0101010101010101 * inttruncfastlog2(fabs(fijd[0]),add,scale),
	   I * (ej - j));
  }
    return;
  case 1 : {
    double g, h, u0, u1;
    ALG_INIT_SCALE_ADD_ONE;
    while (j < ej) {
      u1 = fijd[1];			
      u0 = fijd[0] * ((double) (j + (p >> 1)));
      h = (double) (3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 8, h += 8.) {
	g = fabs(u0+h*u1);
	g += one;
	uint64truncfastlog2(g,add,scale,S,ih);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 2 : {
    double g, h, u0, u1, u2;
    ALG_INIT_SCALE_ADD_ONE;
    while (j < ej) {
      u2 = fijd[2]    ; h  = ((double) (j + (p >> 1)));
      u1 = fijd[1] * h; h *= h;
      u0 = fijd[0] * h;
      h = (double) (3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 8, h += 8.) {
	g = fabs(u0+h*(u1+h*u2));
	g += one;
	uint64truncfastlog2(g,add,scale,S,ih);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 3 : {
    double g, h, u0, u1, u2, u3;
    ALG_INIT_SCALE_ADD_ONE;
    while (j < ej) {
      u3 = fijd[3]    ; h  = ((double) (j + (p >> 1)));
      u2 = fijd[2] * h; g  = h * h;
      u1 = fijd[1] * g; g *= h;
      u0 = fijd[0] * g;
      h = (double) (3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 8, h += 8.) {
	g = fabs(u0+h*(u1+h*(u2+h*u3)));
	g += one;
	uint64truncfastlog2(g,add,scale,S,ih);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 4 : {
    double g, h, u0, u1, u2, u3, u4;
    ALG_INIT_SCALE_ADD_ONE;
    while (j < ej) {
      u4 = fijd[4]    ; h  = ((double) (j + (p >> 1)));
      u3 = fijd[3] * h; g  = h * h;
      u2 = fijd[2] * g; g *= h;
      u1 = fijd[1] * g; g *= h;
      u0 = fijd[0] * g;
      h = (double) (3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 8, h += 8.) {
	g = fabs(u0+h*(u1+h*(u2+h*(u3+h*u4))));
	g += one;
	uint64truncfastlog2(g,add,scale,S,ih);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 5 : {
    double g, h, u0, u1, u2, u3, u4, u5;
    ALG_INIT_SCALE_ADD_ONE;
    while (j < ej) {
      u5 = fijd[5]    ; h  = ((double) (j + (p >> 1)));
      u4 = fijd[4] * h; g  = h * h;
      u3 = fijd[3] * g; g *= h;
      u2 = fijd[2] * g; g *= h;
      u1 = fijd[1] * g; g *= h;
      u0 = fijd[0] * g;
      h = (double) (3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 8, h += 8.) {
	g = fabs(u0+h*(u1+h*(u2+h*(u3+h*(u4+h*u5)))));
	g += one;
	uint64truncfastlog2(g,add,scale,S,ih);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 6 : {
    double g, h, u0, u1, u2, u3, u4, u5, u6;
    ALG_INIT_SCALE_ADD_ONE;
    while (j < ej) {
      u6 = fijd[6]    ; h  = ((double) (j + (p >> 1)));
      u5 = fijd[5] * h; g  = h * h;
      u4 = fijd[4] * g; g *= h;
      u3 = fijd[3] * g; g *= h;
      u2 = fijd[2] * g; g *= h;
      u1 = fijd[1] * g; g *= h;
      u0 = fijd[0] * g;
      h = (double) (3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 8, h += 8.) {
	g = fabs(u0+h*(u1+h*(u2+h*(u3+h*(u4+h*(u5+h*u6))))));
	g += one;
	uint64truncfastlog2(g,add,scale,S,ih);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 7 : {
    double g, h, u0, u1, u2, u3, u4, u5, u6, u7;
    ALG_INIT_SCALE_ADD_ONE;
    while (j < ej) {
      u7 = fijd[7]    ; h  = ((double) (j + (p >> 1)));
      u6 = fijd[6] * h; g  = h * h;
      u5 = fijd[5] * g; g *= h;
      u4 = fijd[4] * g; g *= h;
      u3 = fijd[3] * g; g *= h;
      u2 = fijd[2] * g; g *= h;
      u1 = fijd[1] * g; g *= h;
      u0 = fijd[0] * g;
      h = (double) (3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 8, h += 8.) {
	g = fabs(u0+h*(u1+h*(u2+h*(u3+h*(u4+h*(u5+h*(u6+h*u7)))))));
	g += one;
	uint64truncfastlog2(g,add,scale,S,ih);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 8 : {
    double g, h, u0, u1, u2, u3, u4, u5, u6, u7, u8;
    ALG_INIT_SCALE_ADD_ONE;
    while (j < ej) {
      u8 = fijd[8]    ; h  = ((double) (j + (p >> 1)));
      u7 = fijd[7] * h; g  = h * h;
      u6 = fijd[6] * g; g *= h;
      u5 = fijd[5] * g; g *= h;
      u4 = fijd[4] * g; g *= h;
      u3 = fijd[3] * g; g *= h;
      u2 = fijd[2] * g; g *= h;
      u1 = fijd[1] * g; g *= h;
      u0 = fijd[0] * g;
      h = (double) (3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 8, h += 8.) {
	g = fabs(u0+h*(u1+h*(u2+h*(u3+h*(u4+h*(u5+h*(u6+h*(u7+h*u8))))))));
	g += one;
	uint64truncfastlog2(g,add,scale,S,ih);
      }
      FILL_OTHER_LINES;
    }}
    return;
  case 9 : {
    double g, h, u0, u1, u2, u3, u4, u5, u6, u7, u8, u9;
    ALG_INIT_SCALE_ADD_ONE;
    while (j < ej) {
      u9 = fijd[9]    ; h  = ((double) (j + (p >> 1)));
      u8 = fijd[8] * h; g  = h * h;
      u7 = fijd[7] * g; g *= h;
      u6 = fijd[6] * g; g *= h;
      u5 = fijd[5] * g; g *= h;
      u4 = fijd[4] * g; g *= h;
      u3 = fijd[3] * g; g *= h;
      u2 = fijd[2] * g; g *= h;
      u1 = fijd[1] * g; g *= h;
      u0 = fijd[0] * g;
      h = (double) (3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 8, h += 8.) {
	g = fabs(u0+h*(u1+h*(u2+h*(u3+h*(u4+h*(u5+h*(u6+h*(u7+h*(u8+h*u9)))))))));
	g += one;
	uint64truncfastlog2(g,add,scale,S,ih);
      }
      FILL_OTHER_LINES;
    }}
    return;
  default : {
    double g, h, u[d];
    ALG_INIT_SCALE_ADD_ONE;
    while (j < ej) {
      poly_scale(u, fijd, d, (double) (j + (p >> 1)));
      h = (double) (3 - Idiv2);
      for (ih = -Idiv2; ih < Idiv2; ih += 8, h += 8.) {
	g = u[d]; for (unsigned int k = d; --k != UINT_MAX; g = g*h+u[k]); g = fabs(g);
	g += one;
	uint64truncfastlog2(g,add,scale,S,ih);
      }
      FILL_OTHER_LINES;
    }}
    return;
  } /* End of the big switch of the non SSE version - and end of the function */
#undef ALG_INIT_SCALE_ADD_ONE
#endif /* End of !HAVE_GCC_STYLE_AMD64_INLINE_ASM */
#undef FILL_OTHER_LINES
} /* True end of the function (finally!) */

/* This function is used only to extract the interesting parameters of the
   complex structure si and call the previous function */
inline void init_alg_norms_bucket_region(unsigned char *S,
				  unsigned int j,
				  sieve_info_ptr si)
{
  init_alg_norms_bucket_region_internal (S, j, (size_t) si->I, si->cpoly->alg->deg, si->sides[ALGEBRAIC_SIDE]->scale, si->sides[ALGEBRAIC_SIDE]->fijd);
}

/* return max |g(x)| for x in (0, s) where s can be negative,
   and g(x) = g[d]*x^d + ... + g[1]*x + g[0] */
static double
get_maxnorm_aux (double_poly_srcptr poly, double s)
{
  double_poly_t deriv;
  const int d = poly->deg;

  if (d < 0) {
    return 0;
  } else if (d == 0) {
    return poly->coeff[0];
  }

  double *roots = (double*) malloc (poly->deg * sizeof (double));
  FATAL_ERROR_CHECK(roots == NULL, "malloc failed");

  /* Compute the derivative of polynomial */
  double_poly_init (deriv, d - 1);
  double_poly_derivative (deriv, poly);

  /* Look for extrema of the polynomial, i.e., for roots of the derivative */
  const unsigned int nr_roots = double_poly_compute_roots(roots, deriv, s);

  /* now abscissae of all extrema of poly are 0, roots[0], ..., 
     roots[nr_roots-1], s */
  double gmax = fabs (poly->coeff[0]);
  for (unsigned int k = 0; k <= nr_roots; k++)
    {
      double x = (k < nr_roots) ? roots[k] : s;
      double va = fabs (double_poly_eval (poly, x));
      if (va > gmax)
        gmax = va;
    }
  free (roots);
  double_poly_clear(deriv);
  return gmax;
}

/* Like get_maxnorm_aux(), but for interval [-s, s] */
static double
get_maxnorm_aux_pm (double_poly_srcptr poly, double s)
{
  double norm1 = get_maxnorm_aux(poly, s);
  double norm2 = get_maxnorm_aux(poly, -s);
  return (norm2 > norm1) ? norm2 : norm1;
}

/* returns the maximal norm of |F(x,y)| for -X <= x <= X, 0 <= y <= Y,
   where F(x,y) is a homogeneous polynomial of degree d.
   Let F(x,y) = f(x/y)*y^d, and F(x,y) = rev(F)(y,x).
   Since F is homogeneous, we know M = max |F(x,y)| is attained on the border
   of the rectangle, i.e.:
   (a) either on F(X,y) for -Y <= y <= Y (right-hand-side border, and the
       mirrored image of the left-hand-side border); We want the maximum of
       rev(F)(y,X) = rev(f)(y/X) * X^d for -Y <= j <= Y;
       this is rev(f)(t) * X^d for -Y/X <= t <= Y/X.
   (b) either on F(x,Y) for -X <= x <= X (top border)
       = f(x/Y) * Y^d; this is f(t) * Y^d for -X/Y <= t <= X/J.
   (d) or on F(x,0) for -X <= x <= X (lower border, on the abscissa), but this
       maximum is f[d]*X^d, and is attained in (a).
*/
static double
get_maxnorm_alg (double_poly_srcptr src_poly, const double X, const double Y)
{
  const unsigned int d = src_poly->deg;
  const int debug = 0;
  double norm, max_norm;

  /* Make copy of polynomial as we need to revert the coefficients */
  double_poly_t poly;
  double_poly_init (poly, d);
  double_poly_set (poly, src_poly);

  if (debug)
    double_poly_print(stdout, poly, "# Computing max norm for polynomial ");

  /* (b) determine the maximum of |f(x)| * Y^d for -X/Y <= x <= X/Y */
  max_norm = get_maxnorm_aux_pm (poly, X/Y) * pow(Y, (double)d);

  /* (a) determine the maximum of |g(y)| for -1 <= y <= 1, with g(y) = F(s,y) */
  double_poly_revert(poly);
  norm = get_maxnorm_aux_pm (poly, Y/X) * pow(X, (double)d);
  if (norm > max_norm)
    max_norm = norm;

  double_poly_clear(poly);

  if (debug)
    fprintf(stdout, "# Max norm is %f\n", max_norm);

  return max_norm;
}

/*
   Allocates:
     si->sides[side]->fij
     si->sides[side]->fijd
*/

void sieve_info_init_norm_data(sieve_info_ptr si)
{
  for (int side = 0; side < 2; side++)
    {
      int d = si->cpoly->pols[side]->deg;
      mpz_poly_init (si->sides[side]->fij, d);
      si->sides[side]->fijd = (double *) malloc_aligned((d + 1) * sizeof(double), 16);
      FATAL_ERROR_CHECK(si->sides[side]->fijd == NULL, "malloc failed");
    }

}

void sieve_info_clear_norm_data(sieve_info_ptr si)
{
    for(int side = 0 ; side < 2 ; side++) {
        sieve_side_info_ptr s = si->sides[side];
        mpz_poly_clear (s->fij);
        free(s->fijd);
    }
}

/* return largest possible J by simply bounding the Fij and Gij polynomials

  The image in the a,b-plane of the sieve region might be slanted at an angle
  against the abscissa, thus even though it has the correct skewness, it
  might include larger a,b-coordinates than a non-slanted image would.

  We compute the maximum norm that would occur if we had a perfect lattice
  basis in the sense that its image forms a rectangle A/2 <= a < A/2, 
  0 <= b < B, with A/2/B = skew, and A*B = I*J*q (assuming J=I/2 here).
  Thus we have B = A/2/skew, A*A/2/skew = I*I/2*q, A = I*sqrt(q*skew).
  The optimal maximum norm is then the maximum of |F(a,b)| in this rectangle,
  divided by q on the special-q side.

  Then we reduce J so that the maximum norm in the image of the actual sieve
  regionis no larger than this optimal maximum, times some constant fudge
  factor.
*/
static unsigned int
sieve_info_update_norm_data_Jmax (sieve_info_ptr si)
{
  const int verbose = 0;
  const double fudge_factor = 16.; /* How much bigger a norm than optimal
                                      we're willing to tolerate */
  const double I = (double) (si->I);
  const double q = mpz_get_d(si->doing->p);
  const double skew = si->cpoly->skew;
  const double A = I*sqrt(q*skew);
  const double B = A/2./skew;
  double Jmax = I/2.;

  for (int side = 0; side < 2; side++)
    {
      sieve_side_info_ptr s = si->sides[side];
      mpz_poly_ptr ps = si->cpoly->pols[side];

      /* Compute the best possible maximum norm, i.e., assuming a nice
         rectangular sieve region in the a,b-plane */
      double_poly_t dpoly;
      double_poly_init (dpoly, ps->deg);
      double_poly_set_mpz_poly (dpoly, ps);
      double maxnorm = get_maxnorm_alg (dpoly, A/2., B);
      double_poly_clear (dpoly);
      if (side == si->doing->side)
        maxnorm /= q;
      if (verbose) {
        printf ("Best possible maxnorm for side %d: %g\n", side, maxnorm);
      }

      maxnorm *= fudge_factor;
      if (verbose) {
        printf ("Threshold for acceptable norm for side %d: %g\n", side, maxnorm);
      }

      double_poly_t F;
      F->deg = ps->deg;
      F->coeff = s->fijd;
      
      double v = get_maxnorm_alg (F, I, Jmax);
      if (verbose) {
        printf ("Actual maxnorm for side %d with J=%d: %g\n", side, (int)Jmax, v);
      }
      
      if (v > maxnorm)
        { /* use dichotomy to determine largest Jmax */
          double a, b, c;
          a = 0.0;
          b = Jmax;
          while (trunc (a) != trunc (b))
            {
              c = (a + b) * 0.5;
              v = get_maxnorm_alg (F, I, c);
              if (verbose) {
                printf ("Actual maxnorm for side %d with J=%d: %g\n", side, (int)c, v);
              }
              if (v < maxnorm)
                a = c;
              else
                b = c;
            }
          Jmax = trunc (a) + 1; /* +1 since we don't sieve for j = Jmax */
          if (verbose) {
            printf ("Setting Jmax=%d\n", (int)Jmax);
          }
        }
    }

   if (verbose) {
     printf ("Final maxnorm J=%d\n", (int)Jmax);
   }
  return (unsigned int) Jmax;
}

/* this function initializes the scaling factors and report bounds on the
   rational and algebraic sides */
/*
   Updates:
     si->sides[RATIONAL_SIDE]->logmax
     si->sides[RATIONAL_SIDE]->scale
     si->sides[RATIONAL_SIDE]->cexp2[]
     si->sides[RATIONAL_SIDE]->bound
     
     si->sides[ALGEBRAIC_SIDE]->logmax
     si->sides[ALGEBRAIC_SIDE]->scale
     si->sides[ALGEBRAIC_SIDE]->cexp2[]
     si->sides[ALGEBRAIC_SIDE]->bound
*/

void
sieve_info_update_norm_data (FILE * output, sieve_info_ptr si, int nb_threads)
{
  int64_t H[4] = { si->a0, si->b0, si->a1, si->b1 };

  double step, begin;
  double r, maxlog2;
  double_poly_t poly;
  sieve_side_info_ptr rat = si->sides[RATIONAL_SIDE];
  sieve_side_info_ptr alg = si->sides[ALGEBRAIC_SIDE];

  /* Update floating point version of both polynomials. They will be used in
   * get_maxnorm_alg(). */
  for (int side = 0; side < 2; side++) {
      sieve_side_info_ptr s = si->sides[side];
      mpz_poly_ptr ps = si->cpoly->pols[side];
      mpz_poly_homography (s->fij, ps, H);
      /* On the special-q side, divide all the coefficients of the
         transformed polynomial by q */
      if (si->conf->side == side) {
          for (int i = 0; i <= ps->deg; i++) {
              ASSERT_ALWAYS(mpz_divisible_p(s->fij->coeff[i], si->doing->p));
              mpz_divexact(s->fij->coeff[i], s->fij->coeff[i], si->doing->p);
          }
      }
      for (int k = 0; k <= ps->deg; k++)
          s->fijd[k] = mpz_get_d (s->fij->coeff[k]);
  }

  /************************** rational side **********************************/

  /* Compute the maximum norm of the rational polynomial over the sieve
     region. The polynomial coefficient in fijd are already divided by q
     on the special-q side. */
  poly->deg = si->cpoly->pols[RATIONAL_SIDE]->deg;
  poly->coeff = si->sides[RATIONAL_SIDE]->fijd;
  rat->logmax = log2(get_maxnorm_alg (poly, (double)si->I/2, (double)si->I/2));

  /* we increase artificially 'logmax', to allow larger values of J */
  rat->logmax += 1.;

  /* we know that |G(a,b)| < 2^(rat->logmax) when si->ratq = 0,
     and |G(a,b)/q| < 2^(rat->logmax) when si->ratq <> 0 */

  maxlog2 = rat->logmax;
  fprintf (output, "# Rat. side: log2(maxnorm)=%1.2f logbase=%1.6f",
           maxlog2, exp2 (maxlog2 / ((double) UCHAR_MAX - GUARD)));
  /* we want to map 0 <= x < maxlog2 to GUARD <= y < UCHAR_MAX,
     thus y = GUARD + x * (UCHAR_MAX-GUARD)/maxlog2 */
  rat->scale = ((double) UCHAR_MAX - GUARD) / maxlog2;
  step = 1. / rat->scale;
  begin = -step * GUARD;
  for (unsigned int inc = 0; inc < 257; begin += step) rat->cexp2[inc++] = exp2(begin);
  /* we want to select relations with a cofactor of less than r bits on the
     rational side */
  r = MIN(si->conf->sides[RATIONAL_SIDE]->lambda * (double) si->conf->sides[RATIONAL_SIDE]->lpb, maxlog2 - GUARD / rat->scale);
  rat->bound = (unsigned char) (r * rat->scale + GUARD);
  fprintf (output, " bound=%u\n", rat->bound);
  double max_rlambda = (maxlog2 - GUARD / rat->scale) /
      si->conf->sides[RATIONAL_SIDE]->lpb;
  if (si->conf->sides[RATIONAL_SIDE]->lambda > max_rlambda) {
      fprintf(output, "# Warning, rlambda>%.1f does not make sense (capped to limit)\n", max_rlambda);
  }

  /************************** algebraic side *********************************/

  /* Compute the maximum norm of the algebraic polynomial over the sieve
     region. The polynomial coefficient in fijd are already divided by q
     on the special-q side. */
  poly->deg = si->cpoly->pols[ALGEBRAIC_SIDE]->deg;
  poly->coeff = si->sides[ALGEBRAIC_SIDE]->fijd;
  alg->logmax = log2(get_maxnorm_alg (poly, (double)si->I/2, (double)si->I/2));

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
  step = 1. / alg->scale;
  begin = -step * GUARD;
  for (unsigned int inc = 0; inc < 257; begin += step) alg->cexp2[inc++] = exp2(begin);
  /* we want to report relations with a remaining log2-norm after sieving of
     at most lambda * lpb, which corresponds in the y-range to
     y >= GUARD + lambda * lpb * scale */
  alg->bound = (unsigned char) (r * alg->scale + GUARD);
  fprintf (output, " bound=%u\n", alg->bound);
  double max_alambda = (alg->logmax) / si->conf->sides[ALGEBRAIC_SIDE]->lpb;
  if (si->conf->sides[ALGEBRAIC_SIDE]->lambda > max_alambda) {
      fprintf(output, "# Warning, alambda>%.1f does not make sense (capped to limit)\n", max_alambda);
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
