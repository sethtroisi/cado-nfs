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

#if 0
/* In gcc 4.8.2, fabs(__m128d) is done by X = andpd (X, c_and)
   with const _m128d c_and 0x7FFFFFFFFFFFFFFF7FFFFFFFFFFFFFFF.
   Id for fabs(double) with X = andsd.
   No so good because it needs a read memory operation.
   These 2 SSE operations are much faster. */
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
static inline ffabs_m128d (__m128d *i) {
  __asm__ __volatile__ (
	   "psllq     $0x01,    %0       \n" /* Dont use pabsd! */
	   "psrlq     $0x01,    %0       \n" /* i = fast_abs(i) */
	   : "+&x" (*i));
}
#endif

static inline ffabs_double (double *i) {
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
  ffabs_m128d ((__m128d *) i);
#else
  *i = fabs(*i);
#endif
}
#endif /* 0 */

/* Input: i, double. i > 0 needed! If not, the result has no sense at all.
   Output: o , trunc(o) == trunc(log2(i)) && o <= log2(i) < o + 0.0861.
   Careful: o ~= log2(i) iif add = 0x3FF00000 & scale = 1/0x100000.
   Add & scale are need to compute o'=f(log2(i)) where f is an affine function.
*/
static inline uint8_t lg2 (double i, double add, double scale) {
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
  return (uint8_t) (((double)(*(uint64_t *)tg >> 0x20) - add) * scale);
#endif
}

/* Same than previous with a fabs on i before the computation of the log2.
 */
static inline uint8_t lg2abs (double i, double add, double scale) {
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
  __asm__ __volatile__ (
            "psllq $0x01,  %0    \n"
            "psrlq $0x21,  %0    \n"
	    "cvtdq2pd      %0, %0\n" /* Mandatory in packed double even it's non packed! */
	    : "+&x" (i));             /* Really need + here! i is not modified in C code! */
  return (uint8_t) ((i-add)*scale);
#else
  void *tg = &i;
  return (uint8_t) (((double)((*(uint64_t *)tg << 1) >> 0x21) - add) * scale);
#endif
}

/* Same than previous with the result is duplicated 8 times
   in a "false" double in SSE version, i.e. in a xmm register, or in a 64 bits
   general register in non SSE version, and written at addr[decal].
   The SSE version has 4 interests: no memory use, no General Register use, 
   no xmm -> GR conversion, no * 0x0101010101010101 (3 to 4 cycles) but only
   2 P-instructions (1 cycle).
*/
#ifndef HAVE_GCC_STYLE_AMD64_INLINE_ASM
static inline void w64lg2abs(double i, double add, double scale, uint8_t *addr,
                             ssize_t decal)
{
  void *tg = &i;
  *(uint64_t *)&addr[decal] = 0x0101010101010101 *
    (uint64_t) (((double)((*(uint64_t *)tg << 1) >> 0x21) - add) * scale);
}
#endif

#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
/* This function is for the SSE2 init algebraics.
   I prefer use X86 ASM directly and avoid intrinsics because the trick of
   cvtdq2pd (I insert 2 doubles values and do a cvtdq2pd on them in order to
   compute their log2).
   This function does a double abs of the __m128d i, and computes the double
   log2 of the result.
   Careful: if i == 0, the result is not predictible.
*/
static inline void w128lg2abs(__m128d i, const __m128d add, const __m128d scale, uint8_t *addr, const ssize_t decal) {
  __asm__ __volatile__ (
	   "psllq     $0x01,    %0       \n" /* Dont use pabsd! */
	   "psrlq     $0x01,    %0       \n" /* i = fast_abs(i) */
	   "shufps    $0xED,    %0,    %0\n"
	   "cvtdq2pd  %0,       %0       \n"
	   : "+&x"(i));
  i = _mm_mul_pd(_mm_sub_pd(i, add), scale);
  __asm__ __volatile__ (
	   "cvttpd2dq %0,       %0       \n" /* 0000 0000 000X 000Y */
	   "packssdw  %0,       %0       \n" /* 0000 0000 0000 0X0Y */
	   "punpcklbw %0,       %0       \n" /* 0000 0000 00XX 00YY */
	   "punpcklwd %0,       %0       \n" /* 0000 XXXX 0000 YYYY */
	   "pshufd    $0xA0,    %0,    %0\n" /* XXXX XXXX YYYY YYYY */
	   : "+&x"(i));
  *(__m128d *)&addr[decal] = i; /* addr and decal are 16 bytes aligned: MOVAPD */
}

/* This function is for the SSE2 init algebraics, but for the exact initialization.
   Same than previous but return a SSE2 register with only 16 lowest bits are computed
   as the 2 results (2 8-bits values). 
   CAREFUL! This function returns the computed value in the data i.
   the i value is modified by this function !
*/
static inline __m128i _mm_lg2abs(__m128d *i, const __m128d add, const __m128d scale) {
  __asm__ __volatile__ (
	   "psllq     $0x01,    %0       \n" /* Dont use pabsd! */
	   "psrlq     $0x01,    %0       \n"
	   "shufps    $0xED,    %0,    %0\n"
	   "cvtdq2pd  %0,       %0       \n"
	   : "+&x"(*i));
  *i = _mm_mul_pd(_mm_sub_pd(*i, add), scale);
  __asm__ __volatile__ (
	   "cvttpd2dq %0,       %0       \n" /* 0000 0000 000X 000Y */
	   "packssdw  %0,       %0       \n" /* 0000 0000 0000 0X0Y */
	   "packuswb  %0,       %0       \n" /* 0000 0000 0000 00XY */
	   : "+&x"(*i));
  return *(__m128i *) i;
}
#endif /* HAVE_GCC_STYLE_AMD64_INLINE_ASM */

/* This function computes the roots of the polynome F(i,1) and
   the roots of the polynomes d^2(F(i,1))/d(i)^2, so the roots of
   F(i,1)F"(i,1)-F'(i,1)^2.
   The number of Roots is <= d + 2 * (d - 1), where d is the degree
   of F(i,1).
   After, the special "pseudo" roots 0.0 is add, in order to correct
   the smart initialization in the neigborhood of i = 0.
   These roots are useful to correct the smart initialization of the
   algebraics, because the computed values are false near the roots.

   if deg(f) = 0, nroots = 1 (pseudo root 0.0),
   if deg(f) >= 1, nroots <= 1 (root of f) + 0 (root of f") + 1 (pseudo root),.
   if deg(f) >= 2, nroots <= deg(f) + 2 * deg(f) - 2 + 1 = 3 * deg(f) - 1.
   So, the size of the array roots must be 3 * deg(f) + 1.

   This function is internal. Don't use it. Use the wrapper below.

   Input: F(i,1): degree et coefficients. max_abs_root = Absolute maximum value
   of a root (before and after, the root is useless).
   Output: nroots & roots.
*/
void init_norms_roots_internal (unsigned int degree, double *coeff, double max_abs_root, unsigned int *nroots, double *roots)
{
  const double_poly_t f = {{ degree, coeff }};
  double_poly_t df, ddf, f_ddf, df_df, d2f;
  root_struct Roots[2 * f->deg + 1];
  mpz_t           p[2 * f->deg + 1];
  unsigned int nroots_f, nroots_d2f;
  size_t k;
  for (k = 2 * f->deg + 1; k--; ) {
    mpz_init (p[k]);
    root_struct_init (&(Roots[k]));
  }
  for (k = f->deg + 1; k--; )
    mpz_set_d (p[k], f->coeff[k]);
  
  nroots_f = numberOfRealRoots (p, f->deg, max_abs_root, 0, Roots);
  for (k = nroots_f; k--; )
    roots[k] = rootRefine (&(Roots[k]), p, f->deg);

  double_poly_init (df, MAX(0,((int)f->deg - 1)));
  double_poly_derivative (df, f);
  double_poly_init (df_df, df->deg + df->deg);
  double_poly_product (df_df, df, df);

  double_poly_init (ddf, MAX(0,((int)df->deg - 1)));
  double_poly_derivative (ddf, df);
  double_poly_init (f_ddf, f->deg + ddf->deg);
  double_poly_product (f_ddf, f, ddf);
  
  double_poly_init (d2f, MAX(f_ddf->deg, df_df->deg));
  double_poly_substract (d2f, f_ddf, df_df);
  
  for (k = d2f->deg + 1; k--; )
    mpz_set_d (p[k], d2f->coeff[k]);

  double_poly_clear(df);
  double_poly_clear(ddf);
  double_poly_clear(f_ddf);
  double_poly_clear(df_df);
  double_poly_clear(d2f);

  nroots_d2f = numberOfRealRoots (p, d2f->deg, max_abs_root, 0, Roots);
  for (k = nroots_d2f; k--; )
    roots[k + nroots_f] = rootRefine (&(Roots[k]), p, d2f->deg);

  /* Pseudo root 0.0 */
  roots[nroots_f + nroots_d2f] = 0.0;
  *nroots = nroots_f + nroots_d2f + 1;

  for (k = f->deg * 2 + 1; k--; ) {
    mpz_clear (p[k]);
    root_struct_clear (&(Roots[k]));
  }
}

/* A wrapper for the function above */
void init_norms_roots (sieve_info_ptr si, unsigned int side)
{
  init_norms_roots_internal (si->cpoly->pols[side]->deg, si->sides[side]->fijd, (double) ((si->I + 16) >> 1), &(si->sides[side]->nroots), si->sides[side]->roots);
}

/* Initialize lognorms of the bucket region S[] number J, for F(i,j) with
 * degree = 1.
 * For the moment, nothing clever, wrt discarding (a,b) pairs that are
 * not coprime, except for the line j=0.
 */
/* Internal function, only with simple types, for unit/integration testing */
void init_degree_one_norms_bucket_region_internal (unsigned char *S, uint32_t J, uint32_t I, double scale, double u0, double u1, double *cexp2) {

  /* The computation of the log2 by inttruncfastlog2 needs a value >= 1.
     Here, in the classical rational initialization, the degree of the used
     polynome F(i,j) is hardcoded to 1. So, it's possible to have F(i,j) = 0,
     even if i, j and the coefficients of F(i,j) are integers and != 0.
     So, I add 1.0 on all G values.
     It's not useful to do a fabs(G) here because the code uses always COMPUTE_Y(G)
     with G >= 0. */
#define COMPUTE_Y(G) lg2 ((G) + 1., add, scale)
  /* For internal debug: print all */
  // #define DEBUG_INIT_RAT 1
  double add, d0_init, g, rac, u0J, d0, d1, i;
  size_t ts;
  uint32_t endJ;
  int int_i;
  unsigned int inc;
  uint8_t oy, y;

  ASSERT_ALWAYS (u1 != 0.);

  const int Idiv2 = (int) I >> 1;
  const double Idiv2_double = (double) Idiv2, Idiv2_double_minus_one = Idiv2_double - 1., invu1 = 1./u1;
  
  endJ = LOG_BUCKET_REGION - ctz(I);
  J <<= endJ;
  endJ = (1U << endJ) + J;

#ifdef DEBUG_INIT_RAT
  fprintf (stderr, "Begin: j=%u, u0=%.20e u1=%.20e, scale=%.20e, rac=%.20e\n", J, u0, u1, scale, u0 * j * (-invu1));
#endif

  scale *= 1./0x100000;
  add = 0x3FF00000 - GUARD / scale;
  u0J = u0 * J;
  d0_init = cexp2[((unsigned int)GUARD) - 1U];
  for (; J < endJ ; J++, u0J += u0) {
    int_i = -Idiv2;
    g = u0J + u1 * int_i;
    rac = u0J * (-invu1);
    d0 = d0_init;
    d1 = rac - d0 * rac;
    /* g sign is mandatory at the beginning of the line intialization.
       The sign of g is not significant if fabs(g) * 1ULL<<51 < fabs(u0J)
       In this case, the sign of g is the sign of u1, because g+u1
       is the next g value.
    */
    if (LIKELY(fabs(g) * (1ULL << 51) >= fabs(u0J)))
      if (signbit(g)) {
	g = -g;
	y = COMPUTE_Y(g);
	if (rac >= -Idiv2) goto cas3; else goto cas2;
      }
      else {
	y = COMPUTE_Y(g);
	if (rac >= -Idiv2) goto cas1; else goto cas4;
      }
    else {
      y = GUARD;
      if (signbit(u1)) goto cas2; else goto cas4;
    }
  cas1:
    /* In this case, we exit from the loop when ts == 0, at the exception
       of the first iteration. In this special case, old_i = -Idiv2 and
       int_i = trunc (i), where i=[inverse of the function g](trunc(y)) and
       y=g(old_i).
       So, it's possible if y is very near trunc(y), old_i == int_i, so ts == 0.
       We have to iterate at least one time to avoid this case => this is the
       use of inc here. */
    for (i = rac + cexp2[y] * invu1, inc = 1;; y--) {
      ts = -int_i;
      if (UNLIKELY(i >= Idiv2_double)) {
	ts += Idiv2;
#ifdef DEBUG_INIT_RAT
	fprintf (stderr, "A1.END : i1=%ld i2=%d, ts=%ld, y=%u, rac=%e\n", Idiv2 - ts, Idiv2, ts, y, rac);
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
    g = u0J + u1 * int_i;
    if (UNLIKELY(trunc(rac) >= Idiv2_double_minus_one)) {
      for ( ; int_i < Idiv2; int_i++) {
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
      if (++int_i >= Idiv2) goto nextj;
      oy = y;
      g -= u1;
      y = COMPUTE_Y(g);
    } while (oy != y);
    d0 = 1./d0;
    d1 = rac - d0 * rac;
    y++;
    i = rac - cexp2[(unsigned int)y + 1] * invu1;
    for (;; y++) {
      ts = -int_i;
      if (UNLIKELY(i >= Idiv2_double)) {
	ts += Idiv2;
#ifdef DEBUG_INIT_RAT
	fprintf (stderr, "A4.END : i1=%ld, i2=%d, ts=%ld, y=%u, rac=%e\n", Idiv2 - ts, Idiv2, ts, y, rac);
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
    for (i = rac - cexp2[y] * invu1, inc = 1;; y--) {
      ts = -int_i;
      if (UNLIKELY(i >= Idiv2_double)) {
	ts += Idiv2;
#ifdef DEBUG_INIT_RAT
	fprintf (stderr, "B1.END : i1=%ld, i2=%d, ts=%ld, y=%u, rac=%e\n", Idiv2 - ts, Idiv2, ts, y, rac);
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
    g = -(u0J + u1 * int_i);
    if (UNLIKELY(trunc(rac) >= Idiv2_double_minus_one)) {
      for ( ; int_i < Idiv2; int_i++) {
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
      if (++int_i == Idiv2) goto nextj;
      oy = y;
      g += u1;
      y = COMPUTE_Y(g);
    } while (oy != y);
    d0 = 1./d0;
    d1 = rac - d0 * rac;
    y++;
    i = rac + cexp2[(unsigned int)y + 1] * invu1;
    for (;; y++) {
      ts = -int_i;
      if (UNLIKELY(i >= Idiv2_double)) {
	ts += Idiv2;
#ifdef DEBUG_INIT_RAT
	fprintf (stderr, "B4.END : i1=%ld i2=%d, ts=%ld, y=%u, rac=%e\n", Idiv2 - ts, Idiv2, ts, y, rac);
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
  }
}

static inline void
poly_scale_double (double  *u, const double *t, unsigned int d, const double h)
{
  u[d] = t[d];
  for (double hpow = h; d--; hpow *= h) u[d] = t[d] * hpow;
}

#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
static inline void
poly_scale_m128d (__m128d  *u, const double *t, unsigned int d, const double h)
{
  u[d] = _mm_set1_pd (t[d]);
  for (double hpow = h; d--; hpow *= h) u[d] = _mm_set1_pd (t[d] * hpow);
}
#endif


/* Exact initialisation of F(i,j) with degre >= 2 (not mandatory). Slow.
   Internal function, only with simple types, for unit/integration testing. */
 void init_exact_degree_X_norms_bucket_region_internal (unsigned char *S, uint32_t J, uint32_t I, double scale, unsigned int d, double *fijd)
{ 
  unsigned char *beginS;
  uint32_t beginJ, endJ = LOG_BUCKET_REGION - ctz (I);
  ssize_t ih;
  const ssize_t Idiv2 = (ssize_t) I >> 1;
  J <<= endJ;
  beginJ = J;
  endJ = (1U << endJ) + J;
  S += Idiv2;
  beginS = S;
  scale *= 1./0x100000;
  const double add = 0x3FF00000 - GUARD / scale;

#ifndef HAVE_GCC_STYLE_AMD64_INLINE_ASM /* Light optimization, only log2 */
  
  for (; J < endJ; S += I, J++) {
    double g, h, u[d+1];
    poly_scale_double (u, fijd, d, (double) J);
    for (ih = -Idiv2, h = (double) -Idiv2; ih < Idiv2; h += 1., ++ih) {
      g = u[d]; for (size_t k = d; k--; g = g * h + u[k]);
      S[ih] = lg2abs (g, add, scale);
    }
  }

#else /* HAVE_GCC_STYLE_AMD64_INLINE_ASM : optimized part. Stupid fast code. */

  /****** Some SSE internals macros & functions for this function *********/
#ifndef HAVE_SSSE3
#define _mm_alignr_epi8(a, b, c) \
     (_mm_add_epi16 (_mm_slli_si128 ((a), 16 - (c)), _mm_srli_si128 ((b), (c))))
#endif
  
  /* This macro avoids a stupid & boring C types control */
#define _MM_SHUFFLE_EPI32(A,B,C) __asm__ __volatile__ ("pshufd $" #C ", %1, %0\n":"=x"(A):"x"(B))
  
  /* Initialisation of the data with intrinsics. */
#define ALG_CONST_INIT __m128d _two = _mm_set1_pd(2.), \
    _scale = _mm_set1_pd (scale), _add = _mm_set1_pd (add)
  /*************** End of the macros for this function *****************/

  /* This function is highly expensive in SSE computations,
     so it's seriously optimized.
     It's only a big switch in fact with boring code.
     It's important to declare and initialize the
     variables IN each case for optimization reasons. */
  switch (d) {
  case 2 : {
    ALG_CONST_INIT; __m128i cumul; __m128d h, g, u0, u1, u2;
    for (; J < endJ; S += I, J++) {
      h = _mm_set1_pd (J);
      u0 = _mm_mul_pd (_mm_mul_pd (_mm_load1_pd(fijd    ),h),h);
      u1 =             _mm_mul_pd (_mm_load1_pd(fijd + 1),h);
      u2 =                         _mm_load1_pd(fijd + 2);
      for (ih = -Idiv2, h = _mm_set_pd (1 - Idiv2, -Idiv2); ih < Idiv2; *(__m128i *) &S[ih] = cumul, ih += 16) {
	g = _mm_add_pd(_mm_mul_pd (_mm_add_pd(_mm_mul_pd(h,u2),u1),h),u0);
	h = _mm_add_pd (h, _two); cumul = _mm_slli_si128  (_mm_lg2abs (&g, _add, _scale), 14);
	g = _mm_add_pd(_mm_mul_pd (_mm_add_pd(_mm_mul_pd(h,u2),u1),h),u0);
	h = _mm_add_pd (h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd (_mm_add_pd(_mm_mul_pd(h,u2),u1),h),u0);
	h = _mm_add_pd (h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd (_mm_add_pd(_mm_mul_pd(h,u2),u1),h),u0);
	h = _mm_add_pd (h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd (_mm_add_pd(_mm_mul_pd(h,u2),u1),h),u0);
	h = _mm_add_pd (h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd (_mm_add_pd(_mm_mul_pd(h,u2),u1),h),u0);
	h = _mm_add_pd (h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd (_mm_add_pd(_mm_mul_pd(h,u2),u1),h),u0);
	h = _mm_add_pd (h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd (_mm_add_pd(_mm_mul_pd(h,u2),u1),h),u0);
	h = _mm_add_pd (h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
      }}}
    break;
  case 3 : {
    ALG_CONST_INIT; __m128i cumul; __m128d h, g, u0, u1, u2, u3;
    for (; J < endJ; S += I, J++) {
      h = _mm_set_sd (J); _MM_SHUFFLE_EPI32(g, h, 0x44); g = _mm_mul_pd (g, g);
      u0 = _mm_mul_pd (_mm_mul_sd (_mm_load_pd(fijd    ),h),g);
      u2 =             _mm_mul_sd (_mm_load_pd(fijd + 2),h);
      _MM_SHUFFLE_EPI32 (u1, u0, 0xEE); u0 = _mm_unpacklo_pd (u0, u0);
      _MM_SHUFFLE_EPI32 (u3, u2, 0xEE); u2 = _mm_unpacklo_pd (u2, u2);
      for (ih = -Idiv2, h = _mm_set_pd (1 - Idiv2, -Idiv2); ih < Idiv2; *(__m128i *) &S[ih] = cumul, ih += 16) {
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u3),u2),h),u1),h),u0);
	h = _mm_add_pd (h, _two); cumul = _mm_slli_si128  (_mm_lg2abs (&g, _add, _scale), 14);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u3),u2),h),u1),h),u0);
	h = _mm_add_pd (h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u3),u2),h),u1),h),u0);
	h = _mm_add_pd (h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u3),u2),h),u1),h),u0);
	h = _mm_add_pd (h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u3),u2),h),u1),h),u0);
	h = _mm_add_pd (h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u3),u2),h),u1),h),u0);
	h = _mm_add_pd (h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u3),u2),h),u1),h),u0);
	h = _mm_add_pd (h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u3),u2),h),u1),h),u0);
	h = _mm_add_pd (h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
      }}}
    break;
  case 4 : {
    ALG_CONST_INIT; __m128i cumul; __m128d h, g, u0, u1, u2, u3, u4;
    for (; J < endJ; S += I, J++) {
      h = _mm_set1_pd (J); h = _mm_mul_sd (h, h); _MM_SHUFFLE_EPI32(g, h, 0x44);
      u0 = _mm_mul_pd (_mm_mul_pd (_mm_load_pd  (fijd    ),h),g);
      u2 =             _mm_mul_pd (_mm_load_pd  (fijd + 2),h);
      u4 =                         _mm_load1_pd (fijd + 4);
      _MM_SHUFFLE_EPI32 (u1, u0, 0xEE); u0 = _mm_unpacklo_pd (u0, u0);
      _MM_SHUFFLE_EPI32 (u3, u2, 0xEE); u2 = _mm_unpacklo_pd (u2, u2);
      for (ih = -Idiv2, h = _mm_set_pd (1 - Idiv2, -Idiv2); ih < Idiv2; *(__m128i *) &S[ih] = cumul, ih += 16) {
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u4),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd (h, _two); cumul = _mm_slli_si128  (_mm_lg2abs (&g, _add, _scale), 14);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u4),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd (h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u4),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd (h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u4),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd (h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u4),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd (h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u4),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd (h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u4),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd (h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u4),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd (h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
      }}}
    break;
  case 5 : {
    ALG_CONST_INIT; __m128i cumul; __m128d h, g, u0, u1, u2, u3, u4, u5;
    for (; J < endJ; S += I, J++) {
      h = _mm_set_sd (J); _MM_SHUFFLE_EPI32 (g, h, 0x44); g = _mm_mul_pd (g, g);
      u0 = _mm_mul_pd (_mm_mul_pd (_mm_mul_sd (_mm_load_pd (fijd    ),h),g),g);
      u2 =             _mm_mul_pd (_mm_mul_sd (_mm_load_pd (fijd + 2),h),g);
      u4 =                         _mm_mul_sd (_mm_load_pd (fijd + 4),h);
      _MM_SHUFFLE_EPI32 (u1, u0, 0xEE); u0 = _mm_unpacklo_pd (u0, u0);
      _MM_SHUFFLE_EPI32 (u3, u2, 0xEE); u2 = _mm_unpacklo_pd (u2, u2);
      _MM_SHUFFLE_EPI32 (u5, u4, 0xEE); u4 = _mm_unpacklo_pd (u4, u4);
      for (ih = -Idiv2, h = _mm_set_pd (1 - Idiv2, -Idiv2); ih < Idiv2; *(__m128i *) &S[ih] = cumul, ih += 16) {
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u5),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_slli_si128  (_mm_lg2abs (&g, _add, _scale), 14);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u5),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u5),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u5),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u5),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u5),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u5),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u5),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
      }}}
    break;
  case 6 : {
    ALG_CONST_INIT; __m128i cumul; __m128d h, g, u0, u1, u2, u3, u4, u5, u6;
    for (; J < endJ; S += I, J++) {
      h = _mm_set1_pd (J); h = _mm_mul_sd (h, h); _MM_SHUFFLE_EPI32 (g, h, 0x44);
      u0 = _mm_mul_pd (_mm_mul_pd (_mm_mul_pd (_mm_load_pd  (fijd    ),h),g),g);
      u2 =             _mm_mul_pd (_mm_mul_pd (_mm_load_pd  (fijd + 2),h),g);
      u4 =                         _mm_mul_pd (_mm_load_pd  (fijd + 4),h);
      u6 =                                     _mm_load1_pd (fijd + 6);
      _MM_SHUFFLE_EPI32 (u1, u0, 0xEE); u0 = _mm_unpacklo_pd (u0, u0);
      _MM_SHUFFLE_EPI32 (u3, u2, 0xEE); u2 = _mm_unpacklo_pd (u2, u2);
      _MM_SHUFFLE_EPI32 (u5, u4, 0xEE); u4 = _mm_unpacklo_pd (u4, u4);
      for (ih = -Idiv2, h = _mm_set_pd (1 - Idiv2, -Idiv2); ih < Idiv2; *(__m128i *) &S[ih] = cumul, ih += 16) {
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u6),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_slli_si128  (_mm_lg2abs (&g, _add, _scale), 14);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u6),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u6),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u6),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u6),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u6),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u6),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u6),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
      }}}
    break;
  case 7 : {
    ALG_CONST_INIT; __m128i cumul; __m128d h, g, u0, u1, u2, u3, u4, u5, u6, u7;
    for (; J < endJ; S += I, J++) {
      h = _mm_set_sd (J); _MM_SHUFFLE_EPI32 (g, h, 0x44); g = _mm_mul_pd (g, g);
      u0 = _mm_mul_pd (_mm_mul_sd (_mm_load_pd (fijd    ),h),g);
      u2 =             _mm_mul_sd (_mm_load_pd (fijd + 2),h);
      u4 = _mm_mul_pd (_mm_mul_sd (_mm_load_pd (fijd + 4),h),g);
      u6 =             _mm_mul_sd (_mm_load_pd (fijd + 6),h);
      g = _mm_mul_pd (g, g); u2 = _mm_mul_pd (u2, g); u0 = _mm_mul_pd (u0, g);
      _MM_SHUFFLE_EPI32 (u1, u0, 0xEE); u0 = _mm_unpacklo_pd (u0, u0);
      _MM_SHUFFLE_EPI32 (u3, u2, 0xEE); u2 = _mm_unpacklo_pd (u2, u2);
      _MM_SHUFFLE_EPI32 (u5, u4, 0xEE); u4 = _mm_unpacklo_pd (u4, u4);
      _MM_SHUFFLE_EPI32 (u7, u6, 0xEE); u6 = _mm_unpacklo_pd (u6, u6);
      for (ih = -Idiv2, h = _mm_set_pd (1 - Idiv2, -Idiv2); ih < Idiv2; *(__m128i *) &S[ih] = cumul, ih += 16) {
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u7),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_slli_si128  (_mm_lg2abs (&g, _add, _scale), 14);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u7),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u7),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u7),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u7),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u7),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u7),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u7),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
      }}}
    break;
  case 8 : {
    ALG_CONST_INIT; __m128i cumul; __m128d h, g, u0, u1, u2, u3, u4, u5, u6, u7, u8;
    for (; J < endJ; S += I, J++) {
      h = _mm_set1_pd (J); h = _mm_mul_sd (h, h); _MM_SHUFFLE_EPI32 (g, h, 0x44);
      u0 = _mm_mul_pd (_mm_mul_pd (_mm_load_pd  (fijd    ),h),g);
      u2 =             _mm_mul_pd (_mm_load_pd  (fijd + 2),h);
      u4 = _mm_mul_pd (_mm_mul_pd (_mm_load_pd  (fijd + 4),h),g);
      u6 =             _mm_mul_pd (_mm_load_pd  (fijd + 6),h);
      u8 =                         _mm_load1_pd (fijd + 8);
      g = _mm_mul_pd (g, g); u2 = _mm_mul_pd (u2, g); u0 = _mm_mul_pd (u0, g);
      _MM_SHUFFLE_EPI32 (u1, u0, 0xEE); u0 = _mm_unpacklo_pd (u0, u0);
      _MM_SHUFFLE_EPI32 (u3, u2, 0xEE); u2 = _mm_unpacklo_pd (u2, u2);
      _MM_SHUFFLE_EPI32 (u5, u4, 0xEE); u4 = _mm_unpacklo_pd (u4, u4);
      _MM_SHUFFLE_EPI32 (u7, u6, 0xEE); u6 = _mm_unpacklo_pd (u6, u6);
      for (ih = -Idiv2, h = _mm_set_pd (1 - Idiv2, -Idiv2); ih < Idiv2; *(__m128i *) &S[ih] = cumul, ih += 16) {
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u8),u7),h),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_slli_si128  (_mm_lg2abs (&g, _add, _scale), 14);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u8),u7),h),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u8),u7),h),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u8),u7),h),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u8),u7),h),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u8),u7),h),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u8),u7),h),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u8),u7),h),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
      }}}
    break;
  case 9 : {
    ALG_CONST_INIT; __m128i cumul; __m128d h, g, u0, u1, u2, u3, u4, u5, u6, u7, u8, u9;
    for (; J < endJ; S += I, J++) {
      h = _mm_set_sd (J); _MM_SHUFFLE_EPI32 (g, h, 0x44); g = _mm_mul_pd (g, g);
      u0 =             _mm_mul_sd (_mm_load_pd (fijd    ),h);
      u2 = _mm_mul_pd (_mm_mul_sd (_mm_load_pd (fijd + 2),h),g);
      u4 =             _mm_mul_sd (_mm_load_pd (fijd + 4),h);
      u6 = _mm_mul_pd (_mm_mul_sd (_mm_load_pd (fijd + 6),h),g);
      u8 =             _mm_mul_sd (_mm_load_pd (fijd + 8),h);
      g = _mm_mul_pd (g, g); u4 = _mm_mul_pd (u4, g); u2 = _mm_mul_pd (u2, g);
      g = _mm_mul_pd (g, g); u0 = _mm_mul_pd (u0, g);
      _MM_SHUFFLE_EPI32 (u1, u0, 0xEE); u0 = _mm_unpacklo_pd (u0, u0);
      _MM_SHUFFLE_EPI32 (u3, u2, 0xEE); u2 = _mm_unpacklo_pd (u2, u2);
      _MM_SHUFFLE_EPI32 (u5, u4, 0xEE); u4 = _mm_unpacklo_pd (u4, u4);
      _MM_SHUFFLE_EPI32 (u7, u6, 0xEE); u6 = _mm_unpacklo_pd (u6, u6);
      _MM_SHUFFLE_EPI32 (u9, u8, 0xEE); u8 = _mm_unpacklo_pd (u8, u8);
      for (ih = -Idiv2, h = _mm_set_pd (1 - Idiv2, -Idiv2); ih < Idiv2; *(__m128i *) &S[ih] = cumul, ih += 16) {
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u9),u8),h),u7),h),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_slli_si128  (_mm_lg2abs (&g, _add, _scale), 14);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u9),u8),h),u7),h),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u9),u8),h),u7),h),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u9),u8),h),u7),h),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u9),u8),h),u7),h),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u9),u8),h),u7),h),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u9),u8),h),u7),h),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u9),u8),h),u7),h),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
      }}}
    break;
  default: { /* This default could remplace all (optimized) previous cases. */
    ALG_CONST_INIT; __m128i cumul; __m128d h, g, u[d+1];
    for (; J < endJ; S += I, J++) {
      poly_scale_m128d (u, fijd, d, (double) J);
      for (ih = -Idiv2, h = _mm_set_pd (1 - Idiv2, -Idiv2); ih < Idiv2; *(__m128i *) &S[ih] = cumul, ih += 16) {
	g = u[d]; for (size_t k = d; k--; g = _mm_add_pd(_mm_mul_pd(g,h),u[k]));
	h = _mm_add_pd(h, _two); cumul = _mm_slli_si128  (_mm_lg2abs (&g, _add, _scale), 14);
	g = u[d]; for (size_t k = d; k--; g = _mm_add_pd(_mm_mul_pd(g,h),u[k]));
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = u[d]; for (size_t k = d; k--; g = _mm_add_pd(_mm_mul_pd(g,h),u[k]));
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = u[d]; for (size_t k = d; k--; g = _mm_add_pd(_mm_mul_pd(g,h),u[k]));
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = u[d]; for (size_t k = d; k--; g = _mm_add_pd(_mm_mul_pd(g,h),u[k]));
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = u[d]; for (size_t k = d; k--; g = _mm_add_pd(_mm_mul_pd(g,h),u[k]));
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = u[d]; for (size_t k = d; k--; g = _mm_add_pd(_mm_mul_pd(g,h),u[k]));
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
	g = u[d]; for (size_t k = d; k--; g = _mm_add_pd(_mm_mul_pd(g,h),u[k]));
	h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2abs (&g, _add, _scale), cumul, 2);
      }}}
    break;
  } /* End of the big switch of the SSE version - and end of the function */
#undef _MM_SHUFFLE_EPI32
#undef ALG_CONST_INIT
#endif /* End of HAVE_GCC_STYLE_AMD64_INLINE_ASM */
  /* Special ultra rare case. The correction of log2(F(0,0)) is false, because
     the fast algorithm of log2 is not good in 0.0. */
  if (UNLIKELY(!beginJ)) *beginS = GUARD;
}

/* Smart initialization. Computes only one initialization value F(i+3.5,const j)
   for a bloc of 8 horizontal values F(i,const j)... F(i+7,const j). After this,
   the initialization is corrected around the possible roots of
   d^2(log2(F(i,const j)))/d(i)^2,
   i.e. around the inflexion points of log2(F(i, const j)), and around i = 0.
   NB: in fact, the function is log2(fabs(F(i, const j)))... Don't care ?
   2 versions: SSE version (32bits compatible) & non SSE versions.
   Internal function, only with simple types, for unit/integration testing */
void init_smart_degree_X_norms_bucket_region_internal (unsigned char *S, uint32_t J, uint32_t I, double scale, unsigned int d, double *fijd, unsigned int nroots, double *roots)
{ 
  ASSERT (d >= 2);
  unsigned char *beginS;
  uint32_t beginJ, endJ = LOG_BUCKET_REGION - ctz (I);
  ssize_t ih;
  const ssize_t Idiv2 = (ssize_t) I >> 1;
  J <<= endJ;
  beginJ = J;
  endJ = (1U << endJ) + J;
  S += Idiv2;
  beginS = S;
  scale *= (1./0x100000);
  const double add = ((double) 0x3FF00000) - GUARD / scale;

#ifndef HAVE_GCC_STYLE_AMD64_INLINE_ASM

  for ( ; J < endJ; S += I, ++J) {
    double g, h, u[d];
    poly_scale_double (u, fijd, d, (double) J);
    h = (double) (3.5 - Idiv2);
    for (ih = -Idiv2; ih < Idiv2; ih += 8, h += 8.) {
      g = u[d]; for (size_t k = d; k--; g = g * h + u[k]);
      w64lg2abs (g, add, scale, S, ih);
    }
  }

#else /* HAVE_GCC_STYLE_AMD64_INLINE_ASM, optimized part */

  /*************** Some SSE macros for this function ***********************/
  /* This macro avoids a stupid & boring C types control */
#define _MM_SHUFFLE_EPI32(A,B,C) __asm__ __volatile__ ("pshufd $" #C ", %1, %0\n":"=x"(A):"x"(B))

  /* Initialisation of the data with intrinsics. */
#define ALG_CONST_INIT const __m128d \
	_scale = _mm_set1_pd (scale), _add = _mm_set1_pd(add), _sixteen = _mm_set1_pd(16.)
  /*************** End of the macros for this function *****************/
  
  /* This function is only a big switch in fact */
  switch (d) { /* d >= 2 */
  case 2 : {
    ALG_CONST_INIT; __m128d h, g, u0, u1, u2;
    for ( ; J < endJ; S += I, ++J) {
      h = _mm_set1_pd((double) J);
      u0 = _mm_mul_pd (_mm_mul_pd (_mm_load1_pd(fijd    ),h),h);
      u1 =             _mm_mul_pd (_mm_load1_pd(fijd + 1),h);
      u2 =                         _mm_load1_pd(fijd + 2);
      for (ih = -Idiv2, h = _mm_set_pd (11.5 - Idiv2, 3.5 - Idiv2); ih < Idiv2; ih += 16, h = _mm_add_pd(h, _sixteen)) {
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u2),u1),h),u0);
	w128lg2abs (g, _add, _scale, S, ih);
      }}}	
    break;
  case 3 : {
    ALG_CONST_INIT; __m128d h, g, u0, u1, u2, u3;
    for ( ; J < endJ; S += I, ++J) {
      h = _mm_set_sd((double) J); _MM_SHUFFLE_EPI32(g, h, 0x44); g = _mm_mul_pd (g, g);
      u0 = _mm_mul_pd (_mm_mul_sd (_mm_load_pd(fijd    ),h),g);
      u2 =             _mm_mul_sd (_mm_load_pd(fijd + 2),h);
      _MM_SHUFFLE_EPI32 (u1, u0, 0xEE); u0 = _mm_unpacklo_pd (u0, u0);
      _MM_SHUFFLE_EPI32 (u3, u2, 0xEE); u2 = _mm_unpacklo_pd (u2, u2);
      for (ih = -Idiv2, h = _mm_set_pd (11.5 - Idiv2, 3.5 - Idiv2); ih < Idiv2; ih += 16, h = _mm_add_pd(h, _sixteen)) {
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u3),u2),h),u1),h),u0);
	w128lg2abs (g, _add, _scale, S, ih);
      }}}
    break;
  case 4 : {
    ALG_CONST_INIT; __m128d h, g, u0, u1, u2, u3, u4;
    for ( ; J < endJ; S += I, ++J) {
      h = _mm_set1_pd((double) J); h = _mm_mul_sd (h, h); _MM_SHUFFLE_EPI32(g, h, 0x44);
      u0 = _mm_mul_pd (_mm_mul_pd (_mm_load_pd  (fijd    ),h),g);
      u2 =             _mm_mul_pd (_mm_load_pd  (fijd + 2),h);
      u4 =                         _mm_load1_pd (fijd + 4);
      _MM_SHUFFLE_EPI32 (u1, u0, 0xEE); u0 = _mm_unpacklo_pd (u0, u0);
      _MM_SHUFFLE_EPI32 (u3, u2, 0xEE); u2 = _mm_unpacklo_pd (u2, u2);
      for (ih = -Idiv2, h = _mm_set_pd (11.5 - Idiv2, 3.5 - Idiv2); ih < Idiv2; ih += 16, h = _mm_add_pd(h, _sixteen)) {
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u4),u3),h),u2),h),u1),h),u0);
	w128lg2abs (g, _add, _scale, S, ih);
      }}}
    break;
  case 5 : {
    ALG_CONST_INIT; __m128d h, g, u0, u1, u2, u3, u4, u5;
    for ( ; J < endJ; S += I, ++J) {
      h = _mm_set_sd((double) J); _MM_SHUFFLE_EPI32 (g, h, 0x44); g = _mm_mul_pd (g, g);
      u0 = _mm_mul_pd (_mm_mul_pd (_mm_mul_sd (_mm_load_pd (fijd    ),h),g),g);
      u2 =             _mm_mul_pd (_mm_mul_sd (_mm_load_pd (fijd + 2),h),g);
      u4 =                         _mm_mul_sd (_mm_load_pd (fijd + 4),h);
      _MM_SHUFFLE_EPI32 (u1, u0, 0xEE); u0 = _mm_unpacklo_pd (u0, u0);
      _MM_SHUFFLE_EPI32 (u3, u2, 0xEE); u2 = _mm_unpacklo_pd (u2, u2);
      _MM_SHUFFLE_EPI32 (u5, u4, 0xEE); u4 = _mm_unpacklo_pd (u4, u4);
      for (ih = -Idiv2, h = _mm_set_pd (11.5 - Idiv2, 3.5 - Idiv2); ih < Idiv2; ih += 16, h = _mm_add_pd(h, _sixteen)) {
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u5),u4),h),u3),h),u2),h),u1),h),u0);
	w128lg2abs (g, _add, _scale, S, ih);
      }}}
    break;
  case 6 : {
    ALG_CONST_INIT; __m128d h, g, u0, u1, u2, u3, u4, u5, u6;
    for ( ; J < endJ; S += I, ++J) {
      h = _mm_set1_pd((double) J); h = _mm_mul_sd (h, h); _MM_SHUFFLE_EPI32 (g, h, 0x44);
      u0 = _mm_mul_pd (_mm_mul_pd (_mm_mul_pd (_mm_load_pd  (fijd    ),h),g),g);
      u2 =             _mm_mul_pd (_mm_mul_pd (_mm_load_pd  (fijd + 2),h),g);
      u4 =                         _mm_mul_pd (_mm_load_pd  (fijd + 4),h);
      u6 =                                     _mm_load1_pd (fijd + 6);
      _MM_SHUFFLE_EPI32 (u1, u0, 0xEE); u0 = _mm_unpacklo_pd (u0, u0);
      _MM_SHUFFLE_EPI32 (u3, u2, 0xEE); u2 = _mm_unpacklo_pd (u2, u2);
      _MM_SHUFFLE_EPI32 (u5, u4, 0xEE); u4 = _mm_unpacklo_pd (u4, u4);
      for (ih = -Idiv2, h = _mm_set_pd (11.5 - Idiv2, 3.5 - Idiv2); ih < Idiv2; ih += 16, h = _mm_add_pd(h, _sixteen)) {
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u6),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	w128lg2abs (g, _add, _scale, S, ih);
      }}}
    break;
  case 7 : {
    ALG_CONST_INIT; __m128d h, g, u0, u1, u2, u3, u4, u5, u6, u7;
    for ( ; J < endJ; S += I, ++J) {
      h = _mm_set_sd((double) J); _MM_SHUFFLE_EPI32 (g, h, 0x44); g = _mm_mul_pd (g, g);
      u0 = _mm_mul_pd (_mm_mul_sd (_mm_load_pd (fijd    ),h),g);
      u2 =             _mm_mul_sd (_mm_load_pd (fijd + 2),h);
      u4 = _mm_mul_pd (_mm_mul_sd (_mm_load_pd (fijd + 4),h),g);
      u6 =             _mm_mul_sd (_mm_load_pd (fijd + 6),h);
      g = _mm_mul_pd (g, g); u2 = _mm_mul_pd (u2, g); u0 = _mm_mul_pd (u0, g);
      _MM_SHUFFLE_EPI32 (u1, u0, 0xEE); u0 = _mm_unpacklo_pd (u0, u0);
      _MM_SHUFFLE_EPI32 (u3, u2, 0xEE); u2 = _mm_unpacklo_pd (u2, u2);
      _MM_SHUFFLE_EPI32 (u5, u4, 0xEE); u4 = _mm_unpacklo_pd (u4, u4);
      _MM_SHUFFLE_EPI32 (u7, u6, 0xEE); u6 = _mm_unpacklo_pd (u6, u6);
      for (ih = -Idiv2, h = _mm_set_pd (11.5 - Idiv2, 3.5 - Idiv2); ih < Idiv2; ih += 16, h = _mm_add_pd(h, _sixteen)) {
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u7),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	w128lg2abs (g, _add, _scale, S, ih);
      }}}
    break;
  case 8 : {
    ALG_CONST_INIT; __m128d h, g, u0, u1, u2, u3, u4, u5, u6, u7, u8;
    for ( ; J < endJ; S += I, ++J) {
      h = _mm_set1_pd((double) J); h = _mm_mul_sd (h, h); _MM_SHUFFLE_EPI32 (g, h, 0x44);
      u0 = _mm_mul_pd (_mm_mul_pd (_mm_load_pd  (fijd    ),h),g);
      u2 =             _mm_mul_pd (_mm_load_pd  (fijd + 2),h);
      u4 = _mm_mul_pd (_mm_mul_pd (_mm_load_pd  (fijd + 4),h),g);
      u6 =             _mm_mul_pd (_mm_load_pd  (fijd + 6),h);
      u8 =                         _mm_load1_pd (fijd + 8);
      g = _mm_mul_pd (g, g); u2 = _mm_mul_pd (u2, g); u0 = _mm_mul_pd (u0, g);
      _MM_SHUFFLE_EPI32 (u1, u0, 0xEE); u0 = _mm_unpacklo_pd (u0, u0);
      _MM_SHUFFLE_EPI32 (u3, u2, 0xEE); u2 = _mm_unpacklo_pd (u2, u2);
      _MM_SHUFFLE_EPI32 (u5, u4, 0xEE); u4 = _mm_unpacklo_pd (u4, u4);
      _MM_SHUFFLE_EPI32 (u7, u6, 0xEE); u6 = _mm_unpacklo_pd (u6, u6);
      for (ih = -Idiv2, h = _mm_set_pd (11.5 - Idiv2, 3.5 - Idiv2); ih < Idiv2; ih += 16, h = _mm_add_pd(h, _sixteen)) {
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u8),u7),h),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	w128lg2abs (g, _add, _scale, S, ih);
      }}}
    break;
  case 9 : {
    ALG_CONST_INIT; __m128d h, g, u0, u1, u2, u3, u4, u5, u6, u7, u8, u9;
    for ( ; J < endJ; S += I, ++J) {
      h = _mm_set_sd((double) J); _MM_SHUFFLE_EPI32 (g, h, 0x44); g = _mm_mul_pd (g, g);
      u0 =             _mm_mul_sd (_mm_load_pd (fijd    ),h);
      u2 = _mm_mul_pd (_mm_mul_sd (_mm_load_pd (fijd + 2),h),g);
      u4 =             _mm_mul_sd (_mm_load_pd (fijd + 4),h);
      u6 = _mm_mul_pd (_mm_mul_sd (_mm_load_pd (fijd + 6),h),g);
      u8 =             _mm_mul_sd (_mm_load_pd (fijd + 8),h);
      g = _mm_mul_pd (g, g); u4 = _mm_mul_pd (u4, g); u2 = _mm_mul_pd (u2, g);
      g = _mm_mul_pd (g, g); u0 = _mm_mul_pd (u0, g);
      _MM_SHUFFLE_EPI32 (u1, u0, 0xEE); u0 = _mm_unpacklo_pd (u0, u0);
      _MM_SHUFFLE_EPI32 (u3, u2, 0xEE); u2 = _mm_unpacklo_pd (u2, u2);
      _MM_SHUFFLE_EPI32 (u5, u4, 0xEE); u4 = _mm_unpacklo_pd (u4, u4);
      _MM_SHUFFLE_EPI32 (u7, u6, 0xEE); u6 = _mm_unpacklo_pd (u6, u6);
      _MM_SHUFFLE_EPI32 (u9, u8, 0xEE); u8 = _mm_unpacklo_pd (u8, u8);
      for (ih = -Idiv2, h = _mm_set_pd (11.5 - Idiv2, 3.5 - Idiv2); ih < Idiv2; ih += 16, h = _mm_add_pd(h, _sixteen)) {
	g = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_mul_pd(h,u9),u8),h),u7),h),u6),h),u5),h),u4),h),u3),h),u2),h),u1),h),u0);
	w128lg2abs (g, _add, _scale, S, ih);
      }}}
    break;
  default: { /* This default could remplace all (optimized) previous cases. */
    ALG_CONST_INIT; __m128d h, g, u[d+1];
    for ( ; J < endJ; S += I, ++J) {
      poly_scale_m128d (u, fijd, d, (double) J);
      for (ih = -Idiv2, h = _mm_set_pd (11.5 - Idiv2, 3.5 - Idiv2); ih < Idiv2; ih += 16, h = _mm_add_pd(h, _sixteen)) {
	g = u[d]; for (size_t k = d; k--; g = _mm_add_pd(_mm_mul_pd(g,h),u[k]));
	w128lg2abs (g, _add, _scale, S, ih);
      }}}
    break;
  } /* End of the big switch of the SSE version - and end of the function */
#undef _MM_SHUFFLE_EPI32
#undef ALG_INIT_SCALE_ADD_ONE_SIXTEEN

#endif /* !HAVE_GCC_STYLE_AMD64_INLINE_ASM */
#undef FILL_OTHER_LINES

  /* Now, the value near the roots must be corrected. */
  J = beginJ;
  S = beginS;
  double dIdiv2 = (double) Idiv2;
  for (; J < endJ; S += I, J++) {
    double u[d+1];
    poly_scale_double (u, fijd, d, (double) J);
    for (size_t r = 0; r < nroots; r++) {
      double hl = floor((double)J * roots[r]), hr = hl + 1.;
      ssize_t ih;
      if (hl >= -dIdiv2) {
	unsigned char c1, c2 = 0;
	hl = MIN(hl, dIdiv2 - 1.);
	ih = (ssize_t) hl; 
	do {
	  double g = u[d]; for (size_t k = d; k--; g = g * hl + u[k]);
	  c1 = lg2abs (g,add,scale);
	  --hl;
	  if (UNLIKELY(abs (c1 - S[ih]))) {
	    S[ih] = c1;
	    c2 = 0;
	  } else if (++c2 > 4) break;
	} while (--ih >= -Idiv2);
      }
      if (hr < dIdiv2) {
	unsigned char c1, c2 = 0;
	hr = MAX(hr, -dIdiv2);
	ih = (ssize_t) hr;
	do {
	  double g = u[d]; for (size_t k = d; k--; g = g * hr + u[k]);
	  c1 = lg2abs (g,add,scale);
	  ++hr;
	  if (UNLIKELY(abs (c1 - S[ih]))) {
	    S[ih] = c1;
	    c2 = 0;
	  } else if (++c2 > 4) break;
	} while (++ih < Idiv2);
      }
    }
  }
  /* Special ultra rare case. The correction of log2(F(0,0)) is false, because
     the fast algorithm of log2 is not good in 0.0. */
  if (UNLIKELY(!beginJ)) *beginS = GUARD;
} /* True end of the function (finally!) */

/* This function is used to initialize lognorms (=log2(F(i,j)*scale+GUARD)
   for the bucket_region S[] number J.
   It's a wrapper; except for trivial degree = 0, it extracts the interesting
   parameters of the complex structure si and calls the right function.

   - For degree 0, S[] is initialized by a memset: always exact.
   - A special ultra fast init function is used for degree = 1; it's could be
     considered as exact (the maximal error is always -/+ 1 on S[]).
   - For smart = 0 and others degrees, the exact values F(i,j) are computed with
     a fast log2. The maximal error is always -/+ 1 on S[]. It's obvious slow.

   - For smart != 0 and others degrees, only 1 value each 8 values is computed.
   The value of F(i+3.5,const j) is used to initialize F(i..i+7,const j).
   Then around the neigborhood of the inflexion points of F(i,j), the roots of 
   d^2(F(i,1)/d(i)^2 =  F(i,1)*F"(i,1) - F'^2(i,1), the computed values are
   corrected. It's 4 to 5 times faster than the exact init.

   This smart algo needs :
   -> The roots of d^2(F(i,1)/d(i)^2 must be in si->sides[side]->roots;
   -> si->sides[side]->nroots is the number of roots + 1;
   -> si->roots[si->sides[side]->nroots - 1] must be = 0.0 : it's a "pseudo" root
      in order to correct the neigborhood of F(0, const j).
   Of course, if si->sides[side]->nroots = 0, no correction is done.
*/
void init_norms_bucket_region (unsigned char *S, uint32_t J, sieve_info_ptr si, unsigned int side, unsigned int smart)
{
  unsigned int degree = si->sides[side]->fij->deg;
  switch (degree) {
  case 0 :
    memset (S, (int) (log2(1.+fabs(si->sides[side]->fijd[0])) * si->sides[side]->scale) + GUARD, 1U << LOG_BUCKET_REGION);
    break;
  case 1 :
    init_degree_one_norms_bucket_region_internal (S, J, si->I, si->sides[side]->scale, si->sides[side]->fijd[0], si->sides[side]->fijd[1], si->sides[side]->cexp2);
    break;
  default:
    if (smart)
      init_smart_degree_X_norms_bucket_region_internal (S, J, si->I, si->sides[side]->scale, degree, si->sides[side]->fijd, si->sides[side]->nroots, si->sides[side]->roots);
    else
      init_exact_degree_X_norms_bucket_region_internal (S, J, si->I, si->sides[side]->scale, degree, si->sides[side]->fijd);
    break;
  }
}

/* return max |g(x)| for x in (0, s) where s can be negative,
   and g(x) = g[d]*x^d + ... + g[1]*x + g[0] */
static double
get_maxnorm_aux (double_poly_srcptr poly, double s)
{
  double_poly_t deriv;
  const int d = poly->deg;

  ASSERT_ALWAYS(d >= 0);

  if (d == 0)
    return fabs (poly->coeff[0]);

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
double
get_maxnorm_alg (double_poly_srcptr src_poly, const double X, const double Y)
{
  const unsigned int d = src_poly->deg;
  double norm, max_norm;

  /* Make copy of polynomial as we need to revert the coefficients */
  double_poly_t poly;
  double_poly_init (poly, d);
  double_poly_set (poly, src_poly);

  /* (b) determine the maximum of |f(x)| * Y^d for -X/Y <= x <= X/Y */
  max_norm = get_maxnorm_aux_pm (poly, X/Y) * pow(Y, (double)d);

  /* (a) determine the maximum of |g(y)| for -1 <= y <= 1, with g(y) = F(s,y) */
  double_poly_revert(poly);
  norm = get_maxnorm_aux_pm (poly, Y/X) * pow(X, (double)d);
  if (norm > max_norm)
    max_norm = norm;

  double_poly_clear(poly);

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
      /* Cf init_norms_roots to see the reason of d*3+1 */
      si->sides[side]->roots = (double *) malloc_aligned((d * 3 + 1) * sizeof(double), 16);
      FATAL_ERROR_CHECK(si->sides[side]->roots == NULL, "malloc failed");
      si->sides[side]->nroots = 0;
    }
}

void sieve_info_clear_norm_data(sieve_info_ptr si)
{
    for(int side = 0 ; side < 2 ; side++) {
        sieve_side_info_ptr s = si->sides[side];
        mpz_poly_clear (s->fij);
        free(s->fijd);
        free(s->roots);
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
  region is no larger than this optimal maximum, times some constant fudge
  factor.
*/
static unsigned int
sieve_info_update_norm_data_Jmax (sieve_info_ptr si)
{
  // The following parameter controls the scaling on the norm.
  // Relevant values are between 1.0 and 3.0. A higher value means we
  // select higher values of J, and therefore we find more relations, but
  // this increases the time per relations.
  // The value 2.0 seems to be a good compromise. Setting 1.5 reduces the
  // time per relation and number of relations by about 1.5% on a typical
  // RSA704 benchmark.
  const double fudge_factor = 2.0; 
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
      double maxnorm = get_maxnorm_alg (dpoly, fudge_factor*A/2.,
              fudge_factor*B);
      double_poly_clear (dpoly);
      if (side == si->doing->side)
        maxnorm /= q;

      double_poly_t F;
      F->deg = ps->deg;
      F->coeff = s->fijd;
      
      double v = get_maxnorm_alg (F, I/2, Jmax);
      
      if (v > maxnorm)
        { /* use dichotomy to determine largest Jmax */
          double a, b, c;
          a = 0.0;
          b = Jmax;
          while (trunc (a) != trunc (b))
            {
              c = (a + b) * 0.5;
              v = get_maxnorm_alg (F, I/2, c);
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


/* return 0 if we should discard that special-q, in which case we intend
 * to discard this special-q. For this reason, si->J is then set to an
 * unrounded value, for diagnostic.
 *
 * The current check for discarding is whether we do fill one bucket or
 * not. If we don't even achieve that, we should of course discard.
 *
 * Now for efficiency reasons, the ``minimum reasonable'' number of
 * buckets should be more than that.
 */
int sieve_info_adjust_IJ(sieve_info_ptr si, int nb_threads)/*{{{*/
{
    /* compare skewed max-norms: let u = [a0, b0] and v = [a1, b1],
       and u' = [a0, s*b0], v' = [a1, s*b1] where s is the skewness.
       Assume |u'| <= |v'|.
       We know from Gauss reduction that u' and v' form an angle of at
       least pi/3, thus since their determinant is q*s, we have
       |u'|^2 <= |u'| * |v'| <= q*s/sin(pi/3) = 2*q*s/sqrt(3).
       Define B := sqrt(2/sqrt(3)*q/s), then |a0| <= s*B and |b0| <= B.

       If we take J <= I/2*min(s*B/|a1|, B/|b1|), then for any j <= J:
       |a1|*J <= I/2*s*B and |b1|*J <= I/2*B, thus
       |a| = |a0*i+a1*j| <= s*B*I and |b| <= |b0*i+b1*j| <= B*I.
    */
    const double skewness = si->conf->skewness;
    const int verbose = 0;
    if (verbose) {
        printf("# Called sieve_info_adjust_IJ((a0=%" PRId64 "; b0=%" PRId64
               "; a1=%" PRId64 "; b1=%" PRId64 "), p=%lu, skew=%f, nb_threads=%d)\n",
               si->a0, si->b0, si->a1, si->b1, mpz_get_ui(si->doing->p), skewness, nb_threads);
    }
    double maxab1, maxab0;
    maxab1 = si->b1 * skewness;
    maxab1 = maxab1 * maxab1 + si->a1 * si->a1;
    maxab0 = si->b0 * skewness;
    maxab0 = maxab0 * maxab0 + si->a0 * si->a0;
    if (maxab0 > maxab1) { /* exchange u and v, thus I and J */
        int64_t oa[2] = { si->a0, si->a1 };
        int64_t ob[2] = { si->b0, si->b1 };
        si->a0 = oa[1]; si->a1 = oa[0];
        si->b0 = ob[1]; si->b1 = ob[0];
    }
    maxab1 = MAX(fabs(si->a1), fabs(si->b1) * skewness);
    /* make sure J does not exceed I/2 */
    /* FIXME: We should not have to compute this B a second time. It
     * appears in sieve_info_init_norm_data already */
    double B = sqrt (2.0 * mpz_get_d(si->doing->p) / (skewness * sqrt (3.0)));
    if (maxab1 >= B * skewness)
        si->J = (uint32_t) (B * skewness / maxab1 * (double) (si->I >> 1));
    else
        si->J = si->I >> 1;

    /* Make sure the bucket region size divides the sieve region size, 
       partly covered bucket regions may lead to problems when 
       reconstructing p from half-empty buckets. */
    /* Compute number of i-lines per bucket region, must be integer */
    ASSERT_ALWAYS(LOG_BUCKET_REGION >= si->conf->logI);
    uint32_t i = 1U << (LOG_BUCKET_REGION - si->conf->logI);
    i *= nb_threads;  /* ensures nb of bucket regions divisible by nb_threads */

    /* Bug 15617: if we round up, we are not true to our promises */
    uint32_t nJ = (si->J / i) * i; /* Round down to multiple of i */

    if (verbose) printf("# %s(): Final J=%" PRIu32 "\n", __func__, nJ);
    /* XXX No rounding if we intend to abort */
    if (nJ > 0) si->J = nJ;
    return nJ > 0;
}/*}}}*/


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

  /* Compute the roots of the polynome F(i,1) and the roots of its inflexion points
     d^2(F(i,1))/d(i)^2.
     These roots and 0.0 are need to be corrected the norms initialization on their
     neighboroods in init_smart_degree_X_norms_bucket_region_internal.  */
#ifdef SMART_NORM
  if (si->cpoly->pols[RATIONAL_SIDE]->deg >= 2) init_norms_roots (si, RATIONAL_SIDE);
#endif

  /* Compute the maximum norm of the rational polynomial over the sieve
     region. The polynomial coefficient in fijd are already divided by q
     on the special-q side. */
  poly->deg = si->cpoly->pols[RATIONAL_SIDE]->deg;
  poly->coeff = si->sides[RATIONAL_SIDE]->fijd;
  rat->logmax = log2(get_maxnorm_alg (poly, (double)si->I/2, (double)si->I/2));

  /* we increase artificially 'logmax', to allow larger values of J */
  rat->logmax += 2.0;

  /* we know that |G(a,b)| < 2^(rat->logmax) when si->ratq = 0,
     and |G(a,b)/q| < 2^(rat->logmax) when si->ratq <> 0 */

  maxlog2 = rat->logmax;
  if (output != NULL)
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
  if (output != NULL)
    fprintf (output, " bound=%u\n", rat->bound);
  double max_rlambda = (maxlog2 - GUARD / rat->scale) /
      si->conf->sides[RATIONAL_SIDE]->lpb;
  if (si->conf->sides[RATIONAL_SIDE]->lambda > max_rlambda && output != NULL)
    fprintf (output, "# Warning, rlambda>%.1f does not make sense (capped to limit)\n", max_rlambda);

  /************************** algebraic side *********************************/

  /* Compute the roots of the polynome F(i,1) and the roots of its inflexion points
     d^2(F(i,1))/d(i)^2.
     These roots and 0.0 are need to be corrected the norms initialization on their
     neighboroods in init_smart_degree_X_norms_bucket_region_internal.  */
#ifdef SMART_NORM
  if (si->cpoly->pols[ALGEBRAIC_SIDE]->deg >= 2) init_norms_roots (si, ALGEBRAIC_SIDE);
#endif
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
  maxlog2 = alg->logmax;

  if (output != NULL)
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
  r = MIN(si->conf->sides[ALGEBRAIC_SIDE]->lambda * (double) si->conf->sides[ALGEBRAIC_SIDE]->lpb, maxlog2 - GUARD / alg->scale);
  alg->bound = (unsigned char) (r * alg->scale + GUARD);
  if (output != NULL)
    fprintf (output, " bound=%u\n", alg->bound);
  double max_alambda = (maxlog2 - GUARD / alg->scale) /
      si->conf->sides[ALGEBRAIC_SIDE]->lpb;
  if (si->conf->sides[ALGEBRAIC_SIDE]->lambda > max_alambda && output != NULL)
    fprintf (output, "# Warning, alambda>%.1f does not make sense (capped to limit)\n", max_alambda);

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
