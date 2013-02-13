#define __STDC_LIMIT_MACROS
#include "cado.h"
#include <string.h>
#include <limits.h>
#include <math.h>               /* ceil */
#ifdef SSE_NORM_INIT
#include <emmintrin.h>
#endif

#include "las-config.h"
#include "las-debug.h"
#include "las-norms.h"
#include "utils.h"
#include "mpz_poly.h"
#include "portability.h"

/************************** sieve info stuff *********************************/

/* initialize array C[0..255]: C[i] is zero whenever the log-norm i
   is considered as a potential report, and 127 otherwise.
   The large prime bound is L = 2^l.
*/
static void
sieve_info_init_lognorm (unsigned char *C, unsigned char threshold,
                         unsigned long B MAYBE_UNUSED,
                         unsigned long l MAYBE_UNUSED,
                         double scale MAYBE_UNUSED)
{
  /* for (k = 0; k < 256; k++) C[k] = (k <= threshold) ? 0 : 127; */
  memset (C, 0, threshold + 1);
  memset (C + threshold + 1, 127, 256 - (threshold + 1));

#ifdef COFACTOR_TRICK
  {
    unsigned char k0, k1;
    double lost = 6.0; /* maximal number of bits lost due to prime powers */

    /* for L < R < B^2, a cofactor R cannot be L-smooth, since then it should
       have at least two prime factors in [B,L], and then should be >= B^2.
       Note:
       - the lognorms S[] are GUARD larger than the real values, thus the
         interval [L,R^2] corresponds to [k0+GUARD,k1+GUARD], but to take into
         account the roundoff errors we consider [k0+2*GUARD,k1]
       - 'lost' takes into account the error because we do not sieve
         prime powers; it applies only to the lower bound, and has to be scaled
       - the additional guard value takes into account additional errors
    */
    k0 = (unsigned char) (((double) l + lost) * scale) + 2 * GUARD + 2;
    k1 = (unsigned char) (2.0 * log2 ((double) B) * scale) - 2;
    for (k = k0; k <= k1 && k < 256; k++)
      C[k] = 127;
  }
#endif
}

MAYBE_UNUSED static void
sieve_info_init_lognorm_prob (unsigned char *C, const unsigned long *rels,
    const unsigned long *surv, const double thresh)
{
  const double scale = -126./log(thresh);

  /* Success probability up to thresh get values 0 (probability=1) 
     through 126 (probability=thresh). Probabilities worse than thresh 
     all get 127. */

  for (size_t i = 0; i < 256; i++) {
    if ((double) rels[i] < thresh * (double) surv[i]) {
        C[i] = 127;
    } else {
      const double prob = (double) rels[i] / (double) surv[i];
      C[i] = -log(prob) * scale; /* TODO: How should we round here? */
    }
  }
}

/* {{{ initializing norms */
/* Knowing the norm on the rational side is bounded by 2^(2^k), compute
 * lognorms approximations for k bits of exponent + NORM_BITS-k bits
 * of mantissa.
 * Do the same for the algebraic side (with the corresponding bound for
 * the norms.
 */
void
init_norms (sieve_info_ptr si, int s)
{
      sieve_side_info_ptr sdata = si->sides[s];

      unsigned char *S = sdata->lognorm_table;

      int k = (int) ceil (log2 (sdata->logmax));
      int K = 1 << k;
      ASSERT_ALWAYS(NORM_BITS >= k);
      int l = NORM_BITS - k;
      int L = 1 << l;

      /* extract k bits from the exponent, and l bits from the mantissa */
      double h = 1.0 / (double) L;
      double e,m;
      int i,j;
      for (e = 1.0, i = 0; i < K; i++, e *= 2.0)
      {
          /* e = 2^i for 0 <= i < 2^k */
          for (m = 1.0, j = 0; j < L; j++, m += h)
          {
              /* m = 1 + j/2^l for 0 <= j < 2^l */
              double norm = m * e;
              /* Warning: since sdata->logmax does not usually correspond to
                 a power a two, and we consider full binades here, we have to
                 take care that values > sdata->logmax do not wrap around to 0 */
              norm = log2 (norm);
              if (norm >= sdata->logmax)
                  S[(i << l) + j] = 255;
              else
              {
                  norm = norm * sdata->scale;
                  S[(i << l) + j] = GUARD + (unsigned char) norm;
              }
          }
      }
}

/* }}} */

/* Input: i, double. i > 0 
   Output: o , trunc(o) == trunc(log2(i)) && o <= log2(i) < o + 0.0861.
   Careful: o ~= log2(i) iif add = 0x3FF00000 & scale = 1/0x100000.
   Add & scale are need to compute o'=f(log(i)) where f is an affine function */
inline int inttruncfastlog2(double i, double add, double scale) {
#ifdef HAVE_SSE2
  double dummy;
  int o;

  __asm__ ("movq %2,%1\n"
	   "psrlq $0x20,%1\n"
	   "cvtdq2pd %1, %1\n"
	   "subsd %3, %1\n"
	   "mulsd %4, %1\n"
	   "cvttsd2si %1, %0\n"
	   : "=r" (o), "=x" (dummy) : "x" (i), "x" (add), "x" (scale));
  return (o);
}
#else
/* Same function, but in x86 gcc needs to transfer the input i from a
   xmm register to a classical register. No other way than use memory.
   So this function needs at least 6 to 8 cycles more than the previous,
   which uses ~3 cycles.
*/
/* dummy to avoid gcc warning: dereferencing type-punned pointer */
  void *dummy = &i;
  return (int) ((((double) ((*((uint64_t *) dummy)) >> 0x20)) - add) * scale);
}
#endif

/* Initialize lognorms on the rational side for the bucket_region
 * number N.
 * For the moment, nothing clever, wrt discarding (a,b) pairs that are
 * not coprime, except for the line j=0.
 */
void init_rat_norms_bucket_region(unsigned char *S,
                                 /* no condition array here ! */
                                 const int N,
                                 sieve_info_ptr si)
{
  /* #define DEBUG_INIT_RAT 1 */
  sieve_side_info_ptr rat = si->sides[RATIONAL_SIDE];
  int halfI = (si->I)>>1,
    int_i;
  unsigned int						\
    j0 = N << (LOG_BUCKET_REGION - si->conf->logI),
    j1 = j0 + (1 << (LOG_BUCKET_REGION - si->conf->logI)),
    j = j0,
    inc;
  uint8_t oy, y;
  double						\
    u0 = si->sides[RATIONAL_SIDE]->fijd[0], // gj
    u1 = si->sides[RATIONAL_SIDE]->fijd[1], // gi
    invu1 = 1.0/u1,
    u0j = u0 * j,
    d0_init,
    scale = rat->scale * (1.0/0x100000),
    add = 0x3FF00000 - GUARD / scale,
    g, rac, d0, d1, i;
  size_t ts;

  static double cexp2[257]; /*  cexp2[UCHAR_MAX+1] INCLUDED! It may be used (not in S) */
  static double ratscale = 0;
  
  if (ratscale != rat->scale) {
    double step, i;
    ratscale = rat->scale;
    for (inc = 0, step = 1 / ratscale, i = -step * GUARD;
	 inc <= 256; i += step) cexp2[inc++] = exp2(i);
  }
  d0_init = cexp2[((unsigned int)GUARD) - 1U];
  if (!j) {
    // compute only the norm for i = 1. Everybody else is 255.
    memset(S, 255, halfI<<1);
    S[halfI + 1] = inttruncfastlog2(fabs(u1),add,scale);
    S+= halfI<<1;
    j++;
    u0j += u0;
  }

  for( ; j < j1 ; j++, u0j += u0) {

    int_i = -halfI;
    g = u0j + u1 * int_i;
    rac = u0j * (-invu1);
    d0 = d0_init;
    d1 = rac - d0 * rac;
    if (g > 0) {
      y = inttruncfastlog2 (g, add, scale);
      if (rac > -halfI) goto cas1; else goto cas4;
    }
    else {
      g = -g;
      y = inttruncfastlog2 (g, add, scale);
      if (rac > -halfI) goto cas3; else goto cas2;
    }
  cas1:
    /* In this case, we exit from the loop when ts == 0, at the exception
       of the first iteration. In this special case, old_i = -halfI and
       int_i = trunc (i), where i=[inverse of the function g](trunc(y)) and
       y=g(old_i).
       So, it's possible if y is very near trunc(y), old_i == int_i, so ts == 0.
       We have to iterate at least one time to avoid this case => this is the
       use of inc here. */
    for (i = rac + cexp2[y] * invu1, inc = 1;; y--) {
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
	y = inttruncfastlog2 (g, add, scale);
#ifdef DEBUG_INIT_RAT
	fprintf (stderr, "A2.1 : i=%d, y=%u, rac=%e\n", int_i, y, rac);
#endif
	*S++ = y;
	g += u1;
      }
      goto nextj;
    }
    for (inc = 0; g > 0; g += u1) {
      y = inttruncfastlog2 (g, add, scale);
#ifdef DEBUG_INIT_RAT
      fprintf (stderr, "A2.2 : i=%d, y=%u, rac=%e\n", int_i + inc, y, rac);
#endif
      S[inc++] = y;
    }
    int_i += inc;
    S += inc;
    g = -g;
    y = inttruncfastlog2 (g, add, scale);
  cas2:
    do {
#ifdef DEBUG_INIT_RAT
      fprintf (stderr, "A3 : i=%d, y=%u, rac=%e\n", int_i, y, rac);
#endif
      *S++ = y;
      if (++int_i >= halfI) goto nextj;
      oy = y;
      g -= u1;
      y = inttruncfastlog2 (g, add, scale);
    } while (oy != y);
    d0 = 1.0/d0;
    d1 = rac - d0 * rac;
    y++;
    i = rac - cexp2[(unsigned int)y + 1] * invu1;
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
    for (i = rac - cexp2[y] * invu1, inc = 1;; y--) {
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
	y = inttruncfastlog2 (g, add, scale);
#ifdef DEBUG_INIT_RAT
	fprintf (stderr, "B2.1 : i=%d, y=%u, rac=%e\n", int_i, y, rac);
#endif
	*S++ = y;
	g -= u1;
      }
      goto nextj;
    }
    for (inc = 0; g > 0; g -= u1) {
      y = inttruncfastlog2 (g, add, scale);
#ifdef DEBUG_INIT_RAT
      fprintf (stderr, "B2.2 : i=%d, y=%u, rac=%e\n", int_i + inc, y, rac);
#endif
      S[inc++] = y;
    }
    int_i += inc;
    S += inc;
    g = -g;
    y = inttruncfastlog2 (g, add, scale);
  cas4:
    do {
#ifdef DEBUG_INIT_RAT
      fprintf (stderr, "B3 : i=%d, y=%u, rac=%e\n", int_i, y, rac);
#endif
      *S++ = y;
      if (++int_i == halfI) goto nextj;
      oy = y;
      g += u1;
      y = inttruncfastlog2 (g, add, scale);
    } while (oy != y);
    d0 = 1.0/d0;
    d1 = rac - d0 * rac;
    y++;
    i = rac + cexp2[(unsigned int)y + 1] * invu1;
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
    __asm__("# gcc needs something after a label.\n");
  }
}

/* {{{ some utility stuff */
#ifdef SSE_NORM_INIT
static inline void init_fpoly_v2df(__v2df * F, const double *f, const int deg)
{
    int i;
    for (i = 0; i <= deg; ++i) {
	__v2df tmp = { f[i], f[i] };
	F[i] = tmp;
    }
}

static inline __v2df fpoly_eval_v2df(const __v2df * f, const __v2df x, int d)
{
    __v2df r;
    r = f[d--];
    for (; d >= 0; --d)
	r = r * x + f[d];
    return r;
}
#endif

static inline void fpoly_scale(double * u, const double * t, int d, double h)
{
    u[d] = t[d];        /* This one never changes. */
    double hpow = h;
    u[d-1] = t[d-1] * hpow;
    for (int k = 2; k <= d; k++) {
        hpow *= h;
        u[d - k] = t[d - k] * hpow;
    }
}
/* }}} */

/* Initialize lognorms on the algebraic side for the bucket
 * number N.
 * Only the survivors of the other sieve will be initialized, the
 * others are set to 255. Case GCD(i,j)!=1 also gets 255.
 * return the number of reports (= number of norm initialisations)
 *
 * ``the other side'' is checked via xB[*xS] > 127.
 */
int
init_alg_norms_bucket_region(unsigned char *S,
			     const unsigned char *xS MAYBE_UNUSED, const int N,
			     sieve_info_ptr si)
{
    sieve_side_info_ptr alg = si->sides[ALGEBRAIC_SIDE];

    uint64_t mask = (1 << NORM_BITS) - 1;
    union { double z; uint64_t x; } zx[1];
    int l = NORM_BITS - (int) ceil(log2(alg->logmax));
    int d = si->cpoly->alg->degree;

    const unsigned char * L = alg->lognorm_table;

#ifndef LAZY_NORMS
    /* Used for checking the bound on the other side */
    const unsigned char * xB = si->sides[RATIONAL_SIDE]->Bound;
    const unsigned char * xS0 MAYBE_UNUSED = xS;
#endif

    int32_t I = si->I;
    double halfI = I / 2;
    int logI = si->conf->logI;
    unsigned int j0 = N << (LOG_BUCKET_REGION - logI);
    unsigned int j1 = j0 + (1 << (LOG_BUCKET_REGION - logI));

    int report = 0;

    double u[d+1];

    /* bucket_region is a multiple of I. Let's find the starting j
     * corresponding to N and the last j.
     */
    for (unsigned int j = j0 ; j < j1 ; j++) {
#ifdef LAZY_NORMS
        unsigned int jmask = VERT_NORM_STRIDE - 1;
        if (j > j0 && ((j & jmask) != 0)) {
            // copy norms from the row below.
            unsigned char *old_S;
            old_S = S - si->I;
            for (unsigned int ii = 0; ii < si->I; ++ii) {
                *S++ = *old_S++;
            }
            continue;
        }
#endif
        /* scale by rj^(d-k) the coefficients of fij */
        fpoly_scale(u, alg->fijd, d, j);

#ifdef LAZY_NORMS
        const int normstride=NORM_STRIDE;
        const double fpnormstride=(double)NORM_STRIDE;
        for (int i = -halfI; i < halfI; i += fpnormstride) {
            unsigned char n;
            zx->z = fpoly_eval (u, d, i);
            /* 4607182418800017408 = 1023*2^52 */
            uint64_t y = (zx->x - (uint64_t) 0x3FF0000000000000) >> (52 - l);
            report++;
            n = L[y & mask];
            ASSERT (n > 0);
            for (int ii = 0; ii < normstride; ++ii)
                *S++ = n;
        }
#else  /* full norm computation */


#ifdef SSE_NORM_INIT
	__v2df u_vec[d+1];
	init_fpoly_v2df(u_vec, u, d);
	double previous_ri = 0.0;
	unsigned char * previous_S= NULL;
	int cpt = 0;
#endif
	/* now compute norms */
	for (double ri = -halfI ; ri < halfI; ri += 1.0) {
#ifdef TRACE_K
	    if (trace_on_spot_Nx(N, xS - xS0)) {
		fprintf(stderr, "# init_alg_norms_bucket_region: "
			"xS[%zd] = %u, -log(prob) = %u\n",
			xS - xS0,
			(unsigned int) *xS,
			(unsigned int) xB[*xS]);
	    }
#endif
	    if (xB[*xS] < 127) {
//  The SSE version seems to be slower on the algebraic side...
//  Let's forget it for the moment.
//  TODO: try with single precision.
//#ifndef SSE_NORM_INIT
#if 1
		zx->z = fpoly_eval(u, d, ri);
		/* the constant is 1023*2^52 */
		uint64_t y = (zx->x - UINT64_C(0x3FF0000000000000)) >> (52 - l);
		report++;
		unsigned char n = L[y & mask];
		ASSERT(n > 0);
		*S = n;
#else
		__v2di mask_vec = { mask, mask };
		__v2di cst_vec = {
                    UINT64_C(0x3FF0000000000000),
		    UINT64_C(0x3FF0000000000000),
		};
		__v2di shift_value =
		    { (uint64_t) (52 - l), (uint64_t) (52 - l) };
		report++;
		if (cpt == 0) {
		    previous_ri = ri;
		    previous_S = S;
		    cpt++;
		} else {
		    cpt--;
		    __v2df i_vec = { previous_ri, ri };
		    union {
			__v2df dble;
			__v2di intg;
			struct {
			    uint64_t y0;
			    uint64_t y1;
			} intpair;
		    } fi_vec;
		    fi_vec.dble = fpoly_eval_v2df(u_vec, i_vec, d);
		    fi_vec.intg -= cst_vec;
		    fi_vec.intg = _mm_srl_epi64(fi_vec.intg, shift_value);
		    fi_vec.intg &= mask_vec;
		    *previous_S = L[fi_vec.intpair.y0];
		    *S = L[fi_vec.intpair.y1];
		}
#endif
	    } else {
		*S = UCHAR_MAX;
	    }
	    S++;
	    xS++;
	}
#ifdef SSE_NORM_INIT
	if (cpt == 1) {		// odd number of computations
	    zx->z = fpoly_eval(u, d, previous_ri);
	    uint64_t y = (zx->x - (uint64_t) 0x3FF0000000000000) >> (52 - l);
	    *previous_S = L[y & mask];
	}
#endif
#endif  /* LAZY_NORMS */
    }
    return report;
}

/* }}} */

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
  for (int side = 0; side < 2; side++)
    {
      int d = si->cpoly->pols[side]->degree;
      si->sides[side]->fij = (mpz_t *) malloc((d + 1) * sizeof(mpz_t));
      FATAL_ERROR_CHECK(si->sides[side]->fij == NULL, "malloc failed");
      si->sides[side]->fijd = (double *) malloc((d + 1) * sizeof(double));
      FATAL_ERROR_CHECK(si->sides[side]->fijd == NULL, "malloc failed");
      for (int k = 0; k <= d; k++)
        mpz_init (si->sides[side]->fij[k]);
    }

  double r, maxlog2;
  unsigned char alg_bound, rat_bound;
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

  /* Since |a| <= s*I*B and |b| <= I*B, |G(a,b)| <= (|g[1]|*s+|g[0]|) * I*B */
  r = fabs (mpz_get_d (si->cpoly->rat->f[1])) * si->cpoly->skew
    + fabs (mpz_get_d (si->cpoly->rat->f[0]));
  r *= B * (double) si->I;
  /* if the special-q is on the rational side, divide by it */
  if (qside == RATIONAL_SIDE)
    r /= q0d;
  rat->logmax = maxlog2 = log2 (r);
  /* we know that |G(a,b)| < 2^(rat->logmax) when si->ratq = 0,
     and |G(a,b)/q| < 2^(rat->logmax) when si->ratq <> 0 */

  fprintf (output, "# Rat. side: log2(maxnorm)=%1.2f logbase=%1.6f",
           maxlog2, exp2 (maxlog2 / ((double) UCHAR_MAX - GUARD)));
  /* we want to map 0 <= x < maxlog2 to GUARD <= y < UCHAR_MAX,
     thus y = GUARD + x * (UCHAR_MAX-GUARD)/maxlog2 */
  rat->scale = ((double) UCHAR_MAX - GUARD + 1) / maxlog2 * 0.999999;
  /* fprintf (stderr, "rat->scale, my fomula = %f\n", rat->scale); */
  /* rat->scale = ((double) UCHAR_MAX - GUARD) / maxlog2; */
  /* fprintf (stderr, "rat->scale, old correct fomula = %f\n", rat->scale); */

  /* we want to select relations with a cofactor of less than r bits on the
     rational side */
  r = MIN(si->conf->sides[RATIONAL_SIDE]->lambda * (double) si->conf->sides[RATIONAL_SIDE]->lpb, maxlog2 - GUARD / rat->scale);
  rat_bound = (unsigned char) (r * rat->scale) + GUARD;
  fprintf (output, " bound=%u\n", rat_bound);
  double max_rlambda = (maxlog2 - GUARD / rat->scale) / si->cpoly->rat->lpb;
  if (si->cpoly->rat->lambda > max_rlambda) {
      fprintf(output, "# Warning, rlambda>%.1f does not make sense (capped to limit)\n", max_rlambda);
  }
  sieve_info_init_lognorm (rat->Bound, rat_bound, si->conf->sides[RATIONAL_SIDE]->lim,
                           si->conf->sides[RATIONAL_SIDE]->lpb, rat->scale);

  /************************** algebraic side *********************************/

  alg->logmax = get_maxnorm_alg (si->cpoly, si, q0d, B, qside); /* log2(max norm) */
  /* we know that |F(a,b)/q| < 2^(alg->logmax) when si->ratq = 0,
     and |F(a,b)| < 2^(alg->logmax) when si->ratq <> 0 */

  /* on the algebraic side, we want that the non-reports on the rational
     side, which are set to 255, remain larger than then report bound 'r',
     even if the algebraic norm is totally smooth. For this, we artificially
     increase by 'r' the maximal range */
  r = MIN(si->conf->sides[ALGEBRAIC_SIDE]->lambda * (double) si->conf->sides[ALGEBRAIC_SIDE]->lpb, alg->logmax);
  maxlog2 = alg->logmax + r;

  fprintf (output, "# Alg. side: log2(maxnorm)=%1.2f logbase=%1.6f",
           alg->logmax, exp2 (maxlog2 / ((double) UCHAR_MAX - GUARD)));
  /* we want to map 0 <= x < maxlog2 to GUARD <= y < UCHAR_MAX,
     thus y = GUARD + x * (UCHAR_MAX-GUARD)/maxlog2 */
  alg->scale = ((double) UCHAR_MAX - GUARD) / maxlog2;
  /* we want to report relations with a remaining log2-norm after sieving of
     at most lambda * lpb, which corresponds in the y-range to
     y >= GUARD + lambda * lpb * scale */
  alg_bound = (unsigned char) (r * alg->scale) + GUARD;
  fprintf (output, " bound=%u\n", alg_bound);
  double max_alambda = (alg->logmax) / si->cpoly->alg->lpb;
  if (si->cpoly->alg->lambda > max_alambda) {
      fprintf(output, "# Warning, alambda>%.1f does not make sense (capped to limit)\n", max_alambda);
  }
  sieve_info_init_lognorm (alg->Bound, alg_bound, si->conf->sides[ALGEBRAIC_SIDE]->lim,
                           si->conf->sides[ALGEBRAIC_SIDE]->lpb, alg->scale);
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

void
sieve_info_update_norm_data (sieve_info_ptr si)
{
    int64_t H[4] = { si->a0, si->b0, si->a1, si->b1 };
    /* Update floating point version of algebraic poly (do both, while
     * we're at it...) */
    for (int side = 0; side < 2; side++) {
        sieve_side_info_ptr s = si->sides[side];
        cado_poly_side_ptr ps = si->cpoly->pols[side];
        mp_poly_homography (s->fij, ps->f, ps->degree, H);
        double invq = 1.0;
        if (si->conf->side == side)
            invq /= mpz_get_d(si->doing->p);
        for (int k = 0; k <= ps->degree; k++)
            s->fijd[k] = mpz_get_d (s->fij[k]) * invq;
    }
}
