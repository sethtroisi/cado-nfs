
/* Most of the code in this file was contributed by AF */

/* Same as previous with a fabs on i before the computation of the log2.
 */
static inline double lg2_abs_raw (double i) {
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
  __asm__ __volatile__ (
            "psllq $0x01,  %0    \n"
            "psrlq $0x21,  %0    \n"
            "cvtdq2pd      %0, %0\n" /* Mandatory in packed double even it's non packed! */
            : "+&x" (i));            /* Really need + here! i is not modified in C code! */
  return i;
#else
  void *tg = &i;
  return (double)((*(uint64_t *)tg << 1) >> 0x21);
#endif
}
static inline double lg2_abs (double i, double offset, double scale) {
    return (lg2_abs_raw(i) - offset) * scale;
}

static inline double compute_f (const unsigned int d, const double *u, const double h) {
  size_t k = (size_t) d;
  double f = u[k];
  switch (k) {
  default: do { f = f * h + u[--k]; } while (k > 9);
  case 9: f = f * h + u[8]; no_break();
  case 8: f = f * h + u[7]; no_break();
  case 7: f = f * h + u[6]; no_break();
  case 6: f = f * h + u[5]; no_break();
  case 5: f = f * h + u[4]; no_break();
  case 4: f = f * h + u[3]; no_break();
  case 3: f = f * h + u[2]; no_break();
  case 2: f = f * h + u[1]; no_break();
  case 1: f = f * h + u[0]; no_break();
  case 0: break;
  }
  return f;
}
static inline double compute_f (double_poly_srcptr u, const double h) {
    return compute_f(u->deg, u->coeff, h);
}
// same but another name. We want to get rid of it someday.
static inline void poly_scale_double (double  *u, const double *t, unsigned int d, const double h)
{
    u[d] = t[d];
    for (double hpow = h; d--; hpow *= h) u[d] = t[d] * hpow;
}
#if defined(HAVE_SSE2) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
static inline void
poly_scale_m128d (__m128d  *u, const double *t, unsigned int d, const double h)
{
  u[d] = _mm_set1_pd (t[d]);
  for (double hpow = h; d--; hpow *= h) u[d] = _mm_set1_pd (t[d] * hpow);
}
#endif  /* defined(HAVE_SSE2) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM) */

/***********************************************************************/
/* abandon all hope here. */

/* "smart" lognorms */
/* Details of the smart norms algorithms.
 *
 * This initialization computes first the roots of F, F', F". Each of
 * these defines a line which intercepts (0, 0.).  For each j, the roots
 * of f, f' and f" are computed : there are F, F', F" roots * j.
 *
 * Because of the absolute value, the roots of f have an unstable
 * neighbourhood : f "bounces" on the horizontal axis.
 *
 * The roots of f" have also an unstable neighbourhood (inflexion points of f).
 *
 * The roots of f have a stable neighbourhood.
 *
 * So, the neighbourhoods of f(root(f)) and f(root(f")) are computed
 * until on each side of the root there are SMART_NORM_STABILITY
 * identical values, so until f has a local horizontal stability on the
 * left and on the right of the root.
 * A root has a maximal influence of SMART_NORM_INFLUENCE values on its
 * neighbourhood (on each sides).
 * 
 * These roots define some segments of contiguous values of f(i). Some of
 * them are reduced to a point (f(roots(f')); the length of the others is
 * between SMART_NORM_STABILITY * 2 + 1 and SMART_NORM_MAX_INFLUENCE * 2
 * + 1.

 * 3 artificials roots are inserted: -I/2 and (I/2)-1 as two roots of f'
 * (two one-point segment), and 0.0 as a root of f, because near the
 * neighbourhood of 0, f is very unstable.

 * To compute all missing values of f(i) between ]-I/2,(I/2)-1[, a polygonal
 * approximation is done between all these segments.
 * The polygonal approximation between (i0,f(i0)-(i1,f(i1)) has 2
 * parameters :
 * - The minimal lenght between i0 and i1, SMART_NORM_LENGTH. Below this
 *   value, the values f(i0...i1) are approximated by a line.
 * - The maximal acceptable distance between f((i0+i1)/2) and
 *   (f(i0)+f(i1))/2, SMART_NORM_DISTANCE. Above this value,
 *   (i0,f(i0)-(i1,f(i1)) is approximated by the polygonal approximations
 *   of (i0,f(i0)-((i0+i1)/2,(f(i0)+f(i1))/2) and
 *   ((i0+i1)/2,(f(i0)+f(i1))/2)-(i1,f(i1)).

 * The maximal error of the smart initialization depends on its four
 * parameters, but should be less or egal to the double of
 * SMART_NORM_DISTANCE.  Modifying the parameters of the smart norm
 * algorithm is not a good idea. If you try it, use
 * test_init_norm_bucket_region in the tests part in order to have an
 * idea of the corresponding errors.
 */ 

/* {{{ some utility stuff */
#if defined(HAVE_SSE2) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
/* This function is for the SSE2 init algebraics, but for the exact initialization.
   Same than previous but return a SSE2 register with only 16 lowest bits are computed
   as the 2 results (2 8-bits values).
   CAREFUL! This function returns the computed value in the data i.
   the i value is modified by this function !
*/
static inline __m128i _mm_lg2_abs(__m128d *i, const __m128d offset, const __m128d scale) {
  __asm__ __volatile__ (
           "psllq     $0x01,    %0       \n" /* Dont use pabsd! */
           "psrlq     $0x01,    %0       \n"
           "shufps    $0xED,    %0,    %0\n"
           "cvtdq2pd  %0,       %0       \n"
           : "+&x"(*i));
  *i = _mm_mul_pd(_mm_sub_pd(*i, offset), scale);
  __asm__ __volatile__ (
           "cvttpd2dq %0,       %0       \n" /* 0000 0000 000X 000Y */
           "packssdw  %0,       %0       \n" /* 0000 0000 0000 0X0Y */
           "packuswb  %0,       %0       \n" /* 0000 0000 0000 00XY */
           : "+&x"(*i));
  return *(__m128i *) i;
}
#endif  /* defined(HAVE_SSE2) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM) */
/* }}} */
/* disabled, unused utility functions *//*{{{*/
#if 0
/* unused, untested. */
/* Same than previous with the result is duplicated 8 times
   in a "false" double in SSE version, i.e. in a xmm register, or in a 64 bits
   general register in non SSE version, and written at addr[decal].
   The SSE version has 4 interests: no memory use, no General Register use,
   no xmm -> GR conversion, no * 0x0101010101010101 (3 to 4 cycles) but only
   2 P-instructions (1 cycle).
*/
MAYBE_UNUSED static inline void w64lg2_abs(double i, double offset, double scale, uint8_t *addr, ssize_t decal) {
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
  __asm__ __volatile__ (
            "psllq $0x01,  %0                 \n"
            "psrlq $0x21,  %0                 \n"
            "cvtdq2pd      %0,              %0\n"
            : "+&x" (i));
  i = (i - offset) * scale;
  __asm__ __volatile__ (
            "cvttpd2dq     %0,       %0       \n" /* 0000 0000 0000 000Y */
            "punpcklbw     %0,       %0       \n" /* 0000 0000 0000 00YY */
            "pshuflw    $0x00,       %0,    %0\n" /* 0000 0000 YYYY YYYY */
            : "+&x" (i));
  *(double *)&addr[decal] = i;
#else
  void *tg = &i;
  *(uint64_t *)&addr[decal] = UINT64_C(0x0101010101010101) * (uint64_t) (((double)((*(uint64_t *)tg << 1) >> 0x21) - offset) * scale);
#endif
}
#endif

#if defined(HAVE_SSE2) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)

#if 0
/* unused, untested. */
/* This function is for the SSE2 init algebraics.
   I prefer use X86 ASM directly and avoid intrinsics because the trick of
   cvtdq2pd (I insert 2 doubles values and do a cvtdq2pd on them in order to
   compute their log2).
   This function does a double abs of the __m128d i, and computes the double
   log2 of the result.
   Careful: if i == 0, the result is not predictible.
*/
MAYBE_UNUSED static inline void w128lg2_abs(__m128d i, const __m128d offset, const __m128d scale, uint8_t *addr, const ssize_t decal) {
  __asm__ __volatile__ (
           "psllq     $0x01,    %0       \n" /* Dont use pabsd! */
           "psrlq     $0x01,    %0       \n" /* i = fast_abs(i) */
           "shufps    $0xED,    %0,    %0\n"
           "cvtdq2pd  %0,       %0       \n"
           : "+&x"(i));
  i = _mm_mul_pd(_mm_sub_pd(i, offset), scale);
  __asm__ __volatile__ (
           "cvttpd2dq %0,       %0       \n" /* 0000 0000 000X 000Y */
           "packssdw  %0,       %0       \n" /* 0000 0000 0000 0X0Y */
           "punpcklbw %0,       %0       \n" /* 0000 0000 00XX 00YY */
           "punpcklwd %0,       %0       \n" /* 0000 XXXX 0000 YYYY */
           "pshufd    $0xA0,    %0,    %0\n" /* XXXX XXXX YYYY YYYY */
           : "+&x"(i));
  *(__m128d *)&addr[decal] = i; /* addr and decal are 16 bytes aligned: MOVAPD */
}
#endif

#endif
/*}}}*/

/* This function computes the roots of the polynomial F(i,1) = ln2(f(i))
   and the roots of the first and second derivatives of F.
   These roots are useful to begin the smart initialization of the norms.
   The algorithm is described below

   Inputs :

   F(i,1): degree and coeff(icients)
   max_abs_root = Absolute maximum value of a root (before and after, the root is useless).
   precision = minimal need precision for each root (usually 1/J = 2/I).
   Output:
   roots = the roots of F, F', F", and pseudo root 0.0 added in F roots.

   NEEDED:
   1. roots size must >= 1 + degree  + max (0, degree - 1) + max (0, 2 * degree - 2)
   so degree * 4 + 1.
   2. p & Roots: the biggest polynomial is F", degree = 2 * (degree - 1).
   So, p size is 2 * degree - 1, and Roots size is 2 * degree - 2.
*/

bool operator<(lognorm_oldsmart::smart_norm_root const & a, lognorm_oldsmart::smart_norm_root const & b)
{
  return a.value < b.value;
}

static long lg_page;
/* The norm initialisation code needs the page size; we set it here as this
   function is not thread-safe, to avoid race conditions */
static void set_lg_page()
{
  if (lg_page == 0)
    lg_page = pagesize();
}


/* {{{ init_norms_roots_internal */
/* This works on the integer polynomial obtained by converting
 * the double polynomial obtained by converting the integer polynomial.
 *
 * So we're truncating the input integer coefficients to 53 bits.
 *
 * It's weird, and quite probably needlessly complicated.
 */
void init_norms_roots_internal (cxx_double_poly const & f, double max_abs_root, double precision, std::vector<lognorm_oldsmart::smart_norm_root> & roots)
{
  unsigned int degree = f->deg;
  cxx_mpz_poly fz;
  usp_root_data Roots[degree << 1];
  unsigned int n;

  set_lg_page();
  for (unsigned int k = 0 ; k < (degree << 1) ; k++)
      usp_root_data_init (&(Roots[k]));
  mpz_poly_set_double_poly(fz, f);

  roots.clear();
  roots.reserve(4 * degree + 1);
  /* Pseudo root 0.0 is inserted first as a root of F */
  roots.push_back(lognorm_oldsmart::smart_norm_root());

  if (degree) {
    /* The roots of F are inserted in roots */
    n = numberOfRealRoots (fz->coeff, fz->deg, max_abs_root, 0, Roots);
    for (unsigned int k = 0 ; k < n ; k++) {
        lognorm_oldsmart::smart_norm_root s(0, rootRefine (&(Roots[k]), fz->coeff, degree, precision));
        roots.push_back(s);
    }

    /* Computation of F' */
    cxx_double_poly df;
    double_poly_derivative (df, f);

    /* The roots of F' are inserted in roots */
    mpz_poly_set_double_poly(fz, df);
    if (fz->deg) {
        n = numberOfRealRoots (fz->coeff, fz->deg, max_abs_root, 0, Roots);
        for (unsigned int k = 0 ; k < n ; k++) {
            lognorm_oldsmart::smart_norm_root s(1, rootRefine (&(Roots[k]), fz->coeff, df->deg, precision));
            roots.push_back(s);
        }
    }

    /* Computation of F" */
    /* XXX Hmm. We're computing (f/f')', here...  */
    cxx_double_poly df_df;
    cxx_double_poly ddf;
    cxx_double_poly f_ddf;
    cxx_double_poly d2f;

    double_poly_mul (df_df, df, df);
    double_poly_derivative (ddf, df);
    double_poly_mul (f_ddf, f, ddf);
    double_poly_sub (d2f, f_ddf, df_df);

    /* The roots of F" are inserted in roots */
    mpz_poly_set_double_poly(fz, d2f);
    if (fz->deg) {
        n = numberOfRealRoots (fz->coeff, fz->deg, max_abs_root, 0, Roots);
        for (unsigned int k = 0 ; k < n ; k++) {
            lognorm_oldsmart::smart_norm_root s(2, rootRefine (&(Roots[k]), fz->coeff, d2f->deg, precision));
            roots.push_back(s);
        }
    }

    std::sort(roots.begin(), roots.end());
  }

  for (unsigned int k = 0 ; k < (degree << 1) ; k++)
      usp_root_data_clear (&(Roots[k]));
}
/*}}}*/

/* {{{ lognorm_fill_alg_exact */
/* Exact initialisation of F(i,j) with degre >= 2 (not mandatory). Slow.
   Internal function, only with simple types, for unit/integration testing. */
void lognorm_fill_alg_exact (unsigned char *S, uint32_t N, int logI, double scale, cxx_double_poly const & fijd)
{

#if !defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM) || !defined(HAVE_SSSE3)
    /* Fall back on the reference code */
  lognorm_fill_alg_reference(S, N, logI, scale, fijd);

#else /* HAVE_GCC_STYLE_AMD64_INLINE_ASM && HAVE_SSSE3 : optimized part. Stupid but fast code. */

  unsigned int d = fijd->deg;
  unsigned char *beginS = S;
  uint32_t beginJ, endJ = LOG_BUCKET_REGION - logI;
  scale *= 1./0x100000;
  const double offset = 0x3FF00000 - GUARD / scale;
  uint32_t J = N << endJ;
  beginJ = J;
  endJ = (1U << endJ) + J;
#define BEGIN ".balign 8\n 0:\n add $0x10,%[S]\n"
#define FU(A) "movapd %[u" #A "],%[f]\n"
#define U8 "mulpd %[h],%[f]\n addpd %[u8],%[f]\n"
#define U7 "mulpd %[h],%[f]\n addpd %[u7],%[f]\n"
#define U6 "mulpd %[h],%[f]\n addpd %[u6],%[f]\n"
#define U5 "mulpd %[h],%[f]\n addpd %[u5],%[f]\n"
#define U4 "mulpd %[h],%[f]\n addpd %[u4],%[f]\n"
#define U3 "mulpd %[h],%[f]\n addpd %[u3],%[f]\n"
#define U2 "mulpd %[h],%[f]\n addpd %[u2],%[f]\n"
#define U1 "mulpd %[h],%[f]\n addpd %[u1],%[f]\n"
#define U0 "mulpd %[h],%[f]\n addpd %[u0],%[f]\n"
#define LG2ABS                                                         \
  "psllq $1,%[f]\n"                                                    \
    "psrlq $1,%[f]\n"                                                  \
    "shufps $0xED,%[f],%[f]\n"                                         \
    "cvtdq2pd %[f],%[f]\n"                                             \
    "subpd %[_offset],%[f]\n"                                             \
    "mulpd %[_scale],%[f]\n"                                           \
    "cvttpd2dq %[f],%[f]\n"                                            \
    "packssdw %[f],%[f]\n"                                             \
    "packuswb %[f],%[f]\n"                                             \
    "palignr $2,%[cumul],%[f]\n"                                       \
    "movapd %[f],%[cumul]\n"                                           \
    "addpd %[_two],%[h]\n"
#define END                                     \
  "cmp %[endS],%[S]\n"                          \
    "movapd %[f],-0x10(%[S])\n"                 \
    "jl 0b\n"

  int I = 1 << logI;
  for (; J < endJ; J++) {
    const __m128d _two = _mm_set1_pd(2.), _scale = _mm_set1_pd (scale), _offset = _mm_set1_pd (offset);
    __m128d f, h, u[d+1];
    __m128i cumul;
    unsigned char *endS = S + I;
    poly_scale_m128d (u, fijd->coeff, d, (double) J);
    h = _mm_set_pd ((double) (1 - (int32_t) (I >> 1)), (double) (- (int32_t) (I >> 1)));

    /* These ASM & switch are really ugly. But it's the ONLY way to
       be sure all the line of S is set only with registers. */
    switch (d) {
#if !defined(__ICC) || (__ICC >= 1600)
      /* the Intel compiler icpc fails with "internal error" with the code
         below (version 14.0.3 20140422), version 16.0.0 works */
    case 2:
      __asm__ __volatile__ ( BEGIN
      FU(2) U1 U0 LG2ABS FU(2) U1 U0 LG2ABS
      FU(2) U1 U0 LG2ABS FU(2) U1 U0 LG2ABS
      FU(2) U1 U0 LG2ABS FU(2) U1 U0 LG2ABS
      FU(2) U1 U0 LG2ABS FU(2) U1 U0 LG2ABS END
      : [f]"=&x"(f), [cumul]"=&x"(cumul), [h]"+&x"(h), [S]"+&r"(S)
      : [endS]"r"(endS), [_two]"x"(_two), [_scale]"x"(_scale), [_offset]"x"(_offset),
        [u0]"x"(u[0]), [u1]"x"(u[1]), [u2]"x"(u[2]));
      break;
    case 3:
      __asm__ __volatile__ ( BEGIN
      FU(3) U2 U1 U0 LG2ABS FU(3) U2 U1 U0 LG2ABS
      FU(3) U2 U1 U0 LG2ABS FU(3) U2 U1 U0 LG2ABS
      FU(3) U2 U1 U0 LG2ABS FU(3) U2 U1 U0 LG2ABS
      FU(3) U2 U1 U0 LG2ABS FU(3) U2 U1 U0 LG2ABS END
      : [f]"=&x"(f), [cumul]"=&x"(cumul), [h]"+&x"(h), [S]"+&r"(S)
      : [endS]"r"(endS), [_two]"x"(_two), [_scale]"x"(_scale), [_offset]"x"(_offset),
        [u0]"x"(u[0]), [u1]"x"(u[1]), [u2]"x"(u[2]), [u3]"x"(u[3]));
      break;
    case 4:
      __asm__ __volatile__ ( BEGIN
      FU(4) U3 U2 U1 U0 LG2ABS FU(4) U3 U2 U1 U0 LG2ABS
      FU(4) U3 U2 U1 U0 LG2ABS FU(4) U3 U2 U1 U0 LG2ABS
      FU(4) U3 U2 U1 U0 LG2ABS FU(4) U3 U2 U1 U0 LG2ABS
      FU(4) U3 U2 U1 U0 LG2ABS FU(4) U3 U2 U1 U0 LG2ABS END
      : [f]"=&x"(f), [cumul]"=&x"(cumul), [h]"+&x"(h), [S]"+&r"(S)
      : [endS]"r"(endS), [_two]"x"(_two), [_scale]"x"(_scale), [_offset]"x"(_offset),
        [u0]"x"(u[0]), [u1]"x"(u[1]), [u2]"x"(u[2]), [u3]"x"(u[3]), [u4]"x"(u[4]));
      break;
    case 5:
      __asm__ __volatile__ ( BEGIN
      FU(5) U4 U3 U2 U1 U0 LG2ABS FU(5) U4 U3 U2 U1 U0 LG2ABS
      FU(5) U4 U3 U2 U1 U0 LG2ABS FU(5) U4 U3 U2 U1 U0 LG2ABS
      FU(5) U4 U3 U2 U1 U0 LG2ABS FU(5) U4 U3 U2 U1 U0 LG2ABS
      FU(5) U4 U3 U2 U1 U0 LG2ABS FU(5) U4 U3 U2 U1 U0 LG2ABS END
      : [f]"=&x"(f), [cumul]"=&x"(cumul), [h]"+&x"(h), [S]"+&r"(S)
      : [endS]"r"(endS), [_two]"x"(_two), [_scale]"x"(_scale), [_offset]"x"(_offset),
        [u0]"x"(u[0]), [u1]"x"(u[1]), [u2]"x"(u[2]), [u3]"x"(u[3]), [u4]"x"(u[4]),
        [u5]"x"(u[5]));
      break;
    case 6:
      __asm__ __volatile__ ( BEGIN
      FU(6) U5 U4 U3 U2 U1 U0 LG2ABS FU(6) U5 U4 U3 U2 U1 U0 LG2ABS
      FU(6) U5 U4 U3 U2 U1 U0 LG2ABS FU(6) U5 U4 U3 U2 U1 U0 LG2ABS
      FU(6) U5 U4 U3 U2 U1 U0 LG2ABS FU(6) U5 U4 U3 U2 U1 U0 LG2ABS
      FU(6) U5 U4 U3 U2 U1 U0 LG2ABS FU(6) U5 U4 U3 U2 U1 U0 LG2ABS END
      : [f]"=&x"(f), [cumul]"=&x"(cumul), [h]"+&x"(h), [S]"+&r"(S)
      : [endS]"r"(endS), [_two]"x"(_two), [_scale]"x"(_scale), [_offset]"x"(_offset),
        [u0]"x"(u[0]), [u1]"x"(u[1]), [u2]"x"(u[2]), [u3]"x"(u[3]), [u4]"x"(u[4]),
        [u5]"x"(u[5]), [u6]"x"(u[6]));
      break;
    case 7:
      __asm__ __volatile__ ( BEGIN
      FU(7) U6 U5 U4 U3 U2 U1 U0 LG2ABS FU(7) U6 U5 U4 U3 U2 U1 U0 LG2ABS
      FU(7) U6 U5 U4 U3 U2 U1 U0 LG2ABS FU(7) U6 U5 U4 U3 U2 U1 U0 LG2ABS
      FU(7) U6 U5 U4 U3 U2 U1 U0 LG2ABS FU(7) U6 U5 U4 U3 U2 U1 U0 LG2ABS
      FU(7) U6 U5 U4 U3 U2 U1 U0 LG2ABS FU(7) U6 U5 U4 U3 U2 U1 U0 LG2ABS END
      : [f]"=&x"(f), [cumul]"=&x"(cumul), [h]"+&x"(h), [S]"+&r"(S)
      : [endS]"r"(endS), [_two]"x"(_two), [_scale]"x"(_scale), [_offset]"x"(_offset),
        [u0]"x"(u[0]), [u1]"x"(u[1]), [u2]"x"(u[2]), [u3]"x"(u[3]), [u4]"x"(u[4]),
        [u5]"x"(u[5]), [u6]"x"(u[6]), [u7]"x"(u[7]));
      break;
    case 8:
      __asm__ __volatile__ ( BEGIN
      FU(8) U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS FU(8) U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS
      FU(8) U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS FU(8) U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS
      FU(8) U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS FU(8) U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS
      FU(8) U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS FU(8) U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS END
      : [f]"=&x"(f), [cumul]"=&x"(cumul), [h]"+&x"(h), [S]"+&r"(S)
      : [endS]"r"(endS), [_two]"x"(_two), [_scale]"x"(_scale), [_offset]"x"(_offset),
        [u0]"x"(u[0]), [u1]"x"(u[1]), [u2]"x"(u[2]), [u3]"x"(u[3]), [u4]"x"(u[4]),
        [u5]"x"(u[5]), [u6]"x"(u[6]), [u7]"x"(u[7]), [u8]"x"(u[8]));
      break;
    case 9:
      __asm__ __volatile__ ( BEGIN
      FU(9) U8 U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS FU(9) U8 U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS
      FU(9) U8 U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS FU(9) U8 U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS
      FU(9) U8 U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS FU(9) U8 U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS
      FU(9) U8 U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS FU(9) U8 U7 U6 U5 U4 U3 U2 U1 U0 LG2ABS END
      : [f]"=&x"(f), [cumul]"=&x"(cumul), [h]"+&x"(h), [S]"+&r"(S)
      : [endS]"r"(endS), [_two]"x"(_two), [_scale]"x"(_scale), [_offset]"x"(_offset),
        [u0]"x"(u[0]), [u1]"x"(u[1]), [u2]"x"(u[2]), [u3]"x"(u[3]), [u4]"x"(u[4]),
        [u5]"x"(u[5]), [u6]"x"(u[6]), [u7]"x"(u[7]), [u8]"x"(u[8]), [u9]"x"(u[9]));
      break;
#endif
    default:
      do {
        f = u[d]; for (size_t k = d; k; f = _mm_add_pd(_mm_mul_pd(f,h),u[--k]));
        h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2_abs (&f, _offset, _scale), cumul, 2);
        f = u[d]; for (size_t k = d; k; f = _mm_add_pd(_mm_mul_pd(f,h),u[--k]));
        h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2_abs (&f, _offset, _scale), cumul, 2);
        f = u[d]; for (size_t k = d; k; f = _mm_add_pd(_mm_mul_pd(f,h),u[--k]));
        h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2_abs (&f, _offset, _scale), cumul, 2);
        f = u[d]; for (size_t k = d; k; f = _mm_add_pd(_mm_mul_pd(f,h),u[--k]));
        h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2_abs (&f, _offset, _scale), cumul, 2);
        f = u[d]; for (size_t k = d; k; f = _mm_add_pd(_mm_mul_pd(f,h),u[--k]));
        h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2_abs (&f, _offset, _scale), cumul, 2);
        f = u[d]; for (size_t k = d; k; f = _mm_add_pd(_mm_mul_pd(f,h),u[--k]));
        h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2_abs (&f, _offset, _scale), cumul, 2);
        f = u[d]; for (size_t k = d; k; f = _mm_add_pd(_mm_mul_pd(f,h),u[--k]));
        h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2_abs (&f, _offset, _scale), cumul, 2);
        f = u[d]; for (size_t k = d; k; f = _mm_add_pd(_mm_mul_pd(f,h),u[--k]));
        h = _mm_add_pd(h, _two); cumul = _mm_alignr_epi8 (_mm_lg2_abs (&f, _offset, _scale), cumul, 2);
        *(__m128i *) S = cumul;
        S += 16;
      } while (S < endS);
    } /* End of the switch */
  } /* End of the line */
#undef BEGIN
#undef FU
#undef U0
#undef U1
#undef U2
#undef U3
#undef U4
#undef U5
#undef U6
#undef U7
#undef U8
#undef LG2ABS
#undef END
  /* Special ultra rare case. The correction of log2(F(0,0)) is wrong, because
     the fast algorithm of log2 is not good in 0.0. */
  if (UNLIKELY(!beginJ)) beginS[I>>1] = GUARD;
#endif /* End of HAVE_GCC_STYLE_AMD64_INLINE_ASM */
}
/* }}} */

/* {{{ piecewise linear approximation used by smart norms */
/* Approximates [S[x],S[y][ by the segment (x,fx),(y,fy).
   (y - x) / (fy - fx) should be really larger than one for speed.
   CAREFUL : Need y > x. */
static inline void Fill_S (unsigned char *S, int x, double fx, int y, double fy) {
  int next_f = (int) fx, step_f = (int) fy - next_f;

  if (LIKELY (step_f)) {
    double m = (double) (y - x) / (fy - fx), next_m = (double) (x + 1);
    if (step_f < 0.) {
      m = fabs(m);
      step_f = -1;
      next_m += (fx - floor(fx)) * m;
    }
    else {
      step_f = 1;
      next_m += (ceil(fx) - fx) * m;
    }
    for (;;) {
      int next_x = (int) next_m;
      if (UNLIKELY (next_x >= y)) break;
      memset (S + x, next_f, (size_t) (next_x - x));
      x = next_x;
      next_f += step_f;
      next_m += m;
    }
  }
  memset (S + x, next_f, (size_t) (y - x));
}

/* These segments ((x, F(x)), (y, F(y))) are used in the smart normalization */
typedef struct sg_s {
  int begin, end;
  double f_begin, f_end;
} sg_t;

/* This functions sets all the segments of the contiguous unset values of
 * the line S with F(i, const j) by a polygonal approximation on these
 * segments.  Derecursivate optimal version.
 */
static inline void poly_approx_on_S (unsigned char *S, const unsigned int degree, const double *coeff, const double my_scale MAYBE_UNUSED, const double scale MAYBE_UNUSED, const double offset MAYBE_UNUSED, const unsigned int nsg, const sg_t *sg) {
#define SIZE_STACK 256 /* In fact, on RSA704 benchmark, the max is... 10 */
  typedef struct x_fx_s {
    int x;
    double fx; } x_fx_t;
  size_t current_stack, max_stack = SIZE_STACK;
  x_fx_t begin, current, *stack = (x_fx_t *) malloc_check (sizeof(*stack) *  max_stack);

  /* At least 3 segments: -Idiv2, neighbourhood of 0., Idiv2 - 1 */
  ASSERT(nsg >= 3);
  for (size_t parse_sg = 1; parse_sg < nsg; ++parse_sg) {
    ASSERT (sg[parse_sg - 1].end < sg[parse_sg].begin);
    begin.x = sg[parse_sg - 1].end;
    begin.fx = sg[parse_sg - 1].f_end;
    current.x = sg[parse_sg].begin;
    current.fx = sg[parse_sg].f_begin;
    current_stack = 0;
    for (;;) {
      if (LIKELY((unsigned int) (current.x - begin.x) >= SMART_NORM_LENGTH)) {
        x_fx_t possible_new;
        possible_new.x = (begin.x + current.x) >> 1;
        possible_new.fx = compute_f (degree, coeff, (double) possible_new.x);
        possible_new.fx = lg2_abs (possible_new.fx, offset, scale);
        if (UNLIKELY(fabs (possible_new.fx + possible_new.fx - current.fx - begin.fx)) > (2 * SMART_NORM_DISTANCE)) {
          stack[current_stack++] = current;
          current = possible_new;
          if (UNLIKELY(current_stack == max_stack)) {
            max_stack += (max_stack >> 1);
            if ((stack = (x_fx_t *) realloc (stack, sizeof(*stack) * max_stack)) == NULL) {
              fprintf (stderr, "Error, realloc of %zu bytes failed\n", sizeof(*stack) * max_stack);
              abort ();
            }
          }
          continue;
        }
      }
      Fill_S (S, begin.x, begin.fx, current.x, current.fx);
      if (!current_stack) break;
      begin = current;
      current = stack[--current_stack];
    }
  }
  free (stack);
}
/* }}} */

lognorm_oldsmart::lognorm_oldsmart(siever_config const & sc, cado_poly_srcptr cpoly, int side, qlattice_basis const & Q, int J) : lognorm_base(sc, cpoly, side, Q, J)/*{{{*/
{
    /* See init_degree_one_norms_bucket_region_internal for the explanation of
     * this table. Well, to be honest you won't get an explanation there... */
    for (int inc = 0; inc < 257; inc++)
        cexp2[inc] = exp2((inc-GUARD)/scale);

    int I = 1 << logI;

    init_norms_roots_internal (fijd, (I + 16) >> 1, 1. / (I >> 1), roots);
}/*}}}*/
/* {{{ lognorm_fill_alg_oldsmart */
/* Smart initialization of the normalization. Only useful for degree >= 2.
   Cf las-config.h, SMART_INIT for the algorithm.
   No SSE version: it's completly unreadable, and the gain is not really interesting.
   Internal function, only with simple types, for unit/integration testing */

void lognorm_fill_alg_oldsmart (unsigned char *S, uint32_t N, int logI, double original_scale, cxx_double_poly const & fijd, std::vector<lognorm_oldsmart::smart_norm_root> const & roots)
{
    int I = 1 << logI;
  unsigned int d = fijd->deg;
  ASSERT (d >= 2);
  /* F, F' and F" roots needs stability for their neighbourhood ?
     F roots, sure; F', sure not; F"... maybe. */
  const unsigned char stability_for_derivative[3] =
    { SMART_NORM_STABILITY, 0, SMART_NORM_STABILITY };
  const ssize_t Idiv2 = (ssize_t) (I >> 1);
  sg_t sg[d * 4 + 3]; /* For (F, F', F") roots, -Idiv2, Idiv2 - 1, 0 */
  size_t nsg;

  unsigned char *beginS;
  uint32_t beginJ, endJ = LOG_BUCKET_REGION - logI;
  double scale = original_scale * (1./0x100000);
  const double offset = ((double) 0x3FF00000) - GUARD / scale;
  uint32_t J = N << endJ;
  beginJ = J;
  endJ = (1U << endJ) + J;
  S += Idiv2; /* CAREFUL! Here, *S is the middle of the first line we have to compute! */
  beginS = S;

  for (; J < endJ; S += I, J++) {
    double u[d+1], g, hl, hr;
    unsigned char fg, f1, f2;
    unsigned int cptf2id = 0, cpt;
    ssize_t ih;
    {
      ASSERT_ALWAYS(lg_page != 0);
      for (ssize_t k = I; (k -= lg_page) >= 0; __builtin_prefetch (S + k, 1));
    }
    poly_scale_double (u, fijd->coeff, d, (double) J);

    /* Insertion of point (-Idiv2, F(-Idiv2)) in sg[0] : an artificial one-point segment. */
    g = compute_f (d, u, (double) -Idiv2);
    g = lg2_abs (g, offset, scale);
    sg[0].begin   = sg[0].end = -(int)Idiv2;
    sg[0].f_begin = sg[0].f_end = g;
    nsg = 1;
    S[-Idiv2] = (uint8_t) g;

    for (size_t r = 0; r < roots.size(); r++) {
      ASSERT_ALWAYS(nsg < d * 4 + 3);
      hl = floor((double)J * roots[r].value);
      ih = (ssize_t) hl;
      /* Need stability for this root ? */
      if (UNLIKELY(!stability_for_derivative[roots[r].derivative])) {
        /* No. It's an one-point segment. */
        /* Is it not in the interesting region or in the last segment ? Yes -> next root */
        if (ih <= -Idiv2 || ih >= Idiv2 - 1 || ih <= sg[nsg - 1].end) continue;
        /* Ok, we insert this one-point segment */
        sg[nsg].begin = sg[nsg].end = ih;
        g = compute_f (d, u, hl);
        g = lg2_abs (g, offset, scale);
        sg[nsg].f_begin = sg[nsg].f_end = g;
        S[ih] = (unsigned char) g;
      }
      else {
        /* It's a real (non an one-point) segment. */
        /* Is the root really far after the interesting zone ? Yes -> end of roots */
        if (ih > Idiv2 - 1 + SMART_NORM_INFLUENCE) break;

        /* Is the root really far before the interesting zone or is this segment is
           included (in all sides) in the previous segment ? Yes -> next root */
        if (ih < -Idiv2 - SMART_NORM_INFLUENCE ||
            ih <= sg[nsg - 1].end - SMART_NORM_STABILITY) continue;

        /* OK, we have the right side to compute, and maybe the left side */
        hr = hl;

        /* Special case for left side. Is the root :
           1. Near and before the interesting zone OR
           2. in the previous segment AND [the beginning of the previous segment
           is far enough of the root OR this previous segment is the first, so
              its beginning is the beginning of the interesting zone] ? */
        if (UNLIKELY(ih <= -Idiv2 || (sg[nsg - 1].end >= ih &&
                                      (nsg == 1 || sg[nsg - 1].begin + SMART_NORM_STABILITY <= ih)))) {
          /* OK, this segment and the previous could be fusionned on this left side.
             It's not 100% true in fact, but really very probable. */
          sg[nsg].begin   = sg[nsg - 1].begin;
          sg[nsg].f_begin = sg[nsg - 1].f_begin;
          fg = (ih >= -Idiv2 && ih < Idiv2) ? S[ih] : 0;
        } else {

          /* Here we compute the left side. -Idiv2 < ih <= Idiv2 - 1 + SMART_NORM_INFLUENCE */
          /* First, we compute the root itself */
          g = compute_f (d, u, hl);
          g = lg2_abs (g, offset, scale);
          fg = (unsigned char) g;
          if (LIKELY(ih < Idiv2)) S[ih] = fg;    /* Right guard for the write */

          /* The loop on the left side */
          for (f1 = fg, f2 = 0, cpt = 0; ; ) {
            --ih; hl -= 1.;
            if (UNLIKELY(ih <= -Idiv2)) {        /* Left guard: in fact ==, not <= */
              ASSERT(ih == -Idiv2);
              sg[nsg].begin = sg[0].begin;
              sg[nsg].f_begin = sg[0].f_begin;
              break;
            }
            g = compute_f (d, u, hl);
            g = lg2_abs (g, offset, scale);
            f1 = (unsigned char) g;
            if (LIKELY(ih < Idiv2)) S[ih] = f1;  /* Right guard for the write */
            if (LIKELY(f1 == f2)) {
              if (++cptf2id >= SMART_NORM_STABILITY) goto end_left;
            } else {
              cptf2id = 0;
              f2 = f1;
            }
            if (UNLIKELY(++cpt >= SMART_NORM_INFLUENCE)) {
            end_left:
              sg[nsg].begin = ih;
              sg[nsg].f_begin = g;
              break;
            }
          }
        }
        /* If the left side is totally out of the interesting region ? -> next root */
        if (sg[nsg].begin > Idiv2 - 1) continue;

        /* Now, the right side. */
        ih = (ssize_t) hr; /* -Idiv2-SMART_NORM_INFLUENCE <= ih <= Idiv2-1+SMART_NORM_INFLUENCE */

        /* Immediate right guard for -Idiv2-SMART_NORM_INFLUENCE <= ih < Idiv2-1 */
        if (UNLIKELY(ih >= Idiv2 - 1)) {
          ih = Idiv2 - 1;
          hr = (double) ih;
          g = compute_f (d, u, hr);
          g = lg2_abs (g, offset, scale);
        } else {
          /* The loop on the right side */
          for (f1 = fg, f2 = 0, cpt = 0;;) {
            ++ih; hr += 1.;
            g = compute_f (d, u, hr);
            g = lg2_abs (g, offset, scale);
            f1 = (unsigned char) g;
            if (LIKELY(ih > -Idiv2)) S[ih] = f1; /* Left guard for the write */
            if (UNLIKELY(ih >= Idiv2 - 1)) {     /* Right guard: in fact ==, not >= */
              ASSERT(ih == Idiv2 - 1);
              break;
            }
            if (LIKELY (f1 == f2)) {
              if (++cptf2id >= SMART_NORM_STABILITY) break;
            } else {
              cptf2id = 0;
              f2 = f1;
            }
            if (++cpt >= SMART_NORM_INFLUENCE) break;
          }
          /* If the right side is totally out of the interesting region ? -> next root */
          if (ih <= -Idiv2) continue;
        }
        sg[nsg].end = ih;
        sg[nsg].f_end = g;
        S[ih] = (unsigned char) g;
      }
      /* We have to do the possible fusions */
      while (nsg && sg[nsg - 1].end + 1 >= sg[nsg].begin) {
        if (sg[nsg - 1].begin > sg[nsg].begin) {
          sg[nsg - 1].begin = sg[nsg].begin; sg[nsg - 1].f_begin = sg[nsg].f_begin; }
        if (sg[nsg - 1].end   < sg[nsg].end  ) {
          sg[nsg - 1].end   = sg[nsg].end;   sg[nsg - 1].f_end   = sg[nsg].f_end;   }
        --nsg;
      }
      ++nsg;
    } /* End of the roots */

    /* The last sg end must be Idiv2 - 1 */
    if (LIKELY(sg[nsg - 1].end != Idiv2 - 1)) {
      /* We have to add a segment, except if the last ends at Idiv2 - 2 */
      g = compute_f (d, u, (double) (Idiv2 - 1));
      g = lg2_abs (g, offset, scale);
      S[Idiv2 - 1] = (unsigned char) g;
      if (LIKELY(sg[nsg - 1].end != Idiv2 - 2)) {
        sg[nsg].begin = Idiv2 - 1;
        sg[nsg].f_begin = g;
        ++nsg;
      }
      sg[nsg - 1].end = Idiv2 - 1;
      sg[nsg - 1].f_end = g;
    }

#if 0
    for (size_t k = 0; k < roots.size(); k++)
      fprintf (stderr, "Racine %zu (derivative=%u): %e - %u\n", k, roots[k].derivative, roots[k].value * J, J);
    for (size_t k = 0; k < nsg; k++)
      fprintf (stderr, "Segments %zu: (%d, %e) - (%d, %e)\n", k, sg[k].begin, sg[k].f_begin, sg[k].end, sg[k].f_end);
#endif
    /* Here, we have to set all the unset S values by polygonal approximation on all
       sg[x].end...sg[x+1].begin */
    poly_approx_on_S (S, d, u, original_scale, scale, offset, nsg, sg);
  } /* End of the line J; */

  /* Special ultra rare case. The correction of log2(F(0,0)) could be false,
     because the fast algorithm of log2 is not good in 0.0. */
  if (UNLIKELY(!beginJ)) *beginS = GUARD;
} /* True end of the function (finally!) */
/* }}} */

void lognorm_fill_rat_oldsmart (unsigned char *S, uint32_t N, int logI, double scale, cxx_double_poly const & fijd, const double *cexp2)
{
    double u0 = fijd->coeff[0];
    double u1 = fijd->coeff[1];

  /* The computation of the log2 by inttruncfastlog2 needs a value >= 1.
   * Here, in the classical rational initialization, the degree of the
   * used polynomial F(i,j) is hardcoded to 1. So, it's possible to have
   * F(i,j) = 0, even if i, j and the coefficients of F(i,j) are integers
   * and != 0.  So, I add 1.0 on all G values.  It's not useful to do a
   * fabs(G) here because the code uses always COMPUTE_Y(G) with G >= 0.
   */
#define COMPUTE_Y(G) lg2 ((G) + 1., add, scale)
  /* For internal debug: print all */
  // #define DEBUG_INIT_RAT 1
  double i;
  size_t ts;
  unsigned int inc;
  uint8_t oy, y;
  /* icc-14 gets confused between signbit from math.h and std::signbit ;
   * probably a bug, but fixing is easy. */
  using std::signbit;

  ASSERT_ALWAYS (u1 != 0.);

  LOGNORM_FILL_COMMON_DEFS();
  LOGNORM_COMMON_HANDLE_ORIGIN();

  const int Idiv2 = I >> 1;
  const double Idiv2_double = (double) Idiv2;
  const double Idiv2_double_minus_one = Idiv2_double - 1.;
  const double invu1 = 1./u1;

#ifdef DEBUG_INIT_RAT
  fprintf (stderr, "Begin: j=%u, u0=%.20e u1=%.20e, scale=%.20e, rac=%.20e\n", N, u0, u1, scale, u0 * N * (-invu1));
#endif


  scale *= 1./0x100000;
  double add = 0x3FF00000 - GUARD / scale;
  double u0J = u0 * j0;
  double d0_init = cexp2[((unsigned int)GUARD) - 1U];
  for (uint32_t J = j0 ; J < j1 ; J++, u0J += u0) {
    int int_i = -Idiv2;
    double g = u0J + u1 * int_i;
    double rac = u0J * (-invu1);
    double d0 = d0_init;
    double d1 = rac - d0 * rac;
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
      if (++int_i >= Idiv2) goto nextj;
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

void lognorm_oldsmart::fill(unsigned char * S, int N) const
{
    if (fijd->deg == 1) {
        lognorm_fill_rat_oldsmart (S, N, logI, scale, fijd, cexp2);
    } else {
        lognorm_fill_alg_oldsmart (S, N, logI, scale, fijd, roots);
        /* Adjust, because the "smart" code computes rubbish there
         */
        if (N == 0) {
            int i0 = - (1 << (logI-1));
            int c = S[1-i0];
            memset(S, 255, 1<<logI);
            S[1-i0]=c;
        }
    }
}

