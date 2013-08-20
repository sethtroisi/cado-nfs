#ifndef LAS_ARITH_H_
#define LAS_ARITH_H_

#include "cado.h"
#include <stdint.h>

#include "fb.h"
#include "utils.h"
#include "las-config.h"
#include "utils/misc.h" /* ctzl */

#ifdef __cplusplus
extern "C" {
#endif

// Redc_32 based on 64-bit arithmetic
// Assume:
//   * p is an odd prime < 2^32. FIXME: p < 2^31? (see below)
//   * invp is -1/p mod 2^32.
//   * x is some integer in [0, 2^32*p[
// Compute:
//   * x/2^32 mod p as an integer in [0, p[
static inline uint32_t
redc_u32(const uint64_t x, const uint32_t p, const uint32_t invp)
{
    uint32_t t = x*invp;                            /* t = x * invp mod 2^32 */
    uint32_t y = (x + (uint64_t)t*(uint64_t)p) >> 32; /* x + t*p is bounded
            by 2^32*p-1+(2^32-1)*p < 2*2^32*p: we might want p < 2^31 so that
            there is no overflow */
    if (y >= p)
        y -= p;
    ASSERT(y<p);
    return y;
}

// Signed redc_32 based on 64-bit arithmetic
// Assume:
//   * p is an odd prime < 2^32.
//   * invp is -1/p mod 2^32.
//   * x is some signed integer in ]-2^32*p, 2^32*p[
// Compute:
//   * x/2^32 mod p as an integer in [0, p[
static inline uint32_t
redc_32(const int64_t x, const uint32_t p, const uint32_t invp)
{
    uint32_t t = ((uint64_t)x)*invp;
    uint32_t y = (x + (uint64_t)t*(uint64_t)p) >> 32;
    // might be too large by p, or too small by p.
    if ((int32_t)y < 0)
        y += p;
    else if (y >= p)
        y -= p;
#ifndef NDEBUG
    if (UNLIKELY(y >= p)) {
        fprintf(stderr, "BUG in redc_32. x = %" PRId64
                " p = %u, invp = %u, y = %d\n", x, p, invp, y);
        if (x < 0) {
            fprintf(stderr, "x/2^32 = -%"PRId64, (-x)>>32);
            if (((-x)>>32) < p) {
                fprintf(stderr, ", within bounds\n");
            } else {
                fprintf(stderr, ", OUT OF BOUNDS\n");
            }
        } else {
            fprintf(stderr, "x/2^32 = -%"PRId64, x>>32);
            if ((x>>32) < p) {
                fprintf(stderr, ", within bounds\n");
            } else {
                fprintf(stderr, ", OUT OF BOUNDS\n");
            }
        }
        ASSERT(0);
        /* TODO: Fall back to safer code ? */
    }
#endif
    return y;
}

#ifdef   HAVE_GCC_STYLE_AMD64_INLINE_ASM
#define HAVE_redc_64
/* This does a full mul, and should grok slightly larger bounds than the
 * version above. Presumably, as long as x fits within 64 bits, (63 bits,
 * with sign), things should be ok. (TODO: check).
 */
static inline uint64_t
redc_64(const int64_t x, const uint32_t p, const uint64_t invp)
{
    int64_t t = ((uint64_t)x)*invp;
    uint64_t u;
    /* Need the high part of a 64x64 mul */
    __asm__ (
            "    mulq    %[p]\n"
            "    addq    %[x], %%rax\n"
            "    adcq    $0, %%rdx\n"
            : "=&d" (u)
            : "a" (t), [x] "rm" (x), [p] "r" ((uint64_t) p));
    /* As per the early clobber on rdx, x can't be put in there.
     * Furthermore, since t goes to rax, x doesn't go there either. Thus
     * it is reasonable to assume that it is still unmodified after the
     * asm block */
    u-=x<0;     /* FIXME. Can I get around this ? */
    if ((int64_t)u < 0)
        u += p;
    else if (u >= p)
        u -= p;
    return u;
}
#endif

MAYBE_UNUSED
static inline fbprime_t
invmod_po2 (fbprime_t n)
{
  fbprime_t r;

  ASSERT (n % 2 != 0);

  r = (3 * n) ^ 2;
  r = 2 * r - r * r * n;
  r = 2 * r - r * r * n;
  r = 2 * r - r * r * n;
  if (sizeof (fbprime_t) > 4)
    r = 2 * r - r * r * n;

  return r;
}


NOPROFILE_INLINE int
invmod (uint64_t *pa, uint64_t b)
{
  modulusul_t m;
  residueul_t r;
  int rc;
  modul_initmod_ul (m, b);
  modul_init (r, m);
  modul_set_ul (r, *pa, m); /* With mod reduction */
  if ((rc = modul_inv(r, r, m)))
    *pa = modul_get_ul (r, m);
  modul_clear (r, m);
  modul_clearmod (m);
  return rc;
}


/* TODO: this is a close cousin of modredcul_inv, but the latter does
 * 64-bit redc */

// Compute 2^32/a mod b for b odd,
// and 1/a mod b for b even, by binary xgcd.
// a must be less than b.
// a is modified in place
// return 1 on succes, 0 on failure
NOPROFILE_INLINE int
invmod_redc_32(uint64_t *pa, uint64_t b) {
  uint64_t a, u, v, fix, p = b;
  unsigned int t, lsh;
  a = *pa;
  if (UNLIKELY(!a)) return 0; /* or we get infinite loop */
  if (UNLIKELY(!(b & 1))) return invmod(pa, b); 
  fix = (b+1)>>1;
  ASSERT (a < b);
  u = 1; v = 0; t = 0;

  // make a odd
  lsh = ctzl(a); a >>= lsh; t += lsh;

  // Here a and b are odd, and a < b
  /* T1 & T2 in x86 asm has 8 instructions */
#define T1 b -= a; v += u; lsh = ctzl(b); b >>= lsh; t += lsh; u <<= lsh; if (a >= b) break
#define T2 a -= b; u += v; lsh = ctzl(a); a >>= lsh; t += lsh; v <<= lsh; if (b >= a) break
  do {
    for (;;) { T1; T1; T1; T1; }
    if (a == b) break;
    for (;;) { T2; T2; T2; T2; }
  } while (a != b);
  if (a != 1) return 0;
  
  // Here, the inverse of a is u/2^t mod b.
  /* T3 in x86 asm has 3 instructions; T4 has 5 */
#ifdef __x86_64
#define T3 do {								\
    uint64_t addq;							\
    __asm__ ( "shr $1,%1\n lea (%1,%2),%0\n cmovcq %0, %1\n" :		\
	      "=r"(addq), "+r"(u) : "r"(fix));				\
      } while (0)
#else
#define T3 do { unsigned char sig = (unsigned char) u; u >>= 1; if (sig & 1) u += fix; } while (0)
#endif
#define T4 do { u <<= 1; if (u >= p) u -= p; } while (0)
  
#if 0
  for (; t > 32; --t) T3;
  for (; t < 32; ++t) T4;
#else
  if (t > 32) {
    for (t = t - 33; t > 7; t -= 8) { T3; T3; T3; T3; T3; T3; T3; T3; }
    switch (t) {
    case 7: T3; case 6: T3; case 5: T3; case 4: T3; case 3: T3; case 2: T3; case 1: T3; case 0: T3;
    }
  }
  else if (t < 32) {
    for (t = 31 - t; t > 7; t -= 8) { T4; T4; T4; T4; T4; T4; T4; T4; }
    switch (t) {
    case 7: T4; case 6: T4; case 5: T4; case 4: T4; case 3: T4; case 2: T4; case 1: T4; case 0: T4;
    }
  }
#endif
  *pa = u;
  return 1;
}

/* Only used for together with redc_64 */
NOPROFILE_INLINE int
invmod_redc_64(uint64_t *pa, uint64_t b)
{
    uint64_t a, u, v, fix, p = b;
    int t, lsh;
    a = *pa;

    if (UNLIKELY(*pa == 0)) return 0; /* or we get infinite loop */
    if (UNLIKELY(b % 2UL == 0)) return invmod(pa, b);

    fix = (b+1)>>1;

    ASSERT (a < b);

    u = 1; v = 0; t = 0;

    // make a odd
    lsh = ctzl(a);
    a >>= lsh;
    t += lsh;
    /* v <<= lsh; ??? v is 0 here */

    // Here a and b are odd, and a < b
    do {
        do {
            b -= a; v += u;
            lsh = ctzl(b);
            b >>= lsh;
            t += lsh;
            u <<= lsh;
        } while (a<b);
        if (UNLIKELY(a == b))
            break;
        do {
            a -= b; u += v;
            lsh = ctzl(a);
            a >>= lsh;
            t += lsh;
            v <<= lsh;
        } while (b < a);
    } while (a != b);
    if (a != 1)
        return 0;

    // Here, the inverse of a is u/2^t mod b.
    while (t > 64) {
        unsigned long sig = u & 1UL;
        u >>= 1;
        u += fix & -sig;
        --t;
    }
    while (t < 64) {
        u <<= 1;
        if (u >= p)
            u -= p;
        t ++;
    }
    *pa = u;
    return 1;
}

fbprime_t is_prime_power(fbprime_t q);

#ifdef __cplusplus
}
#endif

#endif	/* LAS_ARITH_H_ */
