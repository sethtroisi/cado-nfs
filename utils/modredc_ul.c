#include "modredc_ul.h"

#include "mod_ul_common.c"

/* Opteron prefers LOOKUP_TRAILING_ZEROS 1, 
   Core2 prefers LOOKUP_TRAILING_ZEROS 0 */
#define LOOKUP_TRAILING_ZEROS 1
#define ctzl(x)         __builtin_ctzl(x)

int
modredcul_inv (residue_t r, const residue_t A, const modulusredcul_t m) 
{
#if LOOKUP_TRAILING_ZEROS
  static const unsigned char trailing_zeros[256] = 
    {8,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     7,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0};
#endif
  unsigned long a, b = m[0].m, u, v;
  int t, lsh;

  ASSERT (A[0] < b);
  ASSERT (b & 1UL);

  if (A[0] == 0UL)
    return 0;

  /* Do two REDC so the result is 1/a*2^w */
  a = modredcul_get_ul (A, m);
  a = modredcul_get_ul (&a, m);

  u = 1UL; v = 0UL; t = 0;

  // make a odd
#if LOOKUP_TRAILING_ZEROS
  do {
    lsh = trailing_zeros [(unsigned char) a];
    a >>= lsh;
    t += lsh;
  } while (lsh == 8);
#else
  lsh = ctzl(a);
  a >>= lsh;
  t = lsh;
#endif
  /* v <<= lsh; ??? v is 0 here */

  // Here a and b are odd, and a < b
  do {
    do {
      b -= a; v += u;
#if LOOKUP_TRAILING_ZEROS
      do {
	lsh = trailing_zeros [(unsigned char) b];
	ASSERT_EXPENSIVE (lsh > 0);
	b >>= lsh;
	t += lsh;
	u <<= lsh;
      } while (lsh == 8);
#else
      lsh = ctzl(b);
      ASSERT_EXPENSIVE (lsh > 0);
      b >>= lsh;
      t += lsh;
      u <<= lsh;
#endif
    } while (b > a); /* ~50% branch taken :( */
    if (a == b)
      break;
    do {
      a -= b; u += v;
#if LOOKUP_TRAILING_ZEROS
      do {
	lsh = trailing_zeros [(unsigned char) a];
	ASSERT_EXPENSIVE (lsh > 0);
	a >>= lsh;
	t += lsh;
	v <<= lsh;
      } while (lsh == 8);
#else
      lsh = ctzl(a);
      ASSERT_EXPENSIVE (lsh > 0);
      a >>= lsh;
      t += lsh;
      v <<= lsh;
#endif
    } while (b < a); /* about 50% branch taken :( */
  } while (a != b);
  
  if (a != 1UL) /* Non-trivial GCD */
    return 0;

  /* TODO: prove that t<128 */
  ASSERT (t < 128);

  /* Here, the inverse of a is u/2^t mod b. To do the division by t,
     we use a variable-width REDC. We want to add a multiple of m to u
     so that the low t bits of the sum are 0 and we can right-shift by t
     with impunity. */
  if (t >= 64)
    {
      unsigned long tlow, thigh;
      tlow = u * m[0].invm; /* tlow <= 2^w-1 */
      modredcul_mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0].m); /* thigh:tlow <= (2^w-1)*m */
      u = thigh + ((u != 0UL) ? 1UL : 0UL);
      /* thigh:tlow + u < (2^w-1)*m + m < 2^w*m. No correction necesary */
      t -= 64;
    }

  if (t > 0)
    {
      unsigned long tlow, thigh;
      /* Necessarily t < 64, so the shift is ok */
      tlow = ((u * m[0].invm) & ((1UL << t) - 1UL)); /* tlow <= 2^t-1 */
      modredcul_mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0].m); /* thigh:tlow <= m*(2^t-1) */
      modredcul_add_ul_2ul (&tlow, &thigh, u); /* thigh:tlow <= m*2^t-1 (since u<m) */
      /* Now the low t bits of tlow are 0 */
      ASSERT_EXPENSIVE ((tlow & ((1UL << t) - 1UL)) == 0UL);
      __asm__ ("#invmod shrq here\n\t"
               "shrdq %1, %0\n": 
               "+r" (tlow), "+r" (thigh) :
               "c" (t)
              ); /* tlow <= (m*2^t-1) / 2^t <= m-1 */
      u = tlow;
      ASSERT_EXPENSIVE ((thigh >> t) == 0UL && u < m[0].m);
    }

  r[0] = u;
  return 1;
}
