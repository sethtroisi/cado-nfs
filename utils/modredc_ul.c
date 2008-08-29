#include "modredc_ul.h"
#include "modredc_ul_default.h"
#include "mod_ul_common.c"

#if defined(__GNUC__) && (__GNUC__ >= 4 || __GNUC__ >= 3 && __GNUC_MINOR__ >= 4)
/* Opteron prefers LOOKUP_TRAILING_ZEROS 1, 
   Core2 prefers LOOKUP_TRAILING_ZEROS 0 */
#ifndef LOOKUP_TRAILING_ZEROS
#define LOOKUP_TRAILING_ZEROS 1
#endif
#ifndef PRECONDITION_T
/* Sometimes faster, sometimes slower, doesn't seem to matter much on 
   average */
#define PRECONDITION_T 0
#endif
#define ctzl(x) __builtin_ctzl(x)
#define clzl(x) __builtin_clzl(x)
#else
/* If we have no ctzl(), we always use the table lookup */
#define LOOKUP_TRAILING_ZEROS 1
#define PRECONDITION_T 0
#endif

#ifdef T_HIST
unsigned int t_hist[256] = {
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};
#endif


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

  /* Let A = x*2^w, so we want the Montgomery representation of 1/x, 
     which is 2^w/x. We start by getting a = x */ 
  a = modredcul_get_ul (A, m);

#if PRECONDITION_T
/* If the inverse exists, the value of t is bounded below by log_2(a+b) - 1.
   Proof.
   
   Assume a is odd, gcd(a, b) = 1.
   
            { 0, a = b (implies a = b = 1)
   t(a,b) = { t(a, b/2)+1, a != b, b even
            { t(b, a), a > b, b odd
            { t(a, b-a) , a < b, b odd
            
   The last case implies a != b-a and b-a even, so can be substituted by
              t(a, (b-a)/2) + 1 , a < b, b odd
   In case 2, a+b  ->  a + b/2 and a + b/2 >= (a+b)/2.
   In case 4, a+b  ->  a + (b-a)/2 and a + (b-a)/2 = (a+b)/2.
   So each time t increases by 1, a+b drops by at most half. The process 
   stops when a = b = 1, i.e. a+b = 2. Hence, t >= log_2(a+b) - 1.

   Before the correction step, the result is 2^t/x, 
   where t >= log_2(a+b)-1. We divide here by 2^(w-ceil(log_2(a+b)-1))
   and init t = -ceil(log_2(a+b)-1), so that the result before the
   correction step is 2^(w+t)/x with t >= 0. This way we can do the 
   correction step via a single REDC of width t. */

  ASSERT (b > 1UL);
  t = clzl (b); /* Since a will change again, we estimate just 
                    log_2(b)-1 <= log_2(a+b)-1. */
  t++; /* Now 1 <= t = w - ceil(log_2(b)-1) <= 63 */
  
  {
    unsigned long tlow, thigh;
    /* Necessarily t < LONG_BIT, so the shift is ok */
    /* Doing a left shift first and then a full REDC needs a modular addition
       at the end due to larger summands and thus is probably slower */
    tlow = ((a * m[0].invm) & ((1UL << t) - 1UL)); /* tlow <= 2^t-1 */
    ularith_mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0].m); /* thigh:tlow <= m*(2^t-1) */
    ularith_add_ul_2ul (&tlow, &thigh, a); /* thigh:tlow <= m*2^t-1 (since u<m) */
    /* Now the low t bits of tlow are 0 */
    ASSERT_EXPENSIVE ((tlow & ((1UL << t) - 1UL)) == 0UL);
    modredcul_shrd (&tlow, thigh, t);
    a = tlow;
    ASSERT_EXPENSIVE ((thigh >> t) == 0UL && a < m[0].m);
  }
  t -= LONG_BIT;
#else
  /* Alternatively, we simply set a = x/2^w and t=0. The result before 
     correction will be 2^(w+t)/x so we have to divide by t, which
     may be >64, so we may have to do a full and a variable width REDC. */
  a = modredcul_get_ul (&a, m);
  /* Now a = x/2^w */
  t = 0;
#endif

  u = 1UL; v = 0UL;

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
  t += lsh;
#endif
  /* v <<= lsh; ??? v is 0 here */

  // Here a and b are odd, and a < b
  do {
    /* Here, a and b are odd, 0 < a < b, u is odd and v is even */
    ASSERT ((u & 1UL) == 1UL && (v & 1UL) == 0UL);
    do {
      b -= a; v += u;
#if LOOKUP_TRAILING_ZEROS
      do {
	lsh = trailing_zeros [(unsigned char) b];
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
    /* Here, a and b are odd, 0 < b =< a, u is even and v is odd */

    if (a == b)
      break;

    /* Here, a and b are odd, 0 < b < a, u is even and v is odd */
    ASSERT ((u & 1UL) == 0UL && (v & 1UL) == 1UL);
    do {
      a -= b; u += v;
#if LOOKUP_TRAILING_ZEROS
      do {
	lsh = trailing_zeros [(unsigned char) a];
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
    /* Here, a and b are odd, 0 < a =< b, u is odd and v is even */
  } while (a != b);
  
  if (a != 1UL) /* Non-trivial GCD */
    return 0;

  ASSERT (t >= 0);

#ifdef T_HIST
  if (t >= 0 && t < 255)
    t_hist[t]++;
  else
    t_hist[255]++;
#endif

  /* Here, u = 2^w * 2^t / x. We want 2^w / x. */

  /* Here, the inverse of a is u/2^t mod b. To do the division by t,
     we use a variable-width REDC. We want to add a multiple of m to u
     so that the low t bits of the sum are 0 and we can right-shift by t
     with impunity. */
  if (t >= LONG_BIT)
    {
      unsigned long tlow, thigh;
      tlow = u * m[0].invm; /* tlow <= 2^w-1 */
      ularith_mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0].m); /* thigh:tlow <= (2^w-1)*m */
      u = thigh + ((u != 0UL) ? 1UL : 0UL);
      /* thigh:tlow + u < (2^w-1)*m + m < 2^w*m. No correction necesary */
      t -= LONG_BIT;
    }

  ASSERT (t < LONG_BIT);
  if (t > 0)
    {
      unsigned long tlow, thigh;
      /* Necessarily t < LONG_BIT, so the shift is ok */
      /* Doing a left shift first and then a full REDC needs a modular addition
	 at the end due to larger summands and thus is probably slower */
      tlow = ((u * m[0].invm) & ((1UL << t) - 1UL)); /* tlow <= 2^t-1 */
      ularith_mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0].m); /* thigh:tlow <= m*(2^t-1) */
      ularith_add_ul_2ul (&tlow, &thigh, u); /* thigh:tlow <= m*2^t-1 (since u<m) */
      /* Now the low t bits of tlow are 0 */
      ASSERT_EXPENSIVE ((tlow & ((1UL << t) - 1UL)) == 0UL);
      ularith_shrd (&tlow, thigh, t);
      u = tlow;
      ASSERT_EXPENSIVE ((thigh >> t) == 0UL && u < m[0].m);
    }

  r[0] = u;
  return 1;
}

