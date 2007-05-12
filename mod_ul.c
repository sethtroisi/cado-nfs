/* Some functions for modular arithmetic with residues and modulus in
   unsigned long variables. Due to inlining, this file must be included
   in the caller's source code with #include */

/* Naming convention: all function start with modul, for 
   MODulus Unsigned Long, followed by underscore, functionality of function
  (add, mul, etc), and possibly underscore and specification of what argument
  types the function takes (_ul, etc).
  There are typedef's that rename all functions to mod_* instead of 
  modul_*, which one day might become an automatic renaming scheme so that
  different modulus sizes can be used simply by #including different mod_*.c
  files, but without changing anything else in the source code. */

#ifndef ASSERT
 #ifdef WANT_ASSERT
  #include <assert.h>
  #define ASSERT(x) assert(x)
 #else
  #define ASSERT(x)
 #endif
#define MOD_UL_ASSERT
#endif

#ifndef GNUC_ATTRIBUTE_UNUSED
#define __GNUC_ATTRIBUTE_UNUSED__ __attribute__ ((unused))
#define MOD_UL_GNUC_ATTRIBUTE_UNUSED
#endif

/* Define the default names and types for arithmetic with these functions */
#define mod_init             modul_init
#define mod_init_set0        modul_init_set0
#define mod_clear            modul_clear
#define mod_set              modul_set
#define mod_set_ul           modul_set_ul
#define mod_set_ul_reduced   modul_set_ul_reduced
#define mod_initmod_ul       modul_initmod_ul
#define mod_clearmod         modul_clearmod
#define mod_get_ul           modul_get_ul
#define mod_cmp              modul_cmp
#define mod_is0              modul_is0
#define mod_add              modul_add
#define mod_add_ul           modul_add_ul
#define mod_sub              modul_sub
#define mod_sub_ul           modul_sub_ul
#define mod_neg              modul_neg
#define mod_mul              modul_mul
#define mod_div2             modul_div2
#define mod_div3             modul_div3
#define mod_inv              modul_inv
#define mod_jacobi           modul_jacobi

typedef unsigned long residueul[1];
typedef unsigned long modulusul[1];
typedef residueul residue;
typedef modulusul modulus;

static inline void
modul_init (residueul r __GNUC_ATTRIBUTE_UNUSED__, 
            modulusul m __GNUC_ATTRIBUTE_UNUSED__)
{
  return;
}

static inline void
modul_init_set0 (residueul r, modulusul m __GNUC_ATTRIBUTE_UNUSED__)
{
  r[0] = 0UL;
  return;
}

static inline void
modul_clear (residueul r __GNUC_ATTRIBUTE_UNUSED__,
              modulusul m __GNUC_ATTRIBUTE_UNUSED__)
{
  return;
}

static inline void
modul_set (residueul r, residueul s, modulusul m __GNUC_ATTRIBUTE_UNUSED__)
{
  r[0] = s[0];
}

static inline void
modul_set_ul (residueul r, unsigned long s, modulusul m)
{
  r[0] = s % m[0];
}

/* Sets the residue to the class represented by the integer s. Assumes that
   s is reduced (mod m), i.e. 0 <= s < m */

static inline void
modul_set_ul_reduced (residueul r, unsigned long s, 
                      modulusul m __GNUC_ATTRIBUTE_UNUSED__)
{
  ASSERT (s < m[0]);
  r[0] = s;
}

static inline void
modul_initmod_ul (modulusul r, unsigned long s)
{
  r[0] = s;
}

static inline void
modul_clearmod (modulusul m __GNUC_ATTRIBUTE_UNUSED__)
{
  return;
}

static inline unsigned long
modul_get_ul (residueul s, modulusul m __GNUC_ATTRIBUTE_UNUSED__)
{
  return s[0];
}

static inline int
modul_cmp (residueul a, residueul b, modulusul m __GNUC_ATTRIBUTE_UNUSED__)
{
  return (a[0] < b[0]) ? -1 : ((a[0] == b[0]) ? 0 : 1);
}

static inline int
modul_is0 (residueul a, modulusul m __GNUC_ATTRIBUTE_UNUSED__)
{
  return (a[0] == 0UL);
}

static inline void
modul_add (residueul r, residueul a, residueul b, modulusul m)
{
  ASSERT(a[0] < m[0] && b[0] < m[0]);
#ifdef MODTRACE
  printf ("addmodul: a = %lu, b = %lu", a[0], b[0]);
#endif
  r[0] = (a[0] >= m[0] - b[0]) ? (a[0] - (m[0] - b[0])) : (a[0] + b[0]);
#ifdef MODTRACE
  printf (", r = %lu\n", r[0]);
#endif
}

static inline void
modul_add_ul (residueul r, residueul a, unsigned long b, modulusul m)
{
  ASSERT(a[0] < m[0] && b < m[0]);
#ifdef MODTRACE
  printf ("addmodul_ul: a = %lu, b = %lu", a[0], b);
#endif
  r[0] = (a[0] >= m[0] - b) ? (a[0] - (m[0] - b)) : (a[0] + b);
#ifdef MODTRACE
  printf (", r = %lu\n", r[0]);
#endif
}

static inline void
modul_sub_ul (residueul r, residueul a, unsigned long b, modulusul m)
{
  ASSERT(a[0] < m[0] && b < m[0]);
#ifdef MODTRACE
  printf ("submodul_ul: a = %lu, b = %lu", a[0], b);
#endif
  r[0] = (a[0] >= b) ? (a[0] - b) : (a[0] + (m[0] - b));
#ifdef MODTRACE
  printf (", r = %lu\n", r[0]);
#endif
}

static inline void
modul_neg (residueul r, residueul a, modulusul m)
{
  if (a[0] == 0)
    r[0] = a[0];
  else
    r[0] = m[0] - a[0];
}

static inline void
modul_sub (residueul r, residueul a, residueul b, modulusul m)
{
  ASSERT(a[0] < m[0] && b[0] < m[0]);
#ifdef MODTRACE
  printf ("submod_ul: a = %lu, b = %lu", a[0], b[0]);
#endif
  r[0] = (a[0] < b[0]) ? (a[0] + (m[0] - b[0])) : (a[0] - b[0]);
#ifdef MODTRACE
  printf (", r = %lu\n", r[0]);
#endif
}

static inline void
modul_mul (residueul r, const residueul a, const residueul b, 
           const modulusul m)
{
  unsigned long _a = a[0], _r;
#if defined(MODTRACE)
  printf ("mulmod_ul: a = %lu, b = %lu", a[0], b[0]);
#endif
#if defined(__x86_64__) && defined(__GNUC__)
  __asm__ ( "mulq %2\n\t"
	    "divq %3"
	    : "=&d" (_r), "+a" (_a)
	    : "rm" (b[0]), "rm" (m[0])
	    : "cc");
#elif defined(__i386__) && defined(__GNUC__)
  __asm__ ( "mull %2\n\t"
	    "divl %3"
	    : "=&d" (_r), "+a" (_a)
	    : "rm" (b[0]), "rm" (m[0])
	    : "cc");
#else
  if (sizeof (unsigned long long) >= 2 * sizeof (unsigned long))
    {
      unsigned long long t;
      t = (unsigned long long) a[0] * (unsigned long long) b[0];
      _r = t % (unsigned long long) m[0];
    }
  else
    {
      /* Write a long product in four multiplies. TBD. */
      abort();
    }
#endif
#if defined(MODTRACE)
  printf (", r = %lu\n", _r);
#endif
  r[0] = _r;
}

static inline void
modul_div2 (residueul r, residueul a, modulusul m)
{
  ASSERT(m[0] % 2 != 0);
  r[0] = (a[0] % 2 == 0) ? (a[0] / 2) : (a[0] / 2 + m[0] / 2 + 1);
}

void
modul_div3 (residueul r, residueul a, modulusul m)
{
  const unsigned long a3 = a[0] % 3;
  const unsigned long m3 = m[0] % 3;
#ifdef WANT_ASSERT
  residueul t;
#endif

  ASSERT(m3 != 0);
  if (a3 == 0)
    r[0] = a[0] / 3;
  else if (a3 + m3 == 3) /* Hence a3 == 1, m3 == 2 or a3 == 3, m3 == 1 */
    r[0] = a[0] / 3 + m[0] / 3 + 1;
  else /* a3 == 1, m3 == 1 or a3 == 2, m3 == 2 */
    r[0] = m[0] - (m[0] - a[0]) / 3;

#ifdef WANT_ASSERT
  modul_add (t, r, r, m);
  modul_add (t, t, r, m);
  assert (t[0] == a[0]);
#endif
}

/* Put 1/s (mod t) in r and return 1 if s is invertible, 
   or return 0 if s is not invertible */
static int
modul_inv (residueul r, residueul s, modulusul t)
{
  long u1, v1;
  unsigned long q, u2, v2;

  ASSERT (t[0] > 0);
  ASSERT (s[0] < t[0]);

  if (s[0] == 0)
    {
      r[0] = 0; /* Not invertible */
      return 0;
    }

  if (s[0] == 1)
    {
      r[0] = 1;
      return 1;
    }

#if 0
  u1 = 1;
  v1 = 0;
  u2 = s;
  v2 = t;

  q = s / t; /* == 0, due to s < t */
  u1 = u1 - (long) q * v1; /* == u1 - 0 == u1 == 1*/
  u2 = u2 - q * v2; /* == u2 - 0 == u2 == s */

  q = v2 / u2; /* == t / s <= t / 2*/
  v1 = v1 - (long) q * u1; /* == v1 - q * 1 == 0 - q == -q */
  v2 = v2 - q * u2; /* == v2 - q * s == t - q * s == t % s */
#endif

  u1 = 1UL;
  u2 = s[0];
  v1 = - (t[0] / s[0]); /* No overflow, since s >= 2 */
  v2 = t[0] % s[0];

  while (v2 != 0)
    {
      /* unroll twice and swap u/v */
      q = u2 / v2;
      u1 = u1 - (long) q * v1;
      u2 = u2 - q * v2;

      if (u2 == 0)
        {
          u1 = v1;
          u2 = v2;
          break;
        }

      q = v2 / u2;
      v1 = v1 - (long) q * u1;
      v2 = v2 - q * u2;
    }

  if (u2 != 1)
    {
      /* printf ("s=%lu t=%lu found %lu\n", s[0], t[0], u2); */
      r[0] = 0; /* non trivial gcd */
      return 0;
    }

  if (u1 < 0)
    u1 = u1 - t[0] * (-u1 / t[0] - 1);

#ifdef WANT_ASSERT
  modul_mul (&u2, (unsigned long *)&u1, s, t);
  ASSERT(u2 == 1);
#endif

  r[0] = (unsigned long) u1;
  return 1;
}

__GNUC_ATTRIBUTE_UNUSED__ static int
modul_jacobi (residueul a_par, modulusul m_par)
{
  unsigned long a, m, s;
  int t = 1;

  ASSERT (a_par[0] < m_par[0]);
  a = a_par[0];
  m = m_par[0];
  
  while (a != 0)
  {
    while (a % 2 == 0)
    {
      a /= 2;
      if (m % 8 == 3 || m % 8 == 5)
        t = -t;
    }
    s = a; a = m; m = s; /* swap */
    if (a % 4 == 3 && m % 4 == 3)
      t = -t;
    a %= m;
  }
  if (m != 1)
    t = 0;
  
#ifdef MODTRACE
  printf ("kronecker(%lu, %lu) == %d\n", a_par[0], m_par[0], t);
#endif
  return t;
}

/* Try to restore the macro defines to the same state they had been in */
#ifdef MOD_UL_ASSERT
#undef ASSERT
#undef MOD_UL_ASSERT
#endif
#ifdef MOD_UL_GNUC_ATTRIBUTE_UNUSED
#undef __GNUC_ATTRIBUTE_UNUSED__
#endif
