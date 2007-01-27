/* Usage: find_primes_hit P B1 sigma - find all primes <= P that are found with
   ECM(B1,sigma) */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "getprime.c"

#if 1
#include "gmp.h"
#include "gmp-impl.h"
#include "longlong.h"
/* return x*y mod n */
static inline mp_limb_t
mulmod (mp_limb_t x, mp_limb_t y, mp_limb_t n)
{
  mp_limb_t h, l, q, r;

  umul_ppmm  (h, l, x, y);
#if 0
  if (h >= n)
    {
      printf ("x=%lu y=%lu h=%lu l=%lu n=%lu\n", x, y, h, l, n);
      abort ();
    }
#endif
  udiv_qrnnd (q, r, h, l, n);
#if 0
  if (r >= n)
    {
      printf ("error in mulmod, r=%lu\n", r);
      exit (1);
    }
#endif
  return r;
}
#else /* 64-bit architecture */
#define    mulmod(x,y,n)   ((x*y)   % n)
#endif

/* return 1/s mod t, or 0 if non-trivial gcd found */
static unsigned long
inverse (unsigned long s, unsigned long t)
{
  long u1, v1;
  unsigned long q, u2, v2;
 
  u1 = 1;
  v1 = 0;
  u2 = s;
  v2 = t;
  
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
      /* printf ("s=%lu t=%lu found %lu\n", s, t, u2); */
      return 0; /* non trivial gcd */
    }
  
  if (u1 < 0)
    u1 = u1 - t * (-u1 / t - 1);

  return (unsigned long) u1;
}

/* safe u-v mod p, assumes 0 <= u, v < p */
#define submod(u,v,p) ((u >= v) ? u - v : (p - (v - u)))
#define addmod(u,v,p) ((u + v >= u) ? ((u + v) % p) : (u + v) - p)
#define div2mod(u,p)  (((u & 1) == 0) ? (u / 2) : (addmod(u,p,p) / 2))

/* adds Q=(x2:z2) and R=(x1:z1) and puts the result in (x3:z3),
     using 6 muls (4 muls and 2 squares), and 6 add/sub.
   One assumes that Q-R=P or R-Q=P where P=(x:z).
     - n : number to factor
   Modifies: x3, z3, u, v, w.
   (x3,z3) may be identical to (x2,z2) and to (x,z)
*/
static void
add3 (ulong *x3, ulong *z3, ulong x2, ulong z2, ulong x1, ulong z1, 
      ulong x, ulong z, ulong n)
{
  unsigned long u, v, w;

  u = submod (x2, z2, n);
  v = addmod (x1, z1, n);
  u = mulmod (u, v, n);
  w = addmod (x2, z2, n);
  v = submod (x1, z1, n);
  v = mulmod (w, v, n);
  w = addmod (u, v, n);
  v = submod (u, v, n);
  w = mulmod (w, w, n);
  v = mulmod (v, v, n);

  u = x; /* save x */
  *x3 = mulmod (w, z, n);
  *z3 = mulmod (u, v, n);
}

/* computes 2P=(x2:z2) from P=(x1:z1), with 5 muls (3 muls and 2 squares)
   and 4 add/sub.
     - p : number to factor
     - b : (a+2)/4 mod n
*/
static void
duplicate (ulong *x2, ulong *z2, ulong x1, ulong z1, ulong p, ulong b)
{
  ulong u, v, w;

  u = addmod (x1, z1, p);
  u = mulmod (u, u, p);
  v = submod (x1, z1, p);
  v = mulmod (v, v, p);
  *x2 = mulmod (u, v, p);
  w = submod (u, v, p);
  u = mulmod (w, b, p);
  u = addmod (u, v, p);
  *z2 = mulmod (w, u, p);
}

/* (x:z) <- e*(x:z) (mod p)
   Assumes e >= 5.
*/
static void
ecm_mul_ui (ulong *x, ulong *z, ulong e, ulong p, ulong b)
{
  unsigned long l, n, x1, z1, x2, z2;

  e --;

  /* compute number of steps needed: we start from (1,2) and go from
     (i,i+1) to (2i,2i+1) or (2i+1,2i+2) */
  for (l = e, n = 0; l > 1; n ++, l /= 2);

  /* start from P1=P, P2=2P */
  x1 = *x;
  z1 = *z;
  duplicate (&x2, &z2, x1, z1, p, b);

  while (n--)
    {
      if ((e >> n) & 1) /* (i,i+1) -> (2i+1,2i+2) */
        {
          /* printf ("(i,i+1) -> (2i+1,2i+2)\n"); */
          add3 (&x1, &z1, x2, z2, x1, z1, *x, *z, p);
          duplicate (&x2, &z2, x2, z2, p, b);
        }
      else /* (i,i+1) -> (2i,2i+1) */
        {
          /* printf ("(i,i+1) -> (2i,2i+1)\n"); */
          add3 (&x2, &z2, x1, z1, x2, z2, *x, *z, p);
          duplicate (&x1, &z1, x1, z1, p, b);
        }
    }
  
  *x = x2;
  *z = z2;
}

static int
ecm (unsigned long p, double B1, unsigned long sigma)
{
  /* small primes */
  static double primes[] = {5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 0};
  unsigned long p0, t, u, v, A, b, x, z, xB, zB;
  double r, *s;

  p0 = p;
#ifdef UDIV_NEEDS_NORMALIZATION
  while (p < (p << 1))
    p <<= 1;
#endif

  /* compute b, x */
  v = mulmod (4, sigma, p);
  t = mulmod (sigma, sigma, p);
  u = submod (t, 5, p);
  t = mulmod (u, u, p);
  x = mulmod (t, u, p);
  t = mulmod (v, v, p);
  z = mulmod (t, v, p);
  t = mulmod (x, v, p);
  b = mulmod (t, 4, p);
  t = mulmod (u, 3, p);
  u = submod (v, u, p);
  v = addmod (t, v, p);
  t = mulmod (u, u, p);
  u = mulmod (t, u, p);
  A = mulmod (u, v, p);
  v = mulmod (b, z, p);
  u = inverse (v, p0);
  if (u == 0) /* non trivial gcd */
    return 1;
  v = mulmod (u, b, p);
  x = mulmod (x, v, p);
  v = mulmod (u, z, p);
  t = mulmod (A, v, p);
  A = submod (t, 2, p);

  /* now start ecm */
  z = 1;
  b = addmod (A, 2, p);
  b = div2mod (b, p);
  b = div2mod (b, p);

  /*  printf ("x=%lu b=%lu\n", x % p0, b % p0); */

  /* prime 2 */
  for (r = 2.0; r <= B1; r *= 2.0)
    duplicate (&x, &z, x, z, p, b);

  /*  printf ("2: x=%lu z=%lu\n", x % p0, z % p0); */

  /* prime 3 */
  for (r = 3.0; r <= B1; r *= 3.0)
    {
      duplicate (&xB, &zB, x, z, p, b);
      add3 (&x, &z, x, z, xB, zB, x, z, p);
    }

  /*  printf ("3: x=%lu z=%lu\n", x % p0, z % p0); */

  /* other primes */
  for (s = primes; s[0] <= B1; s++)
    {
      if (s[0] == 0.0)
        {
          fprintf (stderr, "Error, not enough small primes\n");
          exit (1);
        }
      for (r = s[0]; r <= B1; r *= s[0])
        ecm_mul_ui (&x, &z, (ulong) s[0], p, b);
    }

  u = inverse (z, p0); /* return 0 if non-trivial gcd is found */

  return u == 0;
}

static void
usage (char *s)
{
  printf ("Usage: %s P B1 sigma - find all primes <= P hit by ECM(B1,sigma)\n", s);
  exit (1);
}

int
main (int argc, char *argv[])
{
  unsigned long P, B1, sigma;
  unsigned long primes = 0, found = 0;
  double p;

  if (argc != 4)
    usage (argv[0]);

  P = atoi (argv[1]);
  B1 = atoi (argv[2]);
  sigma = atoi (argv[3]);

  for (p = 2.0; p <= P; p = getprime (p))
    {
      primes ++;
      if (ecm ((unsigned long) p, B1, sigma))
        found ++;
    }
  getprime (0.0);

  printf ("primes=%lu found=%lu (%1.0f%%)\n", primes, found,
	  100.0 * (double) found / (double) primes);
  
  return 0;
}
