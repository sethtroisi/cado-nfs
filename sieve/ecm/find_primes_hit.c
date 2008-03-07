/* Usage: find_primes_hit P B1 sigma - find all primes <= P that are found with
   ECM(B1,sigma) */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "getprime.c"

#undef MODTRACE
#include "mod_ul.h"

#define BRENT12 0
#define MONTY12 1

int verbose = 0;

/* computes 2P=(x2:z2) from P=(x1:z1), with 5 muls (3 muls and 2 squares)
   and 4 add/sub.
     - p : number to factor
     - b : (a+2)/4 mod n
*/
static void
ellM_duplicate (residue_t x2, residue_t z2, residue_t x1, residue_t z1, 
                modulus_t m, residue_t b)
{
  residue_t u, v, w;

  mod_init_noset0 (u, m);
  mod_init_noset0 (v, m);
  mod_init_noset0 (w, m);

  mod_add (u, x1, z1, m);
  mod_mul (u, u, u, m);   /* u = (x1 + z1)^2 */
  mod_sub (v, x1, z1, m);
  mod_mul (v, v, v, m);   /* v = (x1 - z1)^2 */
  mod_mul (x2, u, v, m);  /* x2 = (x1^2 - z1^2)^2 */
  mod_sub (w, u, v, m);   /* w = 4 * x1 * z1 */
  mod_mul (u, w, b, m);   /* u = x1 * z1 * (A + 2) */
  mod_add (u, u, v, m);
  mod_mul (z2, w, u, m);

  mod_clear (w, m);
  mod_clear (v, m);
  mod_clear (u, m);
}


/* For Weierstrass coordinates. Returns 1 if doubling worked normally, 
   0 if the result is point at infinity */

static int
ellW_duplicate (residue_t x3, residue_t y3, residue_t x1, residue_t y1,
	        residue_t a, modulus_t m)
{
  residue_t lambda, u, v;

  mod_init_noset0 (lambda, m);
  mod_init_noset0 (u, m);
  mod_init_noset0 (v, m);

  mod_mul (u, x1, x1, m);
  mod_add (v, u, u, m);
  mod_add (v, v, u, m);
  mod_add (v, v, a, m); /* 3x^2 + a */
  mod_add (u, y1, y1, m);
  if (mod_inv (u, u, m) == 0)    /* 1/(2*y) */
  {
      mod_clear (v, m);
      mod_clear (u, m);
      mod_clear (lambda, m);
      return 0; /* y was 0  =>  result is point at infinity */
  }
  mod_mul (lambda, u, v, m);
  mod_mul (u, lambda, lambda, m);
  mod_sub (u, u, x1, m);
  mod_sub (u, u, x1, m);    /* x3 = u = lambda^2 - 2*x */
  mod_sub (v, x1, u, m);
  mod_mul (v, v, lambda, m);
  mod_sub (y3, v, y1, m);
  mod_set (x3, u, m);
  
  mod_clear (v, m);
  mod_clear (u, m);
  mod_clear (lambda, m);
  return 1;
}


/* adds Q=(x2:z2) and R=(x1:z1) and puts the result in (x3:z3),
     using 6 muls (4 muls and 2 squares), and 6 add/sub.
   One assumes that Q-R=P or R-Q=P where P=(x:z).
     - n : number to factor
   Modifies: x3, z3, u, v, w.
   (x3,z3) may be identical to (x2,z2) and to (x,z)
*/

static void
ellM_add3 (residue_t x3, residue_t z3, residue_t x2, residue_t z2, 
           residue_t x1, residue_t z1, residue_t x, residue_t z, 
           residue_t b __attribute__ ((unused)), modulus_t m)
{
  residue_t u, v, w;

  mod_init_noset0 (u, m);
  mod_init_noset0 (v, m);
  mod_init_noset0 (w, m);

  if (mod_is0 (z1, m)) /* (x1,z1) is at infinity */
    {
      mod_set (x3, x2, m);
      mod_set (z3, z2, m);
      return;
    }
  else if (mod_is0 (z2, m)) /* (x2, z2) is at infinity */
    {
      mod_set (x3, x1, m);
      mod_set (z3, z1, m);
      return;
    }
  else
  {
      mod_sub (u, x2, z2, m);
      mod_add (v, x1, z1, m);
      mod_mul (u, u, v, m);
      mod_add (w, x2, z2, m);
      mod_sub (v, x1, z1, m);
      mod_mul (v, w, v, m);
      mod_add (w, u, v, m);
      mod_sub (v, u, v, m);
      mod_mul (w, w, w, m);
      mod_mul (v, v, v, m);
      mod_set (u, x, m); /* save x */
      mod_mul (x3, w, z, m);
      mod_mul (z3, u, v, m);
  }

  mod_clear (w, m);
  mod_clear (v, m);
  mod_clear (u, m);
}

/* Adds two points (x2, y2) and (x1, y1) on the curve y^2 = x^3 + a*x + b
   in Weierstrass coordinates and puts result in (x3, y3). 
   Returns 1 if the addition worked (i.e. the modular inverse existed) 
   and 0 otherwise (resulting point is point at infinity) */

static int
ellW_add3 (residue_t x3, residue_t y3, residue_t x2, residue_t y2, 
           residue_t x1, residue_t y1, residue_t a, modulus_t m)
{
  residue_t lambda, u, v;
  int r;

  mod_init_noset0 (u, m);
  mod_init_noset0 (v, m);

  mod_sub (u, y2, y1, m);
  mod_sub (v, x2, x1, m);
  r = mod_inv (v, v, m);
  if (r == 0)
  {
      /* Maybe we were trying to add two identical points? If so,
         use the duplicateW() function instead */
      if (mod_equal (x1, x2, m) && mod_equal (y1, y2, m))
	  r = ellW_duplicate (x3, y3, x1, y1, a, m);
      else
	  r = 0; /* No, the points were negatives of each other */
  }
  else
  {
      mod_mul (lambda, u, v, m);
      mod_mul (u, lambda, lambda, m);
      mod_sub (u, u, x1, m);
      mod_sub (u, u, x2, m);    /* x3 = u = lambda^2 - x1 - x2 */
      mod_sub (v, x1, u, m);
      mod_mul (v, v, lambda, m);
      mod_sub (y3, v, y1, m);
      mod_set (x3, u, m);
      r = 1;
  }

  mod_clear (v, m);
  mod_clear (u, m);
  return r;
}


/* (x:z) <- e*(x:z) (mod p)
   Assumes e >= 5.
*/
static void
ellM_mul_ui (residue_t x, residue_t z, unsigned long e, 
	     modulus_t m, residue_t b)
{
  unsigned long l, n;
  residue_t x1, z1, x2, z2;

  mod_init_noset0 (x1, m);
  mod_init_noset0 (z1, m);
  mod_init_noset0 (x2, m);
  mod_init_noset0 (z2, m);

  e --;

  /* compute number of steps needed: we start from (1,2) and go from
     (i,i+1) to (2i,2i+1) or (2i+1,2i+2) */
  for (l = e, n = 0; l > 1; n ++, l /= 2);

  /* start from P1=P, P2=2P */
  mod_set (x1, x, m);
  mod_set (z1, z, m);
  ellM_duplicate (x2, z2, x1, z1, m, b);

  while (n--)
    {
      if ((e >> n) & 1) /* (i,i+1) -> (2i+1,2i+2) */
        {
          /* printf ("(i,i+1) -> (2i+1,2i+2)\n"); */
          ellM_add3 (x1, z1, x2, z2, x1, z1, x, z, b, m);
          ellM_duplicate (x2, z2, x2, z2, m, b);
        }
      else /* (i,i+1) -> (2i,2i+1) */
        {
          /* printf ("(i,i+1) -> (2i,2i+1)\n"); */
          ellM_add3 (x2, z2, x1, z1, x2, z2, x, z, b, m);
          ellM_duplicate (x1, z1, x1, z1, m, b);
        }
    }
  
  mod_set (x, x2, m);
  mod_set (z, z2, m);

  mod_clear (z2, m);
  mod_clear (x2, m);
  mod_clear (z1, m);
  mod_clear (x1, m);
}

static int
ellW_mul_ui (residue_t x, residue_t y, unsigned long e, residue_t a, 
	     modulus_t m)
{
  unsigned long i;
  residue_t xt, yt;
  int tfinite; /* Nonzero iff (xt, yt) is NOT point at infinity */

  if (e == 0)
    return 0; /* signal point at infinity */

  mod_init_noset0 (xt, m);
  mod_init_noset0 (yt, m);

  i = ~(0UL);
  i -= i/2;   /* Now the most significant bit of i is set */
  while ((i & e) == 0)
    i >>= 1;

  mod_set (xt, x, m);
  mod_set (yt, y, m);
  tfinite = 1;
  i >>= 1;

  while (i > 0)
  {
      if (tfinite)
        tfinite = ellW_duplicate (xt, yt, xt, yt, a, m);
      if (e & i)
      {
	  if (tfinite)
	      tfinite = ellW_add3 (xt, yt, x, y, xt, yt, a, m);
	  else
	  {
	      mod_set (xt, x, m);
	      mod_set (yt, y, m);
	      tfinite = 1;
	  }
      }
      i >>= 1;
  }

  if (tfinite)
  {
      mod_set (x, xt, m);
      mod_set (y, yt, m);
  }
  mod_clear (yt, m);
  mod_clear (xt, m);

  return tfinite;
}

/* Produces curve in Montgomery parameterization from sigma value.
   Return 1 if it worked, 0 if a modular inverse failed */

static int
Brent12_curve_from_sigma (residue_t A, residue_t x, residue_t sigma, modulus_t m)
{
  residue_t u, v, t, b, z;
  int r;

  mod_init_noset0 (u, m);
  mod_init_noset0 (v, m);
  mod_init_noset0 (t, m);
  mod_init_noset0 (b, m);
  mod_init_noset0 (z, m);

  /* compute b, x */
  mod_add (v, sigma, sigma, m);
  mod_add (v, v, v, m); /* v = 4*sigma */
  mod_mul (u, sigma, sigma, m);
  mod_set_ul (t, 5UL, m);
  mod_sub (u, u, t, m); /* u = sigma^2 - 5 */
  mod_mul (t, u, u, m);
  mod_mul (x, t, u, m);
  mod_mul (t, v, v, m);
  mod_mul (z, t, v, m);
  mod_mul (t, x, v, m);
  mod_add (b, t, t, m);
  mod_add (b, b, b, m); /* b = 4 * t */
  mod_add (t, u, u, m);
  mod_add (t, t, u, m); /* t = 3 * u */
  mod_sub (u, v, u, m);
  mod_add (v, t, v, m);
  mod_mul (t, u, u, m);
  mod_mul (u, t, u, m);
  mod_mul (A, u, v, m);
  mod_mul (v, b, z, m);

  r = mod_inv (u, v, m);
  if (r) /* non trivial gcd */
  {
      mod_mul (v, u, b, m);
      mod_mul (x, x, v, m);
      mod_mul (v, u, z, m);
      mod_mul (t, A, v, m);
      mod_set_ul (u, 2UL, m);
      mod_sub (A, t, u, m);
  }

  mod_clear (z, m);
  mod_clear (b, m);
  mod_clear (t, m);
  mod_clear (v, m);
  mod_clear (u, m);

  return r;
}

/* Produces curve in Montgomery parameterization from n value, using
   parameters for a torsion 12 curve as in Montgomery's thesis.
   Return 1 if it worked, 0 if a modular inverse failed */

static int
Monty12_curve_from_k (residue_t A, residue_t x, unsigned long n, modulus_t m)
{
  residue_t u, v, u0, v0, a, t2;
  
  /* We want a multiple of the point (-2,4) on the curve Y^2=X^3-12*X */
  mod_init (a, m);
  mod_init (u, m);
  mod_init_noset0 (v, m);
  mod_init (u0, m);
  mod_init (v0, m);

  mod_sub_ul (a, a, 12UL, m);
  mod_sub_ul (u, u, 2UL, m);
  mod_set_ul (v, 4UL, m);
  ellW_mul_ui (u, v, n/2, a, m);
  if (n % 2 == 1)
    ellW_add3 (u, v, u, v, u0, v0, a, m);
  /* Now we have a $u$ so that $u^3-12u$ is a square */
  mod_clear (u0, m);
  mod_clear (v0, m);
  /* printf ("Monty12_curve_from_k: u = %lu\n", mod_get_ul (u)); */
  
  mod_init_noset0 (t2, m);
  mod_div2 (v, u, m);
  mod_mul (t2, v, v, m); /* u^2/4 */
  mod_sub_ul (t2, t2, 3UL, m);
  if (mod_inv (u, u, m) == 0)
  {
    fprintf (stderr, "Monty12_curve_from_k: u = 0\n");
    mod_clear (t2, m);
    mod_clear (v, m);
    mod_clear (u, m);
    mod_clear (a, m);
    return 0;
  }
  mod_mul (t2, t2, u, m); /* t^2 = (u^2/4 - 3)/u = (u^2 - 12)/4u */

  mod_sub_ul (u, t2, 1UL, m);
  mod_add_ul (v, t2, 3UL, m);
  mod_mul (a, u, v, m);
  if (mod_inv (a, a, m) == 0) /* a  = 1/(uv), I want u/v and v/u */
  {
    fprintf (stderr, "Monty12_curve_from_k: (t^2 - 1)(t^2 + 3) = 0\n");
    mod_clear (t2, m);
    mod_clear (v, m);
    mod_clear (u, m);
    mod_clear (a, m);
    return 0;
  }
  mod_mul (u, u, u, m); /* u^2 */
  mod_mul (v, v, v, m); /* v^2 */
  mod_mul (v, v, a, m); /* v^2 * (1/(uv)) = v/u = 1/a*/
  mod_mul (a, a, u, m); /* u^2 * (1/(uv)) = u/v = a*/

  mod_mul (u, a, a, m); /* a^2 */
  mod_add_ul (A, u, 2UL, m); /* a^2 + 2 */
  mod_add (t2, A, A, m);
  mod_add (A, A, t2, m); /* 3*(a^2 + 2) */
  mod_mul (t2, A, a, m);
  mod_set (A, v, m);
  mod_sub (A, A, t2, m); /* 1/a - 3 a (a^2 + 2) */
  mod_div2 (v, v, m); /* v = 1/(2a) */
  mod_mul (t2, v, v, m); /* t2 = 1/(2a)^2 */
  mod_mul (A, A, t2, m);

  mod_add (x, u, u, m);
  mod_add (x, x, u, m); /* 3*a^2 */
  mod_add_ul (x, x, 1UL, m); /* 3*a^2 + 1 */
  mod_div2 (v, v, m); /* v = 1/(4a) */
  mod_mul (x, x, v, m);
  
  mod_clear (t2, m);
  mod_clear (v, m);
  mod_clear (u, m);
  mod_clear (a, m);
  return 1;
}


/* Make a curve of the form y^2 = x^3 + a*x^2 + b with a valid point
   (x, y) from a curve Y^2 = X^3 + A*X^2 + X. The value of b will not
   be computed. 
   x and X may be the same variable. */

static int
curveW_from_Montgomery (residue_t a, residue_t x, residue_t y,
			residue_t X, residue_t A, modulus_t m)
{
  residue_t g, one;
  int r;

  mod_init_noset0 (g, m);
  mod_init_noset0 (one, m);

  mod_set_ul (one, 1UL, m);
  mod_add (g, X, A, m);
  mod_mul (g, g, X, m);
  mod_add_ul (g, g, 1UL, m);
  mod_mul (g, g, X, m); /* G = X^3 + A*X^2 + X */
  /* printf ("curveW_from_Montgomery: Y^2 = %lu\n", g[0]); */

  /* Now (x,1) is on the curve G*Y^2 = X^3 + A*X^2 + X. */
  r = mod_inv (g, g, m);
  if (r != 0)
  {
      mod_set (y, g, m);       /* y = 1/G */
      mod_div3 (a, A, m);
      mod_add (x, X, a, m);
      mod_mul (x, x, g, m); /* x = (X + A/3)/G */
      mod_mul (a, a, A, m);
      mod_sub (a, one, a, m);
      mod_mul (a, a, g, m);
      mod_mul (a, a, g, m); /* a = (1 - (A^2)/3)/G^2 */
  }
  else
    fprintf (stderr, "curveW_from_Montgomery: r = 0\n");

  mod_clear (one, m);
  mod_clear (g, m);

  return r;
}

static int
ecm (unsigned long p_par, double B1, unsigned long sigma_par, 
     int parameterization)
{
  /* small primes */
  static double primes[] = {5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 0};
  residue_t u, A, b, x, z, xB, zB, sigma;
  modulus_t p;
  double r, *s;

  mod_initmod_ul (p, p_par);
  mod_set_ul (sigma, sigma_par, p);

  if (parameterization == BRENT12)
  {
    if (Brent12_curve_from_sigma (A, x, sigma, p) == 0)
      return 1;
  }
  else if (parameterization == MONTY12)
  {
    if (Monty12_curve_from_k (A, x, sigma_par, p) == 0)
      return 1;
  }
  else
  {
    fprintf (stderr, "ecm: Unknown parameterization\n");
    abort();
  }

  mod_add_ul (b, A, 2UL, p);
  mod_div2 (b, b, p);
  mod_div2 (b, b, p);


  /* now start ecm */
  mod_set_ul (z, 1UL, p);

  /* printf ("p=%lu, x=%lu b=%lu\n", p_par, mod_get_ul(x), 
           mod_get_ul(b)); */

  /* prime 2 */
  for (r = 2.0; r <= B1; r *= 2.0)
    ellM_duplicate (x, z, x, z, p, b);

  /* printf ("2: x=%lu z=%lu\n", mod_get_ul(x) , mod_get_ul(z)); */

  /* prime 3 */
  for (r = 3.0; r <= B1; r *= 3.0)
    {
      ellM_duplicate (xB, zB, x, z, p, b);
      ellM_add3 (x, z, x, z, xB, zB, x, z, b, p);
    }

  /* printf ("3: x=%lu z=%lu\n", mod_get_ul(x), mod_get_ul(z)); */

  /* other primes */
  for (s = primes; s[0] <= B1; s++)
    {
      if (s[0] == 0.0)
        {
          fprintf (stderr, "Error, not enough small primes\n");
          exit (1);
        }
      for (r = s[0]; r <= B1; r *= s[0])
        ellM_mul_ui (x, z, (unsigned long) s[0], p, b);
    }

  return (mod_inv (u, z, p) == 0); /* return 1 if non-trivial gcd is found */
}

/* Determine order of a point P on a curve, both defined by the sigma value
   as in ECM. Looks for i in Hasse interval so that i*P = O, has complexity
   O(sqrt(m)). */

unsigned long
ell_pointorder (const unsigned long m_par, const unsigned long sigma_par)
{
  residue_t sigma, A, X, a, xi, yi, x1, y1;
  modulus_t m;
  unsigned long min, max, i, order, p;

  mod_initmod_ul (m, m_par);
  mod_set_ul (sigma, sigma_par, m);

  if (Brent12_curve_from_sigma (A, X, sigma, m) == 0)
    return 0;
  
  /* printf ("Curve parameters in Montgomery form: A = %ld, X = %ld "
     "(mod %ld)\n", A[0], X[0], m[0]); */

  if (curveW_from_Montgomery (a, x1, y1, X, A, m) == 0)
    return 0UL;

  if (verbose >= 2)
    printf ("Finding order of point (%ld, %ld) on curve "
	    "y^2 = x^3 + %ld * x + b (mod %ld)\n", 
	    x1[0], y1[0], a[0], m[0]);

  i = 2 * (unsigned long) sqrt((double) m_par);
  min = m_par - i + 1;
  max = m_par + i + 1;
  mod_set (xi, x1, m);
  mod_set (yi, y1, m);
  if (ellW_mul_ui (xi, yi, min, a, m) == 0)
  {
      i = min;
  }
  else
  {
      for (i = min + 1; i <= max; i++)
      {
	  if (!ellW_add3 (xi, yi, xi, yi, x1, y1, a, m))
	      break;
      }
      
      if (i > max)
      {
	  fprintf (stderr, "ell_order: Error, point at infinity not "
		   "reached with i*(x0, z0), i in [%ld, %ld]\n", min, max);
	  return 0UL;
      }

#ifndef NDEBUG
      /* Check that this is the correct order */
      mod_set (xi, x1, m);
      mod_set (yi, y1, m);
      if (ellW_mul_ui (xi, yi, i, a, m) != 0)
      {
	  fprintf (stderr, "ell_order: Error, %ld*(%ld, %ld) (mod %ld) is "
		   "not the point at infinity\n", 
		   i, mod_get_ul (x1), mod_get_ul (y1), m_par);
	  return 0UL;
      }
#endif
  }
  
  /* Ok, now we have some i so that ord(P) | i. Find ord(P).
     We know that ord(P) > 1 since P is not at infinity */

  order = i;
  for (p = 2; p * p <= order; p++)
  {
      mod_set (xi, x1, m);
      mod_set (yi, y1, m);
      while (order % p == 0 && ellW_mul_ui (xi, yi, order / p, a, m) == 0)
	  order /= p;
  }

  return order;
}


/* Count points on curve using the Jacobi symbol. This has complexity O(m). */

unsigned long 
ellM_curveorderjacobi (residue_t A, residue_t X, modulus_t m)
{
  residue_t t;
  unsigned long order, i;
  int bchar;

  mod_init_noset0 (t, m);

  /* Compute X^3 + A*X^2 + X and see if it is a square */
  mod_set (t, X, m);
  mod_add (t, t, A, m);
  mod_mul (t, t, X, m);
  mod_add_ul (t, t, 1UL, m);
  mod_mul (t, t, X, m);
  bchar = mod_jacobi (t, m);
  ASSERT (bchar != 0);

  order = 2; /* One for (0, 0, 1), one for the point at infinity */
  for (i = 1; i < mod_getmod_ul(m); i++)
    {
      mod_set_ul (X, i, m);
      mod_set (t, X, m);
      mod_add (t, t, A, m);
      mod_mul (t, t, X, m);
      mod_add_ul (t, t, 1UL, m);
      mod_mul (t, t, X, m);
      if (bchar == 1) 
	order = order + 1 + mod_jacobi (t, m);
      else
	order = order + 1 - mod_jacobi (t, m);
    }

  mod_clear (t, m);
  
  return order;
}

unsigned long 
ell_curveorder (const unsigned long m_par, const unsigned long sigma_par,
                int parameterization)
{
  residue_t sigma, A, X;
  modulus_t m;
  unsigned long order;

  mod_initmod_ul (m, m_par);
  mod_set_ul (sigma, sigma_par, m);

  if (parameterization == BRENT12)
  {
    if (Brent12_curve_from_sigma (A, X, sigma, m) == 0)
      return 0UL;
  }
  else if (parameterization == MONTY12)
  {
    if (Monty12_curve_from_k (A, X, sigma_par, m) == 0)
      return 0UL;
  }
  else
  {
    fprintf (stderr, "ell_curveorder: Unknown parameterization\n");
    abort();
  }
  order = ellM_curveorderjacobi (A, X, m);

#ifndef NDEBUG
  ASSERT (parameterization != BRENT12 || order == ell_pointorder (m, sigma));
#endif

  return order;
}


static void
usage (const char *s)
{
  printf ("Usage: %s Pmin Pmax B1 sigma\n", s); 
  printf ("Find all primes >= Pmin, <= Pmax hit by ECM(B1,sigma)\n");
  printf ("Options:\n");
  printf ("-v        Increase verbosity\n");
  printf ("-gnuplot  Make output suitable for gnuplot\n");
  printf ("-c        Print average exponent of 2, 3, 5, 7 in group order\n");
  printf ("-m12      Make curves of rational 12 torsion from Montgomery's thesis\n");
  exit (1);
}

unsigned int
prime_exponent (unsigned long n, unsigned long p)
{
  unsigned int i = 0;
  while (n % p == 0)
  {
    i++;
    n /= p;
  }
  return i;
}

int
main (int argc, char *argv[])
{
  unsigned long Pmin, Pmax, B1, sigma;
  unsigned long primes = 0, found = 0;
  double p;
  unsigned long nr_points, pow2 = 0, pow3 = 0, pow5 = 0, pow7 = 0;
  int do_count = 0, gnuplot = 0;
  int parameterization = BRENT12;
  const char *programname = argv[0];

  while (argc > 1 && argv[1][0] == '-')
    {
	if (strcmp(argv[1], "-v") == 0)
	{
	    verbose++;
	    argv++;
	    argc--;
	}
	else if (strcmp(argv[1], "-gnuplot") == 0)
	{
	    gnuplot = 1;
	    argv++;
	    argc--;
	}
	else if (strcmp(argv[1], "-c") == 0)
	{
	    do_count = 1;
	    argv++;
	    argc--;
	}
	else if (strcmp(argv[1], "-m12") == 0)
	{
	    parameterization = MONTY12;
	    argv++;
	    argc--;
	}
	else 
	  usage (programname);
    }

  if (argc < 5)
    usage (programname);

  Pmin = atol (argv[1]);
  Pmax = atol (argv[2]);
  B1 = atol (argv[3]);
  sigma = atol (argv[4]);

  for (p = 2.0; p < Pmin; p = getprime (p));

  for ( ; p <= Pmax; p = getprime (p))
    {
      primes ++;
      if (verbose)
	{
	  printf ("Testing p = %ld", (unsigned long) p);
	}

      if (ecm ((unsigned long) p, (double) B1, sigma, parameterization))
	{
	  found ++;
	  if (verbose)
	    printf (": smooth");
	}

      if (do_count)
	{
	  nr_points = ell_curveorder ((unsigned long) p, sigma, 
	                              parameterization);
	  if (verbose)
	      printf (", order = %ld", nr_points);
	  if (nr_points > 0)
	    {
	      pow2 += prime_exponent (nr_points, 2UL);
	      pow3 += prime_exponent (nr_points, 3UL);
	      pow5 += prime_exponent (nr_points, 5UL);
	      pow7 += prime_exponent (nr_points, 7UL);
	    }
	}
      if (verbose)
	  printf ("\n");
    }
  getprime (0.0);

  if (do_count)
    {
      printf ("Avg. exponent of 2: %f, 3: %f, 5: %f, 7: %f\n",
	      (double)pow2/(double)primes, (double)pow3/(double)primes, 
	      (double)pow5/(double)primes, (double)pow7/(double)primes);
    }

  if (gnuplot)
      printf ("%lu %f\n", Pmin, (double) found / (double) primes);
  else
      printf ("primes=%lu found=%lu (%1.0f%%)\n", primes, found,
	      100.0 * (double) found / (double) primes);
  
  return 0;
}
