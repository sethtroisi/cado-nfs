/* Usage: find_primes_hit P B1 sigma - find all primes <= P that are found with
   ECM(B1,sigma) */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "getprime.c"

#undef MODTRACE
#include "mod_ul.c"

#define BRENT12 0
#define MONTY12 1

int verbose = 0;

/* computes 2P=(x2:z2) from P=(x1:z1), with 5 muls (3 muls and 2 squares)
   and 4 add/sub.
     - p : number to factor
     - b : (a+2)/4 mod n
*/
static void
ellM_duplicate (residue x2, residue z2, residue x1, residue z1, 
                modulus m, residue b)
{
  residue u, v, w;

  initmod (u);
  initmod (v);
  initmod (w);

  addmod (u, x1, z1, m);
  mulmod (u, u, u, m);   /* u = (x1 + z1)^2 */
  submod (v, x1, z1, m);
  mulmod (v, v, v, m);   /* v = (x1 - z1)^2 */
  mulmod (x2, u, v, m);  /* x2 = (x1^2 - z1^2)^2 */
  submod (w, u, v, m);   /* w = 4 * x1 * z1 */
  mulmod (u, w, b, m);   /* u = x1 * z1 * (A + 2) */
  addmod (u, u, v, m);
  mulmod (z2, w, u, m);

  clearmod (w);
  clearmod (v);
  clearmod (u);
}


/* For Weierstrass coordinates. Returns 1 if doubling worked normally, 
   0 if the result is point at infinity */

static int
ellW_duplicate (residue x3, residue y3, residue x1, residue y1,
	        residue a, modulus m)
{
  residue lambda, u, v;

  initmod (lambda);
  initmod (u);
  initmod (v);

  mulmod (u, x1, x1, m);
  addmod (v, u, u, m);
  addmod (v, v, u, m);
  addmod (v, v, a, m); /* 3x^2 + a */
  addmod (u, y1, y1, m);
  if (invmod (u, u, m) == 0)    /* 1/(2*y) */
  {
      clearmod (v);
      clearmod (u);
      clearmod (lambda);
      return 0; /* y was 0  =>  result is point at infinity */
  }
  mulmod (lambda, u, v, m);
  mulmod (u, lambda, lambda, m);
  submod (u, u, x1, m);
  submod (u, u, x1, m);    /* x3 = u = lambda^2 - 2*x */
  submod (v, x1, u, m);
  mulmod (v, v, lambda, m);
  submod (y3, v, y1, m);
  setmod (x3, u);
  
  clearmod (v);
  clearmod (u);
  clearmod (lambda);
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
ellM_add3 (residue x3, residue z3, residue x2, residue z2, 
           residue x1, residue z1, residue x, residue z, 
           residue b __attribute__ ((unused)), modulus m)
{
  residue u, v, w;

  initmod (u);
  initmod (v);
  initmod (w);

  if (is0mod (z1)) /* (x1,z1) is at infinity */
    {
      setmod (x3, x2);
      setmod (z3, z2);
      return;
    }
  else if (is0mod (z2)) /* (x2, z2) is at infinity */
    {
      setmod (x3, x1);
      setmod (z3, z1);
      return;
    }
  else
  {
      submod (u, x2, z2, m);
      addmod (v, x1, z1, m);
      mulmod (u, u, v, m);
      addmod (w, x2, z2, m);
      submod (v, x1, z1, m);
      mulmod (v, w, v, m);
      addmod (w, u, v, m);
      submod (v, u, v, m);
      mulmod (w, w, w, m);
      mulmod (v, v, v, m);
      setmod (u, x); /* save x */
      mulmod (x3, w, z, m);
      mulmod (z3, u, v, m);
  }

  clearmod (w);
  clearmod (v);
  clearmod (u);
}

/* Adds two points (x2, y2) and (x1, y1) on the curve y^2 = x^3 + a*x + b
   in Weierstrass coordinates and puts result in (x3, y3). 
   Returns 1 if the addition worked (i.e. the modular inverse existed) 
   and 0 otherwise (resulting point is point at infinity) */

static int
ellW_add3 (residue x3, residue y3, residue x2, residue y2, 
           residue x1, residue y1, residue a, modulus m)
{
  residue lambda, u, v;
  int r;

  initmod (u);
  initmod (v);

  submod (u, y2, y1, m);
  submod (v, x2, x1, m);
  r = invmod (v, v, m);
  if (r == 0)
  {
      /* Maybe we were trying to add two identical points? If so,
         use the duplicateW() function instead */
      if (cmpmod (x1, x2) == 0 && cmpmod (y1, y2) == 0)
	  r = ellW_duplicate (x3, y3, x1, y1, a, m);
      else
	  r = 0; /* No, the points were negatives of each other */
  }
  else
  {
      mulmod (lambda, u, v, m);
      mulmod (u, lambda, lambda, m);
      submod (u, u, x1, m);
      submod (u, u, x2, m);    /* x3 = u = lambda^2 - x1 - x2 */
      submod (v, x1, u, m);
      mulmod (v, v, lambda, m);
      submod (y3, v, y1, m);
      setmod (x3, u);
      r = 1;
  }

  clearmod (v);
  clearmod (u);
  return r;
}


/* (x:z) <- e*(x:z) (mod p)
   Assumes e >= 5.
*/
static void
ellM_mul_ui (residue x, residue z, unsigned long e, 
	     modulus m, residue b)
{
  unsigned long l, n;
  residue x1, z1, x2, z2;

  initmod (x1);
  initmod (z1);
  initmod (x2);
  initmod (z2);

  e --;

  /* compute number of steps needed: we start from (1,2) and go from
     (i,i+1) to (2i,2i+1) or (2i+1,2i+2) */
  for (l = e, n = 0; l > 1; n ++, l /= 2);

  /* start from P1=P, P2=2P */
  setmod (x1, x);
  setmod (z1, z);
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
  
  setmod (x, x2);
  setmod (z, z2);

  clearmod (z2);
  clearmod (x2);
  clearmod (z1);
  clearmod (x1);
}

static int
ellW_mul_ui (residue x, residue y, unsigned long e, residue a, 
	     modulus m)
{
  unsigned long i;
  residue xt, yt;
  int tfinite; /* Nonzero iff (xt, yt) is NOT point at infinity */

  if (e == 0)
    return 0; /* signal point at infinity */

  initmod (xt);
  initmod (yt);

  i = ~(0UL);
  i -= i/2;   /* Now the most significant bit of i is set */
  while ((i & e) == 0)
    i >>= 1;

  setmod (xt, x);
  setmod (yt, y);
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
	      setmod (xt, x);
	      setmod (yt, y);
	      tfinite = 1;
	  }
      }
      i >>= 1;
  }

  if (tfinite)
  {
      setmod (x, xt);
      setmod (y, yt);
  }
  clearmod (yt);
  clearmod (xt);

  return tfinite;
}

/* Produces curve in Montgomery parameterization from sigma value.
   Return 1 if it worked, 0 if a modular inverse failed */

static int
Brent12_curve_from_sigma (residue A, residue x, residue sigma, modulus m)
{
  residue u, v, t, b, z;
  int r;

  initmod (u);
  initmod (v);
  initmod (t);
  initmod (b);
  initmod (z);

  /* compute b, x */
  addmod (v, sigma, sigma, m);
  addmod (v, v, v, m); /* v = 4*sigma */
  mulmod (u, sigma, sigma, m);
  setmod_ul (t, 5UL, m);
  submod (u, u, t, m); /* u = sigma^2 - 5 */
  mulmod (t, u, u, m);
  mulmod (x, t, u, m);
  mulmod (t, v, v, m);
  mulmod (z, t, v, m);
  mulmod (t, x, v, m);
  addmod (b, t, t, m);
  addmod (b, b, b, m); /* b = 4 * t */
  addmod (t, u, u, m);
  addmod (t, t, u, m); /* t = 3 * u */
  submod (u, v, u, m);
  addmod (v, t, v, m);
  mulmod (t, u, u, m);
  mulmod (u, t, u, m);
  mulmod (A, u, v, m);
  mulmod (v, b, z, m);

  r = invmod (u, v, m);
  if (r) /* non trivial gcd */
  {
      mulmod (v, u, b, m);
      mulmod (x, x, v, m);
      mulmod (v, u, z, m);
      mulmod (t, A, v, m);
      setmod_ul (u, 2UL, m);
      submod (A, t, u, m);
  }

  clearmod (z);
  clearmod (b);
  clearmod (t);
  clearmod (v);
  clearmod (u);

  return r;
}

/* Produces curve in Montgomery parameterization from n value, using
   parameters for a torsion 12 curve as in Montgomery's thesis.
   Return 1 if it worked, 0 if a modular inverse failed */

static int
Monty12_curve_from_k (residue A, residue x, unsigned long n, modulus m)
{
  residue u, v, u0, v0, a, t2;
  
  /* We want a multiple of the point (-2,4) on the curve Y^2=X^3-12*X */
  initmod_set0 (a);
  initmod_set0 (u);
  initmod (v);
  initmod_set0 (u0);
  initmod_set0 (v0);

  submod_ul (a, a, 12UL, m);
  submod_ul (u, u, 2UL, m);
  setmod_ul (v, 4UL, m);
  ellW_mul_ui (u, v, n/2, a, m);
  if (n % 2 == 1)
    ellW_add3 (u, v, u, v, u0, v0, a, m);
  /* Now we have a $u$ so that $u^3-12u$ is a square */
  clearmod (u0);
  clearmod (v0);
  /* printf ("Monty12_curve_from_k: u = %lu\n", getmod_ul (u)); */
  
  initmod (t2);
  div2mod (v, u, m);
  mulmod (t2, v, v, m); /* u^2/4 */
  submod_ul (t2, t2, 3UL, m);
  if (invmod (u, u, m) == 0)
  {
    fprintf (stderr, "Monty12_curve_from_k: u = 0\n");
    clearmod (t2);
    clearmod (v);
    clearmod (u);
    clearmod (a);
    return 0;
  }
  mulmod (t2, t2, u, m); /* t^2 = (u^2/4 - 3)/u = (u^2 - 12)/4u */

  submod_ul (u, t2, 1UL, m);
  addmod_ul (v, t2, 3UL, m);
  mulmod (a, u, v, m);
  if (invmod (a, a, m) == 0) /* a  = 1/(uv), I want u/v and v/u */
  {
    fprintf (stderr, "Monty12_curve_from_k: (t^2 - 1)(t^2 + 3) = 0\n");
    clearmod (t2);
    clearmod (v);
    clearmod (u);
    clearmod (a);
    return 0;
  }
  mulmod (u, u, u, m); /* u^2 */
  mulmod (v, v, v, m); /* v^2 */
  mulmod (v, v, a, m); /* v^2 * (1/(uv)) = v/u = 1/a*/
  mulmod (a, a, u, m); /* u^2 * (1/(uv)) = u/v = a*/

  mulmod (u, a, a, m); /* a^2 */
  addmod_ul (A, u, 2UL, m); /* a^2 + 2 */
  addmod (t2, A, A, m);
  addmod (A, A, t2, m); /* 3*(a^2 + 2) */
  mulmod (t2, A, a, m);
  setmod (A, v);
  submod (A, A, t2, m); /* 1/a - 3 a (a^2 + 2) */
  div2mod (v, v, m); /* v = 1/(2a) */
  mulmod (t2, v, v, m); /* t2 = 1/(2a)^2 */
  mulmod (A, A, t2, m);

  addmod (x, u, u, m);
  addmod (x, x, u, m); /* 3*a^2 */
  addmod_ul (x, x, 1UL, m); /* 3*a^2 + 1 */
  div2mod (v, v, m); /* v = 1/(4a) */
  mulmod (x, x, v, m);
  
  clearmod (t2);
  clearmod (v);
  clearmod (u);
  clearmod (a);
  return 1;
}


/* Make a curve of the form y^2 = x^3 + a*x^2 + b with a valid point
   (x, y) from a curve Y^2 = X^3 + A*X^2 + X. The value of b will not
   be computed. 
   x and X may be the same variable. */

static int
curveW_from_Montgomery (residue a, residue x, residue y,
			residue X, residue A, modulus m)
{
  residue g, one;
  int r;

  initmod (g);
  initmod (one);

  setmod_ul (one, 1UL, m);
  addmod (g, X, A, m);
  mulmod (g, g, X, m);
  addmod_ul (g, g, 1UL, m);
  mulmod (g, g, X, m); /* G = X^3 + A*X^2 + X */
  /* printf ("curveW_from_Montgomery: Y^2 = %lu\n", g[0]); */

  /* Now (x,1) is on the curve G*Y^2 = X^3 + A*X^2 + X. */
  r = invmod (g, g, m);
  if (r != 0)
  {
      setmod (y, g);       /* y = 1/G */
      div3mod (a, A, m);
      addmod (x, X, a, m);
      mulmod (x, x, g, m); /* x = (X + A/3)/G */
      mulmod (a, a, A, m);
      submod (a, one, a, m);
      mulmod (a, a, g, m);
      mulmod (a, a, g, m); /* a = (1 - (A^2)/3)/G^2 */
  }
  else
    fprintf (stderr, "curveW_from_Montgomery: r = 0\n");

  clearmod (one);
  clearmod (g);

  return r;
}

static int
ecm (unsigned long p_par, double B1, unsigned long sigma_par, 
     int parameterization)
{
  /* small primes */
  static double primes[] = {5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 0};
  residue u, A, b, x, z, xB, zB, sigma;
  modulus p;
  double r, *s;

  setmodulus_ul (p, p_par);
  setmod_ul (sigma, sigma_par, p);

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

  addmod_ul (b, A, 2UL, p);
  div2mod (b, b, p);
  div2mod (b, b, p);


  /* now start ecm */
  setmod_ul (z, 1UL, p);

  /* printf ("p=%lu, x=%lu b=%lu\n", p_par, getmod_ul(x), 
           getmod_ul(b)); */

  /* prime 2 */
  for (r = 2.0; r <= B1; r *= 2.0)
    ellM_duplicate (x, z, x, z, p, b);

  /* printf ("2: x=%lu z=%lu\n", getmod_ul(x) , getmod_ul(z)); */

  /* prime 3 */
  for (r = 3.0; r <= B1; r *= 3.0)
    {
      ellM_duplicate (xB, zB, x, z, p, b);
      ellM_add3 (x, z, x, z, xB, zB, x, z, b, p);
    }

  /* printf ("3: x=%lu z=%lu\n", getmod_ul(x), getmod_ul(z)); */

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

  return (invmod (u, z, p) == 0); /* return 1 if non-trivial gcd is found */
}

/* Determine order of a point P on a curve, both defined by the sigma value
   as in ECM. Looks for i in Hasse interval so that i*P = O, has complexity
   O(sqrt(m)). */

unsigned long
ell_pointorder (const unsigned long m_par, const unsigned long sigma_par)
{
  residue sigma, A, X, a, xi, yi, x1, y1;
  modulus m;
  unsigned long min, max, i, order, p;

  setmodulus_ul (m, m_par);
  setmod_ul (sigma, sigma_par, m);

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
  setmod (xi, x1);
  setmod (yi, y1);
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
      setmod (xi, x1);
      setmod (yi, y1);
      if (ellW_mul_ui (xi, yi, i, a, m) != 0)
      {
	  fprintf (stderr, "ell_order: Error, %ld*(%ld, %ld) (mod %ld) is "
		   "not the point at infinity\n", 
		   i, getmod_ul (x1), getmod_ul (y1), m_par);
	  return 0UL;
      }
#endif
  }
  
  /* Ok, now we have some i so that ord(P) | i. Find ord(P).
     We know that ord(P) > 1 since P is not at infinity */

  order = i;
  for (p = 2; p * p <= order; p++)
  {
      setmod (xi, x1);
      setmod (yi, y1);
      while (order % p == 0 && ellW_mul_ui (xi, yi, order / p, a, m) == 0)
	  order /= p;
  }

  return order;
}


/* Count points on curve using the Jacobi symbol. This has complexity O(m). */

unsigned long 
ellM_curveorderjacobi (residue A, residue X, modulus m)
{
  residue t;
  unsigned long order, i;
  int bchar;

  initmod (t);

  /* Compute X^3 + A*X^2 + X and see if it is a square */
  setmod (t, X);
  addmod (t, t, A, m);
  mulmod (t, t, X, m);
  addmod_ul (t, t, 1UL, m);
  mulmod (t, t, X, m);
  bchar = jacobimod (t, m);
  ASSERT (bchar != 0);

  order = 2; /* One for (0, 0, 1), one for the point at infinity */
  for (i = 1; i < getmod_ul(m); i++)
    {
      setmod_ul (X, i, m);
      setmod (t, X);
      addmod (t, t, A, m);
      mulmod (t, t, X, m);
      addmod_ul (t, t, 1UL, m);
      mulmod (t, t, X, m);
      if (bchar == 1) 
	order = order + 1 + jacobimod (t, m);
      else
	order = order + 1 - jacobimod (t, m);
    }

  clearmod (t);
  
  return order;
}

unsigned long 
ell_curveorder (const unsigned long m_par, const unsigned long sigma_par,
                int parameterization)
{
  residue sigma, A, X;
  modulus m;
  unsigned long order;

  setmodulus_ul (m, m_par);
  setmod_ul (sigma, sigma_par, m);

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
