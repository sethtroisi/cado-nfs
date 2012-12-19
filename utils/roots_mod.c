/* Routines to solve x^k = a (mod p).

   Reference:
   [1] Adleman-Manders-Miller Root Extraction Method Revisited, Zhengjun Cao,
   Qian Sha, Xiao Fan, 2011, http://arxiv.org/abs/1111.4877.
*/

#include "cado.h"
#include <stdlib.h>
#include <stdint.h>
#include <gmp.h>

#include "gmp_aux.h"
#include "rootfinder.h"
#include "modredc_ul.h"
#include "modredc_ul_default.h"

static int 
mod_roots (residue_t *, residue_t, int, modulus_t);

/* For i < 50, isprime_table[i] == 1 iff i is prime */
static unsigned char isprime_table[] = {0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 
  0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 
  1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0};

/* return 1/a mod b */
static uint64_t
invert_mod (uint64_t a, uint64_t b)
{
  modulus_t bb;
  residue_t aa;
  uint64_t ret;

  mod_initmod_ul (bb, b);
  mod_init (aa, bb);
  mod_set_ul (aa, a, bb);
  mod_inv (aa, aa, bb);
  ret = mod_get_ul (aa, bb);
  mod_clearmod (bb);
  return ret;
}

/* Roots of x^d = a (mod p) for d even and a 64-bit word */
static int
roots2 (residue_t *rr, residue_t aa, int d, modulus_t pp)
{
  uint64_t q, s, n, i, j, k = 0, l;
  residue_t hh, delta, dd, bb, zz;
  const uint64_t p =  mod_getmod_ul(pp);

  if (mod_jacobi (aa, pp) != 1)
    return 0;

  /* find the roots of x^(d/2) = a (mod p) */
  n = mod_roots (rr + d / 2, aa, d / 2, pp);

  if (n == 0)
    return n;

  /* write p-1 = q*2^s with q odd */
  for (q = p-1, s = 0; (q&1) == 0; q/=2, s++);

  if (s == 1) /* p = 3 (mod 4) */
    {
      /* solutions are +/-a^((p+1)/4) */
      for (i = 0; i < n; i++)
        {
          /* If d = 2 (mod 4), then every root of x^(d/2) = a (mod p) gives
             two roots of x^d = a (mod p) since a is a quadratic residue.
             However if d = 0 (mod 4), then a root of x^(d/2) = a (mod p) can
             give no root of x^d = a (mod p), thus we have to check the
             Legendre symbol again.
             Consider for example x^4 = 3 (mod 11): x^2 = 3 mod 11 has
             two roots (5 and 6), then x^2-5 mod 11 has two roots (4 and 7)
             but x^2-6 mod 11 has no root. */
          
          if ((d & 3) == 0 && mod_jacobi (rr[d/2+i], pp) != 1)
            continue;

          mod_pow_ul (rr[2*k], rr[d/2+i], (p + 1) >> 2, pp);
          mod_neg (rr[2*k+1], rr[2*k], pp);
          k ++;
        }
      return 2*k;
    }

  /* case p = 1 (mod 4). Uses Tonelli-Shanks, more precisely
     Algorithm from Table 1 in [1] */
  
  mod_init (zz, pp);
  mod_set1 (zz, pp);
  i = 1;
  do {
    /* zz is equal to i (mod pp) */
    mod_add1 (zz, zz, pp);
    i++;
  } while (!isprime_table[i] || mod_jacobi (zz, pp) != -1);
  mod_init (hh, pp);
  mod_init (aa, pp);
  mod_init (delta, pp);
  mod_init (dd, pp);
  mod_init (bb, pp);
  for (i = 0; i < n; i++)
    {
      /* For d divisible by 4, we have to check the Legendre symbol again.
         Consider for example x^4 = 10 (mod 13): x^2 = 10 (mod 13) has two
         roots (6 and 7) but none of them has a root mod 13. */
      if ((d & 3) == 0 && mod_jacobi (rr[d/2+i], pp) != 1)
        continue;

      mod_pow_ul (aa, zz, q, pp);
      mod_pow_ul (bb, rr[d/2+i], q, pp);
      mod_set1 (hh, pp);
      for (j = 1; j < s; j++)
        {
          mod_set (dd, bb, pp);
          for (l = 0; l < s - 1 - j; l++)
            mod_sqr (dd, dd, pp);
          if (!mod_is1 (dd, pp))
            {
              mod_mul (hh, hh, aa, pp);
              mod_sqr (aa, aa, pp);
              mod_mul (bb, bb, aa, pp);
            }
          else
            mod_sqr (aa, aa, pp);
        }
      mod_pow_ul (delta, rr[d/2+i], (q + 1) >> 1, pp);
      mod_mul (hh, hh, delta, pp);
      mod_set (rr[2*k], hh, pp);
      mod_neg (rr[2*k+1], rr[2*k], pp);
      k++;
    }

  mod_clear (hh, pp);
  mod_clear (aa, pp);
  mod_clear (delta, pp);
  mod_clear (dd, pp);
  mod_clear (bb, pp);
  mod_clearmod (pp);

  return 2*k;
}

/* Put in rr a root of x^3 = ddelta (mod pp), assuming one exists */
/* rr and ddelta overlapping is permissible. */
static void 
one_cubic_root (residue_t rr, residue_t ddelta, modulus_t pp)
{
  residue_t rho, a, aprime, b, h, d;
  uint64_t i, j, s, t, l;
  const uint64_t p = mod_getmod_ul(pp);

  /* when p = 2 (mod 3), then 1/3 = (2p-1)/3 mod (p-1), thus a cubic root
     is delta^((2p-1)/3) mod p. We rewrite exponent as (p+1)/3*2-1 to 
     avoid overflow */

  if ((p % 3) == 2) {
    mod_pow_ul (rr, ddelta, (p + 1) / 3 * 2 - 1, pp);
    return;
  }

  /* now p = 1 (mod 3), use Algorithm from Table 3 of [1] */
  s = (p - 1) / 3;
  t = 1;
  while ((s % 3) == 0)
    s /= 3, t++;
  mod_init (rho, pp);
  mod_init (a, pp);
  mod_init (aprime, pp);
  mod_init (b, pp);
  mod_init (h, pp);
  mod_init (d, pp);
  
  mod_set1 (rho, pp);
  for (i = 2; i < p; i++)
    {
      mod_add1 (rho, rho, pp);
      /* rho is equal to i (mod pp) */
      if (!isprime_table[i])
        continue;
      mod_pow_ul (a, rho, s, pp); /* a = rho^s */
      mod_set (aprime, a, pp);
      for (j = 0; j < t - 1; j++)
        mod_pow_ul (aprime, aprime, 3, pp);
      /* aprime = rho^(3^(t-1)*s) = rho^((p-1)/3) */
      if (!mod_is1 (aprime, pp))
        break;
    }
  mod_pow_ul (b, ddelta, s, pp);
  mod_set1 (h, pp);
  for (i = 1; i < t ; i++)
    {
      mod_set (d, b, pp);
      for (j = 0; j < t - 1 - i; j++)
        mod_pow_ul (d, d, 3, pp);
      if (mod_is1 (d, pp))
        {
          mod_pow_ul (a, a, 3, pp);
        }
      else if (mod_equal (d, aprime, pp))
        {
          mod_sqr (d, a, pp);
          mod_mul (h, h, d, pp);
          mod_mul (a, d, a, pp);
          mod_sqr (d, a, pp);
          mod_mul (b, b, d, pp);
        }
      else
        {
          mod_mul (h, h, a, pp);
          mod_pow_ul (a, a, 3, pp);
          mod_mul (b, b, a, pp);
        }
    }
  l = (s + 1) / 3;
  mod_set (d, ddelta, pp);
  mod_pow_ul (d, d, l, pp);
  mod_mul (h, h, d, pp);
  if (s == 3 * l + 1)
    mod_pow_ul (h, h, p - 2, pp);
  mod_set (rr, h, pp);

  mod_clear (rho, pp);
  mod_clear (a, pp);
  mod_clear (aprime, pp);
  mod_clear (b, pp);
  mod_clear (h, pp);
  mod_clear (d, pp);
}

/* Return a r-th root of delta (mod p), assuming one exists,
   using the algorithm from Table 4 in reference [1]. */
static uint64_t
one_rth_root (uint64_t r, uint64_t delta, uint64_t p)
{
  uint64_t rho, i, j, s, t, alpha, a, c, h, b, d;

  /* when p-1 is not divisible by r, there is exactly one root */
  rho = (p - 1) % r;
  if (rho != 0)
    {
      for (i = 1; (i * (p - 1) + 1) % r != 0; i++);
      return uint64_pow_mod (delta, (i * (p - 1) + 1) / r, p);
    }

  /* now p-1 is divisible by r */
  
  for (s = (p - 1) / r, t = 1; (s % r) == 0; s /= r, t++);
  for (rho = 2; rho < p; rho++)
    {
      c = uint64_pow_mod (rho, s, p);
      for (a = c, i = 0; i < t - 1; i++)
        a = uint64_pow_mod (a, r, p);
      if (a != 1)
        break;
    }
  alpha = invert_mod (r, s);
  b = uint64_pow_mod (delta, r * alpha - 1, p);
  h = 1;
  for (i = 1; i < t; i++)
    {
      d = b;
      for (j = 0; j < t - 1- i; j++)
        d = uint64_pow_mod (d, r, p);
      if (d == 1)
        c = uint64_pow_mod (c, r, p);
      else
        {
          
        }
    }
}

/* return 1 iff a is a cube mod p */
static int
is_cube (residue_t aa, modulus_t pp)
{
  if ((mod_getmod_ul(pp) % 3) == 1)
    {
      residue_t cc;
      int r;
      mod_init_noset0 (cc, pp);
      mod_pow_ul (cc, aa, (mod_getmod_ul(pp) - 1) / 3, pp);
      r = mod_is1(cc, pp);
      mod_clear (cc, pp);
      return r;
    }
  else /* p = 2 (mod 3) */
    return 1;
}


/* Roots of x^d = a (mod p) for d divisible by 3 and a 64-bit word */
static int
roots3 (residue_t *rr, residue_t aa, int d, modulus_t pp)
{
  uint64_t i, n;
  const uint64_t p = mod_getmod_ul(pp);

  if (!is_cube (aa, pp))
    return 0;

  /* find the roots of x^(d/3) = a (mod p) */
  n = mod_roots (rr + 2 * (d / 3), aa, d / 3, pp);

  if (n == 0)
    return n;

  if ((p % 3) == 1)
    {
      residue_t zeta, t;

      mod_init (zeta, pp);
      mod_init (t, pp);
      mod_set1 (t, pp);
      i = 1;

      do {
        mod_add1 (t, t, pp);
        i++;
        if (!isprime_table[i])
          continue;
        mod_pow_ul (zeta, t, (p - 1) / 3, pp);
      } while (mod_is1 (zeta, pp));

      /* zeta is a cubic root of 1 */
      for (i = 0; i < n; i++)
        {
          one_cubic_root (rr[3*i], rr[2 * (d / 3) + i], pp);
          mod_mul (rr[3*i+1], rr[3*i], zeta, pp);
          mod_mul (rr[3*i+2], rr[3*i+1], zeta, pp);
        }
      mod_clear (zeta, pp);
      mod_clear (t, pp);
      return 3*n;
    }
  else /* p = 2 (mod 3): exactly one root each */
    {
      for (i = 0; i < n; i++) {
        one_cubic_root (rr[i], rr[2 * (d / 3) + i], pp);
      }
      return n;
    }
}


static int 
mod_roots (residue_t *rr, residue_t aa, int d, modulus_t pp)
{
  if (d == 1)
    {
      mod_set (rr[0], aa, pp);
      return 1;
    }
  else if (d % 2 == 0) /* d is even */
    {
      return roots2 (rr, aa, d, pp);
    }
  else if (d % 3 == 0)
    {
      return roots3 (rr, aa, d, pp);
    }
  else
    abort ();
}


/* sort the roots r[0], ..., r[n-1] in increasing order */
static void
sort_roots (uint64_t *r, int n)
{
  int i, j;
  uint64_t t;

  for (i = 1; i < n; i++)
    {
      t = r[i];
      for (j = i; j > 0 && r[j-1] > t; j--)
        r[j] = r[j-1];
      r[j] = t;
    }
}


/* put in r[0], r[1], ... the roots of x^d = a (mod p),
   and return the number of roots.
   Assumes 0 <= a < p.
*/
int
roots_mod_uint64 (uint64_t *r, uint64_t a, int d, uint64_t p)
{
  mpz_t *f;
  int n = -1, nn = -1, i;
  uint64_t r2[10];
  const int do_both = 0; /* Compute both ways to test? */

  if (d == 1)
    {
      r[0] = a;
      return 1;
    }

  if (sizeof (unsigned long) == 8)
    {
      modulus_t pp;
      residue_t aa, rr[10];
      
      mod_initmod_ul (pp, p);
      mod_init (aa, pp);
      mod_set_ul (aa, a, pp);

      n = mod_roots (rr, aa, d, pp);

      for (i = 0; i < n; i++)
        r[i] = mod_get_ul (rr[i], pp);
      sort_roots (r, n);

      mod_clear (aa, pp);
      mod_clearmod (pp);
    }

  if (n == -1 || do_both) {
    f = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
    for (i = 0; i <= d; i++)
      mpz_init (f[i]);
    mpz_set_ui (f[d], 1);
    mpz_set_uint64 (f[0], p - a);
    nn = poly_roots_uint64 (r2, f, d, p);
    for (i = 0; i <= d; i++)
      mpz_clear (f[i]);
    free (f);
    sort_roots (r2, nn);
  }

  if (n != -1 && nn != -1) {
    /* If we ran both, compare to verify results */
    ASSERT_ALWAYS(n == nn);
    for (i = 0; i < n; i++) {
      ASSERT_ALWAYS(r[i] == r2[i]);
    }
  } else if (n == -1) {
    /* We ran poly_roots_uint64(). Copy result */
    n = nn;
    for (i = 0; i < n; i++) {
      r[i] = r2[i];
    }
  }

  return n;
}
