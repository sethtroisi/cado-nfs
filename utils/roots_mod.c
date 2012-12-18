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

/* this seems faster than calling mpz_legendre (which needs converting both
   a and p) or mpz_ui_kronecker (which needs converting p) */
static int
legendre (unsigned long a, unsigned long p)
{
  mpz_t aa;
  int ret;

  mpz_init_set_ui (aa, a);
  ret = mpz_kronecker_ui (aa, p);
  mpz_clear (aa);
  return ret;
}

/* Return b^e mod m. Assumes a 64-bit word. */
uint64_t
uint64_pow_mod (uint64_t b, uint64_t e, uint64_t m)
{
  modulusredcul_t mm;
  residueredcul_t rr, bb;
  uint64_t r;

  modredcul_initmod_ul (mm, m);
  modredcul_init (rr, mm);
  modredcul_init (bb, mm);
  modredcul_set_ul (bb, b, mm);
  modredcul_pow_ul (rr, bb, e, mm);
  r = modredcul_get_ul (rr, mm);
  modredcul_clear (rr, mm);
  modredcul_clear (bb, mm);
  modredcul_clearmod (mm);
  return r;
}

/* Roots of x^d = a (mod p) for d even and a 64-bit word */
static int
roots2 (uint64_t *r, uint64_t a, int d, uint64_t p)
{
  uint64_t q, s, z, n, i, j, k = 0, l;
  modulusredcul_t pp;
  residueredcul_t hh, aa, delta, dd, bb;

  if (legendre (a, p) != 1)
    return 0;

  /* find the roots of x^(d/2) = a (mod p) */
  n = roots_mod_uint64 (r + d / 2, a, d / 2, p);

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
          if ((d & 3) == 0 && legendre (r[d/2+i], p) != 1)
            continue;
          r[2*k] = uint64_pow_mod (r[d/2+i], (p + 1) >> 2, p);
          r[2*k+1] = p - r[2*k];
          k ++;
        }
      return 2*n;
    }

  /* case p = 1 (mod 4). Uses Tonelli-Shanks, more precisely
     Algorithm from Table 1 in [1] */
  for (z = 2; legendre (z, p) != -1; z++);
  modredcul_initmod_ul (pp, p);
  modredcul_init (hh, pp);
  modredcul_init (aa, pp);
  modredcul_init (delta, pp);
  modredcul_init (dd, pp);
  modredcul_init (bb, pp);
  for (i = 0; i < n; i++)
    {
      /* For d divisible by 4, we have to check the Legendre symbol again.
         Consider for example x^4 = 10 (mod 13): x^2 = 10 (mod 13) has two
         roots (6 and 7) but none of them has a root mod 13. */
      if ((d & 3) == 0 && legendre (r[d/2+i], p) != 1)
        continue;

      modredcul_set_ul (delta, r[d/2+i], pp);
      modredcul_set_ul (aa, z, pp);
      modredcul_pow_ul (aa, aa, q, pp);
      modredcul_pow_ul (bb, delta, q, pp);
      modredcul_set1 (hh, pp);
      for (j = 1; j < s; j++)
        {
          modredcul_set (dd, bb, pp);
          for (l = 0; l < s - 1 - j; l++)
            modredcul_sqr (dd, dd, pp);
          if (modredcul_is1 (dd, pp) == 0)
            {
              modredcul_mul (hh, hh, aa, pp);
              modredcul_sqr (aa, aa, pp);
              modredcul_mul (bb, bb, aa, pp);
            }
          else
            modredcul_sqr (aa, aa, pp);
        }
      modredcul_pow_ul (delta, delta, (q + 1) >> 1, pp);
      modredcul_mul (hh, hh, delta, pp);
      r[2*k] = modredcul_get_ul (hh, pp);
      r[2*k+1] = p - r[2*k];
      k++;
    }

  modredcul_clear (hh, pp);
  modredcul_clear (aa, pp);
  modredcul_clear (delta, pp);
  modredcul_clear (dd, pp);
  modredcul_clear (bb, pp);
  modredcul_clearmod (pp);

  return 2*k;
}

/* return a root of x^3 = delta (mod p), assuming one exists */
static uint64_t
one_cubic_root (uint64_t delta, uint64_t p)
{
  modulusredcul_t pp;
  residueredcul_t rho, a, aprime, b, h, d;
  uint64_t i, j, s, t, l, r;

  /* when p = 2 (mod 3), then 1/3 = (2p-1)/3 mod (p-1), thus a cubic root
     is delta^((2p-1)/3) mod p */

  if ((p % 3) == 2)
    return uint64_pow_mod (delta, (2 * p - 1) / 3, p);

  /* now p = 1 (mod 3), use Algorithm from Table 3 of [1] */
  s = (p - 1) / 3;
  t = 1;
  while ((s % 3) == 0)
    s /= 3, t++;
  modredcul_initmod_ul (pp, p);
  modredcul_init (rho, pp);
  modredcul_init (a, pp);
  modredcul_init (aprime, pp);
  modredcul_init (b, pp);
  modredcul_init (h, pp);
  modredcul_init (d, pp);
  for (i = 2; i < p; i++)
    {
      modredcul_set_ul (rho, i, pp);
      modredcul_pow_ul (a, rho, s, pp); /* a = rho^s */
      modredcul_set (aprime, a, pp);
      for (j = 0; j < t - 1; j++)
        modredcul_pow_ul (aprime, aprime, 3, pp);
      /* aprime = rho^(3^(t-1)*s) = rho^((p-1)/3) */
      if (modredcul_is1 (aprime, pp) == 0)
        break;
    }
  modredcul_set_ul (b, delta, pp);
  modredcul_pow_ul (b, b, s, pp);
  modredcul_set1 (h, pp);
  for (i = 1; i < t ; i++)
    {
      modredcul_set (d, b, pp);
      for (j = 0; j < t - 1 - i; j++)
        modredcul_pow_ul (d, d, 3, pp);
      if (modredcul_is1 (d, pp))
        {
          modredcul_pow_ul (a, a, 3, pp);
        }
      else if (modredcul_intequal (d, aprime))
        {
          modredcul_sqr (d, a, pp);
          modredcul_mul (h, h, d, pp);
          modredcul_mul (a, d, a, pp);
          modredcul_sqr (d, a, pp);
          modredcul_mul (b, b, d, pp);
        }
      else
        {
          modredcul_mul (h, h, a, pp);
          modredcul_pow_ul (a, a, 3, pp);
          modredcul_mul (b, b, a, pp);
        }
    }
  l = (s + 1) / 3;
  modredcul_set_ul (d, delta, pp);
  modredcul_pow_ul (d, d, l, pp);
  modredcul_mul (h, h, d, pp);
  if (s == 3 * l + 1)
    modredcul_pow_ul (h, h, p - 2, pp);
  r = modredcul_get_ul (h, pp);

  modredcul_clear (rho, pp);
  modredcul_clear (a, pp);
  modredcul_clear (aprime, pp);
  modredcul_clear (b, pp);
  modredcul_clear (h, pp);
  modredcul_clear (d, pp);
  modredcul_clearmod (pp);

  return r;
}

/* return 1 iff a is a cube mod p */
static int
is_cube (uint64_t a, uint64_t p)
{
  uint64_t c;

  if ((p % 3) == 1)
    {
      c = uint64_pow_mod (a, (p - 1) / 3, p);
      return c == 1;
    }
  else /* p = 2 (mod 3) */
    return 1;
}

/* Roots of x^d = a (mod p) for d divisible by 3 and a 64-bit word */
static int
roots3 (uint64_t *r, uint64_t a, int d, uint64_t p)
{
  uint64_t zeta = 0, i, n;

  if (is_cube (a, p) == 0)
    return 0;

  /* find the roots of x^(d/3) = a (mod p) */
  n = roots_mod_uint64 (r + 2 * (d / 3), a, d / 3, p);

  if (n == 0)
    return n;

  if ((p % 3) == 1)
    {
      modulusredcul_t pp;
      residueredcul_t z, t;

      for (i = 2; i < p; i++)
        if ((zeta = uint64_pow_mod (i, (p - 1) / 3, p)) != 1)
          break;
      /* zeta is a cubic root of 1 */
      modredcul_initmod_ul (pp, p);
      modredcul_init (z, pp);
      modredcul_init (t, pp);
      modredcul_set_ul (z, zeta, pp);
      for (i = 0; i < n; i++)
        {
          r[3*i] = one_cubic_root (r[2 * (d / 3) + i], p);
          modredcul_set_ul (t, r[3*i], pp);
          modredcul_mul (t, t, z, pp);
          r[3*i+1] = modredcul_get_ul (t, pp);
          modredcul_mul (t, t, z, pp);
          r[3*i+2] = modredcul_get_ul (t, pp);
        }
      modredcul_clear (z, pp);
      modredcul_clear (t, pp);
      modredcul_clearmod (pp);
      return 3*n;
    }
  else /* p = 2 (mod 3): only one root */
    {
      for (i = 0; i < n; i++)
        r[i] = one_cubic_root (r[2 * (d / 3) + i], p);
      return n;
    }
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
  int n, i;

  if (d == 1)
    {
      r[0] = a;
      return 1;
    }

  if (sizeof (unsigned long) == 8)
    {
      if ((d & 1) == 0) /* d is even */
        {
          n = roots2 (r, a, d, p);
          sort_roots (r, n);
          goto sort_and_exit;
        }

      if ((d % 3) == 0)
        {
          n = roots3 (r, a, d, p);
          sort_roots (r, n);
          goto sort_and_exit;
        }
    }

  f = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
  for (i = 0; i <= d; i++)
    mpz_init (f[i]);
  mpz_set_ui (f[d], 1);
  mpz_set_uint64 (f[0], p - a);
  n = poly_roots_uint64 (r, f, d, p);
  for (i = 0; i <= d; i++)
    mpz_clear (f[i]);
  free (f);

 sort_and_exit:
  sort_roots (r, n);
  return n;
}
