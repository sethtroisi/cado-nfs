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
  modulus_t mm;
  residue_t rr, bb;
  uint64_t r;

  mod_initmod_ul (mm, m);
  mod_init (rr, mm);
  mod_init (bb, mm);
  mod_set_ul (bb, b, mm);
  mod_pow_ul (rr, bb, e, mm);
  r = mod_get_ul (rr, mm);
  mod_clear (rr, mm);
  mod_clear (bb, mm);
  mod_clearmod (mm);
  return r;
}

/* Roots of x^d = a (mod p) for d even and a 64-bit word */
static int
roots2 (uint64_t *r, uint64_t a, int d, uint64_t p)
{
  uint64_t q, s, z, n, i, j, k = 0, l;
  modulus_t pp;
  residue_t hh, aa, delta, dd, bb;

  if (legendre (a, p) != 1)
    return 0;

  /* find the roots of x^(d/2) = a (mod p) */
  n = roots_mod_uint64 (r + d / 2, a, d / 2, p);

  if (n == 0)
    return n;

  /* write p-1 = q*2^s with q odd */
  for (q = p-1, s = 0; (q&1) == 0; q/=2, s++);

  mod_initmod_ul (pp, p);
  if (s == 1) /* p = 3 (mod 4) */
    {
      /* solutions are +/-a^((p+1)/4) */
      for (i = 0; i < n; i++)
        {
          residue_t rr, bb;
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

          // r[2*k] = uint64_pow_mod (r[d/2+i], (p + 1) >> 2, p);
          mod_init (rr, pp);
          mod_init (bb, pp);
          mod_set_ul (bb, r[d/2+i], pp);
          mod_pow_ul (rr, bb, (p + 1) >> 2, pp);
          r[2*k] = mod_get_ul (rr, pp);
          mod_clear (rr, pp);
          mod_clear (bb, pp);
          r[2*k+1] = p - r[2*k];
          k ++;
        }
      mod_clearmod (pp);
      return 2*n;
    }

  /* case p = 1 (mod 4). Uses Tonelli-Shanks, more precisely
     Algorithm from Table 1 in [1] */
  for (z = 2; legendre (z, p) != -1; z++);
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
      if ((d & 3) == 0 && legendre (r[d/2+i], p) != 1)
        continue;

      mod_set_ul (delta, r[d/2+i], pp);
      mod_set_ul (aa, z, pp);
      mod_pow_ul (aa, aa, q, pp);
      mod_pow_ul (bb, delta, q, pp);
      mod_set1 (hh, pp);
      for (j = 1; j < s; j++)
        {
          mod_set (dd, bb, pp);
          for (l = 0; l < s - 1 - j; l++)
            mod_sqr (dd, dd, pp);
          if (mod_is1 (dd, pp) == 0)
            {
              mod_mul (hh, hh, aa, pp);
              mod_sqr (aa, aa, pp);
              mod_mul (bb, bb, aa, pp);
            }
          else
            mod_sqr (aa, aa, pp);
        }
      mod_pow_ul (delta, delta, (q + 1) >> 1, pp);
      mod_mul (hh, hh, delta, pp);
      r[2*k] = mod_get_ul (hh, pp);
      r[2*k+1] = p - r[2*k];
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
  residue_t rho, a, aprime, b, h, d, one;
  uint64_t i, j, s, t, l;

  /* when p = 2 (mod 3), then 1/3 = (2p-1)/3 mod (p-1), thus a cubic root
     is delta^((2p-1)/3) mod p. We rewrite exponent as (p+1)/3*2-1 to 
     avoid overflow */

  if ((mod_getmod_ul(pp) % 3) == 2) {
    mod_pow_ul (rr, ddelta, (mod_getmod_ul(pp) + 1) / 3 * 2 - 1, pp);
    return;
  }

  /* now p = 1 (mod 3), use Algorithm from Table 3 of [1] */
  s = (mod_getmod_ul(pp) - 1) / 3;
  t = 1;
  while ((s % 3) == 0)
    s /= 3, t++;
  mod_init (rho, pp);
  mod_init (a, pp);
  mod_init (aprime, pp);
  mod_init (b, pp);
  mod_init (h, pp);
  mod_init (d, pp);
  mod_init (one, pp);
  mod_set1 (one, pp);

  mod_add(rho, one, one, pp);
  for (i = 2; i < mod_getmod_ul(pp); i++)
    {
      mod_pow_ul (a, rho, s, pp); /* a = rho^s */
      mod_set (aprime, a, pp);
      for (j = 0; j < t - 1; j++)
        mod_pow_ul (aprime, aprime, 3, pp);
      /* aprime = rho^(3^(t-1)*s) = rho^((p-1)/3) */
      if (mod_is1 (aprime, pp) == 0)
        break;
      mod_add (rho, rho, one, pp);
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
    mod_pow_ul (h, h, mod_getmod_ul(pp) - 2, pp);
  mod_set (rr, h, pp);

  mod_clear (rho, pp);
  mod_clear (a, pp);
  mod_clear (aprime, pp);
  mod_clear (b, pp);
  mod_clear (h, pp);
  mod_clear (d, pp);
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
roots3 (uint64_t *r, uint64_t a, int d, uint64_t p)
{
  uint64_t i, n;
  modulus_t pp;
  residue_t aa;

  mod_initmod_ul (pp, p);
  mod_init (aa, pp);
  mod_set_ul (aa, a, pp);

  if (is_cube (aa, pp) == 0) {
    mod_clear (aa, pp);
    mod_clearmod (pp);
    return 0;
  }

  /* find the roots of x^(d/3) = a (mod p) */
  n = roots_mod_uint64 (r + 2 * (d / 3), a, d / 3, p);

  if (n == 0) {
    mod_clear (aa, pp);
    mod_clearmod (pp);
    return n;
  }

  if ((p % 3) == 1)
    {
      residue_t zeta, t, one;

      mod_init (zeta, pp);
      mod_init (t, pp);
      mod_init (one, pp);
      mod_set1 (one, pp);
      mod_add (t, one, one, pp);

      for (i = 2; i < p; i++) {
        mod_pow_ul (zeta, t, (p - 1) / 3, pp);
        if (!mod_is1 (zeta, pp))
          break;
        mod_add (t, t, one, pp);
      }
      /* zeta is a cubic root of 1 */
      for (i = 0; i < n; i++)
        {
          mod_set_ul (t, r[2 * (d / 3) + i], pp);
          one_cubic_root (t, t, pp);
          r[3*i] = mod_get_ul (t, pp);
          mod_mul (t, t, zeta, pp);
          r[3*i+1] = mod_get_ul (t, pp);
          mod_mul (t, t, zeta, pp);
          r[3*i+2] = mod_get_ul (t, pp);
        }
      mod_clear (zeta, pp);
      mod_clear (t, pp);
      mod_clear (one, pp);
      mod_clear (aa, pp);
      mod_clearmod (pp);
      return 3*n;
    }
  else /* p = 2 (mod 3): only one root */
    {
      residue_t t;
      mod_init (t, pp);
      for (i = 0; i < n; i++) {
        mod_set_ul(t, r[2 * (d / 3) + i], pp);
        one_cubic_root (t, t, pp);
        r[i] = mod_get_ul (t, pp);
      }
      mod_clear (t, pp);
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
      if ((d & 1) == 0) /* d is even */
        {
          n = roots2 (r, a, d, p);
          sort_roots (r, n);
        }
      else if ((d % 3) == 0)
        {
          n = roots3 (r, a, d, p);
          sort_roots (r, n);
        }
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

 exit:
  
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
