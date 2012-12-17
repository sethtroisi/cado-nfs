/* Routines to solve x^k = a (mod p) */

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
  uint64_t q, s, z, m, n, i, j, k = 0;
  modulusredcul_t pp;
  residueredcul_t cc0, aa, RR, tt, cc, uu, bb;

  if (legendre (a, p) != 1)
    return 0;

  /* find the roots of x^(d/2) = a (mod p) */
  n = roots_mod_uint64 (r + d / 2, a, d / 2, p);

  /* then use Tonelli-Shanks */

  /* write p-1 = q*2^s with q odd */
  for (q = p-1, s = 0; (q&1) == 0; q/=2, s++);

  if (s == 1) /* p = 3 (mod 4) */
    {
      /* solutions are +/-a^((p+1)/4) */
      for (i = 0; i < n; i++)
        {
          if (legendre (r[d/2+i], p) != 1)
            continue;
          r[2*k] = uint64_pow_mod (r[d/2+i], (p + 1) >> 2, p);
          r[2*k+1] = p - r[2*k];
          k ++;
        }
      return 2*n;
    }

  /* case p = 1 (mod 4) */
  for (z = 2; legendre (z, p) != -1; z++);
  modredcul_initmod_ul (pp, p);
  modredcul_init (cc0, pp);
  modredcul_init (aa, pp);
  modredcul_init (RR, pp);
  modredcul_init (tt, pp);
  modredcul_init (cc, pp);
  modredcul_init (uu, pp);
  modredcul_init (bb, pp);
  modredcul_set_ul (cc0, z, pp);
  modredcul_pow_ul (cc0, cc0, q, pp);
  for (i = 0; i < n; i++)
    {
      if (legendre (r[d/2+i], p) != 1)
        continue;

      modredcul_set_ul (aa, r[d/2+i], pp);
      modredcul_pow_ul (RR, aa, (q + 1) >> 1, pp);
      modredcul_pow_ul (tt, aa, q, pp);
      m = s;
      modredcul_set (cc, cc0, pp);
      do
        {
          if (modredcul_is1 (tt, pp))
            {
              r[2*k] = modredcul_get_ul (RR, pp);
              r[2*k+1] = p - r[2*k];
              k ++;
              break;
            }
          for (j = 0, modredcul_set (uu, tt, pp); modredcul_is1 (uu, pp) == 0;
               modredcul_sqr (uu, uu, pp), j++);
          ASSERT_ALWAYS(j < m);
          for (m = m - j - 1, modredcul_set (bb, cc, pp); m > 0;
               modredcul_sqr (bb, bb, pp), m--);
          modredcul_mul (RR, RR, bb, pp);
          modredcul_sqr (cc, bb, pp);
          modredcul_mul (tt, tt, cc, pp);
          m = j;
        }
      while (1);
    }

  modredcul_clear (cc0, pp);
  modredcul_clear (aa, pp);
  modredcul_clear (RR, pp);
  modredcul_clear (tt, pp);
  modredcul_clear (cc, pp);
  modredcul_clear (uu, pp);
  modredcul_clear (bb, pp);
  modredcul_clearmod (pp);

  return 2*k;
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

  if ((d & 1) == 0 && sizeof (unsigned long) == 8) /* d is even */
    return roots2 (r, a, d, p);

  f = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
  for (i = 0; i <= d; i++)
    mpz_init (f[i]);
  mpz_set_ui (f[d], 1);
  mpz_set_uint64 (f[0], p - a);
  n = poly_roots_uint64 (r, f, d, p);
  for (i = 0; i <= d; i++)
    mpz_clear (f[i]);
  free (f);
  return n;
}
