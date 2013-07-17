/* auxiliary routines on GMP data-types */

#include "cado.h"
#include <stdint.h>
#include "gmp_aux.h"
#include "macros.h"

/* Set z to q. Warning: on 32-bit machines, we cannot use mpz_set_ui! */
void
mpz_set_uint64 (mpz_t z, uint64_t q)
{
  if (sizeof (unsigned long) == 8)
    mpz_set_ui (z, (unsigned long) q);
  else
    {
      ASSERT_ALWAYS (sizeof (unsigned long) == 4);
      mpz_set_ui (z, (unsigned long) (q >> 32));
      mpz_mul_2exp (z, z, 32);
      /* The & here should be optimized into a direct cast to a 32-bit
       * register in most cases (TODO: check) */
      mpz_add_ui (z, z, (unsigned long) (q & 4294967295UL));
    }
}

void
mpz_set_int64 (mpz_t z, int64_t q)
{
  if (sizeof (long) == 8)
    mpz_set_si (z, (long) q);
  else if (q >= 0)
    mpz_set_uint64 (z, q);
  else
    {
      mpz_set_uint64 (z, -q);
      mpz_neg (z, z);
    }
}

uint64_t
mpz_get_uint64 (mpz_srcptr z)
{
    uint64_t q;

    if (sizeof (unsigned long) == 8)
        q = mpz_get_ui (z);
    else
    {
        ASSERT_ALWAYS (sizeof (unsigned long) == 4);
        ASSERT_ALWAYS (sizeof (mp_limb_t) == 4);
        ASSERT_ALWAYS (GMP_LIMB_BITS == 32);
        q = mpz_get_ui (z); /* get the low word of z */
        q += ((uint64_t) mpz_getlimbn(z,1)) << 32;
    }
    return q;
}

int64_t
mpz_get_int64 (mpz_srcptr z)
{
    return mpz_get_uint64(z) * (int64_t) mpz_sgn(z);
}

int mpz_fits_int64_p(mpz_srcptr z)
{
    int l = mpz_sizeinbase(z, 2);
    if (l <= 63) return 1;
    /* Also accept -2^63, which is INT64_MIN */
    if (mpz_sgn(z) < 0 && l == 64  && mpz_scan1(z, 0) == 63) return 1;
    return 0;
}

int mpz_fits_uint64_p(mpz_srcptr z)
{
    ASSERT_ALWAYS(mpz_sgn(z) >= 0);
    return mpz_sizeinbase(z, 2) <= 64;
}

void
mpz_mul_uint64 (mpz_t a, mpz_srcptr b, uint64_t c)
{
  if (sizeof (unsigned long) == 8)
    mpz_mul_ui (a, b, (unsigned long) c);
  else
    {
      mpz_t d;
      mpz_init (d);
      mpz_set_uint64 (d, c);
      mpz_mul (a, b, d);
      mpz_clear (d);
    }
}

void
mpz_mul_int64 (mpz_t a, mpz_srcptr b, int64_t c)
{
  if (sizeof (long) == 8)
    mpz_mul_si (a, b, (long) c);
  else
    {
      mpz_t d;
      mpz_init (d);
      mpz_set_int64 (d, c);
      mpz_mul (a, b, d);
      mpz_clear (d);
    }
}

void
mpz_addmul_int64 (mpz_t a, mpz_srcptr b, int64_t c)
{
  if (sizeof (long) == 8)
    mpz_addmul_si (a, b, (long) c);
  else
    {
      mpz_t d;
      mpz_init (d);
      mpz_set_int64 (d, c);
      mpz_addmul (a, b, d);
      mpz_clear (d);
    }
}

/* returns the smallest prime > q */
uint64_t
uint64_nextprime (uint64_t q)
{
  mpz_t p;

  mpz_init (p);
  mpz_set_uint64 (p, q);
  mpz_nextprime (p, p);
  q = mpz_get_uint64 (p);
  mpz_clear (p);
  return q;
}

/* same as above, for an unsigned long */
unsigned long
ulong_nextprime (unsigned long q)
{
  mpz_t p;

  mpz_init (p);
  mpz_set_ui (p, q);
  mpz_nextprime (p, p);
  ASSERT_ALWAYS (mpz_fits_ulong_p (p));
  q = mpz_get_ui (p);
  mpz_clear (p);
  return q;
}

#define REPS 1 /* number of Miller-Rabin tests in isprime */

int
isprime (unsigned long p)
{
  mpz_t P;
  int res;
  
  mpz_init_set_ui (P, p);
  res = mpz_probab_prime_p (P, REPS);
  mpz_clear (P);
  return res;
}

/* return the number of bits of p, counting from the least significant end */
int nbits (uintmax_t p)
{
  int k;

  for (k = 0; p != 0; p >>= 1, k ++);
  return k;
}

/* q <- n/d rounded to nearest, assuming d <> 0
   r <- n - q*d
*/
void
mpz_ndiv_qr (mpz_t q, mpz_t r, mpz_t n, mpz_t d)
{
  int s;

  ASSERT (mpz_cmp_ui (d, 0) != 0);
  mpz_fdiv_qr (q, r, n, d); /* round towards -inf, r has same sign as d */
  mpz_mul_2exp (r, r, 1);
  s = mpz_cmpabs (r, d);
  mpz_div_2exp (r, r, 1);
  if (s > 0) /* |r| > |d|/2 */
    {
      mpz_add_ui (q, q, 1);
      mpz_sub (r, r, d);
    }
}

void
mpz_ndiv_q (mpz_t q, mpz_t n, mpz_t d)
{
  mpz_t r;

  mpz_init (r);
  mpz_ndiv_qr (q, r, n, d);
  mpz_clear (r);
}


