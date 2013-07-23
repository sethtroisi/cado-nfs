#include "cado.h"
#include "las-config.h"
#include "las-types.h"
#include "las-coordinates.h"
#include "gmp_aux.h"

/*  Conversions between different representations for sieve locations:
 *
 * Coordinate systems for identifying an (a,b) pair.
 *
 * (a, b): This is the one from textbooks. We always have b>=0 (only free
 *         relations have b=0, and we don't see them here).
 * (i, j): For a given special q, this is a point in the q-lattice. Given
 *         the lattice basis given by (a0 b0 a1 b1), this corresponds to
 *         the (a,b) pair equal to i*(a0,b0)+j*(a1,b1). By construction
 *         this should lead to one of the norms having si->doing->p as a factor.
 *         i is within [-I/2, I/2[, and j is within [1, J[
 * (N, x): bucket number N, location x. N is within [0,si->nb_buckets[
 *         and x within [0,bucket_region[ ; we have:
 *         N*bucket_region+x == (I/2+i)+j*I
 *
 * There are some change of coordinate functions at the end of this file.
 * The coordinate system (N, x) is almost always referred to as just x,
 * where x is a 32-bit value equal to N*bucket_region+x. Thus from this
 * wide x, it is possible to recover both N and (the short) x.
 */
void xToIJ(int *i, unsigned int *j, const unsigned int X, sieve_info_srcptr si)
{
    *i = (X % (si->I)) - (si->I >> 1);
    *j = X / si->I;
}

void NxToIJ(int *i, unsigned int *j, const unsigned int N, const unsigned int x, sieve_info_srcptr si)
{
    unsigned int X = x + (N << LOG_BUCKET_REGION);
    return xToIJ(i, j, X, si);
}

void IJTox(unsigned int * x, int i, unsigned int j, sieve_info_srcptr si)
{
    *x = i + (si->I)*j + (si->I>>1);
}

void IJToNx(unsigned int *N, unsigned int * x, int i, unsigned int j, sieve_info_srcptr si)
{
    IJTox(x, i, j, si);
    *N = *x >> LOG_BUCKET_REGION;
    *x &= ((1 << LOG_BUCKET_REGION) - 1);
}

void IJToAB(int64_t *a, uint64_t *b, const int i, const unsigned int j, 
       sieve_info_srcptr si)
{
    int64_t s, t;
    s = (int64_t)i * (int64_t) si->a0 + (int64_t)j * (int64_t)si->a1;
    t = (int64_t)i * (int64_t) si->b0 + (int64_t)j * (int64_t)si->b1;
    if (t >= 0)
      {
        *a = s;
        *b = t;
      }
    else
      {
        *a = -s;
        *b = -t;
      }
}

int ABToIJ(int *i, unsigned int *j, const int64_t a, const uint64_t b, sieve_info_srcptr si)
{
    /* Both a,b and the coordinates of the lattice basis can be quite
     * large. However the result should be small.
     */
    mpz_t za,zb,ii,jj,a0,b0,a1,b1;
    mpz_init(ii);
    mpz_init(jj);
    mpz_init(za); mpz_init(zb);
    mpz_init(a0); mpz_init(b0);
    mpz_init(a1); mpz_init(b1);
    mpz_set_int64(za, a); mpz_set_uint64(zb, b);
    mpz_set_int64(a0, si->a0);
    mpz_set_int64(b0, si->b0);
    mpz_set_int64(a1, si->a1);
    mpz_set_int64(b1, si->b1);
    int ok = 1;
    mpz_mul(ii, za, b1); mpz_submul(ii, zb, a1);
    mpz_mul(jj, zb, a0); mpz_submul(jj, za, b0);
    /*
    int64_t ii =   a * (int64_t) si->b1 - b * (int64_t)si->a1;
    int64_t jj = - a * (int64_t) si->b0 + b * (int64_t)si->a0;
    */
    if (!mpz_divisible_p(ii, si->doing->p)) ok = 0;
    if (!mpz_divisible_p(jj, si->doing->p)) ok = 0;
    mpz_divexact(ii, ii, si->doing->p);
    mpz_divexact(jj, jj, si->doing->p);
    if (mpz_sgn(jj) < 0 || (mpz_sgn(jj) == 0 && mpz_sgn(ii) < 0)) {
        mpz_neg(ii, ii);
        mpz_neg(jj, jj);
    }
    *i = mpz_get_si(ii);
    *j = mpz_get_ui(jj);
    mpz_clear(a0);
    mpz_clear(b0);
    mpz_clear(a1);
    mpz_clear(b1);
    mpz_clear(za);
    mpz_clear(zb);
    mpz_clear(ii);
    mpz_clear(jj);
    return ok;
}
int ABTox(unsigned int *x, const int64_t a, const uint64_t b, sieve_info_srcptr si)
{
    int i;
    unsigned int j;
    if (!ABToIJ(&i, &j, a, b, si)) return 0;
    IJTox(x, a, b, si);
    return 1;
}

int ABToNx(unsigned int * N, unsigned int *x, const int64_t a, const uint64_t b, sieve_info_srcptr si)
{
    int i;
    unsigned int j;
    if (!ABToIJ(&i, &j, a, b, si)) return 0;
    IJToNx(N, x, a, b, si);
    return 1;
}
/*  */

