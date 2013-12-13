#include "cado.h"
#include <stdlib.h> /* for malloc and free */
#include "implicit_mpz_poly.h"
#include "gmp_aux.h"

static void mp_poly_evalz(mpz_t r, mpz_t * poly, int deg, mpz_srcptr a)
{
    int i;

    mpz_set(r, poly[deg]);
    for (i = deg - 1; i >= 0; i--) {
	mpz_mul(r, r, a);
	mpz_add(r, r, poly[i]);
    }
}

// Evaluate the derivative of poly at a.
static void
mp_poly_eval_diffz(mpz_t r, mpz_t * poly, int deg, mpz_srcptr a)
{
    int i;

    mpz_mul_ui(r, poly[deg], (unsigned long) deg);
    for (i = deg - 1; i >= 1; i--) {
	mpz_mul(r, r, a);
	mpz_addmul_ui(r, poly[i], (unsigned long) i);
    }
}

static void mp_poly_eval(mpz_t r, mpz_t * poly, int deg, long a)
{
    int i;

    mpz_set(r, poly[deg]);
    for (i = deg - 1; i >= 0; i--) {
	mpz_mul_si(r, r, a);
	mpz_add(r, r, poly[i]);
    }
}

// Evaluate the derivative of poly at a.
static void mp_poly_eval_diff(mpz_t r, mpz_t * poly, int deg, long a)
{
    int i;

    mpz_mul_ui(r, poly[deg], (unsigned long) deg);
    for (i = deg - 1; i >= 1; i--) {
	mpz_mul_si(r, r, a);
	mpz_addmul_ui(r, poly[i], (unsigned long) i);
    }
}

// assuming that pk is odd, that r is a root mod p^l, and that pk is a
// power of p equal to at most p^(2l), compute a lift of r mod pk.
int lift_root(mpz_t * f, int d, const unsigned long pk, unsigned long *r)
{
    unsigned long res;
    mpz_t aux, aux2, mp_p;

    mpz_init(aux);
    mpz_init(aux2);
    mpz_init_set_ui(mp_p, pk);
    mp_poly_eval_diff(aux, f, d, *r);
    mpz_mod_ui(aux, aux, pk);
    if (!mpz_invert(aux, aux, mp_p)) {
	return 0;
    }
    mp_poly_eval(aux2, f, d, *r);
    mpz_mod(aux2, aux2, mp_p);
    mpz_mul(aux2, aux2, aux);
    mpz_neg(aux2, aux2);
    mpz_add_ui(aux2, aux2, *r);
    mpz_mod(aux2, aux2, mp_p);

    res = mpz_get_ui(aux2);
    mpz_clear(aux);
    mpz_clear(aux2);
    mpz_clear(mp_p);

    *r = res;
    return 1;
}

int lift_rootz(mpz_t * f, int d, mpz_t pk, mpz_t r)
{
    mpz_t aux, aux2;

    mpz_init(aux);
    mpz_init(aux2);
    mp_poly_eval_diffz(aux, f, d, r);
    mpz_mod(aux, aux, pk);
    if (!mpz_invert(aux, aux, pk)) {
	return 0;
    }
    mp_poly_evalz(aux2, f, d, r);
    mpz_mod(aux2, aux2, pk);
    mpz_mul(aux2, aux2, aux);
    mpz_sub(r, r, aux2);
    mpz_mod(r, r, pk);
    mpz_clear(aux);
    mpz_clear(aux2);
    return 1;
}

/* return 0 if f and g (both of degree d) are equal, non-zero otherwise */
int
mp_poly_cmp (mpz_t *f, mpz_t *g, int d)
{
  int i;

  for (i = 0; i <= d; i++)
    if (mpz_cmp (f[i], g[i]))
      return 1; /* f and g differ */
  return 0;
}

/* v <- |f(i,j)|, where f is of degree d */
void
mp_poly_homogeneous_eval_siui (mpz_t v, mpz_t *f, const unsigned int d,
			       const int64_t i, const uint64_t j)
{
  unsigned int k;
  mpz_t jpow;

  mpz_init_set_ui (jpow, 1);
  mpz_set (v, f[d]);
  for (k = d; k-- > 0;)
    {
      mpz_mul_int64 (v, v, i);
      mpz_mul_uint64 (jpow, jpow, j);
      mpz_addmul (v, f[k], jpow);
    }
  mpz_abs (v, v); /* avoids problems with negative norms */
  mpz_clear (jpow);
}

/*  Homographic transform on polynomials */
/* Put in fij[] the coefficients of f'(i) = F(a0*i+a1, b0*i+b1).
   Assumes the coefficients of fij[] are initialized.
   Put in fijd[] a double-precision approximation of fij[]/q.
*/
void
mp_poly_homography (mpz_t *fij, mpz_t *f, const int d, int64_t H[4])
{
  int k, l;
  mpz_t *g; /* will contain the coefficients of (b0*i+b1)^l */
  mpz_t f0;

  for (k = 0; k <= d; k++)
    mpz_set (fij[k], f[k]);

  g = malloc ((d + 1) * sizeof (mpz_t));
  for (k = 0; k <= d; k++)
    mpz_init (g[k]);
  mpz_init (f0);

  /* Let h(x) = quo(f(x), x), then F(x,y) = H(x,y)*x + f0*y^d, thus
     F(a0*i+a1, b0*i+b1) = H(a0*i+a1, b0*i+b1)*(a0*i+a1) + f0*(b0*i+b1)^d.
     We use that formula recursively. */

  mpz_set_ui (g[0], 1); /* g = 1 */

  for (k = d - 1; k >= 0; k--)
    {
      /* invariant: we have already translated coefficients of degree > k,
         in f[k+1..d], and g = (b0*i+b1)^(d - (k+1)), with coefficients in
         g[0..d - (k+1)]:
         f[k] <- a1*f[k+1]
         ...
         f[l] <- a0*f[l]+a1*f[l+1] for k < l < d
         ...
         f[d] <- a0*f[d] */
      mpz_swap (f0, fij[k]); /* save the new constant coefficient */
      mpz_mul_si (fij[k], fij[k + 1], H[2]);
      for (l = k + 1; l < d; l++)
        {
          mpz_mul_si (fij[l], fij[l], H[0]);
          mpz_addmul_si (fij[l], fij[l + 1], H[2]);
        }
      mpz_mul_si (fij[d], fij[d], H[0]);

      /* now compute (b0*i+b1)^(d-k) from the previous (b0*i+b1)^(d-k-1):
         g[d-k] = b0*g[d-k-1]
         ...
         g[l] = b1*g[l]+b0*g[l-1] for 0 < l < d-k
         ...
         g[0] = b1*g[0]
      */
      mpz_mul_si (g[d - k], g[d - k - 1], H[1]);
      for (l = d - k - 1; l > 0; l--)
        {
          mpz_mul_si (g[l], g[l], H[3]);
          mpz_addmul_si (g[l], g[l-1], H[1]);
        }
      mpz_mul_si (g[0], g[0], H[3]);

      /* now g has degree d-k, and we add f0*g */
      for (l = k; l <= d; l++)
        mpz_addmul (fij[l], g[l - k], f0);
    }

  mpz_clear (f0);
  for (k = 0; k <= d; k++)
    mpz_clear (g[k]);
  free (g);
}

/* put in c the content of f */
void
mp_poly_content (mpz_t c, mpz_t *f, const int d)
{
  int i;

  mpz_set (c, f[0]);
  for (i = 1; i <= d; i++)
    mpz_gcd (c, c, f[i]);
  mpz_abs (c, c);
}
