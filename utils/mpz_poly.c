/**
   These files (mpz_poly.*) implement arithmetics of polynomials whose
   coefficients are in multiprecision integers (using mpz_t from GNU MP).
   We use them in sqrt/algsqrt.c to represent rings of integers.

   Please see the file README_POLYNOMIALS for more details and
   comparisons to other poly formats.
*/


#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "utils.h"
#include "portability.h"

#ifndef max
#define max(a,b) ((a)<(b) ? (b) : (a))
#endif


/* -----------------
   Static functions
   ----------------- */


/* normalize h so that h->coeff[deg(h)] <> 0 mod p */
static void
mpz_poly_normalize_modp (mpz_poly_t h, const mpz_t p)
{
  int dh = h->deg;
  while (dh >= 0 && mpz_divisible_p (h->coeff[dh],p))
    dh --;
  h->deg = dh;
}

/* Return f=g*h, where g has degree r, and h has degree s. */
static int
mpz_poly_mul_basecase (mpz_t *f, mpz_t *g, int r, mpz_t *h, int s) {
  int i, j;
  for (i = 0; i <= r + s; i++)
    mpz_set_ui (f[i], 0);
  for (i = 0; i <= r; ++i)
    for (j = 0; j <= s; ++j)
      mpz_addmul(f[i+j], g[i], h[j]);
  return r + s;
}

/* Given f[0]...f[t] that contain respectively f(0), ..., f(t),
   put in f[0]...f[t] the coefficients of f. Assumes t <= MAX_T. */
static void
mpz_poly_mul_tc_interpolate (mpz_t *f, int t) {
#define MAX_T 13
  uint64_t M[MAX_T+1][MAX_T+1], g, h;
  int i, j, k;

  ASSERT_ALWAYS (t <= MAX_T); /* Ensures that all M[i][j] fit in uint64_t,
                                 and similarly for all intermediate
                                 computations on M[i][j]. This avoids the
                                 use of mpz_t to store the M[i][j]. */

  /* initialize M[i][j] = i^j */
  for (i = 0; i <= t; i++)
    for (j = 0; j <= t; j++)
      M[i][j] = (j == 0) ? 1 : i * M[i][j-1];

  /* forward Gauss: zero the under-diagonal coefficients while going down */
  for (i = 1; i <= t; i++)
    for (j = 0; j < i; j++)
      if (M[i][j] != 0)
      {
        g = gcd_uint64 (M[i][j], M[j][j]);
        h = M[i][j] / g;
        g = M[j][j] / g;
        /* f[i] <- g*f[i] - h*f[j] */
        mpz_mul_uint64 (f[i], f[i], g);
        mpz_submul_uint64 (f[i], f[j], h);
        for (k = j; k <= t; k++)
          M[i][k] = g * M[i][k] - h * M[j][k];
      }

  /* now zero upper-diagonal coefficients while going up */
  for (i = t; i >= 0; i--)
  {
    for (j = i + 1; j <= t; j++)
      /* f[i] = f[i] - M[i][j] * f[j] */
      mpz_submul_ui (f[i], f[j], M[i][j]);
    ASSERT (mpz_divisible_uint64_p (f[i], M[i][i]));
    mpz_divexact_uint64 (f[i], f[i], M[i][i]);
  }
}

/* Generic Toom-Cook implementation: stores in f[0..r+s] the coefficients
   of g*h, where g has degree r and h has degree s, and their coefficients
   are in g[0..r] and h[0..s].
   Assumes f differs from g and h, and f[0..r+s] are already allocated.
   Assumes r + s <= MAX_T (MAX_T = 13 ensures all the matrix coefficients
   in the inversion fit into uint64_t, and those used in mpz_mul_ui calls
   fit into uint32_t).
   Returns the degree of f.
*/
static int
mpz_poly_mul_tc (mpz_t *f, mpz_t *g, int r, mpz_t *h, int s)
{
  int t = r + s; /* product has t+1 coefficients */
  int i, j;
  mpz_t tmp;

  if ((r == -1) || (s == -1)) /* g or h is 0 */
    return -1;

  ASSERT_ALWAYS (t <= MAX_T);

  mpz_init (tmp);

  /* first store g(i)*h(i) in f[i] for 0 <= i <= t */
  for (i = 0; i <= t; i++)
  {
    /* f[i] <- g(i) */
    mpz_set (f[i], g[r]);
    for (j = r - 1; j >= 0; j--)
    {
      mpz_mul_ui (f[i], f[i], i);
      mpz_add (f[i], f[i], g[j]);
    }
    /* tmp <- h(i) */
    mpz_set (tmp, h[s]);
    for (j = s - 1; j >= 0; j--)
    {
      mpz_mul_ui (tmp, tmp, i);
      mpz_add (tmp, tmp, h[j]);
    }
    /* f[i] <- g(i)*h(i) */
    mpz_mul (f[i], f[i], tmp);
  }

  mpz_poly_mul_tc_interpolate (f, t);

  mpz_clear (tmp);
  return t;
}

/* Same as mpz_poly_mul_tc for the squaring: store in f[0..2r] the coefficients
   of g^2, where g has degree r, and their coefficients are in g[0..r].
   Assumes f differs from g, and f[0..2r] are already allocated.
   Assumes 2r <= MAX_T (MAX_T = 17 ensures all the matrix coefficients
   in the inversion fit into uint64_t).
   Returns the degree of f.
*/
static int
mpz_poly_sqr_tc (mpz_t *f, mpz_t *g, int r)
{
  int t = 2 * r; /* product has t+1 coefficients */
  int i, j;

  if (r == -1) /* g is 0 */
    return -1;

  ASSERT_ALWAYS (t <= MAX_T);

  /* first store g(i)^2 in f[i] for 0 <= i <= t */
  for (i = 0; i <= t; i++)
  {
    /* f[i] <- g(i) */
    mpz_set (f[i], g[r]);
    for (j = r - 1; j >= 0; j--)
    {
      mpz_mul_ui (f[i], f[i], i);
      mpz_add (f[i], f[i], g[j]);
    }
    mpz_mul (f[i], f[i], f[i]);
  }

  mpz_poly_mul_tc_interpolate (f, t);

  return t;
}

/* Pseudo-reduce a plain polynomial p modulo a non-monic polynomial F.
   The result is of type polymodF_t P, and satisfies:
   P->p = lc(F)^P->v * p mod F.
   WARNING: this function destroys its input p !!! */
static void
mpz_poly_reducemodF(polymodF_t P, mpz_poly_t p, const mpz_poly_t F) {
  int v = 0;

  if (p->deg < F->deg) {
    mpz_poly_copy(P->p, p);
    P->v = 0;
    return;
  }

  const int d = F->deg;

  while (p->deg >= d) {
    const int k = p->deg;
    int i;

    /* We compute F[d]*p - p[k]*F. In case F[d] divides p[k], we can simply
       compute p - p[k]/F[d]*F. However this will happen rarely with
       Kleinjung's polynomial selection, since lc(F) is large. */

    /* FIXME: in msieve, Jason Papadopoulos reduces by F[d]^d*F(x/F[d])
       instead of F(x). This might avoid one of the for-loops below. */

    // temporary hack: account for the possibility that we're indeed
    // using f_hat instead of f.
    if (mpz_cmp_ui(F->coeff[d], 1) != 0) {
      v++; /* we consider p/F[d]^v */
      for (i = 0; i < k; ++i)
        mpz_mul (p->coeff[i], p->coeff[i], F->coeff[d]);
    }

    for (i = 0; i < d; ++i)
      mpz_submul (p->coeff[k-d+i], p->coeff[k], F->coeff[i]);

    mpz_poly_cleandeg (p, k-1);
  }

  mpz_poly_copy(P->p, p);
  P->v = v;
}


/* ----------------
   Public functions
   ---------------- */

int
mpz_poly_normalized_p (const mpz_poly_t f) {
  return (f->deg == -1) || mpz_cmp_ui (f->coeff[f->deg], 0) != 0;
}

/* Allocate a polynomial that can contain 'd+1' coefficients and set to zero.
   We allow d=-1.
 */
void mpz_poly_init(mpz_poly_t f, int d) {
  int i;
  f->alloc = d+1;
  f->deg = -1;
  if (d == -1)
    f->coeff = (mpz_t *) NULL;
  else
    {
      f->coeff = (mpz_t *) malloc ((d+1)*sizeof(mpz_t));
      ASSERT (f->coeff != NULL);
    }
  for (i = 0; i <= d; ++i)
    mpz_init(f->coeff[i]);
}


/* realloc f to (at least) nc coefficients */
void
mpz_poly_realloc (mpz_poly_t f, int nc)
{
  int i;
  if (f->alloc < nc)
    {
      f->coeff = (mpz_t*) realloc (f->coeff, nc * sizeof (mpz_t));
      if (f->coeff == NULL)
        {
          fprintf (stderr, "Error, not enough memory\n");
          exit (1);
        }
      for (i = f->alloc; i < nc; i++)
        mpz_init (f->coeff[i]);
      f->alloc = nc;
    }
}

/* Free polynomial f in mpz_poly_t. */
void mpz_poly_clear(mpz_poly_t f) {
  int i;
  for (i = 0; i < f->alloc; ++i)
    mpz_clear(f->coeff[i]);
  free(f->coeff);
}

/* Print coefficients of f. */
void
mpz_poly_fprintf (FILE *fp, const mpz_poly_t f)
{
  int i;

  if (f->deg == -1)
    {
      fprintf (fp, "0\n");
      return;
    }
  else if (f->deg == 0)
    {
      gmp_fprintf (fp, "%Zd\n", f->coeff[0]);
      return;
    }
  gmp_fprintf (fp, "%Zd + ", f->coeff[0]);
  for (i = 1; i < f->deg; ++i)
    gmp_fprintf (fp, "%Zd*x^%d + ", f->coeff[i], i);
  gmp_fprintf (fp, "%Zd*x^%d\n", f->coeff[f->deg], f->deg);
}

/* Set polynomial degree */
void mpz_poly_set_deg(mpz_poly_t f, int deg)
{
  f->deg = deg;
}


/* Find polynomial degree. */
void mpz_poly_cleandeg(mpz_poly_t f, int deg) {
  ASSERT(deg >= -1);
  while ((deg >= 0) && (mpz_cmp_ui(f->coeff[deg], 0)==0))
    deg--;
  f->deg = deg;
}

/* Sets f to the polynomial of degree d, of coefficients
   given by coeffs. */
void
mpz_poly_set (mpz_poly_t f, mpz_t * coeffs, int d)
{
  int i;
  for (i=d; i>=0; --i)
    mpz_poly_setcoeff(f,i,coeffs[i]);
  mpz_poly_cleandeg(f,d);
}

/* Set a zero polynomial. */
void mpz_poly_set_zero(mpz_poly_t f) {
  f->deg = -1;
}

/* Set mpz_t coefficient for the i-th term. */
void mpz_poly_setcoeff(mpz_poly_t f, int i, const mpz_t z) {
  int j;
  if (i >= f->alloc) {
    f->coeff = (mpz_t *)realloc(f->coeff, (i+1)*sizeof(mpz_t));
    ASSERT (f->coeff != NULL);
    for (j = f->alloc; j <= i; ++j)
      mpz_init(f->coeff[j]);
    f->alloc = i+1;
  }
  mpz_set(f->coeff[i], z);
  if (i >= f->deg)
    mpz_poly_cleandeg(f, i);
}

/* Set int64 coefficient for the i-th term. */
void mpz_poly_setcoeff_int64(mpz_poly_t f, int i, int64_t z) {
  int j;
  if (i >= f->alloc) {
    f->coeff = (mpz_t *)realloc(f->coeff, (i+1)*sizeof(mpz_t));
    ASSERT (f->coeff != NULL);
    for (j = f->alloc; j <= i; ++j)
      mpz_init(f->coeff[j]);
    f->alloc = i+1;
  }
  mpz_set_int64(f->coeff[i], z);
  if (i >= f->deg)
    mpz_poly_cleandeg(f, i);
}


/* Set signed int coefficient for the i-th term. */
void mpz_poly_setcoeff_si(mpz_poly_t f, int i, int z) {
  int j;
  if (i >= f->alloc) {
    f->coeff = (mpz_t *)realloc(f->coeff, (i+1)*sizeof(mpz_t));
    ASSERT (f->coeff != NULL);
    for (j = f->alloc; j <= i; ++j)
      mpz_init(f->coeff[j]);
    f->alloc = i+1;
  }
  mpz_set_si(f->coeff[i], z);
  if (i >= f->deg)
    mpz_poly_cleandeg(f, i);
}

/* Get coefficient for the i-th term. */
void mpz_poly_getcoeff(mpz_t res, int i, const mpz_poly_t f) {
  // The code below will work anyway,
  // this assert is better called a warning.
  ASSERT_ALWAYS( f->deg == -1 ||  f->deg>=i );
  if(i > f->deg)
    mpz_set_ui(res,0);
  else
    mpz_set(res,f->coeff[i]);
}

/* Copy f to g, where g must be initialized (but there is no need it has
   enough allocated coefficients, since mpz_poly_setcoeff reallocates if
   needed). */
void mpz_poly_copy(mpz_poly_t g, const mpz_poly_t f) {
  int i;
  g->deg = f->deg;
  for (i = f->deg; i >= 0; --i)
    mpz_poly_setcoeff (g, i, f->coeff[i]);
}


/* -------------------------------------------------------------------------- */
/* return 0 if f and g are equal, non-zero otherwise */
int mpz_poly_cmp (mpz_poly_t f, mpz_poly_t g)
{
  int i;

  if (f->deg != g->deg)
    return 1;

  if (f->deg == -1)
    return 0;

  for (i = 0; i <= f->deg; i++)
    if (mpz_cmp (f->coeff[i], g->coeff[i]))
      return 1; /* f and g differ */
  return 0;
}
/* -------------------------------------------------------------------------- */


/* Set f=g+h.
   Note: f can be the same as g or h;
         g can be the same as h. */
void mpz_poly_add(mpz_poly_t f, const mpz_poly_t g, const mpz_poly_t h) {

  if(f==h || f==g) {
    mpz_poly_t aux;
    mpz_poly_init(aux,-1);
    mpz_poly_add(aux,g,h);
    mpz_poly_copy(f,aux);
    mpz_poly_clear(aux);
    return;
  }
  int i, maxdeg;
  mpz_t z;
  mpz_init(z);

  maxdeg = max(g->deg, h->deg);
  f->deg = maxdeg;
  for (i = maxdeg; i >= 0; --i) {
    if (i <= g->deg)
      mpz_set(z, g->coeff[i]);
    else
      mpz_set_ui(z, 0);
    if (i <= h->deg) {
      //ASSERT_ALWAYS(h->alloc >= h->deg+1);
      mpz_add(z, z, h->coeff[i]);
    }
    mpz_poly_setcoeff(f, i, z);
  }
  mpz_clear(z);
  mpz_poly_cleandeg(f, maxdeg);
}

/* Set f=g-h.
   Note: f can be the same as g */
void mpz_poly_sub(mpz_poly_t f, const mpz_poly_t g, const mpz_poly_t h) {
  int i, maxdeg;
  mpz_t z;
  mpz_init(z);
  maxdeg = max(g->deg, h->deg);
  f->deg = maxdeg;
  for (i = maxdeg; i >= 0; --i) {
    if (i <= g->deg)
      mpz_set(z, g->coeff[i]);
    else
      mpz_set_ui(z, 0);
    if (i <= h->deg)
      mpz_sub(z, z, h->coeff[i]);
    mpz_poly_setcoeff(f, i, z);
  }
  mpz_clear(z);
  mpz_poly_cleandeg(f, maxdeg);
}

/* Set f=f-a where a is unsigned long. */
void
mpz_poly_sub_ui (mpz_poly_t f, unsigned long a)
{
  mpz_t aux;
  mpz_init(aux);
  mpz_sub_ui(aux, f->coeff[0], a);
  mpz_poly_setcoeff(f, 0, aux);
  mpz_clear(aux);
}

/* Set f=g-h (mod m). */
void
mpz_poly_sub_mod_mpz (mpz_poly_t f, const mpz_poly_t g, const mpz_poly_t h, const mpz_t m)
{
  int i, maxdeg;
  int gdeg = g->deg; /* save in case f=g */
  int hdeg = h->deg; /* save in case f=h */
  mpz_t z;
  mpz_init(z);
  maxdeg = max(gdeg, hdeg);
  f->deg = maxdeg;
  for (i = maxdeg; i >= 0; --i) {
    if (i <= gdeg)
      mpz_set(z, g->coeff[i]);
    else
      mpz_set_ui(z, 0);
    if (i <= hdeg)
      mpz_sub(z, z, h->coeff[i]);
    mpz_poly_setcoeff(f, i, z);
  }
  mpz_clear(z);
  mpz_poly_cleandeg(f, maxdeg);
  for (i = 0; i <= f->deg; ++i) {
    mpz_mod(f->coeff[i], f->coeff[i], m);
  }
  mpz_poly_cleandeg(f, f->deg);
}

/* Set f=g*h. Note: f might equal g or h. */
void
mpz_poly_mul (mpz_poly_t f, const mpz_poly_t g, const mpz_poly_t h) {
  int i, maxdeg;
  mpz_poly_t prd;

  if (f == h || f == g)
    {
      mpz_poly_t aux;
      mpz_poly_init (aux, -1);
      mpz_poly_mul (aux, g, h);
      mpz_poly_copy (f,aux);
      mpz_poly_clear (aux);
      return;
    }

  if ((g->deg == -1) || (h->deg == -1)) {
    f->deg = -1;
    return;
  }

  if (g->deg == 0) {
    mpz_poly_mul_mpz (f, h, g->coeff[0]);
    return;
  }

  if (h->deg == 0) {
    mpz_poly_mul_mpz (f, g, h->coeff[0]);
    return;
  }

  ASSERT_ALWAYS(mpz_cmp_ui (g->coeff[g->deg], 0) != 0);
  ASSERT_ALWAYS(mpz_cmp_ui (h->coeff[h->deg], 0) != 0);
  ASSERT_ALWAYS(f != g);
  ASSERT_ALWAYS(f != h);

  maxdeg = g->deg + h->deg;
  mpz_poly_init(prd, maxdeg);

#if 1
  if (g->deg + h->deg <= MAX_T) {
    prd->deg = mpz_poly_mul_tc (prd->coeff, g->coeff, g->deg, h->coeff, h->deg);
  } else {
    /* naive product */
    /* currently we have to resort to this for larger degree, because
     * the generic toom implementation is bounded in degree.
     */
    prd->deg = mpz_poly_mul_basecase (prd->coeff, g->coeff, g->deg, h->coeff, h->deg);
  }
#else /* segmentation, this code has problem with huge runs, for example
         degree 5 with lifting to 631516975 bits */
  {
    mpz_t G, H;
    size_t sg, sh, s;
    mpz_init (G);
    mpz_init (H);
    for (sg = 0, i = 0; i <= g->deg; i++)
    {
      s = mpz_sizeinbase (g->coeff[i], 2);
      if (s > sg)
        sg = s;
    }
    for (sh = 0, i = 0; i <= h->deg; i++)
    {
      s = mpz_sizeinbase (h->coeff[i], 2);
      if (s > sh)
        sh = s;
    }
    /* the +1 accounts for a possible sign */
    for (s = sg + sh + 1, i = h->deg; i > 1; i = (i + 1) / 2, s++);
    mpz_set (G, g->coeff[g->deg]);
    for (i = g->deg - 1; i >= 0; i--)
    {
      mpz_mul_2exp (G, G, s);
      mpz_add (G, G, g->coeff[i]);
    }
    /* sanity check: G should have sizeinbase(lc(g))+d*s bits (or -1) */
    size_t size_g = mpz_sizeinbase (g->coeff[g->deg], 2) + g->deg * s;
    ASSERT_ALWAYS(mpz_sizeinbase (G, 2) == size_g ||
                  mpz_sizeinbase (G, 2) == size_g - 1);
    mpz_set (H, h->coeff[h->deg]);
    for (i = h->deg - 1; i >= 0; i--)
    {
      mpz_mul_2exp (H, H, s);
      mpz_add (H, H, h->coeff[i]);
    }
    /* sanity check: H should have sizeinbase(lc(h))+d*s bits (or -1) */
    size_t size_h = mpz_sizeinbase (h->coeff[h->deg], 2) + h->deg * s;
    ASSERT_ALWAYS(mpz_sizeinbase (H, 2) == size_h ||
                  mpz_sizeinbase (H, 2) == size_h - 1);
    size_g = mpz_sizeinbase (G, 2);
    size_h = mpz_sizeinbase (H, 2);
    /* sanity check: we verify that the product agrees both mod B and B-1 */
    mp_limb_t g0 = mpz_getlimbn (G, 0), h0 = mpz_getlimbn (H, 0);
    mp_limb_t g1 = mod_base_minus_1 (G), h1 = mod_base_minus_1 (H);
    mpz_mul (G, G, H);
    ASSERT_ALWAYS (mpz_getlimbn (G, 0) == g0 * h0);
    mpz_set_ui (H, g1);
    mpz_mul_ui (H, H, h1);
    ASSERT_ALWAYS (mod_base_minus_1 (G) == mod_base_minus_1 (H));
    ASSERT_ALWAYS(mpz_sizeinbase (G, 2) == size_g + size_h ||
                  mpz_sizeinbase (G, 2) == size_g + size_h - 1);
    for (i = 0; i < g->deg + h->deg; i++)
    {
      mpz_fdiv_r_2exp (prd->coeff[i], G, s);
      if (mpz_sizeinbase (prd->coeff[i], 2) == s)
      {
        mpz_cdiv_r_2exp (prd->coeff[i], G, s);
        mpz_cdiv_q_2exp (G, G, s);
      }
      else
        mpz_fdiv_q_2exp (G, G, s);
      ASSERT_ALWAYS(mpz_sizeinbase (prd->coeff[i], 2) < s);
    }
    mpz_set (prd->coeff[i], G);
    mpz_clear (G);
    mpz_clear (H);
  }
#endif

  for (i = maxdeg; i >= 0; --i)
    mpz_poly_setcoeff(f, i, prd->coeff[i]);
  mpz_poly_cleandeg(f, maxdeg);
  ASSERT_ALWAYS(mpz_cmp_ui (f->coeff[f->deg], 0) != 0);

  mpz_poly_clear(prd);
}

/* Set Q=a*P, where a is an mpz_t */
void
mpz_poly_mul_mpz (mpz_poly_t Q, const mpz_poly_t P, const mpz_t a)
{
  int i;
  mpz_t aux;

  mpz_init (aux);
  Q->deg = P->deg;
  for (i = 0; i <= P->deg; ++i)
    {
      mpz_mul (aux, P->coeff[i], a);
      mpz_poly_setcoeff (Q, i, aux);
    }
  mpz_clear (aux);
}

/* h=rem(h, f) mod p, f not necessarily monic. */
static void
mpz_poly_div_r (mpz_poly_t h, const mpz_poly_t f, const mpz_t p)
{
  int i, d = f->deg, dh = h->deg;
  mpz_t tmp, inv;

  mpz_init_set_ui (inv, 1);
  mpz_init_set_ui (tmp, 1);

  mpz_set (tmp, f->coeff[d]);
  if (mpz_cmp_ui(tmp, 1) != 0)
    /* inv is 1/f[d] mod p */
    mpz_invert (inv, tmp, p);

  while (dh >= d)
  {
    /* subtract h[dh]/f[d]*x^(dh-d)*f from h */
    if (mpz_cmp_ui(inv, 1) != 0) {
      mpz_mul (h->coeff[dh], h->coeff[dh], inv);
      mpz_fdiv_r (h->coeff[dh], h->coeff[dh], p);
    }

    for (i = 0; i < d; i++) {
      mpz_mul (tmp, h->coeff[dh], f->coeff[i]);
      mpz_fdiv_r (tmp, tmp, p);
      mpz_sub (h->coeff[dh - d + i], h->coeff[dh - d + i], tmp);
      mpz_fdiv_r (h->coeff[dh - d + i], h->coeff[dh - d + i], p);
    }

    do {
      dh --;
    }
    while (dh >= 0 && mpz_divisible_p (h->coeff[dh],p));

    h->deg = dh;
  }

  mpz_clear (inv);
  mpz_clear (tmp);
}

/* 
   computes q, r such that f = q*g + r mod p, with deg(r) < deg(g) and p in mpz_t
   q and r must be allocated!                                   
*/
void mpz_poly_div_qr (mpz_poly_t q, mpz_poly_t r, const mpz_poly_t f, const mpz_poly_t g, const mpz_t p)
{
  int k, j, df = f->deg, dg = g->deg, dq = df - dg;
  mpz_t tmp, lg, invlg;

  ASSERT_ALWAYS(df >= dg);
  ASSERT_ALWAYS(q->alloc >= dq+1);

  mpz_init(tmp);
  mpz_init(lg);
  mpz_init_set_ui(invlg, 1);

  mpz_poly_copy(r, f);
  q->deg = dq;

  mpz_set (lg, g->coeff[dg]);
  mpz_mod (lg, lg, p);
  /* invlg = 1/g[dg] mod p */
  if (mpz_cmp_ui(lg, 1) != 0)
    mpz_invert(invlg, lg, p);

  for (k = df-dg ; k >=0 ; k--) {
    mpz_mul(q->coeff[k], r->coeff[k+dg], invlg);
    mpz_mod(q->coeff[k], q->coeff[k], p);
    for (j = dg+k ; j >= k ; j--) {
      mpz_mul(tmp, q->coeff[k], g->coeff[j-k]);
      mpz_sub(r->coeff[j], r->coeff[j], tmp);
      mpz_mod(r->coeff[j], r->coeff[j], p);
    }
  }
  mpz_poly_cleandeg(r, r->deg);

  mpz_clear(invlg);
  mpz_clear(lg);
  mpz_clear(tmp);
}

/* q=divexact(h, f) mod p, f not necessarily monic. Clobbers h. */
static void
mpz_poly_divexact (mpz_poly_t q, mpz_poly_t h, const mpz_poly_t f, const mpz_t p) {
  int i, d = f->deg, dh = h->deg;
  mpz_t t, aux;

  mpz_init (t);
  mpz_init (aux);
  ASSERT (d >= 0);
  ASSERT (dh >= 0);
  ASSERT (dh >= d);

  mpz_poly_realloc (q, dh + 1 - d);
  q->deg = dh - d;
  /* t is 1/f[d] mod p */
  mpz_set (aux, f->coeff[d]);
  mpz_fdiv_r (aux, aux, p);
  if (mpz_cmp_ui(aux, 1) != 0)
    mpz_invert (t, aux, p);

  while (dh >= d) {

    /* subtract h[dh]/f[d]*x^(dh-d)*f from h */
    if (mpz_cmp_ui(t, 1) != 0) {
      mpz_mul (h->coeff[dh], h->coeff[dh], t);
      mpz_fdiv_r (h->coeff[dh], h->coeff[dh], p);
    }
    mpz_set (q->coeff[dh-d], h->coeff[dh]);
    mpz_fdiv_r (q->coeff[dh-d], q->coeff[dh - d], p);

    /* we only need to update the coefficients of degree >= d of h,
       i.e., we want i >= 2d - dh */
    for (i = (2 * d > dh) ? 2 * d - dh : 0; i < d; i++) {
      mpz_mul (aux, h->coeff[dh], f->coeff[i]);
      mpz_sub (h->coeff[dh - d + i], h->coeff[dh - d + i], aux);
      mpz_fdiv_r (h->coeff[dh - d + i], h->coeff[dh - d + i], p);
    }
    dh --;
  }
  mpz_poly_normalize_modp (q, p);

  mpz_clear (t);
  mpz_clear (aux);
}

/* Set f=g/2 (mod m). */
void mpz_poly_div_2_mod_mpz (mpz_poly_t f, const mpz_poly_t g, const mpz_t m)
{
  int i;
  mpz_t aux;

  mpz_init (aux);
  for (i = g->deg; i >= 0; --i)
  {
    if (mpz_scan1 (g->coeff[i], 0) == 0)
    {
      mpz_add (aux, g->coeff[i], m);
      mpz_div_2exp (aux, aux, 1);
    }
    else
      mpz_div_2exp (aux, g->coeff[i], 1);
    mpz_mod(aux, aux, m);
    mpz_poly_setcoeff(f, i, aux);
  }
  mpz_clear (aux);
}

/* Set res=f(x) */
void mpz_poly_eval(mpz_t res, const mpz_poly_t f, const mpz_t x) {
  int i, d;
  d = f->deg;
  if (d == -1) {
    mpz_set_ui(res, 0);
    return;
  }
  mpz_set(res, f->coeff[d]);
  for (i = d-1; i>=0; --i) {
    mpz_mul(res, res, x);
    mpz_add(res, res, f->coeff[i]);
  }
}

/* Set res=f(x) (mod m) */
void mpz_poly_eval_mod_mpz(mpz_t res, const mpz_poly_t f, const mpz_t x,
                       const mpz_t m)
{
    mpz_poly_eval_mod_mpz_barrett(res, f, x, m, NULL);
}

void mpz_poly_eval_mod_mpz_barrett(mpz_t res, const mpz_poly_t f, const mpz_t x,
                       const mpz_t m, const mpz_t mx) {
  int i, d;
  d = f->deg;
  if (d == -1) {
    mpz_set_ui(res, 0);
    return;
  }
  mpz_mod(res, f->coeff[d], m);
  for (i = d-1; i>=0; --i) {
    mpz_mul(res, res, x);
    mpz_add(res, res, f->coeff[i]);
    barrett_mod(res, res, m, mx);
  }
}

/* This evaluates several polynomials at the same point w. It is possible
 * to do fewer modular reductions in this case.
 *
 * When k==1, we use mpz_poly_eval_mod_mpz instead, since it's faster
 * then.
 */
void mpz_poly_eval_several_mod_mpz_barrett(mpz_ptr * r, mpz_poly_srcptr * f, int k, const mpz_t x,
                       const mpz_t m, const mpz_t mx)
{
    int i;

    if (k == 1) {
        mpz_poly_eval_mod_mpz_barrett(r[0], f[0], x, m, mx);
        return;
    }

    mpz_t w;
    mpz_init(w);
    int maxdeg = -1;
    for(int j = 0 ; j < k ; j++) {
        if (f[j]->deg >= 0)
            mpz_set(r[j],f[j]->coeff[0]);
        else
            mpz_set_ui(r[j], 0);
        if (f[j]->deg > maxdeg)
            maxdeg = f[j]->deg;
    }

    mpz_set(w, x);
    for(int j = 0 ; j < k ; j++) {
        if (f[j]->deg >= 1)
            mpz_addmul(r[j], w, f[j]->coeff[1]);
    }
    for(i = 2 ; i <= maxdeg ; i++) {
        mpz_mul(w, w, x);
        barrett_mod(w, w, m, mx);
        for(int j = 0 ; j < k ; j++)
            if (f[j]->deg >= i)
                mpz_addmul(r[j], w, f[j]->coeff[i]);
    }
    for(int j = 0 ; j < k ; j++) {
        barrett_mod(r[j], r[j], m, mx);
    }
    mpz_clear(w);
}

/* Set Q=P1*P2 (mod F). Warning: Q might equal P1 (or P2). */
void polymodF_mul (polymodF_t Q, const polymodF_t P1, const polymodF_t P2,
                   const mpz_poly_t F) {
  mpz_poly_t prd;
  int v;

  mpz_poly_init (prd, P1->p->deg + P2->p->deg);

  ASSERT_ALWAYS(mpz_poly_normalized_p (P1->p));
  ASSERT_ALWAYS(mpz_poly_normalized_p (P2->p));

  mpz_poly_mul (prd, P1->p, P2->p);
  v = P1->v + P2->v;

  mpz_poly_reducemodF (Q, prd, F);
  Q->v += v;
  mpz_poly_clear (prd);
}

/* Set Q = P (mod m) */
void mpz_poly_reduce_mod_mpz (mpz_poly_t Q, const mpz_poly_t P, const mpz_t m)
{
  int i;
  mpz_t aux;

  mpz_init(aux);
  Q->deg = P->deg;
  for (i = 0; i <= P->deg; ++i) {
    mpz_mod(aux, P->coeff[i], m);
    mpz_poly_setcoeff(Q, i, aux);
  }
  mpz_clear(aux);
  mpz_poly_cleandeg(Q,Q->deg);
}

/* Set Q = P/lc(P) (mod m). Q and P might be identical. */
static void
mpz_poly_reduce_makemonic_mod_mpz (mpz_poly_t Q, const mpz_poly_t P, const mpz_t m)
{
  int i;
  mpz_t aux, aux2;
  mpz_init(aux);
  mpz_init(aux2);
  i = P->deg;
  do {
    mpz_mod(aux, P->coeff[i], m);
    i--;
  } while ((i>=0) && (mpz_cmp_ui(aux, 0) == 0));
  /* either i+1 >= 0 and P[i+1] <> 0 (mod m),
     or i = -1 and P is identically zero modulo m */
  if (i>=0) {
    Q->deg = i+1;
    mpz_invert(aux2, aux, m);
    for (i = 0; i < Q->deg; ++i) {
      mpz_mul(aux, aux2, P->coeff[i]);
      mpz_mod(aux, aux, m);
      mpz_poly_setcoeff(Q, i, aux);
    }
    /* we can directly set the leading coefficient to 1 */
    mpz_set_ui (Q->coeff[Q->deg], 1);
  } else {
    if (mpz_cmp_ui(aux, 0) == 0)
      Q->deg = -1;
    else {
      Q->deg = 0;
      mpz_set_ui(aux, 1);
      mpz_poly_setcoeff(Q, 0, aux);
    }
  }
  mpz_clear(aux);
  mpz_clear(aux2);
}

/* Reduce R[d]*x^d + ... + R[0] mod f[df]*x^df + ... + f[0] modulo m.
   Return the degree of the remainder. */
int
mpz_poly_mod_f_mod_mpz (mpz_t *R, int d, mpz_t *f, int df, const mpz_t m,
                    const mpz_t invm)
{
  int i;
  mpz_t aux, c;

  mpz_init (aux);
  mpz_init (c);
  mpz_invert (aux, m, f[df]); /* aux = 1/m mod lc(f). We could precompute it
                                 but it should not cost much. */
  // FIXME: write a subquadratic variant
  while (d >= df)
  {
    /* First reduce the leading coefficient of R mod m. However this is
       quite expensive, since it costs O(D(n)) where m has size n, whereas
       if we leave R[d] to size 2n, we have an overhead of O(n) only. */
    // barrett_mod (R[d], R[d], m, invm);
    /* Here m is large (typically several million bits) and lc(f) is small
       (typically one word). We first add to R lambda * m * x^(d-df)
       --- which is zero mod m ---
       such that the new coefficient of degree d
       is divisible by lc(f), i.e., lambda = -R[d]/m mod lc(f).
       Then if c = (R[d] + lambda * m) / lc(f), we simply subtract
       c * x^(d-df) * f.
    */
    mpz_mod (c, R[d], f[df]); /* R[d] mod lc(f) */
    mpz_mul (c, c, aux);
    mpz_mod (c, c, f[df]);    /* R[d]/m mod lc(f) */
    mpz_submul (R[d], m, c);  /* R[d] - m * (R[d] / m mod lc(f)) */
    ASSERT (mpz_divisible_p (R[d], f[df]));
    mpz_divexact (c, R[d], f[df]);
    for (i = d - 1; i >= d - df; --i)
      mpz_submul (R[i], c, f[df-d+i]);
    d--;
  }

  /* reduce lower coefficients */
  for (i = 0; i <= d; ++i)
    barrett_mod (R[i], R[i], m, invm);
  mpz_clear (aux);
  mpz_clear (c);

  return d;
}

/*  Reduce frac (= num / denom) mod F mod m ,
    i.e. compute num * denom^-1 mod F mod m .
    The return value is in num, denom is set to constant polynomial 1
    */
void
mpz_poly_reduce_frac_mod_f_mod_mpz (mpz_poly_t num, mpz_poly_t denom,
                                    const mpz_poly_t F, const mpz_t m)
{
  if (denom->deg == 0)
  {
    mpz_t inv;
    mpz_init (inv);
    mpz_poly_getcoeff (inv, 0, denom); /* inv <- denom[0] */
    mpz_invert (inv, inv, m);          /* inv <- denom[0]^-1 */
    mpz_poly_mul_mpz (num, num, inv);  /* num <- num * inv */
    mpz_poly_reduce_mod_mpz (num, num, m); /* num <- num * inv mod m */
    mpz_clear (inv);
  }
  else
  {
    int d;
    mpz_poly_t g, U, V;
    mpz_poly_init (g, 0);
    mpz_poly_init (U, 0);
    mpz_poly_init (V, 0);
    mpz_poly_xgcd_mpz (g, F, denom, U, V, m);
    mpz_poly_mul (num, num, V);
    d=mpz_poly_mod_f_mod_mpz (num->coeff, num->deg, F->coeff, F->deg, m, NULL);
    mpz_poly_cleandeg(num, d);
    mpz_poly_clear (g);
    mpz_poly_clear (U);
    mpz_poly_clear (V);
  }
  mpz_poly_set_zero (denom);
  mpz_poly_setcoeff_si (denom, 0, 1);
}


// Q = P1*P2 mod f, mod m
// f is the original algebraic polynomial (non monic but small coefficients)
void
mpz_poly_mul_mod_f_mod_mpz (mpz_poly_t Q, const mpz_poly_t P1, const mpz_poly_t P2,
                        const mpz_poly_t f, const mpz_t m, const mpz_t invm)
{
  int d1 = P1->deg;
  int d2 = P2->deg;
  int d = d1+d2;
  int df = f->deg;
  mpz_poly_t R;

  mpz_poly_init(R, d);

  mpz_poly_mul_tc (R->coeff, P1->coeff, d1, P2->coeff, d2);

  // reduce mod f
  d = mpz_poly_mod_f_mod_mpz (R->coeff, d, f->coeff, df, m, invm);

  mpz_poly_cleandeg(R, d);
  mpz_poly_copy(Q, R);
  mpz_poly_clear(R);
}

// Q = P^2 mod f, mod m
// f is the original algebraic polynomial (non monic but small coefficients)
void
mpz_poly_sqr_mod_f_mod_mpz (mpz_poly_t Q, const mpz_poly_t P, const mpz_poly_t f,
                        const mpz_t m, const mpz_t invm)
{
  int d1 = P->deg;
  int d = d1 + d1;
  int df = f->deg;
  mpz_poly_t R;

  mpz_poly_init(R, d);

  /* Fast squaring in 2d1+1 squares, i.e., 2d-1 squares.
     For d=5, this gives 9 squares. */
  mpz_poly_sqr_tc (R->coeff, P->coeff, d1);

  // reduce mod f
  d = mpz_poly_mod_f_mod_mpz (R->coeff, d, f->coeff, df, m, invm);

  mpz_poly_cleandeg(R, d);
  mpz_poly_copy(Q, R);
  mpz_poly_clear(R);
}

/* Affects the derivative of f to df. Assumes df different from f.
   Assumes df has been initialized with degree at least f->deg-1. */
void mpz_poly_derivative(mpz_poly_t df, const mpz_poly_t f) {
  int n;
  if (f->deg <= 0)
    {
      df->deg = -1;
      return;
    }
  // at this point, f->deg >=1
  df->deg = f->deg-1;
  for(n=0; n<=f->deg-1; n++)
    mpz_mul_si(df->coeff[n],f->coeff[n+1],n+1);
}

/* Q = P^a mod f, mod p (f is the algebraic polynomial, non monic) */
void
mpz_poly_power_mod_f_mod_ui (mpz_poly_t Q, const mpz_poly_t P, const mpz_poly_t f,
                         const mpz_t a, unsigned long p)
{
  mpz_t m;
  int k = mpz_sizeinbase(a, 2);
  mpz_poly_t R;

  if (mpz_cmp_ui(a, 0) == 0) {
    Q->deg = 0;
    mpz_set_ui(Q->coeff[0], 1);
    return;
  }

  mpz_init_set_ui(m, p);
  mpz_poly_init(R, 2*f->deg);

  // Initialize R to P
  mpz_poly_copy(R, P);

  // Horner
  for (k -= 2; k >= 0; k--)
  {
    mpz_poly_sqr_mod_f_mod_mpz(R, R, f, m, NULL);  // R <- R^2
    if (mpz_tstbit(a, k))
      mpz_poly_mul_mod_f_mod_mpz(R, R, P, f, m, NULL);  // R <- R*P
  }

  mpz_poly_copy(Q, R);
  mpz_clear (m);
  mpz_poly_clear(R);
}

/* Q = P^a mod f, mod p. Note, p is mpz_t */
void
mpz_poly_power_mod_f_mod_mpz (mpz_poly_t Q, const mpz_poly_t P, const mpz_poly_t f,
                          const mpz_t a, const mpz_t p)
{
  int k = mpz_sizeinbase(a, 2);
  mpz_poly_t R;

  if (mpz_cmp_ui(a, 0) == 0) {
    Q->deg = 0;
    mpz_set_ui(Q->coeff[0], 1);
    return;
  }

  mpz_poly_init(R, 2*f->deg);

  // Initialize R to P
  mpz_poly_copy(R, P);

  // Horner
  for (k -= 2; k >= 0; k--)
  {
    mpz_poly_sqr_mod_f_mod_mpz(R, R, f, p, NULL);  // R <- R^2
    if (mpz_tstbit(a, k))
      mpz_poly_mul_mod_f_mod_mpz(R, R, P, f, p, NULL);  // R <- R*P
  }

  mpz_poly_copy(Q, R);
  mpz_poly_clear(R);
}

/* store in invm the value of floor(B^(2k)/m), where m has k limbs,
   and B is the limb base */
void barrett_init (mpz_ptr invm, mpz_srcptr m)
{
  size_t k = mpz_size (m);

  mpz_set_ui (invm, 1);
  mpz_mul_2exp (invm, invm, 2 * k * GMP_NUMB_BITS);
  mpz_tdiv_q (invm, invm, m);
}

/* a <- b mod m */
void barrett_mod (mpz_ptr a, mpz_srcptr b, mpz_srcptr m, mpz_srcptr invm)
{
  size_t k = mpz_size (m), sizeb, l;
  mpz_t c;

  sizeb = mpz_size (b);
  if (invm == NULL || sizeb <= k || k == 1)
  {
    mpz_mod (a, b, m);
    return;
  }

  /* now sizeb > k */
  l = sizeb - k;
  if (l > k)
    l = k;
  /* l > 0 is the number of limbs to reduce */
  mpz_init (c);
  mpz_div_2exp (c, b, (sizeb - l) * GMP_NUMB_BITS); /* c has l limbs */
  mpz_mul (c, c, invm);                             /* c has k+l limbs */
  mpz_div_2exp (c, c, k * GMP_NUMB_BITS);           /* c has l limbs */
  mpz_mul (c, c, m);                                /* c has k+l limbs */
  mpz_mul_2exp (c, c, (sizeb - k - l) * GMP_NUMB_BITS);
  mpz_sub (a, b, c);
  while (mpz_size (a) > k)
  {
    sizeb = mpz_size (a);
    l = sizeb - k;
    if (l > k)
      l = k;
    mpz_div_2exp (c, a, (sizeb - l) * GMP_NUMB_BITS);
    mpz_mul (c, c, invm);
    mpz_div_2exp (c, c, k * GMP_NUMB_BITS);
    mpz_mul (c, c, m);
    mpz_mul_2exp (c, c, (sizeb - k - l) * GMP_NUMB_BITS);
    mpz_sub (a, a, c);
  }
  mpz_clear (c);
  mpz_mod (a, a, m);
}

/* Q = P^a mod f, mod p. Note, p is mpz_t */
/* Same as mpz_poly_power_mod_f_mod_mpz but use barrett for reduction mod p */
void
mpz_poly_power_mod_f_mod_mpz_Barrett (mpz_poly_t Q, const mpz_poly_t P,
                                      const mpz_poly_t f, const mpz_t a,
                                      const mpz_t p, const mpz_t invp)
{
  int k = mpz_sizeinbase(a, 2);
  mpz_poly_t R;

  if (mpz_cmp_ui(a, 0) == 0) {
    Q->deg = 0;
    mpz_set_ui(Q->coeff[0], 1);
    return;
  }

  mpz_poly_init(R, 2*f->deg);

  // Initialize R to P
  mpz_poly_copy(R, P);

  // Horner
  for (k -= 2; k >= 0; k--)
  {
    mpz_poly_sqr_mod_f_mod_mpz(R, R, f, p, invp);  // R <- R^2
    if (mpz_tstbit(a, k))
      mpz_poly_mul_mod_f_mod_mpz(R, R, P, f, p, invp);  // R <- R*P
  }

  mpz_poly_copy(Q, R);
  mpz_poly_clear(R);
}


/* Return a list of polynomials P[0], P[1], ..., P[l] such that
   P0 = Q[l-1] + p^K[1]*P[l]
   Q[l-1] = Q[l-2] + p^K[2]*P[l-1]
   ...
   Q[2] = Q[1] + p^K[l-1]*P[2]
   Q[1] = Q[0] + p*P[1]         where K[l] = 1
   Q[0] = P[0]
   ...
   With all coefficients of P[i] smaller than p^(K[l-i]-K[l-(i-1)]).
   Assume K[l]=1, K[l-1]=2, ..., K[i] = 2*K[i+1] or 2*K[i+1]-2.
   P0 = P[0] + p*P[1] + p^2*P[2] + ... + p^K[1]*P[l] < p^K[0]
   The end of the list is P[l+1]=0.
*/
mpz_poly_t*
mpz_poly_base_modp_init (const mpz_poly_t P0, int p, int *K, int l)
{
  mpz_poly_t *P;
  int k, i, j;
  mpz_t *pk;

  ASSERT_ALWAYS (K[l] == 1);

  /* initialize pk[i] = p^K[l-i] for 0 <= i < l */
  pk = (mpz_t*) malloc (l * sizeof (mpz_t));
  mpz_init_set_ui (pk[0], p);
  for (i = 1; i < l; i++)
  {
    mpz_init (pk[i]);
    mpz_mul (pk[i], pk[i-1], pk[i-1]);
    if (K[l-i] & 1)
    {
      ASSERT_ALWAYS (K[l-i] == 2 * K[l-i+1] - 1);
      mpz_div_ui (pk[i], pk[i], p);
    }
    else
      ASSERT_ALWAYS (K[l-i] == 2 * K[l-i+1]);
  }

  /* now decompose P0: we need P[0], P[1] for factor p, P[2] for p^2,
     ..., P[l] for p^K[1], and one for the end of list,
     thus l+2 polynomials */
  P = (mpz_poly_t*) malloc ((l + 2) * sizeof(mpz_poly_t));
  for (i = 0; i < l + 2; i++)
    mpz_poly_init (P[i], P0->deg);
  /* P[l+1] is initialized to 0 by mpz_poly_init */

  /* initialize P[k], and put remainder in P[k-1] */
  for (k = l, i = 0; i <= P0->deg; i++)
    mpz_tdiv_qr (P[k]->coeff[i], P[k-1]->coeff[i], P0->coeff[i], pk[k-1]);
  mpz_poly_cleandeg (P[k], P0->deg);

  /* now go down */
  for (j = l-1; j >= 1; j--)
  {
    /* reduce P[j] into P[j] and P[j-1] */
    for (i = 0; i <= P0->deg; i++)
      mpz_tdiv_qr (P[j]->coeff[i], P[j-1]->coeff[i], P[j]->coeff[i],
                   pk[j-1]);
    mpz_poly_cleandeg (P[j], P0->deg);
  }
  mpz_poly_cleandeg (P[0], P0->deg);

  for (i = 0; i < l; i++)
    mpz_clear (pk[i]);
  free (pk);

  return P;
}

/* a <- a + pk*P[k] */
void
mpz_poly_base_modp_lift (mpz_poly_t a, mpz_poly_t *P, int k, mpz_t pk)
{
  int i;

  /* first check P[k] exists and is not zero */
  for (i = 0; i <= k; i++)
    if (P[i]->deg == -1)
      return;

  for (i = 0; i <= P[k]->deg; i++)
    mpz_addmul (a->coeff[i], P[k]->coeff[i], pk);

  mpz_poly_cleandeg (a, (a->deg >= P[k]->deg) ? a->deg : P[k]->deg);
}

void
mpz_poly_base_modp_clear (mpz_poly_t *P)
{
  mpz_poly_t *t = P;

  while (1)
  {
    mpz_poly_clear (t[0]);
    if (t[0]->deg == -1)
      break;
    t ++;
  }
  free (P);
}

/* return the maximal size of the coefficients of f in base b */
size_t
mpz_poly_sizeinbase (mpz_poly_t f, int d, int b)
{
  size_t S = 0, s;
  int i;

  for (i = 0; i <= d; i++)
  {
    s = mpz_sizeinbase (f->coeff[i], b);
    if (s > S)
      S = s;
  }
  return S;
}

/* swap f and g */
void
mpz_poly_swap (mpz_poly_t f, mpz_poly_t g)
{
  int i;
  mpz_t *t;

  i = f->alloc;
  f->alloc = g->alloc;
  g->alloc = i;
  i = f->deg;
  f->deg = g->deg;
  g->deg = i;
  t = f->coeff;
  f->coeff = g->coeff;
  g->coeff = t;
}

/* f=gcd(f,g) mod p, with p in mpz_t */
static void
mpz_poly_gcd_mpz (mpz_poly_t f, mpz_poly_t g, const mpz_t p)
{
  while (g->deg >= 0)
    {
      mpz_poly_div_r (f, g, p);
      /* now deg(f) < deg(g): swap f and g */
      mpz_poly_swap (f, g);
    }
}

/* computes d = gcd(f,g) = u*f + v*g mod p, with p in mpz_t */
void
mpz_poly_xgcd_mpz (mpz_poly_t d, const mpz_poly_t f, const mpz_poly_t g, mpz_poly_t u, mpz_poly_t v, const mpz_t p)
{
  mpz_poly_t q, tmp;
  mpz_poly_t gg;

  mpz_poly_init(gg, g->alloc);
  mpz_poly_copy(d, f);
  mpz_poly_copy(gg, g);

  mpz_poly_t uu, vv;
  mpz_poly_init (uu, 0);
  mpz_poly_init (vv, 0);

  u->deg = 0;
  mpz_poly_setcoeff_si(u, 0, 1);
  mpz_poly_set_zero(uu);

  vv->deg = 0;
  mpz_poly_setcoeff_si(vv, 0, 1);
  mpz_poly_set_zero(v);

  mpz_poly_init(q, d->deg);
  mpz_poly_init(tmp, d->deg + gg->deg);
    
  while (gg->deg >= 0)
    {

      /* q, r := f div g mod p */
      mpz_poly_div_qr (q, d, d, gg, p);

      /* u := u - q * uu mod p */
      mpz_poly_mul(tmp, q, uu);
      mpz_poly_sub_mod_mpz(u, u, tmp, p);
      mpz_poly_swap (u, uu);

      /* v := v - q * vv mod p */
      mpz_poly_mul(tmp, q, vv);
      mpz_poly_sub_mod_mpz(v, v, tmp, p);
      mpz_poly_swap (v, vv);
 
      /* now deg(f) < deg(g): swap f and g */
      mpz_poly_swap (d, gg);
    }

  /* make monic */
  mpz_t inv;
  if (mpz_cmp_ui(d->coeff[d->deg], 1) != 0)
    {
      mpz_init(inv);
      mpz_invert(inv, d->coeff[0], p);
      mpz_poly_mul_mpz(d, d, inv);
      mpz_poly_reduce_mod_mpz(d, d, p);
      mpz_poly_mul_mpz(u, u, inv);
      mpz_poly_reduce_mod_mpz(u, u, p);
      mpz_poly_mul_mpz(v, v, inv);
      mpz_poly_reduce_mod_mpz(v, v, p);
      mpz_clear(inv);
    }

  mpz_poly_clear(gg);
  mpz_poly_clear(uu);
  mpz_poly_clear(vv);
  mpz_poly_clear(q);
  mpz_poly_clear(tmp);
}


/* Assuming f is a (squarefree) product of linear factors mod p, splits it
   and put the corresponding roots mod p in r[]. Return number of roots
   which should be degree of f. Assumes p is odd, and deg(f) >= 1. */
static int
mpz_poly_cantor_zassenhaus (mpz_t *r, mpz_poly_t f, const mpz_t p, int depth)
{
  mpz_t a, aux;
  mpz_poly_t q, h, ff;
  int d = f->deg, dq, n, m;

  mpz_init (a);
  mpz_init (aux);

  /* linear polynomial */
  if (d == 1) {
    mpz_neg (aux, f->coeff[1]);
    mpz_invert (a, aux, p);
    mpz_mul (r[0], a, f->coeff[0]);
    mpz_fdiv_r (r[0], r[0], p);
    n = 1;
    goto clear_a;
  }

  /* if f has degree d, then q,h may have up to degree 2d-1 in the
     powering algorithm */
  mpz_poly_init (q, 2 * d - 1);
  mpz_poly_init (h, 2 * d - 1);
  mpz_poly_init (ff, d);

  /* random polynomial by a */
  mpz_set_ui (a, lrand48());
  mpz_fdiv_r (a, a, p);
  for (;;)
  {
    /* q=x+a */
    mpz_set_ui (aux, 1);
    mpz_poly_setcoeff(q, 1, aux);
    mpz_poly_setcoeff(q, 0, a);

    /* h=(x+a)^((p-1)/2) mod (f, p) */
    mpz_sub_ui (aux, p, 1);
    mpz_divexact_ui (aux, aux, 2);
    mpz_poly_power_mod_f_mod_mpz (h, q, f, aux, p);
    mpz_poly_sub_ui (h, 1);
    
    /* q = gcd(f,h) */
    mpz_poly_copy(q, f);
    mpz_poly_gcd_mpz (q, h, p);
    dq = q->deg;
    ASSERT (dq >= 0);

    /* recursion-split */
    if (0 < dq && dq < d)
    {
      n = mpz_poly_cantor_zassenhaus (r, q, p, depth+1);
      ASSERT (n == dq);

      mpz_poly_copy (ff, f);
      /* mpz_poly_divexact clobbers its 2nd arg */
      mpz_poly_divexact (h, ff, q, p);
      m = mpz_poly_cantor_zassenhaus (r + n, h, p, depth + 1);
      ASSERT (m == h->deg);
      n += m;
      break;
    }

    mpz_add_ui (a, a, 1);
    if (mpz_cmp (a, p) > 0)
      mpz_sub (a, a, p);
  }

  mpz_poly_clear (q);
  mpz_poly_clear (h);
  mpz_poly_clear (ff);

clear_a:
  mpz_clear (a);
  mpz_clear (aux);
  return n;
}

typedef int (*sortfunc_t) (const void *, const void *);

static int mpz_poly_coeff_cmp(const mpz_t *a, const mpz_t *b) {
  if (mpz_cmp(*a, *b) < 0) return -1;
  if (mpz_cmp(*a, *b) > 0) return 1;
  return 0;
}

/* Solve f(x)=0 (mod p), where p is a prime. Return nr. */
int
mpz_poly_roots_mpz (mpz_t *r, mpz_t *f, int d, const mpz_t p)
{
  int nr;
  mpz_t tmp;
  mpz_poly_t mpz_poly_fp, mpz_poly_f, g, h;

  mpz_init (tmp);
  mpz_poly_init (mpz_poly_fp, d);
  mpz_poly_init (mpz_poly_f, d);
  mpz_poly_init (g, 2*d-1);
  mpz_poly_init (h, 2*d-1);

  /* reduce f to monic and modulo p */
  mpz_poly_set (mpz_poly_f, f, d);
  mpz_poly_reduce_makemonic_mod_mpz (mpz_poly_fp, mpz_poly_f, p);
  /* h=x^p-x (mod mpz_poly_fp) */
  mpz_set_ui (tmp, 1UL);
  mpz_poly_setcoeff (g, 1, tmp);
  mpz_poly_power_mod_f_mod_mpz (h, g, mpz_poly_fp, p, p);
  mpz_poly_sub (h, h, g);
  /* g = gcd (mpz_poly_fp, h) */
  mpz_poly_gcd_mpz (mpz_poly_fp, h, p);
  /* mpz_poly_fp contains gcd(x^p-x, f) */
  nr = mpz_poly_fp->deg;
  ASSERT (nr >= 0);

  /* If r is NULL, we only return the number of roots. */
  if (r != NULL && nr > 0)
  {
    int n MAYBE_UNUSED = mpz_poly_cantor_zassenhaus (r, mpz_poly_fp, p, 0);
    ASSERT (n == nr);
  }

  mpz_poly_clear(mpz_poly_fp);
  mpz_poly_clear(mpz_poly_f);
  mpz_poly_clear(g);
  mpz_poly_clear(h);
  mpz_clear (tmp);

  /* Sort the roots */
  if (r && nr)
    qsort(r, nr, sizeof(mpz_t), (sortfunc_t) &mpz_poly_coeff_cmp);
  
  return nr;
}

/*  Homographic transform on polynomials */
/* Put in fij[] the coefficients of f'(i) = F(a0*i+a1, b0*i+b1).
   Assumes the coefficients of fij[] are initialized.
*/
void
mpz_poly_homography (mpz_poly_t Fij, mpz_poly_t F, int64_t H[4])
{
  int k, l;
  mpz_t *g; /* will contain the coefficients of (b0*i+b1)^l */
  mpz_t f0;
  mpz_t *f = F->coeff;
  int d = F->deg;
  mpz_t *fij = Fij->coeff;

  mpz_poly_realloc (Fij, d + 1);

  for (k = 0; k <= d; k++)
    mpz_set (fij[k], f[k]);

  mpz_poly_set_deg(Fij, d);

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

/* v <- |f(i,j)|, where f is homogeneous of degree d */
void mpz_poly_homogeneous_eval_siui (mpz_t v, mpz_poly_srcptr f, const int64_t i, const uint64_t j)
{
  unsigned int k;
  mpz_t jpow;

  mpz_init_set_ui (jpow, 1);
  mpz_set (v, f->coeff[f->deg]);
  for (k = f->deg; k-- > 0;)
    {
      mpz_mul_int64 (v, v, i);
      mpz_mul_uint64 (jpow, jpow, j);
      mpz_addmul (v, f->coeff[k], jpow);
    }
  mpz_abs (v, v); /* avoids problems with negative norms */
  mpz_clear (jpow);
}
