/**
   These files (mpz_poly.*) implement arithmetics of polynomials whose
   coefficients are in multiprecision integers (using mpz_t GNUMP).
   We use them in sqrt/algsqrt.c to represent rings of integers
   and in fast_rootsieve.c, to represent the involved polynomials.

   Please see the file README_POLYNOMIALS for more details and
   comparisons to other poly formats.
*/


#include "cado.h"
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "portability.h"

#ifndef max
#define max(a,b) ((a)<(b) ? (b) : (a))
#endif


/* -----------------
   Static functions
   ----------------- */


static int
mpz_poly_normalized_p (const mpz_poly_t f) {
  return (f->deg == -1) || mpz_cmp_ui (f->coeff[f->deg], 0) != 0;
}

/* normalize h so that h->coeff[deg(h)] <> 0 */
static void
mpz_poly_normalize (mpz_poly_t h, const mpz_t p)
{
  int dh = h->deg;
  while (dh >= 0 && mpz_divisible_p(h->coeff[dh],p))
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
        mpz_mul_ui (f[i], f[i], g);
        mpz_submul_ui (f[i], f[j], h);
        for (k = j; k <= t; k++)
          M[i][k] = g * M[i][k] - h * M[j][k];
      }

  /* now zero upper-diagonal coefficients while going up */
  for (i = t; i >= 0; i--)
  {
    for (j = i + 1; j <= t; j++)
      /* f[i] = f[i] - M[i][j] * f[j] */
      mpz_submul_ui (f[i], f[j], M[i][j]);
    ASSERT (mpz_divisible_ui_p (f[i], M[i][i]));
    mpz_divexact_ui (f[i], f[i], M[i][i]);
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

/* Reduce a plain polynomial p modulo a non-monic polynomial F.
   The result is of type polymodF_t P.
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

    cleandeg (p, k-1);
  }

  mpz_poly_copy(P->p, p);
  P->v = v;
}


/* ----------------
   Public functions
   ---------------- */

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


/* Allocate a polynomial that can contain 'nc' coefficients and set to zero.
   We allow nbcoeff=0.
 */
/* void mpz_poly_alloc(mpz_poly_t f, int nc) { */
/*   int i; */
/*   f->alloc = nc; */
/*   f->deg = -1; */
/*   if (nc == 0) */
/*     f->coeff = (mpz_t *) NULL; */
/*   else */
/*     { */
/*       f->coeff = (mpz_t *) malloc (nc*sizeof(mpz_t)); */
/*       ASSERT (f->coeff != NULL); */
/*     } */
/*   for (i = 0; i < nc; ++i) */
/*     mpz_init(f->coeff[i]); */
/* } */

/* realloc f to (at least) nc coefficients */
void mpz_poly_realloc (mpz_poly_t f, int nc) {
  if (f->alloc < nc) {
    f->coeff = (mpz_t*) realloc (f->coeff, nc * sizeof (mpz_t));
    if (f->coeff == NULL) {
      fprintf (stderr, "Error, not enough memory\n");
      exit (1);
    }
    f->alloc = nc;
  }
}

/* Free polynomial f in mpz_poly_t. */
void mpz_poly_free(mpz_poly_t f) {
  int i;
  for (i = 0; i < f->alloc; ++i)
    mpz_clear(f->coeff[i]);
  free(f->coeff);
}

/* Print coefficients of f. */
void mpz_poly_print(const mpz_poly_t f) {
  int i;
  if (f->deg == -1) {
    printf("0\n");
    return;
  } else if (f->deg == 0) {
    gmp_printf("%Zd\n", f->coeff[0]);
    return;
  }
  gmp_printf("%Zd + ", f->coeff[0]);
  for (i = 1; i < f->deg; ++i)
    gmp_printf("%Zd*x^%d + ", f->coeff[i], i);
  gmp_printf("%Zd*x^%d\n", f->coeff[f->deg], f->deg);
}

/* Find polynomial degree. */
void cleandeg(mpz_poly_t f, int deg) {
  ASSERT(deg >= -1);
  while ((deg >= 0) && (mpz_cmp_ui(f->coeff[deg], 0)==0))
    deg--;
  f->deg = deg;
}

/* Sets f to the polynomial of degree d, of coefficients
   given by coeffs. */
void mpz_poly_set(mpz_poly_t f, mpz_t * coeffs, int d) {
  int i;
  for (i=d; i>=0; --i)
    mpz_poly_setcoeff(f,i,coeffs[i]);
  cleandeg(f,d);
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
    cleandeg(f, i);
}

/* Set uint64 coefficient for the i-th term. */
void mpz_poly_setcoeff_uint64(mpz_poly_t f, int i, uint64_t z) {
  int j;
  if (i >= f->alloc) {
    f->coeff = (mpz_t *)realloc(f->coeff, (i+1)*sizeof(mpz_t));
    ASSERT (f->coeff != NULL);
    for (j = f->alloc; j <= i; ++j)
      mpz_init(f->coeff[j]);
    f->alloc = i+1;
  }
  mpz_set_uint64(f->coeff[i], z);
  if (i >= f->deg)
    cleandeg(f, i);
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
    cleandeg(f, i);
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
    cleandeg(f, i);
}

/* Set coefficient for the i-th term with string input. */
void mpz_poly_setcoeff_str(mpz_poly_t f, int i,char *str, int base) {
  mpz_t z;
  mpz_init_set_str(z,str,base);
  mpz_poly_setcoeff(f,i,z);
  mpz_clear(z);
}

/* Get coefficients from f. Array coeffs[] needs length d+1. */
void mpz_poly_get(mpz_poly_t f, mpz_t * coeffs, int d) {
  int i;
  for (i=d; i>=0; --i)
    mpz_poly_getcoeff(coeffs[i],i,f);
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

/* Copy f to g, where g must be allocated. */
void mpz_poly_copy(mpz_poly_t g, const mpz_poly_t f) {
  int i;
  g->deg = f->deg;
  for (i = f->deg; i >= 0; --i)
    mpz_poly_setcoeff(g, i, f->coeff[i]);
}

/* Set f=g+h, which is used in fast_rootsieve.
   Note: f can be the same as g or h;
         g can be the same as h. */
void mpz_poly_add(mpz_poly_t f, const mpz_poly_t g, const mpz_poly_t h) {

  if(f==h || f==g) {
    mpz_poly_t aux;
    mpz_poly_init(aux,-1);
    mpz_poly_add(aux,g,h);
    mpz_poly_copy(f,aux);
    mpz_poly_free(aux);
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
  cleandeg(f, maxdeg);
}

/* Set f=g-h, which is used in fast_rootsieve.
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
  cleandeg(f, maxdeg);
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
  cleandeg(f, maxdeg);
  for (i = 0; i <= f->deg; ++i) {
    mpz_mod(f->coeff[i], f->coeff[i], m);
  }
  cleandeg(f, f->deg);
}

/* Set f=g*h. Note: f differs from g and h. */
void mpz_poly_mul(mpz_poly_t f, const mpz_poly_t g, const mpz_poly_t h) {
  int i, maxdeg;
  mpz_poly_t prd;

  if(f==h || f==g) {
    mpz_poly_t aux;
    mpz_poly_init(aux,-1);
    mpz_poly_mul(aux,g,h);
    mpz_poly_copy(f,aux);
    mpz_poly_free(aux);
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
  cleandeg(f, maxdeg);
  ASSERT_ALWAYS(mpz_cmp_ui (f->coeff[f->deg], 0) != 0);

  mpz_poly_free(prd);
}

/* Set f=g*a where a is unsigned long. */
void mpz_poly_mul_ui(mpz_poly_t f, const mpz_poly_t g, unsigned long a) {
  int i;
  mpz_t aux;
  mpz_init(aux);
  cleandeg(f,g->deg);
  for (i = g->deg; i >= 0; --i) {
    mpz_mul_ui(aux, g->coeff[i], a);
    mpz_poly_setcoeff(f, i, aux);
  }
  mpz_clear(aux);
}

/* Set Q=a*P, where a is an mpz_t */
void
mpz_poly_mul_mpz(mpz_poly_t Q, const mpz_poly_t P, const mpz_t a) {
  int i;
  mpz_t aux;
  mpz_init(aux);
  Q->deg = P->deg;
  for (i = 0; i <= P->deg; ++i) {
    mpz_mul(aux, P->coeff[i], a);
    mpz_poly_setcoeff(Q, i, aux);
  }
  mpz_clear(aux);
}

/* h=rem(h, f) mod p, f not necessarily monic. */
void
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
  cleandeg(r, r->deg);

  mpz_clear(invlg);
  mpz_clear(lg);
  mpz_clear(tmp);
}

/* q=divexact(h, f) mod p, f not necessarily monic. Clobbers h. */
void
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
  mpz_poly_normalize (q, p);

  mpz_clear (t);
  mpz_clear (aux);
}

/* Set f=g/a. Note: it calls mpz_divexact_ui. */
void mpz_poly_div_ui(mpz_poly_t f, const mpz_poly_t g, unsigned long a) {
  int i;
  mpz_t aux;
  mpz_init(aux);
  for (i = g->deg ; i>=0; --i) {
    mpz_divexact_ui(aux,g->coeff[i],a);
    mpz_poly_setcoeff(f,i,aux);
  }
  mpz_clear(aux);
}

/* Set f=g/a (mod m) where a, m are unsigned long. */
void mpz_poly_div_ui_mod_ui(mpz_poly_t f, const mpz_poly_t g, const unsigned long a,
                        const unsigned long m) {
  int i;
  mpz_t aux;
  mpz_init(aux);
  for (i = g->deg ; i>=0; --i) {
    mpz_divexact_ui(aux,g->coeff[i],a);
    mpz_mod_ui(aux,aux,m);
    mpz_poly_setcoeff(f,i,aux);
  }
  mpz_clear(aux);
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
                       const mpz_t m) {
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
    mpz_mod(res, res, m);
  }
}

/* Set Q=P1*P2 (mod F). Warning: Q might equal P1. */
void polymodF_mul (polymodF_t Q, const polymodF_t P1, const polymodF_t P2,
                   const mpz_poly_t F) {
  mpz_poly_t prd;
  int v;

  mpz_poly_init(prd, P1->p->deg+P2->p->deg);

  ASSERT_ALWAYS(mpz_poly_normalized_p (P1->p));
  ASSERT_ALWAYS(mpz_poly_normalized_p (P2->p));

  mpz_poly_mul(prd, P1->p, P2->p);
  v = P1->v + P2->v;

  mpz_poly_reducemodF(Q, prd, F);
  Q->v += v;
  mpz_poly_free(prd);
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
  cleandeg(Q,Q->deg);
}

/* Set Q = P/lc(P) (mod m). Q and P might be identical. */
void
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

#if 0 /* quadratic product */
  int j;
  // product mod m
  for (i = 0; i <= d1; ++i)
    for (j = 0; j <= d2; ++j)
      mpz_addmul(R->coeff[i+j], P1->coeff[i], P2->coeff[j]);
#else /* fast product */
  mpz_poly_mul_tc (R->coeff, P1->coeff, d1, P2->coeff, d2);
#endif

  // reduce mod f
  d = mpz_poly_mod_f_mod_mpz (R->coeff, d, f->coeff, df, m, invm);

  cleandeg(R, d);
  mpz_poly_copy(Q, R);
  mpz_poly_free(R);
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

#if 0 /* Naive squaring in d1*(d1+1)/2 products and d1+1 squares, i.e.,
         d*(d-1)/2 products and d squares for an algebraic polynomial of
         degree d. For d=5 (d1=4), this gives 10 products and 5 squares. */
  int j;
  // product mod m
  /* first accumulate triangular terms */
  for (i = 0; i < d1; ++i)
    for (j = i + 1; j <= d1; ++j)
      mpz_addmul(R->coeff[i+j], P->coeff[i], P->coeff[j]);
  /* then double triangular terms */
  for (i = 1; i < d; i++)
    mpz_mul_2exp (R->coeff[i], R->coeff[i], 1);
  /* finally accumulate diagonal terms */
  for (i = 0; i <= d1; ++i)
    mpz_addmul(R->coeff[2*i], P->coeff[i], P->coeff[i]);
#else
  /* Fast squaring in 2d1+1 squares, i.e., 2d-1 squares.
     For d=5, this gives 9 squares. */
  mpz_poly_sqr_tc (R->coeff, P->coeff, d1);
#endif

  // reduce mod f
  d = mpz_poly_mod_f_mod_mpz (R->coeff, d, f->coeff, df, m, invm);

  cleandeg(R, d);
  mpz_poly_copy(Q, R);
  mpz_poly_free(R);
}

/* Affects the derivative of f to df. Assumes df different from f.
   Assumes df has been initialized with degree at least f->deg-1. */
void mpz_poly_derivative(mpz_poly_t df, const mpz_poly_t f) {
  int n;
  if (mpz_poly_is_constant(f)) {
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
  mpz_poly_free(R);
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
  mpz_poly_free(R);
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
  cleandeg (P[k], P0->deg);

  /* now go down */
  for (j = l-1; j >= 1; j--)
  {
    /* reduce P[j] into P[j] and P[j-1] */
    for (i = 0; i <= P0->deg; i++)
      mpz_tdiv_qr (P[j]->coeff[i], P[j-1]->coeff[i], P[j]->coeff[i],
                   pk[j-1]);
    cleandeg (P[j], P0->deg);
  }
  cleandeg (P[0], P0->deg);

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

  cleandeg (a, (a->deg >= P[k]->deg) ? a->deg : P[k]->deg);
}

void
mpz_poly_base_modp_clear (mpz_poly_t *P)
{
  mpz_poly_t *t = P;

  while (1)
  {
    mpz_poly_free (t[0]);
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

int mpz_poly_is_constant(const mpz_poly_t f) {
  return (f->deg <= 0);
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
void
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

  mpz_poly_init(d, f->alloc);
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

  mpz_poly_free(gg);
  mpz_poly_free(uu);
  mpz_poly_free(vv);
  mpz_poly_free(q);
  mpz_poly_free(tmp);
}


/* Assuming f is a (squarefree) product of linear factors mod p, splits it
   and put the corresponding roots mod p in r[]. Return number of roots
   which should be degree of f. Assumes p is odd, and deg(f) >= 1. */
int
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

  mpz_poly_free (q);
  mpz_poly_free (h);
  mpz_poly_free (ff);

clear_a:
  mpz_clear (a);
  mpz_clear (aux);
  return n;
}

typedef int (*sortfunc_t) (const void *, const void *);

int mpz_poly_coeff_cmp(const mpz_t *a, const mpz_t *b) {
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

  mpz_poly_free(mpz_poly_fp);
  mpz_poly_free(mpz_poly_f);
  mpz_poly_free(g);
  mpz_poly_free(h);
  mpz_clear (tmp);

  /* Sort the roots */
  if (r && nr)
    qsort(r, nr, sizeof(mpz_t), (sortfunc_t) &mpz_poly_coeff_cmp);
  
  return nr;
}

#if 0
/* returns z mod (B-1) where B is 2^GMP_NUMB_BITS */
static mp_limb_t
mod_base_minus_1 (mpz_t z)
{
  size_t i;
  mp_limb_t r = (mp_limb_t) 0, s;
  for (i = 0; i <= mpz_size (z); i++)
  {
    s = mpz_getlimbn (z, i);
    r += s;
    /* If there is a carry, it should be added back to r.
       In that case, since r, s <= B-1, r+s <=2B-2, thus
       the remainder without the carry is <= B-2, and the
       following cannot give a carry. */
    r += (r < s);
  }
  /* if r = B-1, reduce to 0 */
  if (r + 1 == (mp_limb_t) 0)
    r = (mp_limb_t) 0;
  return r;
}
#endif


