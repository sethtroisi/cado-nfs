#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cado.h"
#include "poly.h"
#include "utils/utils.h"

/* store in invm the value of floor(B^(2k)/m), where m has k limbs,
   and B is the limb base */
void
barrett_init (mpz_t invm, const mpz_t m)
{
  size_t k = mpz_size (m);

  mpz_set_ui (invm, 1);
  mpz_mul_2exp (invm, invm, 2 * k * GMP_NUMB_BITS);
  mpz_tdiv_q (invm, invm, m);
}

/* a <- b mod m */
static void
barrett_mod (mpz_t a, mpz_t b, const mpz_t m, const mpz_t invm)
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

// allocate a polynomial that can contain d+1 coeffs, and set to zero.
void poly_alloc(poly_t f, int d) {
  int i;
  f->alloc = d+1;
  f->deg = -1;
  f->coeff = (mpz_t *)malloc((d+1)*sizeof(mpz_t));
  ASSERT (f->coeff != NULL);
  for (i = 0; i <= d; ++i)
    mpz_init(f->coeff[i]);
}

void poly_free(poly_t f) {
  int i;
  for (i = 0; i < f->alloc; ++i)
    mpz_clear(f->coeff[i]);
  free(f->coeff);
}

void poly_print(const poly_t f) {
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


void cleandeg(poly_t f, int deg) {
  while ((deg >= 0) && (mpz_cmp_ui(f->coeff[deg], 0)==0))
    deg--;
  f->deg = deg;
}


void poly_setcoeff(poly_t f, int i, const mpz_t z) {
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

// g must be allocated. 
void poly_copy(poly_t g, const poly_t f) {
  int i;
  g->deg = f->deg;
  for (i = f->deg; i >= 0; --i)
    poly_setcoeff(g, i, f->coeff[i]);
}

#ifndef max
#define max(a,b) ((a)<(b) ? (b) : (a))
#endif

#if 0 /* unused */
void poly_add(poly_t f, const poly_t g, const poly_t h) {
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
      mpz_add(z, z, h->coeff[i]);
    poly_setcoeff(f, i, z);
  }
  mpz_clear(z);
  cleandeg(f, maxdeg);
}
#endif

#if 0 /* unused */
void poly_sub(poly_t f, const poly_t g, const poly_t h) {
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
    poly_setcoeff(f, i, z);
  }
  mpz_clear(z);
  cleandeg(f, maxdeg);
}
#endif

/* f <- g - h mod m */
void
poly_sub_mod_mpz (poly_t f, const poly_t g, const poly_t h, const mpz_t m)
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
    poly_setcoeff(f, i, z);
  }
  mpz_clear(z);
  cleandeg(f, maxdeg);
  for (i = 0; i <= f->deg; ++i) {
    mpz_mod(f->coeff[i], f->coeff[i], m);
  }
  cleandeg(f, f->deg);
}

#if 0 /* unused */
void poly_mul_ui(poly_t f, const poly_t g, unsigned long a) {
  int i;
  mpz_t aux;
  mpz_init(aux);
  for (i = g->deg; i >= 0; --i) {
    mpz_mul_ui(aux, g->coeff[i], a);
    poly_setcoeff(f, i, aux);
  }
  mpz_clear(aux);
}
#endif

/* f <- f - a */
void
poly_sub_ui (poly_t f, unsigned long a)
{
  mpz_t aux;
  mpz_init(aux);
  mpz_sub_ui(aux, f->coeff[0], a);
  poly_setcoeff(f, 0, aux);
  mpz_clear(aux);
}

/* f <- g/2 mod m */
void
poly_div_2_mod_mpz (poly_t f, const poly_t g, const mpz_t m)
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
      poly_setcoeff(f, i, aux);
  }
      mpz_clear (aux);
}

// compute res := f(x) mod m
void poly_eval_mod_mpz(mpz_t res, const poly_t f, const mpz_t x, const mpz_t m) 
{
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

static void poly_mul(poly_t f, const poly_t g, const poly_t h) {
  int i, j, maxdeg;
  poly_t prd;

  if ((g->deg == -1) || (h->deg == -1)) {
    f->deg = -1;
    return;
  }

  maxdeg = g->deg+h->deg;
  poly_alloc(prd, maxdeg);
  for (i = maxdeg; i >= 0; --i)
    mpz_set_ui(prd->coeff[i], 0);
  
  for (i = 0; i <= g->deg; ++i)
    for (j = 0; j <= h->deg; ++j)
      mpz_addmul(prd->coeff[i+j], g->coeff[i], h->coeff[j]);

  for (i = maxdeg; i >= 0; --i) 
    poly_setcoeff(f, i, prd->coeff[i]);
  cleandeg(f, maxdeg);

  poly_free(prd);
}


// Reduce a plain polynomial p modulo a non-monic polynomial F.
// The result is of type polymodF_t P.
// WARNING: this function destroys its input p !!!
static void
poly_reducemodF(polymodF_t P, poly_t p, const poly_t F) {
  int v = 0;

  if (p->deg < F->deg) {
    poly_copy(P->p, p);
    P->v = 0;
    return;
  }

  const int d = F->deg;

  while (p->deg >= d) {
    const int k = p->deg;
    int i;

    /* We compute F[d]*p - p[k]*F. In case F[d] divides p[k],
       we can simply compute p - p[k]/F[d]*F */

    if (mpz_divisible_p (p->coeff[k], F->coeff[d]))
      mpz_divexact (p->coeff[k], p->coeff[k], F->coeff[d]);
    else
      {
        v++; /* we consider p/F[d]^v */
        for (i = 0; i < k; ++i)
          mpz_mul(p->coeff[i], p->coeff[i], F->coeff[d]);
      }

    for (i = 0; i < d; ++i) 
      mpz_submul(p->coeff[k-d+i], p->coeff[k], F->coeff[i]);

    cleandeg(p, k-1);
  }
  
  poly_copy(P->p, p);
  P->v = v;
}

/* Q <- P1 * P2 mod F

   Warning: Q might equal P1.
 */
void
polymodF_mul (polymodF_t Q, const polymodF_t P1, const polymodF_t P2,
              const poly_t F) {
  poly_t prd;
  int v;

  poly_alloc(prd, P1->p->deg+P2->p->deg);

  poly_mul(prd, P1->p, P2->p);
  v = P1->v + P2->v;

  poly_reducemodF(Q, prd, F);
  Q->v += v;
  poly_free(prd);
}

// Q = a*P, where a is an mpz_t
void
poly_mul_mpz(poly_t Q, const poly_t P, const mpz_t a) {
  int i;
  mpz_t aux;
  mpz_init(aux);
  Q->deg = P->deg;
  for (i = 0; i <= P->deg; ++i) {
    mpz_mul(aux, P->coeff[i], a);
    poly_setcoeff(Q, i, aux);
  }
  mpz_clear(aux);
}

#if 0 /* unused */
/* Q <- P mod m, where m = p^k */
void     
poly_reduce_mod_mpz (poly_t Q, const poly_t P, const mpz_t m)
{
  int i;
  mpz_t aux;

  mpz_init(aux);
  Q->deg = P->deg;
  for (i = 0; i <= P->deg; ++i) {
    mpz_mod(aux, P->coeff[i], m);
    poly_setcoeff(Q, i, aux);
  }
  mpz_clear(aux);
}
#endif

/* return a list of polynomials P[0], P[1], ..., P[k] such that
   P0 = Q[k-1] + p^(2^(k-1))*P[k]
   Q[k-1] = Q[k-2] + p^(2^(k-2))*P[k-1]
   ...
   Q[2] = Q[1] + p^2*P[2]
   Q[1] = Q[0] + p*P[1]
   Q[0] = P[0]
   ...
   With all coefficients of P[i] smaller than p^(2^(i-1)).

   P0 = P[0] + p*P[1] + p^2*P[2] + ... + p^(2^(k-1))*P[k] < p^(2^k)

   The end of the list is P[k+1]=0.
*/
poly_t*
poly_base_modp_init (const poly_t P0, int p)
{
  poly_t *P;
  int k, i, j;
  mpz_t *pk;
  size_t l;

  /* first compute p^(2^k) until the coefficients of P0 are smaller than
     p^(2^k) */
  pk = (mpz_t*) malloc (sizeof (mpz_t));
  k = 1;
  mpz_init_set_ui (pk[0], p); /* pk[k] = p^(2^k) */
  while (1)
    {
      /* Check if all coefficients of P0 are smaller than pk[k-1]^2.
         We use an approximate test to avoid squaring pk[k-1]. */
      l = mpz_sizeinbase (pk[k-1], 2);
      /* 2^(2l-2) <= pk[k-1] */
      for (i = 0; i <= P0->deg; i++)
        if (mpz_sizeinbase (P0->coeff[i], 2) > 2 * l - 2)
          break;
      if (i > P0->deg) /* all coeffs are smaller than pk[k-1]^2 */
        break;
      k ++;
      pk = (mpz_t*) realloc (pk, k * sizeof (mpz_t));
      mpz_init (pk[k-1]);
      mpz_mul (pk[k-1], pk[k-2], pk[k-2]);
    }

  /* now decompose P0: we need P[0], P[1] for factor p, P[2] for p^2,
     ..., P[k] for p^(2^(k-1)), and one for the end of list,
     thus k+2 polynomials */
  P = (poly_t*) malloc ((k + 2) * sizeof(poly_t));
  for (i = 0; i < k + 2; i++)
    poly_alloc (P[i], P0->deg);
  /* P[k+1] is initialized to 0 by poly_alloc */

  /* initialize P[k], and put remainder in P[k-1] */
  for (i = 0; i <= P0->deg; i++)
    mpz_tdiv_qr (P[k]->coeff[i], P[k-1]->coeff[i], P0->coeff[i], pk[k-1]);
  cleandeg (P[k], P0->deg);

  /* now go down */
  for (j = k-1; j >= 1; j--)
    {
      /* reduce P[j] into P[j] and P[j-1] */
      for (i = 0; i <= P0->deg; i++)
        mpz_tdiv_qr (P[j]->coeff[i], P[j-1]->coeff[i], P[j]->coeff[i],
                     pk[j-1]);
      cleandeg (P[j], P0->deg);
    }
  cleandeg (P[0], P0->deg);

  for (i = 0; i < k; i++)
    mpz_clear (pk[i]);
  free (pk);

  return P;
}

/* a <- a + pk*P[k] */
void
poly_base_modp_lift (poly_t a, poly_t *P, int k, mpz_t pk)
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
poly_base_modp_clear (poly_t *P)
{
  poly_t *t = P;

  while (1)
    {
      poly_free (t[0]);
      if (t[0]->deg == -1)
        break;
      t ++;
    }
  free (P);
}

void
poly_reduce_makemonic_mod_mpz(poly_t Q, const poly_t P, const mpz_t m) {
  int i;
  mpz_t aux, aux2;
  mpz_init(aux);
  mpz_init(aux2);
  i = P->deg;
  do {
    mpz_mod(aux, P->coeff[i], m);
    i--;
  } while ((i>=0) && (mpz_cmp_ui(aux, 0) == 0));
  if (i>=0) {
    Q->deg = i+1;
    mpz_invert(aux2, aux, m);
    for (i = 0; i <= Q->deg; ++i) {
      mpz_mul(aux, aux2, P->coeff[i]);
      mpz_mod(aux, aux, m);
      poly_setcoeff(Q, i, aux);
    }
  } else {
    if (mpz_cmp_ui(aux, 0) == 0)
      Q->deg = -1;
    else {
      Q->deg = 0;
      mpz_set_ui(aux, 1);
      poly_setcoeff(Q, 0, aux);
    }
  }
  mpz_clear(aux);
  mpz_clear(aux2);
}

// Q = P1*P2 mod F, mod m (assume F monic)
void
poly_mul_mod_f_mod_mpz(poly_t Q, const poly_t P1, const poly_t P2,
                       const poly_t F, const mpz_t m, const mpz_t invm)
{
  int d1 = P1->deg;
  int d2 = P2->deg;
  int d = d1+d2;
  int dF = F->deg;
  poly_t R;
  int i, j;
  mpz_t aux;

  poly_alloc(R, d);
  mpz_init(aux);

  // product mod m
  for (i = 0; i <= d1; ++i)
    for (j = 0; j <= d2; ++j)
      mpz_addmul(R->coeff[i+j], P1->coeff[i], P2->coeff[j]);

  // reduce mod F
  while (d >= dF) {
    barrett_mod (R->coeff[d], R->coeff[d], m, invm);
    for (i = d-1; i >= d-dF; --i) {
      mpz_submul(R->coeff[i], R->coeff[d], F->coeff[dF-d+i]);
    }
    d--;
  }

  /* reduce lower coefficients */
  for (i = 0; i <= d; ++i)
    barrett_mod(R->coeff[i], R->coeff[i], m, invm);

  cleandeg(R, d);
  poly_copy(Q, R);
  poly_free(R);
  mpz_clear(aux);
}

// Q = P^2 mod F, mod m (assume F monic)
void
poly_sqr_mod_f_mod_mpz (poly_t Q, const poly_t P, const poly_t F,
                        const mpz_t m, const mpz_t invm)
{
  int d1 = P->deg;
  int d = d1 + d1;
  int dF = F->deg;
  poly_t R;
  int i, j;
  mpz_t aux;

  poly_alloc(R, d);
  mpz_init(aux);

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

  // reduce mod F
  while (d >= dF) {
    barrett_mod (R->coeff[d], R->coeff[d], m, invm);
    for (i = d-1; i >= d-dF; --i) {
      mpz_submul(R->coeff[i], R->coeff[d], F->coeff[dF-d+i]);
    }
    d--;
  }

  /* reduce lower coefficients */
  for (i = 0; i <= d; ++i)
    barrett_mod (R->coeff[i], R->coeff[i], m, invm);

  cleandeg(R, d);
  poly_copy(Q, R);
  poly_free(R);
  mpz_clear(aux);
}

// Q = P^a mod F, mod p
void 
poly_power_mod_f_mod_ui(poly_t Q, const poly_t P, const poly_t F,
    const mpz_t a, unsigned long p)
{
  mpz_t m;
  int k = mpz_sizeinbase(a, 2);
  poly_t R;

  if (mpz_cmp_ui(a, 0) == 0) {
    Q->deg = 0;
    mpz_set_ui(Q->coeff[0], 1);
    return;
  }

  mpz_init_set_ui(m, p);
  poly_alloc(R, 2*F->deg);
  
  // Initialize R to P
  poly_copy(R, P);

  // Horner
  for (k -= 2; k >= 0; k--)
  {
    poly_sqr_mod_f_mod_mpz(R, R, F, m, NULL);  // R <- R^2
    if (mpz_tstbit(a, k))
      poly_mul_mod_f_mod_mpz(R, R, P, F, m, NULL);  // R <- R*P
  }

  poly_copy(Q, R);
  mpz_clear (m);
  poly_free(R);
}
