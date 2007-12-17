#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include "poly.h"
#include <assert.h>

// allocate a polynomial that can contain d+1 coeffs, and set to zero.
void poly_alloc(poly_t f, int d) {
  int i;
  f->alloc = d+1;
  f->deg = -1;
  f->coeff = (mpz_t *)malloc((d+1)*sizeof(mpz_t));
  assert (f->coeff != NULL);
  for (i = 0; i <= d; ++i)
    mpz_init(f->coeff[i]);
}

void poly_free(poly_t f) {
  int i;
  for (i = 0; i < f->alloc; ++i)
    mpz_clear(f->coeff[i]);
  free(f->coeff);
}



void cleandeg(poly_t f, int deg) {
  while ((deg >= 0) && (mpz_cmp_ui(f->coeff[deg], 0)==0))
    deg--;
  f->deg = deg;
}


void poly_setcoeff(poly_t f, int i, mpz_t z) {
  int j;
  if (i >= f->alloc) {
    f->coeff = (mpz_t *)realloc(f->coeff, (i+1)*sizeof(mpz_t));
    assert (f->coeff != NULL);
    for (j = f->alloc; j <= i; ++j)
      mpz_init(f->coeff[j]);
    f->alloc = i+1;
  }
  mpz_set(f->coeff[i], z);
  if (i >= f->deg) 
    cleandeg(f, i);
}

// g must be allocated. 
void poly_copy(poly_t g, poly_t f) {
  int i;
  g->deg = f->deg;
  for (i = f->deg; i >= 0; --i)
    poly_setcoeff(g, i, f->coeff[i]);
}

#ifndef max
#define max(a,b) ((a)<(b) ? (b) : (a))
#endif

void poly_add(poly_t f, poly_t g, poly_t h) {
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

void poly_sub(poly_t f, poly_t g, poly_t h) {
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

void poly_mul(poly_t f, poly_t g, poly_t h) {
  int i, j, maxdeg;
  mpz_t z;
  poly_t prd;

  if ((g->deg == -1) || (h->deg == -1)) {
    f->deg = -1;
    return;
  }

  mpz_init(z);

  maxdeg = g->deg+h->deg;
  poly_alloc(prd, maxdeg);
  for (i = maxdeg; i >= 0; --i)
    mpz_set_ui(prd->coeff[i], 0);
  
  for (i = 0; i <= g->deg; ++i)
    for (j = 0; j <= h->deg; ++j)
      mpz_addmul(prd->coeff[i+j], g->coeff[i], h->coeff[j]);

  mpz_clear(z);

  for (i = maxdeg; i >= 0; --i) 
    poly_setcoeff(f, i, prd->coeff[i]);
  cleandeg(f, maxdeg);

  poly_free(prd);
}


// Reduce a plain polynomial p modulo a non-monic polynomial F.
// The result is of type polymodF_t P.
// WARNING: this function destroys its input p !!!
void
poly_reducemodF(polymodF_t P, poly_t p, poly_t F) {
  int v = 0;

  if (p->deg < F->deg) {
    poly_copy(P->p, p);
    P->v = 0;
    return;
  }

  while (p->deg >= F->deg) {
    const int k = p->deg;
    const int d = F->deg;
    int i;

    v++;
    for (i = 0; i < k; ++i)
      mpz_mul(p->coeff[i], p->coeff[i], F->coeff[d]);

    for (i = 0; i < d; ++i) 
      mpz_submul(p->coeff[k-d+i], p->coeff[k], F->coeff[i]);

    cleandeg(p, k-1);
  }

  poly_copy(P->p, p);
  P->v = v;
}

void
polymodF_mul(polymodF_t Q, polymodF_t P1, polymodF_t P2, poly_t F) {
  poly_t prd;
  int v;
  poly_alloc(prd, P1->p->deg+P2->p->deg);

  poly_mul(prd, P1->p, P2->p);
  v = P1->v + P2->v;

  poly_reducemodF(Q, prd, F);
  Q->v += v;
  poly_free(prd);
}

