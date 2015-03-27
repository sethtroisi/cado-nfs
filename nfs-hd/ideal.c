#include "ideal.h"
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <math.h>
#include "macros.h"
#include <inttypes.h>

/* Ideal */

void ideal_init(ideal_ptr ideal)
{
  mpz_poly_init(ideal->h, -1);
  ideal->r = 0;
}

void ideal_set_part(ideal_ptr ideal, uint64_t r, mpz_poly_srcptr h)
{
  mpz_poly_set(ideal->h, h);
  ideal->r = r;
}

void ideal_clear(ideal_ptr ideal)
{
  mpz_poly_clear(ideal->h);
  ideal->r = 0;
}

void ideal_fprintf(FILE * file, ideal_srcptr ideal)
{
  fprintf(file, "r: %" PRIu64 ", h: ", ideal->r);
  mpz_poly_fprintf(file, ideal->h);
}

/* Ideal with h of degree 1 */

void ideal_1_init(ideal_1_ptr ideal)
{
  ideal_init(ideal->ideal);
  ideal->log = 0;
}

/*
  If you want to line sieve, you can need to have Tr only with the coefficient
   modulo r. Enable this with LINESIEVE.
*/
void ideal_1_set_part(ideal_1_ptr ideal, uint64_t r, mpz_poly_srcptr h,
    unsigned int t)
{
  ASSERT(h->deg == 1);
  ASSERT(mpz_cmp_ui(mpz_poly_lc_const(h), 1) == 0);

  mpz_poly_set(ideal->ideal->h, h);
  //if r == 0, Tr is not created.
  if ((ideal->ideal->r) == 0) {
    ideal->Tr = (mpz_t * ) malloc(sizeof(mpz_t) * (t - 1));
    for (unsigned int i = 0; i < t - 1; i++) {
      mpz_init(ideal->Tr[i]);
    }
  }
  mpz_set(ideal->Tr[0], h->coeff[0]);

  for (unsigned int i = 1; i < t - 1; i++) {
    mpz_set_si(ideal->Tr[i], -1);
    mpz_mul(ideal->Tr[i], ideal->Tr[i], ideal->Tr[i - 1]);
    mpz_mul(ideal->Tr[i], ideal->Tr[i], ideal->Tr[0]);

#ifdef LINESIEVE
    mpz_mod_ui(ideal->Tr[i], ideal->Tr[i], r);
    ASSERT(mpz_cmp(ideal->Tr[i], r) < 0);
#endif // LINESIEVE
  }

  ideal->ideal->r = r;
  ideal->log = (unsigned char)(log2((double) r));
}

void ideal_1_clear(ideal_1_ptr ideal, unsigned int t)
{
  ASSERT(ideal->ideal->h->deg == 1 && ideal->ideal->r != 0);

  for (unsigned int i = 0; i < t - 1; i++) {
    mpz_clear(ideal->Tr[i]);
  }
  free(ideal->Tr);

  ideal_clear(ideal->ideal);
  ideal->log = 0;
}

void ideal_1_set(ideal_1_ptr ideal_new, ideal_1_srcptr ideal_old,
    unsigned int t)
{
  mpz_poly_set(ideal_new->ideal->h, ideal_old->ideal->h);
  ideal_new->ideal->r = ideal_old->ideal->r;
  ideal_new->Tr = (mpz_t * ) malloc(sizeof(mpz_t) * (t - 1));
  for (unsigned int i = 0; i < t - 1; i++) {
    mpz_init(ideal_new->Tr[i]);
    mpz_set(ideal_new->Tr[i], ideal_old->Tr[i]);
  }
  ideal_new->log = ideal_old->log;
}

void ideal_1_fprintf(FILE * file, ideal_1_srcptr ideal, unsigned int t)
{
  ASSERT(ideal->ideal->h->deg == 1);

  ideal_fprintf(file, ideal->ideal);
  fprintf(file, "log: %u\n", ideal->log);
  fprintf(file, "Tr: [");
  for (unsigned int i = 0; i < t - 2; i++) {
    gmp_fprintf(file, "%Zd, ", ideal->Tr[i]);
  }
  gmp_fprintf(file, "%Zd]\n", ideal->Tr[t - 2]);
}

/* Ideal with h of degree u > 1 */

void ideal_u_init(ideal_u_ptr ideal)
{
  ideal_init(ideal->ideal);
  ideal->log = 0;
}

/*
  If you want to line sieve, you can need to have Tr only with the coefficient
  modulo r.
*/
void ideal_u_set_part(ideal_u_ptr ideal, uint64_t r, mpz_poly_srcptr h,
                      unsigned int t)
{
  ASSERT(h->deg > 1);
  ASSERT(mpz_cmp_ui(mpz_poly_lc_const(h), 1) == 0);

  mpz_poly_set(ideal->ideal->h, h);
  //if r == 0, Tr is not created.
  if ((ideal->ideal->r) == 0) {
    ideal->Tr = (mpz_t ** ) malloc(sizeof(mpz_t * ) * (unsigned int)(h->deg));
    for (int row = 0; row < h->deg; row++) {
      ideal->Tr[row] = (mpz_t * ) malloc(sizeof(mpz_t) * (t - (unsigned
                                                               int)h->deg));
    }
    for (int row = 0; row < h->deg; row++) {
      for (int col = 0; col < ((int)t - h->deg); col++) {
        mpz_init(ideal->Tr[row][col]);
      }
    }
  }

#ifdef NDEBUG
  else {
    ASSERT(h->deg == ideal->h->deg);
  }
#endif // NDEBUG

  //Initialize the elements of Tr.
  for (int col = 0; col < (int) t - h->deg; col++) {
    for (int row = col; row < h->deg; row++) {
      mpz_set(ideal->Tr[row][col], h->coeff[row-col]);

#ifdef LINESIEVE
      mpz_mod_ui(ideal->Tr[row][col], ideal->Tr[row][col], r);
      ASSERT(mpz_cmp(ideal->Tr[row][col], r) < 0);
#endif // LINESIEVE
    }
  }

  //Perform the linear combination.
  for (int col = 1; col < (int) t - h->deg; col++) {
    for (int row = 0; row < h->deg; row++) {
      for (int k = 0; k < h->deg; k++) {
        if (col - h->deg + k >= 0) {
          mpz_submul(ideal->Tr[row][col], h->coeff[k],
                     ideal->Tr[row][col - h->deg + k]);

#ifdef LINESIEVE
      mpz_mod_ui(ideal->Tr[row][col], ideal->Tr[row][col], r);
      ASSERT(mpz_cmp(ideal->Tr[row][col], r) < 0);
#endif // LINESIEVE
        }
      }
    }
  }

  ideal->ideal->r = r;
  ideal->log = (unsigned char)(log2((double) r) * h->deg);
}


void ideal_u_clear(ideal_u_ptr ideal, unsigned int t)
{
  ASSERT(ideal->ideal->h->deg > 1 && ideal->ideal->r != 0);

  for (int row = 0; row < ideal->ideal->h->deg; row++) {
    for (int col = 0; col < ((int)t - ideal->ideal->h->deg); col++) {
      mpz_clear(ideal->Tr[row][col]);
    }
  }
  for (int row = 0; row < ideal->ideal->h->deg; row++) {
    free(ideal->Tr[row]);
  }
  free(ideal->Tr);

  ideal_clear(ideal->ideal);
  ideal->log = 0;
}

void ideal_u_fprintf(FILE * file, ideal_u_srcptr ideal, unsigned int t)
{
  ASSERT(ideal->ideal->h->deg > 1);

  ideal_fprintf(file, ideal->ideal);
  fprintf(file, "log: %u\n", ideal->log);

  fprintf(file, "Tr: [");
  for (int row = 0; row < ideal->ideal->h->deg - 1; row++) {
    fprintf(file, "[");
    for (int col = 0; col < (int)t - ideal->ideal->h->deg - 1; col++) {
      gmp_fprintf(file, "%Zd, ", ideal->Tr[row][col]);
    }
    gmp_fprintf(file, "%Zd],\n", ideal->Tr[row]
                [t - (unsigned int)ideal->ideal->h->deg - 1]);
  }
  fprintf(file, "[");
  for (int col = 0; col < (int)t - ideal->ideal->h->deg - 1; col++) {
    gmp_fprintf(file, "%Zd, ", ideal->Tr[ideal->ideal->h->deg - 1][col]);
  }
  gmp_fprintf(file, "%Zd]]\n", ideal->Tr[ideal->ideal->h->deg - 1]
              [t - (unsigned int)ideal->ideal->h->deg - 1]);
}

/* Ideal projective root and Tr */

void ideal_pr_init(ideal_pr_ptr ideal)
{
  ideal_init(ideal->ideal);
  ideal->log = 0;
}


void ideal_pr_set_part(ideal_pr_ptr ideal, uint64_t r, unsigned int t)
{
  mpz_poly_t h;
  mpz_poly_init(h, -1);
  mpz_poly_setcoeff_si(h, 1, 1);
  mpz_poly_set(ideal->ideal->h, h);
  mpz_poly_clear(h);

  if ((ideal->ideal->r) == 0) {
    ideal->Tr = (mpz_t * ) malloc(sizeof(mpz_t) * (t - 1));
    for (unsigned int i = 0; i < t - 1; i++) {
      mpz_init(ideal->Tr[i]);
    }
  }
  for (unsigned int i = 0; i < t - 1; i++) {
    mpz_set_ui(ideal->Tr[i], 0);

  }
  ideal->ideal->r = r;
  ideal->log = (unsigned char)(log2((double) r));
}

void ideal_pr_clear(ideal_pr_ptr ideal, unsigned int t)
{
  ASSERT(ideal->ideal->r != 0);

  for (unsigned int i = 0; i < t - 1; i++) {
    mpz_clear(ideal->Tr[i]);
  }
  free(ideal->Tr);

  ideal_clear(ideal->ideal);
  ideal->log = 0;
}

void ideal_pr_fprintf(FILE * file, ideal_pr_srcptr ideal, unsigned int t)
{
  fprintf(file, "Projective root, ");
  ideal_fprintf(file, ideal->ideal);
  fprintf(file, "log: %u\n", ideal->log);
  fprintf(file, "Tr: [");
  for (unsigned int i = 0; i < t - 2; i++) {
    gmp_fprintf(file, "%Zd, ", ideal->Tr[i]);
  }
  gmp_fprintf(file, "%Zd]\n", ideal->Tr[t - 2]);
}
