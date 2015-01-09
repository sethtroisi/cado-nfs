#include <stdint.h>
#include <gmp.h>
#include "macros.h"
#include "getprime.h"
#include "makefb.h"
#include "utils_int64.h"

/*
  mpz_poly_factor does not work in characteristic 2, so we do it the naive way.
*/

/*
  Add 1 to f. If the constant term is equal to 1, set this term to 0 and
  propagate the addition of 1 to the next coefficient, and so on.

  f: the polynomial on which the addition is computed, the modifications are
  made on f.
*/
static void mpz_poly_add_one_in_F2(mpz_poly_ptr f)
{
  ASSERT(f->deg >= 1);

  int i = 0;
  while (mpz_cmp_ui(f->coeff[i], 1) == 0) {
    mpz_poly_setcoeff_si(f, i, 0);
    i++;
    if (i > f->deg) {
      break;
    }
  }
  mpz_poly_setcoeff_si(f, i, 1);
}

/*
  Factorise a polynomial in characteristic 2.

  list: list of factors of f.
  f: the polynomial we want to factorise.
*/
static void mpz_poly_factor2(mpz_poly_factor_list_ptr list, mpz_poly_srcptr f)
{
  mpz_t p;
  mpz_init(p);
  mpz_set_ui(p, 2);

  mpz_poly_t fcopy;
  mpz_poly_init(fcopy, f->deg);
  mpz_poly_set(fcopy, f);

  mpz_poly_factor_list_flush(list);

  if (mpz_poly_is_irreducible(f, p)) {
    mpz_poly_factor_list_push(list, f, 1);
  } else {
    mpz_poly_t tmp;
    mpz_poly_init(tmp, 1);
    mpz_poly_setcoeff_int64(tmp, 1, 1);

    while (tmp->deg <= fcopy->deg) {
      if (mpz_poly_is_irreducible(tmp, p)) {
        mpz_poly_t q;
        mpz_poly_init(q, 0);
        mpz_poly_t r;
        mpz_poly_init(r, 0);
        mpz_poly_div_qr(q, r, fcopy, tmp, p);
        unsigned int m = 0;
        while (r->deg == -1) {
          m++;
          mpz_poly_set(fcopy, q);
          if (fcopy->deg == 0 || fcopy->deg == -1) {
            break;
          }
          mpz_poly_div_qr(q, r, fcopy, tmp, p);
        }
        if (m != 0) {
          mpz_poly_factor_list_push(list, tmp, m);
        }
        mpz_poly_clear(q);
        mpz_poly_clear(r);
      }
      mpz_poly_add_one_in_F2(tmp);
    }
    mpz_poly_clear(tmp);
  }

#ifndef NDBEBUG
  mpz_poly_cleandeg(fcopy, -1);
  mpz_t coeff;
  mpz_init(coeff);
  for (int i = 0; i <= f->deg; i++) {
    mpz_poly_getcoeff(coeff, i, f);
    mpz_mod(coeff, coeff, p);
    mpz_poly_setcoeff(fcopy, i, coeff);
  }

  mpz_poly_t fmul;
  mpz_poly_init(fmul, -1);
  mpz_poly_set(fmul, list->factors[0]->f);
  for (int j = 1; j < list->factors[0]->m; j++) {
    mpz_poly_mul(fmul, fmul, list->factors[0]->f);
  }
  for (int i = 1; i < list->size ; i++) {
    for (int j = 0; j < list->factors[i]->m; j++) {
      mpz_poly_mul(fmul, fmul, list->factors[i]->f);
    }
  }
  for (int i = 0; i <= fmul->deg; i++) {
    mpz_poly_getcoeff(coeff, i, fmul);
    mpz_mod(coeff, coeff, p);
    mpz_poly_setcoeff(fmul, i, coeff);
  }
  mpz_clear(coeff);

  ASSERT(mpz_poly_cmp(fcopy, fmul) == 0);

  mpz_poly_clear(fmul);
#endif // NDBEBUG

  mpz_poly_clear(fcopy);
  mpz_clear(p);
}

/*
  Set an ideal_1 at an index. Do not forget to define LINESIEVE if you
   want to set an ideal mod r.

  fb: the factor base.
  index: index in the factor base.
  r: the r of the ideal (r, h).
  h: the h of the ideal (r, h).
  fbb: factor base bound for this side.
*/
static void add_ideal_1(factor_base_ptr fb, uint64_t * index, uint64_t r,
                        mpz_poly_srcptr h, uint64_t fbb, unsigned int t)
{
  ASSERT(h->deg == 1);

  if (r <= fbb) {
    factor_base_set_ideal_1_part(fb, * index, r, h, t);
    * index = * index + 1;
  }
}

/*
  Set an ideal_1 at an index. Do not forget to define LINESIEVE if you
   want to set an ideal mod r.

  fb: the factor base.
  index: index in the factor base.
  r: the r of the ideal (r, h).
  h: the h of the ideal (r, h).
  fbb: factor base bound for this side.
  lpb: large prime bound.
*/
static void add_ideal_u(factor_base_ptr fb, uint64_t * index, uint64_t r,
                        mpz_poly_srcptr h, uint64_t fbb, mpz_t lpb,
                        unsigned int t)
{
  ASSERT(h->deg > 1);

  if (mpz_cmp_ui(lpb, pow_uint64_t(r, h->deg)) >= 0 && r <= fbb) {
    factor_base_set_ideal_u_part(fb, * index, r, h, t);
    * index = * index + 1;
  }
}

/*
  Set an ideal_pr at an index. Do not forget to define LINESIEVE if you
   want to set an ideal mod r.

  fb: the factor base.
  index: index in the factor base.
  r: the r of the ideal (r, h).
  fbb: factor base bound for this side.
*/
static void add_ideal_pr(factor_base_ptr fb, uint64_t * index, uint64_t r,
                         uint64_t fbb, unsigned int t)
{
  if (r <= fbb) {
    factor_base_set_ideal_pr(fb, * index, r, t);
    * index = * index + 1;
  }
}

void makefb(factor_base_t * fb, mpz_poly_t * f, uint64_t * fbb, unsigned int t,
            mpz_t * lpb, unsigned int V)
{

  ASSERT(V >= 2);
  ASSERT(t >= 2);

#ifndef NDEBUG
  for (unsigned int i = 0; i < V; i++) {
    ASSERT(f[i]->deg >= 1);
    ASSERT(fbb[i] > 2);
    ASSERT(mpz_cmp_ui(lpb[i], fbb[i]) >= 0);
  }
#endif // NDEBUG

  uint64_t q = 2;
  uint64_t * index1 = malloc(sizeof(uint64_t) * V);
  uint64_t * indexu = malloc(sizeof(uint64_t) * V);
  uint64_t * indexpr = malloc(sizeof(uint64_t) * V);
  gmp_randstate_t state;
  mpz_t a;
  mpz_poly_factor_list l;

  mpz_t zero;
  mpz_init(zero);
  mpz_set_ui(zero, 0);

  mpz_t lc;
  mpz_init(lc);

  for (unsigned int k = 0; k < V; k++) {
    index1[k] = 0;
    indexu[k] = 0;
    indexpr[k] = 0;
  }

  gmp_randinit_default(state);
  mpz_init(a);
  mpz_poly_factor_list_init(l);

  mpz_set_ui(a, 2);
  for (unsigned int k = 0; k < V; k++) {
    mpz_set(lc, mpz_poly_lc_const(f[k]));
    if (mpz_congruent_p(lc, zero, a) != 0) {
      add_ideal_pr(fb[k], indexpr + k, q, fbb[k], t);
    }
    mpz_poly_factor2(l, f[k]);
    for (int i = 0; i < l->size ; i++) {
      if (l->factors[i]->f->deg == 1) {
        add_ideal_1(fb[k], index1 + k, q, l->factors[i]->f, fbb[k], t);
      } else if (l->factors[i]->f->deg < (int)t) {
        add_ideal_u(fb[k], indexu + k, q, l->factors[i]->f, fbb[k], lpb[k], t);
      }
    }
  }

  q = getprime(q);
  uint64_t qmax = fbb[0];
  for (unsigned int k = 1; k < V; k++) {
    qmax = MAX(qmax, fbb[k]);
  }

  for ( ; q <= qmax; q = getprime(q)) {
    mpz_set_si(a, q);
    for (unsigned int k = 0; k < V; k++) {
      mpz_set(lc, mpz_poly_lc_const(f[k]));
      if (mpz_congruent_p(lc, zero, a) != 0) {
        add_ideal_pr(fb[k], indexpr + k, q, fbb[k], t);
      }
      if (q <= fbb[k]) {
        mpz_poly_factor(l, f[k], a, state);
        for (int i = 0; i < l->size ; i++) {
          if (l->factors[i]->f->deg == 1) {
            add_ideal_1(fb[k], index1 + k, q, l->factors[i]->f, fbb[k], t);
          } else if (l->factors[i]->f->deg < (int)t) {
            add_ideal_u(fb[k], indexu + k, q, l->factors[i]->f, fbb[k], lpb[k],
                        t);
          }
        }
      }
    }
  }

  mpz_poly_factor_list_clear(l);

  for (unsigned int k = 0; k < V; k++) {
    factor_base_realloc(fb[k], index1[k], indexu[k], indexpr[k]);
  }

  mpz_clear(lc);
  mpz_clear(zero);
  free(index1);
  free(indexu);
  free(indexpr);
  gmp_randclear(state);
  mpz_clear(a);
  getprime(0);
}
