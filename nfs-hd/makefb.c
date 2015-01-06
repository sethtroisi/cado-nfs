#include <stdint.h>
#include <gmp.h>
#include "macros.h"
#include "getprime.h"
#include "makefb.h"
#include "utils_int64.h"


// mpz_poly_factor does not work in char 2, so we do it the naive way.

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
        mpz_poly_div_qr (q, r, fcopy, tmp, p);
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
  mpz_poly_clear(fcopy);
  mpz_clear(p);
}

void add_ideal_1(factor_base_ptr fb, uint64_t * index, uint64_t r,
                 mpz_poly_srcptr h, uint64_t fbb, unsigned int t)
{
  ASSERT(h->deg == 1);

  if (r <= fbb) {
    factor_base_set_ideal_1_part(fb, * index, r, h, t);
    * index = * index + 1;
  }
}

void add_ideal_u(factor_base_ptr fb, uint64_t * index, uint64_t r,
                 mpz_poly_srcptr h, uint64_t fbb, mpz_t lpb, unsigned int t)
{
  ASSERT(h->deg > 1);

  if (mpz_cmp_ui(lpb, pow_uint64_t(r, h->deg)) >= 0 && r <= fbb) {
    factor_base_set_ideal_u_part(fb, * index, r, h, t);
    * index = * index + 1;
  }
}

void makefb(factor_base_t * fb, mpz_poly_t * f, uint64_t * fbb, unsigned int t,
            mpz_t * lpb, unsigned int V)
{
  uint64_t q = 2;
  uint64_t * index1 = malloc(sizeof(uint64_t) * V);
  uint64_t * indexu = malloc(sizeof(uint64_t) * V);
  gmp_randstate_t state;
  mpz_t a;
  mpz_poly_factor_list l;

  for (unsigned int k = 0; k < V; k++) {
    factor_base_init(fb[k], fbb[k], fbb[k]);
    index1[k] = 0;
    indexu[k] = 0;
  }

  gmp_randinit_default(state);
  mpz_init(a);
  mpz_poly_factor_list_init(l);

  for (unsigned int k = 0; k < V; k++) {
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
    for (unsigned int k = 0; k < V; k++) {
      mpz_set_si(a, q);
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
    factor_base_realloc(fb[k], index1[k], indexu[k]);
  }

  free(index1);
  free(indexu);
  gmp_randclear(state);
  mpz_clear(a);
  getprime(0);
}
