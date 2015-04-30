#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <time.h>
#include "cado.h"
#include "utils.h"
#include "functions.h"

MAYBE_UNUSED static void extended_euclidean_algorithm(mpz_ptr r0, mpz_ptr t0,
    mpz_srcptr f, mpz_srcptr g)
{
  mpz_t s0;
  mpz_init(s0);
  mpz_t r1;
  mpz_init(r1);
  mpz_t s1;
  mpz_init(s1);
  mpz_t t1;
  mpz_init(t1);
  mpz_t r2;
  mpz_init(r2);
  mpz_t s2;
  mpz_init(s2);
  mpz_t t2;
  mpz_init(t2);
  mpz_t q;
  mpz_init(q);

  mpz_set(r0, f);
  mpz_set_si(s0, 1);
  mpz_set_si(t0, 0);
  mpz_set(r1, g);
  mpz_set_si(s1, 0);
  mpz_set_si(t1, 1);

  unsigned int i = 1;
  while (mpz_cmp_ui(r1, 0) != 0) {
    mpz_fdiv_qr(q, r2, r0, r1);
    mpz_set(r2, r0);
    mpz_submul(r2, q, r1);
    mpz_set(s2, s0);
    mpz_submul(s2, q, s1);
    mpz_set(t2, t0);
    mpz_submul(t2, q, t1);

#ifndef NDEBUG
    mpz_t tmp;
    mpz_init(tmp);
    mpz_set(tmp, f);
    mpz_mul(tmp, tmp, s1);
    mpz_addmul(tmp, t1, g);
    ASSERT(mpz_cmp(tmp, r1) == 0);
    mpz_clear(tmp);
#endif

    mpz_set(r0, r1);
    mpz_set(s0, s1);
    mpz_set(t0, t1);
    mpz_set(r1, r2);
    mpz_set(s1, s2);
    mpz_set(t1, t2);
    i++;
  }

  mpz_clear(s0);
  mpz_clear(r1);
  mpz_clear(s1);
  mpz_clear(t1);
  mpz_clear(r2);
  mpz_clear(s2);
  mpz_clear(t2);
  mpz_clear(q);
}

#ifdef EEA_BOUND
static void extended_euclidean_algorithm_bound(mpz_ptr r0,
    mpz_ptr t0, mpz_srcptr f, mpz_srcptr g, mpz_srcptr bound)
{
  mpz_t s0;
  mpz_init(s0);
  mpz_t r1;
  mpz_init(r1);
  mpz_t s1;
  mpz_init(s1);
  mpz_t t1;
  mpz_init(t1);
  mpz_t r2;
  mpz_init(r2);
  mpz_t s2;
  mpz_init(s2);
  mpz_t t2;
  mpz_init(t2);
  mpz_t q;
  mpz_init(q);

  mpz_set(r0, f);
  mpz_set_si(s0, 1);
  mpz_set_si(t0, 0);
  mpz_set(r1, g);
  mpz_set_si(s1, 0);
  mpz_set_si(t1, 1);

  unsigned int i = 1;
  while (mpz_cmp_ui(r1, 0) != 0 && mpz_cmp(r1, bound) <= 0) {
    mpz_fdiv_qr(q, r2, r0, r1);
    mpz_set(r2, r0);
    mpz_submul(r2, q, r1);
    mpz_set(s2, s0);
    mpz_submul(s2, q, s1);
    mpz_set(t2, t0);
    mpz_submul(t2, q, t1);

#ifndef NDEBUG
    mpz_t tmp;
    mpz_init(tmp);
    mpz_set(tmp, f);
    mpz_mul(tmp, tmp, s1);
    mpz_addmul(tmp, t1, g);
    ASSERT(mpz_cmp(tmp, r1) == 0);
    mpz_clear(tmp);
#endif

    mpz_set(r0, r1);
    mpz_set(s0, s1);
    mpz_set(t0, t1);
    mpz_set(r1, r2);
    mpz_set(s1, s2);
    mpz_set(t1, t2);
    i++;
  }

  mpz_clear(s0);
  mpz_clear(r1);
  mpz_clear(s1);
  mpz_clear(t1);
  mpz_clear(r2);
  mpz_clear(s2);
  mpz_clear(t2);
  mpz_clear(q);
}
#else
static void extended_euclidean_algorithm_stop_2(mpz_ptr r2, mpz_ptr t2,
    mpz_srcptr f, mpz_srcptr g)
{
  mpz_t r0;
  mpz_init(r0);
  mpz_t s0;
  mpz_init(s0);
  mpz_t t0;
  mpz_init(t0);
  mpz_t r1;
  mpz_init(r1);
  mpz_t s1;
  mpz_init(s1);
  mpz_t t1;
  mpz_init(t1);
  mpz_t s2;
  mpz_init(s2);
  mpz_t q;
  mpz_init(q);

  mpz_set(r0, f);
  mpz_set_si(s0, 1);
  mpz_set_si(t0, 0);
  mpz_set(r1, g);
  mpz_set_si(s1, 0);
  mpz_set_si(t1, 1);

  mpz_fdiv_qr(q, r2, r0, r1);
  mpz_set(r2, r0);
  mpz_submul(r2, q, r1);
  mpz_set(s2, s0);
  mpz_submul(s2, q, s1);
  mpz_set(t2, t0);
  mpz_submul(t2, q, t1);

  mpz_clear(r0);
  mpz_clear(s0);
  mpz_clear(t0);
  mpz_clear(r1);
  mpz_clear(s1);
  mpz_clear(t1);
  mpz_clear(s2);
  mpz_clear(q);
}
#endif

static double find_epsilon(mpz_srcptr q, mpz_srcptr p, unsigned int t)
{
  double q_d = mpz_get_d(q);
  double p_d = mpz_get_d(p);
  return log(q_d) / (2 * (double)(t - 1) * log(p_d));
}

static void random_mpz_poly(mpz_poly_ptr g, int coeff0, int coeff1, int degree)
{
  for (int i = 0; i <= degree; i++) {
    mpz_poly_setcoeff_si(g, i, rand() % (coeff1 - coeff0) + coeff0);
  }
}

void gen_poly_special_q(mpz_poly_ptr f0, mpz_poly_ptr f1, mpz_poly_ptr g,
    mpz_ptr a, mpz_ptr b, mpz_ptr c, mpz_srcptr p, mpz_poly_srcptr h,
    int coeff0, int coeff1, mpz_srcptr q, unsigned int t)
{
  ASSERT(g->alloc > 0);
  srand(time(NULL));

  double epsilon = find_epsilon(q, p, t);
  mpz_set_d(c, pow(mpz_get_d(p), 0.5 - epsilon) + 1);

  mpz_t coeff;
  mpz_init(coeff);
  factor_t factor;
  for (int i = coeff0; i < coeff1; i++) {
    mpz_mul(coeff, c, h->coeff[h->deg]);
    if (i < 0) {
      mpz_add_ui(coeff, coeff, (unsigned int)(-i));
    } else if (i > 0) {
      mpz_add_ui(coeff, coeff, (unsigned int)(i));
    }
    gmp_brute_force_factorize(factor, coeff);
    factor_fprintf(stdout, factor);
  }
  mpz_clear(coeff);

  while (1) {
    random_mpz_poly(g, coeff0, coeff1, h->deg);
    mpz_poly_set(f0, h);
    mpz_poly_mul_mpz(f0, f0, c);
    mpz_poly_add(f0, f0, g);
    if (mpz_poly_is_irreducible(f0, p)) {
#ifdef EEA_BOUND 
      extended_euclidean_algorithm_bound(a, b, p, c, c);
#else
      extended_euclidean_algorithm_stop_2(a, b, p, c);
#endif
      mpz_poly_mul_mpz(f1, g, b);
      mpz_poly_t tmp;
      mpz_poly_init(tmp, -1);
      mpz_poly_mul_mpz(tmp, h, a);
      mpz_poly_add(f1, f1, tmp);
      mpz_poly_clear(tmp);
      break;
    }
  }
}

void function_special_q(mpz_poly_ptr f0, mpz_poly_ptr f1, mpz_poly_ptr g,
    mpz_ptr a, mpz_ptr b, mpz_ptr c, mpz_srcptr p, mpz_poly_srcptr h,
    int coeff0, int coeff1, mpz_srcptr q, unsigned int t)
{
  gen_poly_special_q(f0, f1, g, a, b, c, p, h, coeff0, coeff1, q, t);
  mpz_t p2;
  mpz_init(p2);
  mpz_fdiv_q_ui(p2, p, 2);
  mpz_t tmp;
  mpz_init(tmp);
  for (int i = 0; i <= f1->deg; i++) {
    mpz_poly_getcoeff(tmp, i, f1);
    if (mpz_cmp(tmp, p2) > 0) {
      mpz_sub(tmp, tmp, p);
      mpz_poly_setcoeff(f1, i, tmp);
    }
  }
  mpz_clear(p2);
  mpz_clear(tmp);
}
