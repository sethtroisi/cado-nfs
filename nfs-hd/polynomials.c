#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include "cado.h"
#include "utils.h"
#include "polynomials.h"
#include "utils_cofactorisation.h"
#include "polyselect/auxiliary.h"

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
#else // EEA_BOUND
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
#endif // EEA_BOUND

static double find_epsilon(unsigned int q, mpz_srcptr p, unsigned int t)
{
  double q_d = pow(2.0, (double)q);
  double p_d = mpz_get_d(p);

  return log(q_d) / (2 * (double)(t - 1) * log(p_d));
}

//Generate a random number between [offset, offset + length[
static inline void sub_rand_mpz(mpz_ptr rand_Z, mpz_srcptr length,
    mpz_srcptr offset, gmp_randstate_t state)
{
  mpz_urandomm(rand_Z, state, length);
  mpz_add(rand_Z, rand_Z, offset);
}

static void rand_mpz(mpz_ptr rand_Z, mpz_srcptr min, mpz_srcptr max,
    gmp_randstate_t state)
{
  mpz_t length;
  mpz_init(length);
  mpz_sub(length, max, min);
  sub_rand_mpz(rand_Z, length, min, state);
  mpz_clear(length);
}

static void random_mpz_poly(mpz_poly_ptr g, mpz_srcptr min,
    mpz_srcptr max, int degree, int lc, gmp_randstate_t state)
{
  mpz_t rand_Z;
  mpz_init(rand_Z);
  for (int i = 0; i < degree; i++) {
    rand_mpz(rand_Z, min, max, state);
    mpz_poly_setcoeff(g, i, rand_Z);
  }
  if (lc == 0) {
    while(mpz_cmp_ui(rand_Z, 0) == 0) {
      rand_mpz(rand_Z, min, max, state);
    }
    mpz_poly_setcoeff(g, degree, rand_Z);
  } else {
    ASSERT(lc != 0);

    mpz_poly_setcoeff_si(g, degree, lc);
  }
  mpz_clear(rand_Z);
}

static void random_mpz_poly_constraint(mpz_poly_ptr g, mpz_srcptr min,
    mpz_srcptr max, int degree, gmp_randstate_t state,
    mpz_poly_srcptr h)
{
  ASSERT(degree >= 1);
  ASSERT(degree <= h->deg);

  mpz_t rand_Z;
  mpz_init(rand_Z);
  mpz_t zero;
  mpz_init(zero);

  ASSERT(mpz_cmp_ui(zero, 0) == 0);

  for (int i = 0; i < degree; i++) {
    if (mpz_cmp_ui(h->coeff[i], 0) < 0) {
      rand_mpz(rand_Z, zero, max, state);

      ASSERT(mpz_cmp_ui(rand_Z, 0) >= 0);
    } else {
      rand_mpz(rand_Z, min, zero, state);

      ASSERT(mpz_cmp_ui(rand_Z, 0) < 0);
    }
    mpz_poly_setcoeff(g, i, rand_Z);
  }

  if (mpz_cmp_ui(h->coeff[degree], 0) < 0) {
    rand_mpz(rand_Z, zero, max, state);

    ASSERT(mpz_cmp_ui(rand_Z, 0) >= 0);
  } else {
    rand_mpz(rand_Z, min, zero, state);

    ASSERT(mpz_cmp_ui(rand_Z, 0) < 0);
  }
  while(mpz_cmp_ui(rand_Z, 0) == 0) {
    if (mpz_cmp_ui(h->coeff[degree], 0) < 0) {
      rand_mpz(rand_Z, min, zero, state);
    } else {
      rand_mpz(rand_Z, zero, max, state);
    }
  }
  mpz_poly_setcoeff(g, degree, rand_Z);
  mpz_clear(rand_Z);
}

#ifndef ALPHA
static double double_mpz_poly_infinty_norm(mpz_poly_srcptr a)
{
  ASSERT(a->deg >= 1);

  double max = abs(mpz_get_d(a->coeff[0]));
  for (int i = 1; i <= a->deg; i++) {
    max = MAX(max, abs(mpz_get_d(a->coeff[i])));
  }

  return max;
}
#endif // ALPHA

static double gen_poly_special_q(mpz_poly_ptr f0, mpz_poly_ptr f1,
    mpz_poly_ptr g, mpz_ptr a, mpz_ptr b, mpz_ptr c, mpz_srcptr p,
    mpz_poly_ptr h, int * coeff, unsigned int q, unsigned int t,
    unsigned int nb_times, double * weight, int c_tol, gmp_randstate_t state,
    unsigned int h_set)
{
  ASSERT(g->alloc > 0);

  mpz_t * coeff_Z = (mpz_t *) malloc(sizeof(mpz_t) * 2);
  for (unsigned int i = 0; i < 2; i++) {
    mpz_init(coeff_Z[i]);
    mpz_set_si(coeff_Z[i], coeff[i]);
  }

  mpz_t * c_tol_Z = (mpz_t *) malloc(sizeof(mpz_t) * 2);
  for (unsigned int i = 0; i < 2; i++) {
    mpz_init(c_tol_Z[i]);
  }
  mpz_set_si(c_tol_Z[0], -c_tol);
  mpz_set_si(c_tol_Z[1], c_tol + 1);

  double epsilon = find_epsilon(q, p, t);
  mpz_t c_root;
  mpz_init(c_root);
  mpz_set_d(c_root, pow(mpz_get_d(p), 0.5 - epsilon) + 1.0);

  mpz_t rand_Z;
  mpz_init(rand_Z);
  mpz_t c_tmp;
  mpz_init(c_tmp);

  //TODO: maximize the number of projective roots.
  /*int lc = coeff0;*/
  /*mpz_t smooth;*/
  /*mpz_init(smooth);*/
  /*mpz_t coeff;*/
  /*mpz_init(coeff);*/
  /*for (int i = coeff0; i < coeff1; i++) {*/
    /*factor_t factor;*/
    /*mpz_mul(coeff, c, h->coeff[h->deg]);*/
    /*if (i < 0) {*/
      /*mpz_sub_ui(coeff, coeff, (unsigned int)(-i));*/
    /*} else if (i > 0) {*/
      /*mpz_add_ui(coeff, coeff, (unsigned int)(i));*/
    /*}*/
    /*mpz_abs(coeff, coeff);*/
    /*gmp_brute_force_factorize(factor, coeff);*/
    /*if (mpz_cmp_ui(smooth, 0) == 0) {*/
      /*mpz_set(smooth, factor->factorization[factor->number - 1]);*/
      /*lc = i;*/
    /*} else if (mpz_cmp(smooth, factor->factorization[factor->number - 1]) > 0) {*/
      /*mpz_set(smooth, factor->factorization[factor->number - 1]);*/
      /*lc = i;*/
    /*}*/
    /*factor_clear(factor);*/
  /*}*/
  /*mpz_clear(smooth);*/
  /*mpz_clear(coeff);*/

  double alpha0 = 0.0;
  double alpha1 = 0.0;
  double sum_alpha = DBL_MAX;
  mpz_poly_t f0_tmp;
  mpz_poly_init(f0_tmp, -1);
  mpz_poly_t f1_tmp;
  mpz_poly_init(f1_tmp, -1);

  unsigned int k = 0;
  while (k < nb_times) {
    rand_mpz(rand_Z, c_tol_Z[0], c_tol_Z[1], state);
    mpz_add(c_tmp, c_root, rand_Z);
    if (h_set) {
      random_mpz_poly_constraint(g, coeff_Z[0], coeff_Z[1],
          (rand() % h->deg) + 1, state, h);
    } else {
      random_mpz_poly(h, coeff_Z[0], coeff_Z[1], h->deg, 0, state);
      random_mpz_poly_constraint(g, coeff_Z[0], coeff_Z[1],
          (rand() % h->deg) + 1, state, h);
    }
    mpz_poly_set(f0_tmp, h);
    mpz_poly_mul_mpz(f0_tmp, f0_tmp, c_tmp);
    mpz_poly_add(f0_tmp, f0_tmp, g);
#ifdef ALPHA
    alpha0 = get_alpha(f0_tmp, ALPHA_BOUND);
#else // ALPHA
    alpha0 = log2(double_mpz_poly_infinty_norm(f0_tmp));
#endif // ALPHA
    if (mpz_poly_is_irreducible(f0_tmp, p)) {
      k++;
#ifdef EEA_BOUND
      extended_euclidean_algorithm_bound(a, b, p, c_tmp, c_tmp);
#else // EEA_BOUND
      extended_euclidean_algorithm_stop_2(a, b, p, c_tmp);
#endif // EEA_BOUND

      mpz_poly_mul_mpz(f1_tmp, g, b);
      mpz_poly_t tmp;
      mpz_poly_init(tmp, -1);
      mpz_poly_mul_mpz(tmp, h, a);
      mpz_poly_add(f1_tmp, f1_tmp, tmp);
#ifdef ALPHA
      alpha1 = get_alpha(f1_tmp, ALPHA_BOUND);
#else // ALPHA
      alpha1 = log2(double_mpz_poly_infinty_norm(f1_tmp));
#endif // ALPHA
      if (weight[0] * alpha0 + weight[1] * alpha1 < sum_alpha) {
        mpz_poly_set(f0, f0_tmp);
        mpz_poly_set(f1, f1_tmp);
        sum_alpha = weight[0] * alpha0 + weight[1] * alpha1;
        mpz_set(c, c_tmp);
      }
      mpz_poly_clear(tmp);
    }
  }
  mpz_poly_clear(f0_tmp);
  mpz_poly_clear(f1_tmp);

  mpz_clear(rand_Z);
  mpz_clear(c_tmp);
  mpz_clear(c_root);

  for (unsigned int i = 0; i < 2; i++) {
    mpz_clear(c_tol_Z[i]);
  }
  free(c_tol_Z);

  for (unsigned int i = 0; i < 2; i++) {
    mpz_clear(coeff_Z[i]);
  }
  free(coeff_Z);

  return epsilon;
}

double function_special_q(mpz_poly_ptr f0, mpz_poly_ptr f1, mpz_poly_ptr g,
    mpz_ptr a, mpz_ptr b, mpz_ptr c, mpz_srcptr p, mpz_poly_ptr h,
    int * coeff, unsigned int q, unsigned int t,
    unsigned int nb_times, double * weight, int c_tol, gmp_randstate_t state,
    unsigned int h_set)
{
  double epsilon = gen_poly_special_q(f0, f1, g, a, b, c, p, h, coeff,
      q, t, nb_times, weight, c_tol, state, h_set);
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

  return epsilon;
}

void function_classical(mpz_poly_ptr f0, mpz_poly_ptr f1, mpz_srcptr p,
    unsigned int n, int * coeff, unsigned int nb_times,
    double * weight, gmp_randstate_t state)
{
  mpz_t * coeff_Z = (mpz_t *) malloc(sizeof(mpz_t) * 2);
  for (unsigned int i = 0; i < 2; i++) {
    mpz_init(coeff_Z[i]);
    mpz_set_si(coeff_Z[i], coeff[i]);
  }

  mpz_poly_t f0_tmp;
  mpz_poly_init(f0_tmp, -1);
  mpz_poly_t f1_tmp;
  mpz_poly_init(f1_tmp, -1);

  double sum_alpha = DBL_MAX;
  double alpha0 = 0.0;
  double alpha1 = 0.0;
  unsigned int k = 0;
  while (k < nb_times) {
    random_mpz_poly(f0_tmp, coeff_Z[0], coeff_Z[1], n, 1, state);
    if (mpz_poly_is_irreducible(f0_tmp, p)) {
      k++;
      mpz_poly_set(f1_tmp, f0_tmp);
      if (mpz_cmp_ui(f1_tmp->coeff[0], 0) < 0) {
        mpz_add(f1_tmp->coeff[0], f1_tmp->coeff[0], p);
      } else {
        mpz_sub(f1_tmp->coeff[0], f1_tmp->coeff[0], p);
      }
      alpha0 = get_alpha(f0_tmp, ALPHA_BOUND);
      alpha1 = get_alpha(f1_tmp, ALPHA_BOUND);
      if (weight[0] * alpha0 + weight[1] * alpha1 < sum_alpha) {
        mpz_poly_set(f0, f0_tmp);
        mpz_poly_set(f1, f1_tmp);
        sum_alpha = weight[0] * alpha0 + weight[1] * alpha1;
      }
    }
  }

  mpz_poly_clear(f0_tmp);
  mpz_poly_clear(f1_tmp);

  for (unsigned int i = 0; i < 2; i++) {
    mpz_clear(coeff_Z[i]);
  }
  free(coeff_Z);
}
