/**
   These files (mpz_vector.*) implement vectors whose coefficients are
   multiprecision integers (using mpz_t from GNU MP).
   We use them in twoquadractics and twocubics.
*/

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "portability.h"
#include "utils.h"


void
mpz_vector_init (mpz_vector_t v, unsigned int d)
{
  ASSERT_ALWAYS (d > 0);
  v->dim = d;
  v->c = (mpz_t *) malloc (d * sizeof(mpz_t));
  ASSERT_ALWAYS (v->c != NULL);
  for (unsigned int i = 0; i < d; i++)
    mpz_init (v->c[i]);
}

void
mpz_vector_clear(mpz_vector_t v)
{
  for (unsigned int i = 0; i < v->dim; i++)
    mpz_clear (v->c[i]);
  free (v->c);
}

void
mpz_vector_swap (mpz_vector_t v1, mpz_vector_t v2)
{
  ASSERT_ALWAYS (v1->dim == v2->dim);
  mpz_t *tmp = v1->c;
  v1->c = v2->c;
  v2->c = tmp;
}

void
mpz_vector_set (mpz_vector_t v, mpz_vector_t s)
{
  ASSERT_ALWAYS (v->dim == s->dim);
  for (unsigned int i = 0; i < v->dim; i++)
    mpz_set (v->c[i], s->c[i]);
}

void
mpz_vector_setcoordinate (mpz_vector_t v, unsigned int i, mpz_t z)
{
  ASSERT_ALWAYS (i < v->dim);
  mpz_set(v->c[i], z);
}

void
mpz_vector_setcoordinate_ui (mpz_vector_t v, unsigned int i, unsigned int z)
{
  ASSERT_ALWAYS (i < v->dim);
  mpz_set_ui (v->c[i], z);
}

int
mpz_vector_is_coordinate_zero (mpz_vector_t v, unsigned int i)
{
  ASSERT_ALWAYS (i < v->dim);
  return (mpz_sgn (v->c[i]) == 0);
}

void
mpz_vector_get_mpz_poly (mpz_poly_t p, mpz_vector_t v)
{
  for (unsigned int i = 0; i < v->dim; i++)
    mpz_poly_setcoeff (p, i, v->c[i]);
}

void
mpz_vector_dot_product (mpz_t res, mpz_vector_t a, mpz_vector_t b)
{
  ASSERT_ALWAYS (a->dim == b->dim);
  mpz_mul(res, a->c[0], b->c[0]);
  for (unsigned int i = 1; i < a->dim; i++)
    mpz_addmul(res, a->c[i], b->c[i]);
}

void
mpz_vector_norm (mpz_t res, mpz_vector_t a)
{
  mpz_vector_dot_product (res, a, a);
}

void
mpz_vector_skew_dot_product (mpz_t res, mpz_vector_t a, mpz_vector_t b,
                             mpz_t skewness)
{
  ASSERT_ALWAYS (a->dim == b->dim);
  mpz_t tmp, s, s2;
  mpz_init (tmp);
  mpz_init (s2);
  mpz_init_set_ui (s, 1);
  mpz_mul (s2, skewness, skewness);

  mpz_mul(res, a->c[0], b->c[0]);
  for (unsigned int i = 1; i < a->dim; i++)
  {
    mpz_mul(tmp, a->c[i], b->c[i]);
    mpz_mul(s, s, s2);
    mpz_addmul(res, tmp, s);
  }

  mpz_clear(tmp);
  mpz_clear(s2);
  mpz_clear(s);
}

void
mpz_vector_skew_norm (mpz_t res, mpz_vector_t a, mpz_t skewness)
{
  mpz_vector_skew_dot_product (res, a, a, skewness);
}


/* compute r <- r - q*v */
void
mpz_vector_submul (mpz_vector_t r, mpz_t q, mpz_vector_t v)
{
  ASSERT_ALWAYS (r->dim == v->dim);
  for (unsigned int i = 0; i < r->dim; i++)
    mpz_submul(r->c[i], q, v->c[i]);
}

/* Lagrange Algorithm.
   Reduce lattice of rank 2.
   Reduce the lattice spanned by u and v and set a and b to the reduced basis.
   The obtained basis is always the two shortest vector.
   See, Nguyen And Vallee, "The LLL algorithm: Survey and applications" p.41
   If skewness parameter is lesser or equal to 1, then the non-skew version of
   mpz_vector_norm and mpz_vector_dot_product is used.
   */
void
mpz_vector_Lagrange (mpz_vector_t a, mpz_vector_t b,
                     mpz_vector_t u, mpz_vector_t v, mpz_t skewness)
{
  mpz_t norm_a, norm_b, q, r, tmp;
  mpz_init (norm_a);
  mpz_init (norm_b);
  mpz_init (q);
  mpz_init (r);
  mpz_init (tmp);

  mpz_vector_set (a, u);
  mpz_vector_set (b, v);

  int skew = (mpz_cmp_ui (skewness, 1) > 0);

  if (skew)
  {
    mpz_vector_skew_norm (norm_a, a, skewness);
    mpz_vector_skew_norm (norm_b, b, skewness);
  }
  else
  {
    mpz_vector_norm (norm_a, a);
    mpz_vector_norm (norm_b, b);
  }

  if (mpz_cmp (norm_a, norm_b) < 0)
  {
    mpz_vector_swap (a, b);
    mpz_swap (norm_a, norm_b);
  }
  while (mpz_cmp (norm_a, norm_b) > 0)
  {
    if (skew)
      mpz_vector_skew_dot_product (tmp, a, b, skewness);
    else
      mpz_vector_dot_product (tmp, a, b);
    mpz_ndiv_q (q, tmp, norm_b);
    mpz_vector_submul (a, q, b);
    mpz_vector_swap (a, b);
    mpz_set (norm_a, norm_b);
    if (skew)
      mpz_vector_skew_norm (norm_b, b, skewness);
    else
      mpz_vector_norm (norm_b, b);
  }

  mpz_clear (norm_a);
  mpz_clear (norm_b);
  mpz_clear (q);
  mpz_clear (r);
  mpz_clear (tmp);
}


