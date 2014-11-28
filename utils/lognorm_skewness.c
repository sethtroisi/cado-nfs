/* Code for inner product and norm computation for the 
     the canonical L2 norm (*_L2 functions)
     the circular L2 norm (*_cir_L2 functions
     the square L2 norm (*_sq_L2 functions)
   It exists in 3 versions: mpz_poly, mpz_double, mpz_vector.
   The code for mpz_vector just call mpz_poly with the right parameters
*/


#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <gmp.h>
#include "portability.h"
#include "params.h"
#include "mpz_poly.h"
#include "mpz_vector.h"
#include "double_poly.h"


/*******************************************************************/
/* Functions for mpz_poly and mpz_vector: inner products and norms
   are return as a pair (res_num, res_den), the real value is the
   fraction res_num/res_den. If mpz_t corresponding to res_num or
   res_den is set to NULL, the value is not computed.              */
/*******************************************************************/


/* Compute inner product for the classical L2 norm. */
void
mpz_poly_innerproduct_L2 (mpz_t res_num, mpz_t res_den, mpz_poly_srcptr f,
                          mpz_poly_srcptr g)
{
  FATAL_ERROR_CHECK (f->deg != g->deg, "f->deg should be equal to g->deg");

  mpz_t *a = f->coeff;
  mpz_t *b = g->coeff;

  if (res_num != NULL) // We want the numerator
  {
    mpz_mul (res_num, a[0], b[0]);
    for (int i = 1; i <= f->deg; i++)
      mpz_addmul (res_num, a[i], b[i]);
  }
  if (res_den != NULL) // We want the denominator
  {
    mpz_set_ui (res_den, 1);
  }
}

/* Compute the square of the norm for the classical L2 norm. */
void
mpz_poly_norm_L2 (mpz_t res_num, mpz_t res_den, mpz_poly_srcptr f)
{
  mpz_poly_innerproduct_L2 (res_num, res_den, f, f);
}

/* Compute inner product for the classical L2 norm, with skewness. */
void
mpz_poly_skew_innerproduct_L2 (mpz_t res_num, mpz_t res_den,
                               mpz_poly_srcptr f, mpz_poly_srcptr g,
                               mpz_t skew)
{
  FATAL_ERROR_CHECK (f->deg != g->deg, "f->deg should be equal to g->deg");

  mpz_t *a = f->coeff;
  mpz_t *b = g->coeff;

  mpz_t skew2, s, tmp;
  mpz_init (skew2);
  mpz_init_set_ui (s, 1);
  mpz_init (tmp);

  mpz_mul (skew2, skew, skew);

  if (res_num != NULL) // We want the numerator
  {
    mpz_mul (res_num, a[0], b[0]);
    for (int i = 1; i <= f->deg; i++)
    {
      mpz_mul (s, s, skew2);            /* s = skew^(2*i) */
      mpz_mul (tmp, a[i], b[i]);
      mpz_addmul (res_num, tmp, s);
    }
  }
  if (res_den != NULL) // We want the denominator
  {
    mpz_pow_ui (res_den, skew, f->deg);
  }

  mpz_clear (tmp);
  mpz_clear (skew2);
  mpz_clear (s);
}

/* Compute the square of the norm for the classical L2 norm, with skewness. */
void
mpz_poly_skew_norm_L2 (mpz_t res_num, mpz_t res_den, mpz_poly_srcptr f, mpz_t skew)
{
  mpz_poly_skew_innerproduct_L2 (res_num, res_den, f, f, skew);
}

/* Compute inner product for the classical L2 norm. */
void
mpz_vector_innerproduct_L2 (mpz_t res_num, mpz_t res_den, mpz_vector_srcptr u,
                            mpz_vector_srcptr v)
{
  mpz_poly_t f,g;
  f->coeff = u->c;
  f->deg = u->dim-1;
  g->coeff = v->c;
  g->deg = v->dim-1;
  mpz_poly_innerproduct_L2 (res_num, res_den, f, g);
}

/* Compute the square of the norm for the classical L2 norm. */
void
mpz_vector_norm_L2 (mpz_t res_num, mpz_t res_den, mpz_vector_srcptr u)
{
  mpz_poly_t f;
  f->coeff = u->c;
  f->deg = u->dim-1;
  mpz_poly_norm_L2 (res_num, res_den, f);
}

/* Compute inner product for the classical L2 norm, with skewness. */
void
mpz_vector_skew_innerproduct_L2 (mpz_t res_num, mpz_t res_den, mpz_vector_srcptr u,
                            mpz_vector_srcptr v ,mpz_t skew)
{
  mpz_poly_t f,g;
  f->coeff = u->c;
  f->deg = u->dim-1;
  g->coeff = v->c;
  g->deg = v->dim-1;
  mpz_poly_skew_innerproduct_L2 (res_num, res_den, f, g, skew);
}

/* Compute the square of the norm for the classical L2 norm, with skewness. */
void
mpz_vector_skew_norm_L2 (mpz_t res_num, mpz_t res_den, mpz_vector_srcptr u,
                            mpz_t skew)
{
  mpz_poly_t f;
  f->coeff = u->c;
  f->deg = u->dim-1;
  mpz_poly_skew_norm_L2 (res_num, res_den, f, skew);
}


/* Compute inner product for the circular L2 norm. */
void
mpz_poly_innerproduct_cir_L2 (mpz_t res_num, mpz_t res_den, mpz_poly_srcptr f,
                              mpz_poly_srcptr g)
{
  FATAL_ERROR_CHECK (f->deg != g->deg, "f->deg should be equal to g->deg");

  mpz_t *a = f->coeff;
  mpz_t *b = g->coeff;

  mpz_t tmp;
  mpz_init (tmp);

  if (f->deg == 0)
  {
    if (res_num != NULL) // We want the numerator
    {
      mpz_mul (res_num, a[0], b[0]);
    }
    if (res_den != NULL) // We want the denominator
    {
       mpz_set_ui (res_den, 1);
    }
  }
  else if (f->deg == 1)
  {
    if (res_num != NULL) // We want the numerator
    {
      mpz_mul (res_num, a[0], b[0]);

      mpz_addmul (res_num, a[1], b[1]);
    }
    if (res_den != NULL) // We want the denominator
    {
       mpz_set_ui (res_den, 4);
    }
  }
  else if (f->deg == 2)
  {
    if (res_num != NULL) // We want the numerator
    {
      mpz_mul_ui (res_num, a[0], 3);
      mpz_mul (res_num, res_num, b[0]);

      mpz_addmul (res_num, a[0], b[2]);
      mpz_addmul (res_num, a[1], b[1]);
      mpz_addmul (res_num, a[2], b[0]);

      mpz_mul_ui (tmp, a[2], 3);
      mpz_addmul (res_num, tmp, b[2]);
    }
    if (res_den != NULL) // We want the denominator
    {
       mpz_set_ui (res_den, 24);
    }
  }
  else if (f->deg == 3)
  {
    if (res_num != NULL) // We want the numerator
    {
      mpz_mul_ui (res_num, a[0], 5);
      mpz_mul (res_num, res_num, b[0]);

      mpz_addmul (res_num, a[0], b[2]);
      mpz_addmul (res_num, a[1], b[1]);
      mpz_addmul (res_num, a[2], b[0]);

      mpz_addmul (res_num, a[1], b[3]);
      mpz_addmul (res_num, a[2], b[2]);
      mpz_addmul (res_num, a[3], b[1]);

      mpz_mul_ui (tmp, a[3], 5);
      mpz_addmul (res_num, tmp, b[3]);
    }
    if (res_den != NULL) // We want the denominator
    {
       mpz_set_ui (res_den, 64);
    }
  }
  else if (f->deg == 4)
  {
    if (res_num != NULL) // We want the numerator
    {
      mpz_mul_ui (res_num, a[0], 35);
      mpz_mul (res_num, res_num, b[0]);

      mpz_mul (tmp, a[0], b[2]);
      mpz_addmul (tmp, a[1], b[1]);
      mpz_addmul (tmp, a[2], b[0]);
      mpz_addmul_ui (res_num, tmp, 5);

      mpz_mul (tmp, a[0], b[4]);
      mpz_addmul (tmp, a[1], b[3]);
      mpz_addmul (tmp, a[2], b[2]);
      mpz_addmul (tmp, a[3], b[1]);
      mpz_addmul (tmp, a[4], b[0]);
      mpz_addmul_ui (res_num, tmp, 3);

      mpz_mul (tmp, a[2], b[4]);
      mpz_addmul (tmp, a[3], b[3]);
      mpz_addmul (tmp, a[4], b[2]);
      mpz_addmul_ui (res_num, tmp, 5);

      mpz_mul_ui (tmp, a[4], 35);
      mpz_addmul (res_num, tmp, b[4]);
    }
    if (res_den != NULL) // We want the denominator
    {
       mpz_set_ui (res_den, 640);
    }
  }
  else if (f->deg == 5)
  {
    if (res_num != NULL) // We want the numerator
    {
      mpz_mul_ui (res_num, a[0], 63);
      mpz_mul (res_num, res_num, b[0]);

      mpz_mul (tmp, a[0], b[2]);
      mpz_addmul (tmp, a[1], b[1]);
      mpz_addmul (tmp, a[2], b[0]);
      mpz_addmul_ui (res_num, tmp, 7);

      mpz_mul (tmp, a[0], b[4]);
      mpz_addmul (tmp, a[1], b[3]);
      mpz_addmul (tmp, a[2], b[2]);
      mpz_addmul (tmp, a[3], b[1]);
      mpz_addmul (tmp, a[4], b[0]);
      mpz_addmul_ui (res_num, tmp, 3);

      mpz_mul (tmp, a[1], b[5]);
      mpz_addmul (tmp, a[2], b[4]);
      mpz_addmul (tmp, a[3], b[3]);
      mpz_addmul (tmp, a[4], b[2]);
      mpz_addmul (tmp, a[5], b[1]);
      mpz_addmul_ui (res_num, tmp, 3);

      mpz_mul (tmp, a[3], b[5]);
      mpz_addmul (tmp, a[4], b[4]);
      mpz_addmul (tmp, a[5], b[3]);
      mpz_addmul_ui (res_num, tmp, 7);

      mpz_mul_ui (tmp, a[5], 63);
      mpz_addmul (res_num, tmp, b[5]);
    }
    if (res_den != NULL) // We want the denominator
    {
       mpz_set_ui (res_den, 1536);
    }
  }
  else if (f->deg == 6)
  {
    if (res_num != NULL) // We want the numerator
    {
      mpz_mul_ui (res_num, a[0], 231);
      mpz_mul (res_num, res_num, b[0]);

      mpz_mul (tmp, a[0], b[2]);
      mpz_addmul (tmp, a[1], b[1]);
      mpz_addmul (tmp, a[2], b[0]);
      mpz_addmul_ui (res_num, tmp, 21);

      mpz_mul (tmp, a[0], b[4]);
      mpz_addmul (tmp, a[1], b[3]);
      mpz_addmul (tmp, a[2], b[2]);
      mpz_addmul (tmp, a[3], b[1]);
      mpz_addmul (tmp, a[4], b[0]);
      mpz_addmul_ui (res_num, tmp, 7);

      mpz_mul (tmp, a[0], b[6]);
      mpz_addmul (tmp, a[1], b[5]);
      mpz_addmul (tmp, a[2], b[4]);
      mpz_addmul (tmp, a[3], b[3]);
      mpz_addmul (tmp, a[4], b[2]);
      mpz_addmul (tmp, a[5], b[1]);
      mpz_addmul (tmp, a[6], b[0]);
      mpz_addmul_ui (res_num, tmp, 5);

      mpz_mul (tmp, a[2], b[6]);
      mpz_addmul (tmp, a[3], b[5]);
      mpz_addmul (tmp, a[4], b[4]);
      mpz_addmul (tmp, a[5], b[3]);
      mpz_addmul (tmp, a[6], b[2]);
      mpz_addmul_ui (res_num, tmp, 7);

      mpz_mul (tmp, a[4], b[6]);
      mpz_addmul (tmp, a[5], b[5]);
      mpz_addmul (tmp, a[6], b[4]);
      mpz_addmul_ui (res_num, tmp, 21);

      mpz_mul_ui (tmp, a[6], 231);
      mpz_addmul (res_num, tmp, b[6]);
    }
    if (res_den != NULL) // We want the denominator
    {
       mpz_set_ui (res_den, 7168);
    }
  }
  else if (f->deg == 7)
  {
    if (res_num != NULL) // We want the numerator
    {
      mpz_mul_ui (res_num, a[0], 429);
      mpz_mul (res_num, res_num, b[0]);

      mpz_mul (tmp, a[0], b[2]);
      mpz_addmul (tmp, a[1], b[1]);
      mpz_addmul (tmp, a[2], b[0]);
      mpz_addmul_ui (res_num, tmp, 33);

      mpz_mul (tmp, a[0], b[4]);
      mpz_addmul (tmp, a[1], b[3]);
      mpz_addmul (tmp, a[2], b[2]);
      mpz_addmul (tmp, a[3], b[1]);
      mpz_addmul (tmp, a[4], b[0]);
      mpz_addmul_ui (res_num, tmp, 9);

      mpz_mul (tmp, a[0], b[6]);
      mpz_addmul (tmp, a[1], b[5]);
      mpz_addmul (tmp, a[2], b[4]);
      mpz_addmul (tmp, a[3], b[3]);
      mpz_addmul (tmp, a[4], b[2]);
      mpz_addmul (tmp, a[5], b[1]);
      mpz_addmul (tmp, a[6], b[0]);
      mpz_addmul_ui (res_num, tmp, 5);

      mpz_mul (tmp, a[1], b[7]);
      mpz_addmul (tmp, a[2], b[6]);
      mpz_addmul (tmp, a[3], b[5]);
      mpz_addmul (tmp, a[4], b[4]);
      mpz_addmul (tmp, a[5], b[3]);
      mpz_addmul (tmp, a[6], b[2]);
      mpz_addmul (tmp, a[7], b[1]);
      mpz_addmul_ui (res_num, tmp, 5);

      mpz_mul (tmp, a[3], b[7]);
      mpz_addmul (tmp, a[4], b[6]);
      mpz_addmul (tmp, a[5], b[5]);
      mpz_addmul (tmp, a[6], b[4]);
      mpz_addmul (tmp, a[7], b[3]);
      mpz_addmul_ui (res_num, tmp, 9);

      mpz_mul (tmp, a[5], b[7]);
      mpz_addmul (tmp, a[6], b[6]);
      mpz_addmul (tmp, a[7], b[5]);
      mpz_addmul_ui (res_num, tmp, 33);

      mpz_mul_ui (tmp, a[7], 429);
      mpz_addmul (res_num, tmp, b[7]);
    }
    if (res_den != NULL) // We want the denominator
    {
       mpz_set_ui (res_den, 16384);
    }
  }
  else
  {
    fprintf (stderr, "Not implemented yet, deg =  %d\n", f->deg);
    abort();
  }
  mpz_clear (tmp);
}

/* Compute the square of the norm for the circular L2 norm. */
void
mpz_poly_norm_cir_L2 (mpz_t res_num, mpz_t res_den, mpz_poly_srcptr f)
{
  mpz_poly_innerproduct_cir_L2 (res_num, res_den, f, f);
}

/* Compute inner product for the circular L2 norm, with skewness. */
void
mpz_poly_skew_innerproduct_cir_L2 (mpz_t res_num, mpz_t res_den,
                                   mpz_poly_srcptr f, mpz_poly_srcptr g,
                                   mpz_t skew)
{
  FATAL_ERROR_CHECK (f->deg != g->deg, "f->deg should be equal to g->deg");

  mpz_t *a = f->coeff;
  mpz_t *b = g->coeff;

  mpz_t skew2, s, tmp;
  mpz_init (skew2);
  mpz_init_set_ui (s, 1);
  mpz_init (tmp);

  mpz_mul (skew2, skew, skew);

  if (f->deg == 0)
  {
    if (res_num != NULL) // We want the numerator
    {
      /* s = 1 = skew^0 */
      mpz_mul (res_num, a[0], b[0]);
    }
    if (res_den != NULL) // We want the denominator
    {
      mpz_pow_ui (res_den, skew, f->deg);
      mpz_mul_ui (res_den, res_den, 1);
    }
  }
  else if (f->deg == 1)
  {
    if (res_num != NULL) // We want the numerator
    {
      /* s = 1 = skew^0 */
      mpz_mul (res_num, a[0], b[0]);

      mpz_mul (s, s, skew2);        /* s = skew^2 */
      mpz_mul (tmp, a[1], b[1]);
      mpz_addmul (res_num, tmp, s);
    }
    if (res_den != NULL) // We want the denominator
    {
      mpz_pow_ui (res_den, skew, f->deg);
      mpz_mul_ui (res_den, res_den, 4);
    }
  }
  else if (f->deg == 2)
  {
    if (res_num != NULL) // We want the numerator
    {
      /* s = 1 = skew^0 */
      mpz_mul_ui (res_num, a[0], 3);
      mpz_mul (res_num, res_num, b[0]);

      mpz_mul (s, s, skew2);        /* s = skew^2 */
      mpz_mul (tmp, a[0], b[2]);
      mpz_addmul (tmp, a[1], b[1]);
      mpz_addmul (tmp, a[2], b[0]);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^4 */
      mpz_mul (tmp, a[2], b[2]);
      mpz_mul_ui (tmp, tmp, 3);
      mpz_addmul (res_num, tmp, s);
    }
    if (res_den != NULL) // We want the denominator
    {
      mpz_pow_ui (res_den, skew, f->deg);
      mpz_mul_ui (res_den, res_den, 24);
    }
  }
  else if (f->deg == 3)
  {
    if (res_num != NULL) // We want the numerator
    {
      /* s = 1 = skew^0 */
      mpz_mul_ui (res_num, a[0], 5);
      mpz_mul (res_num, res_num, b[0]);

      mpz_mul (s, s, skew2);        /* s = skew^2 */
      mpz_mul (tmp, a[0], b[2]);
      mpz_addmul (tmp, a[1], b[1]);
      mpz_addmul (tmp, a[2], b[0]);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^4 */
      mpz_mul (tmp, a[1], b[3]);
      mpz_addmul (tmp, a[2], b[2]);
      mpz_addmul (tmp, a[3], b[1]);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^6 */
      mpz_mul (tmp, a[3], b[3]);
      mpz_mul_ui (tmp, tmp, 5);
      mpz_addmul (res_num, tmp, s);
    }
    if (res_den != NULL) // We want the denominator
    {
      mpz_pow_ui (res_den, skew, f->deg);
      mpz_mul_ui (res_den, res_den, 64);
    }
  }
  else if (f->deg == 4)
  {
    if (res_num != NULL) // We want the numerator
    {
      /* s = 1 = skew^0 */
      mpz_mul_ui (res_num, a[0], 35);
      mpz_mul (res_num, res_num, b[0]);

      mpz_mul (s, s, skew2);        /* s = skew^2 */
      mpz_mul (tmp, a[0], b[2]);
      mpz_addmul (tmp, a[1], b[1]);
      mpz_addmul (tmp, a[2], b[0]);
      mpz_mul_ui (tmp, tmp, 5);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^4 */
      mpz_mul (tmp, a[0], b[4]);
      mpz_addmul (tmp, a[1], b[3]);
      mpz_addmul (tmp, a[2], b[2]);
      mpz_addmul (tmp, a[3], b[1]);
      mpz_addmul (tmp, a[4], b[0]);
      mpz_mul_ui (tmp, tmp, 3);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^6 */
      mpz_mul (tmp, a[2], b[4]);
      mpz_addmul (tmp, a[3], b[3]);
      mpz_addmul (tmp, a[4], b[2]);
      mpz_mul_ui (tmp, tmp, 5);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^8 */
      mpz_mul (tmp, a[4], b[4]);
      mpz_mul_ui (tmp, tmp, 35);
      mpz_addmul (res_num, tmp, s);
    }
    if (res_den != NULL) // We want the denominator
    {
      mpz_pow_ui (res_den, skew, f->deg);
      mpz_mul_ui (res_den, res_den, 640);
    }
  }
  else if (f->deg == 5)
  {
    if (res_num != NULL) // We want the numerator
    {
      /* s = 1 = skew^0 */
      mpz_mul_ui (res_num, a[0], 63);
      mpz_mul (res_num, res_num, b[0]);

      mpz_mul (s, s, skew2);        /* s = skew^2 */
      mpz_mul (tmp, a[0], b[2]);
      mpz_addmul (tmp, a[1], b[1]);
      mpz_addmul (tmp, a[2], b[0]);
      mpz_mul_ui (tmp, tmp, 7);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^4 */
      mpz_mul (tmp, a[0], b[4]);
      mpz_addmul (tmp, a[1], b[3]);
      mpz_addmul (tmp, a[2], b[2]);
      mpz_addmul (tmp, a[3], b[1]);
      mpz_addmul (tmp, a[4], b[0]);
      mpz_mul_ui (tmp, tmp, 3);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^6 */
      mpz_mul (tmp, a[1], b[5]);
      mpz_addmul (tmp, a[2], b[4]);
      mpz_addmul (tmp, a[3], b[3]);
      mpz_addmul (tmp, a[4], b[2]);
      mpz_addmul (tmp, a[5], b[1]);
      mpz_mul_ui (tmp, tmp, 3);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^8 */
      mpz_mul (tmp, a[3], b[5]);
      mpz_addmul (tmp, a[4], b[4]);
      mpz_addmul (tmp, a[5], b[3]);
      mpz_mul_ui (tmp, tmp, 7);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^10 */
      mpz_mul (tmp, a[5], b[5]);
      mpz_mul_ui (tmp, tmp, 63);
      mpz_addmul (res_num, tmp, s);
    }
    if (res_den != NULL) // We want the denominator
    {
      mpz_pow_ui (res_den, skew, f->deg);
      mpz_mul_ui (res_den, res_den, 1536);
    }
  }
  else if (f->deg == 6)
  {
    if (res_num != NULL) // We want the numerator
    {
      /* s = 1 = skew^0 */
      mpz_mul_ui (res_num, a[0], 231);
      mpz_mul (res_num, res_num, b[0]);

      mpz_mul (s, s, skew2);        /* s = skew^2 */
      mpz_mul (tmp, a[0], b[2]);
      mpz_addmul (tmp, a[1], b[1]);
      mpz_addmul (tmp, a[2], b[0]);
      mpz_mul_ui (tmp, tmp, 21);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^4 */
      mpz_mul (tmp, a[0], b[4]);
      mpz_addmul (tmp, a[1], b[3]);
      mpz_addmul (tmp, a[2], b[2]);
      mpz_addmul (tmp, a[3], b[1]);
      mpz_addmul (tmp, a[4], b[0]);
      mpz_mul_ui (tmp, tmp, 7);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^6 */
      mpz_mul (tmp, a[0], b[6]);
      mpz_addmul (tmp, a[1], b[5]);
      mpz_addmul (tmp, a[2], b[4]);
      mpz_addmul (tmp, a[3], b[3]);
      mpz_addmul (tmp, a[4], b[2]);
      mpz_addmul (tmp, a[5], b[1]);
      mpz_addmul (tmp, a[6], b[0]);
      mpz_mul_ui (tmp, tmp, 5);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^8 */
      mpz_mul (tmp, a[2], b[6]);
      mpz_addmul (tmp, a[3], b[5]);
      mpz_addmul (tmp, a[4], b[4]);
      mpz_addmul (tmp, a[5], b[3]);
      mpz_addmul (tmp, a[6], b[2]);
      mpz_mul_ui (tmp, tmp, 7);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^10 */
      mpz_mul (tmp, a[4], b[6]);
      mpz_addmul (tmp, a[5], b[5]);
      mpz_addmul (tmp, a[6], b[4]);
      mpz_mul_ui (tmp, tmp, 21);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^12 */
      mpz_mul (tmp, a[6], b[6]);
      mpz_mul_ui (tmp, tmp, 231);
      mpz_addmul (res_num, tmp, s);
    }
    if (res_den != NULL) // We want the denominator
    {
      mpz_pow_ui (res_den, skew, f->deg);
      mpz_mul_ui (res_den, res_den, 7168);
    }
  }
  else if (f->deg == 7)
  {
    if (res_num != NULL) // We want the numerator
    {
      /* s = 1 = skew^0 */
      mpz_mul_ui (res_num, a[0], 429);
      mpz_mul (res_num, res_num, b[0]);

      mpz_mul (s, s, skew2);        /* s = skew^2 */
      mpz_mul (tmp, a[0], b[2]);
      mpz_addmul (tmp, a[1], b[1]);
      mpz_addmul (tmp, a[2], b[0]);
      mpz_mul_ui (tmp, tmp, 33);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^4 */
      mpz_mul (tmp, a[0], b[4]);
      mpz_addmul (tmp, a[1], b[3]);
      mpz_addmul (tmp, a[2], b[2]);
      mpz_addmul (tmp, a[3], b[1]);
      mpz_addmul (tmp, a[4], b[0]);
      mpz_mul_ui (tmp, tmp, 9);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^6 */
      mpz_mul (tmp, a[0], b[6]);
      mpz_addmul (tmp, a[1], b[5]);
      mpz_addmul (tmp, a[2], b[4]);
      mpz_addmul (tmp, a[3], b[3]);
      mpz_addmul (tmp, a[4], b[2]);
      mpz_addmul (tmp, a[5], b[1]);
      mpz_addmul (tmp, a[6], b[0]);
      mpz_mul_ui (tmp, tmp, 5);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^8 */
      mpz_mul (tmp, a[1], b[7]);
      mpz_addmul (tmp, a[2], b[6]);
      mpz_addmul (tmp, a[3], b[5]);
      mpz_addmul (tmp, a[4], b[4]);
      mpz_addmul (tmp, a[5], b[3]);
      mpz_addmul (tmp, a[6], b[2]);
      mpz_addmul (tmp, a[7], b[1]);
      mpz_mul_ui (tmp, tmp, 5);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^10 */
      mpz_mul (tmp, a[3], b[7]);
      mpz_addmul (tmp, a[4], b[6]);
      mpz_addmul (tmp, a[5], b[5]);
      mpz_addmul (tmp, a[6], b[4]);
      mpz_addmul (tmp, a[7], b[3]);
      mpz_mul_ui (tmp, tmp, 9);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^12 */
      mpz_mul (tmp, a[5], b[7]);
      mpz_addmul (tmp, a[6], b[6]);
      mpz_addmul (tmp, a[7], b[5]);
      mpz_mul_ui (tmp, tmp, 33);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^14 */
      mpz_mul (tmp, a[7], b[7]);
      mpz_mul_ui (tmp, tmp, 429);
      mpz_addmul (res_num, tmp, s);
    }
    if (res_den != NULL) // We want the denominator
    {
      mpz_pow_ui (res_den, skew, f->deg);
      mpz_mul_ui (res_den, res_den, 16384);
    }
  }
  else
  {
    fprintf (stderr, "Not implemented yet, deg =  %d\n", f->deg);
    abort();
  }

  mpz_clear (tmp);
  mpz_clear (skew2);
  mpz_clear (s);
}

/* Compute the square of the norm for the circular L2 norm, with skewness. */
void
mpz_poly_skew_norm_cir_L2 (mpz_t res_num, mpz_t res_den, mpz_poly_srcptr f, mpz_t skew)
{
  mpz_poly_skew_innerproduct_cir_L2 (res_num, res_den, f, f, skew);
}

/* Compute inner product for the circular L2 norm. */
void
mpz_vector_innerproduct_cir_L2 (mpz_t res_num, mpz_t res_den, mpz_vector_srcptr u,
                                mpz_vector_srcptr v)
{
  mpz_poly_t f,g;
  f->coeff = u->c;
  f->deg = u->dim-1;
  g->coeff = v->c;
  g->deg = v->dim-1;
  mpz_poly_innerproduct_cir_L2 (res_num, res_den, f, g);
}

/* Compute the square of the norm for the circular L2 norm. */
void
mpz_vector_norm_cir_L2 (mpz_t res_num, mpz_t res_den, mpz_vector_srcptr u)
{
  mpz_poly_t f;
  f->coeff = u->c;
  f->deg = u->dim-1;
  mpz_poly_norm_cir_L2 (res_num, res_den, f);
}

/* Compute inner product for the circular L2 norm, with skewness. */
void
mpz_vector_skew_innerproduct_cir_L2 (mpz_t res_num, mpz_t res_den, mpz_vector_srcptr u,
                                mpz_vector_srcptr v ,mpz_t skew)
{
  mpz_poly_t f,g;
  f->coeff = u->c;
  f->deg = u->dim-1;
  g->coeff = v->c;
  g->deg = v->dim-1;
  mpz_poly_skew_innerproduct_cir_L2 (res_num, res_den, f, g, skew);
}

/* Compute the square of the norm for the circular L2 norm, with skewness. */
void
mpz_vector_skew_norm_cir_L2 (mpz_t res_num, mpz_t res_den, mpz_vector_srcptr u,
                                mpz_t skew)
{
  mpz_poly_t f;
  f->coeff = u->c;
  f->deg = u->dim-1;
  mpz_poly_skew_norm_cir_L2 (res_num, res_den, f, skew);
}


/* Compute inner product for the square L2 norm. */
void
mpz_poly_innerproduct_sq_L2 (mpz_t res_num, mpz_t res_den, mpz_poly_srcptr f,
                             mpz_poly_srcptr g)
{
  FATAL_ERROR_CHECK (f->deg != g->deg, "f->deg should be equal to g->deg");

  mpz_t *a = f->coeff;
  mpz_t *b = g->coeff;

  mpz_t tmp;
  mpz_init (tmp);

  if (f->deg == 0)
  {
    if (res_num != NULL) // We want the numerator
    {
      mpz_mul (res_num, a[0], b[0]);
    }
    if (res_den != NULL) // We want the denominator
    {
       mpz_set_ui (res_den, 1);
    }
  }
  else if (f->deg == 1)
  {
    if (res_num != NULL) // We want the numerator
    {
      mpz_mul (res_num, a[0], b[0]);

      mpz_addmul (res_num, a[1], b[1]);
    }
    if (res_den != NULL) // We want the denominator
    {
       mpz_set_ui (res_den, 3);
    }
  }
  else if (f->deg == 2)
  {
    if (res_num != NULL) // We want the numerator
    {
      mpz_mul_ui (res_num, a[0], 9);
      mpz_mul (res_num, res_num, b[0]);

      mpz_mul (tmp, a[0], b[2]);
      mpz_addmul (tmp, a[1], b[1]);
      mpz_addmul (tmp, a[2], b[0]);
      mpz_addmul_ui (res_num, tmp, 5);

      mpz_mul_ui (tmp, a[2], 9);
      mpz_addmul (res_num, tmp, b[2]);
    }
    if (res_den != NULL) // We want the denominator
    {
       mpz_set_ui (res_den, 45);
    }
  }
  else if (f->deg == 3)
  {
    if (res_num != NULL) // We want the numerator
    {
      mpz_mul_ui (res_num, a[0], 15);
      mpz_mul (res_num, res_num, b[0]);

      mpz_mul (tmp, a[0], b[2]);
      mpz_addmul (tmp, a[1], b[1]);
      mpz_addmul (tmp, a[2], b[0]);
      mpz_addmul_ui (res_num, tmp, 7);

      mpz_mul (tmp, a[1], b[3]);
      mpz_addmul (tmp, a[2], b[2]);
      mpz_addmul (tmp, a[3], b[1]);
      mpz_addmul_ui (res_num, tmp, 7);

      mpz_mul_ui (tmp, a[3], 15);
      mpz_addmul (res_num, tmp, b[3]);
    }
    if (res_den != NULL) // We want the denominator
    {
       mpz_set_ui (res_den, 105);
    }
  }
  else if (f->deg == 4)
  {
    if (res_num != NULL) // We want the numerator
    {
      mpz_mul_ui (res_num, a[0], 175);
      mpz_mul (res_num, res_num, b[0]);

      mpz_mul (tmp, a[0], b[2]);
      mpz_addmul (tmp, a[1], b[1]);
      mpz_addmul (tmp, a[2], b[0]);
      mpz_addmul_ui (res_num, tmp, 75);

      mpz_mul (tmp, a[0], b[4]);
      mpz_addmul (tmp, a[1], b[3]);
      mpz_addmul (tmp, a[2], b[2]);
      mpz_addmul (tmp, a[3], b[1]);
      mpz_addmul (tmp, a[4], b[0]);
      mpz_addmul_ui (res_num, tmp, 63);

      mpz_mul (tmp, a[2], b[4]);
      mpz_addmul (tmp, a[3], b[3]);
      mpz_addmul (tmp, a[4], b[2]);
      mpz_addmul_ui (res_num, tmp, 75);

      mpz_mul_ui (tmp, a[4], 175);
      mpz_addmul (res_num, tmp, b[4]);
    }
    if (res_den != NULL) // We want the denominator
    {
       mpz_set_ui (res_den, 1575);
    }
  }
  else if (f->deg == 5)
  {
    if (res_num != NULL) // We want the numerator
    {
      mpz_mul_ui (res_num, a[0], 945);
      mpz_mul (res_num, res_num, b[0]);

      mpz_mul (tmp, a[0], b[2]);
      mpz_addmul (tmp, a[1], b[1]);
      mpz_addmul (tmp, a[2], b[0]);
      mpz_addmul_ui (res_num, tmp, 385);

      mpz_mul (tmp, a[0], b[4]);
      mpz_addmul (tmp, a[1], b[3]);
      mpz_addmul (tmp, a[2], b[2]);
      mpz_addmul (tmp, a[3], b[1]);
      mpz_addmul (tmp, a[4], b[0]);
      mpz_addmul_ui (res_num, tmp, 297);

      mpz_mul (tmp, a[1], b[5]);
      mpz_addmul (tmp, a[2], b[4]);
      mpz_addmul (tmp, a[3], b[3]);
      mpz_addmul (tmp, a[4], b[2]);
      mpz_addmul (tmp, a[5], b[1]);
      mpz_addmul_ui (res_num, tmp, 297);

      mpz_mul (tmp, a[3], b[5]);
      mpz_addmul (tmp, a[4], b[4]);
      mpz_addmul (tmp, a[5], b[3]);
      mpz_addmul_ui (res_num, tmp, 385);

      mpz_mul_ui (tmp, a[5], 945);
      mpz_addmul (res_num, tmp, b[5]);
    }
    if (res_den != NULL) // We want the denominator
    {
       mpz_set_ui (res_den, 10395);
    }
  }
  else if (f->deg == 6)
  {
    if (res_num != NULL) // We want the numerator
    {
      mpz_mul_ui (res_num, a[0], 24255);
      mpz_mul (res_num, res_num, b[0]);

      mpz_mul (tmp, a[0], b[2]);
      mpz_addmul (tmp, a[1], b[1]);
      mpz_addmul (tmp, a[2], b[0]);
      mpz_addmul_ui (res_num, tmp, 9555);

      mpz_mul (tmp, a[0], b[4]);
      mpz_addmul (tmp, a[1], b[3]);
      mpz_addmul (tmp, a[2], b[2]);
      mpz_addmul (tmp, a[3], b[1]);
      mpz_addmul (tmp, a[4], b[0]);
      mpz_addmul_ui (res_num, tmp, 7007);

      mpz_mul (tmp, a[0], b[6]);
      mpz_addmul (tmp, a[1], b[5]);
      mpz_addmul (tmp, a[2], b[4]);
      mpz_addmul (tmp, a[3], b[3]);
      mpz_addmul (tmp, a[4], b[2]);
      mpz_addmul (tmp, a[5], b[1]);
      mpz_addmul (tmp, a[6], b[0]);
      mpz_addmul_ui (res_num, tmp, 6435);

      mpz_mul (tmp, a[2], b[6]);
      mpz_addmul (tmp, a[3], b[5]);
      mpz_addmul (tmp, a[4], b[4]);
      mpz_addmul (tmp, a[5], b[3]);
      mpz_addmul (tmp, a[6], b[2]);
      mpz_addmul_ui (res_num, tmp, 7007);

      mpz_mul (tmp, a[4], b[6]);
      mpz_addmul (tmp, a[5], b[5]);
      mpz_addmul (tmp, a[6], b[4]);
      mpz_addmul_ui (res_num, tmp, 9555);

      mpz_mul_ui (tmp, a[6], 24255);
      mpz_addmul (res_num, tmp, b[6]);
    }
    if (res_den != NULL) // We want the denominator
    {
       mpz_set_ui (res_den, 315315);
    }
  }
  else if (f->deg == 7)
  {
    if (res_num != NULL) // We want the numerator
    {
      mpz_mul_ui (res_num, a[0], 3003);
      mpz_mul (res_num, res_num, b[0]);

      mpz_mul (tmp, a[0], b[2]);
      mpz_addmul (tmp, a[1], b[1]);
      mpz_addmul (tmp, a[2], b[0]);
      mpz_addmul_ui (res_num, tmp, 1155);

      mpz_mul (tmp, a[0], b[4]);
      mpz_addmul (tmp, a[1], b[3]);
      mpz_addmul (tmp, a[2], b[2]);
      mpz_addmul (tmp, a[3], b[1]);
      mpz_addmul (tmp, a[4], b[0]);
      mpz_addmul_ui (res_num, tmp, 819);

      mpz_mul (tmp, a[0], b[6]);
      mpz_addmul (tmp, a[1], b[5]);
      mpz_addmul (tmp, a[2], b[4]);
      mpz_addmul (tmp, a[3], b[3]);
      mpz_addmul (tmp, a[4], b[2]);
      mpz_addmul (tmp, a[5], b[1]);
      mpz_addmul (tmp, a[6], b[0]);
      mpz_addmul_ui (res_num, tmp, 715);

      mpz_mul (tmp, a[1], b[7]);
      mpz_addmul (tmp, a[2], b[6]);
      mpz_addmul (tmp, a[3], b[5]);
      mpz_addmul (tmp, a[4], b[4]);
      mpz_addmul (tmp, a[5], b[3]);
      mpz_addmul (tmp, a[6], b[2]);
      mpz_addmul (tmp, a[7], b[1]);
      mpz_addmul_ui (res_num, tmp, 715);

      mpz_mul (tmp, a[3], b[7]);
      mpz_addmul (tmp, a[4], b[6]);
      mpz_addmul (tmp, a[5], b[5]);
      mpz_addmul (tmp, a[6], b[4]);
      mpz_addmul (tmp, a[7], b[3]);
      mpz_addmul_ui (res_num, tmp, 819);

      mpz_mul (tmp, a[5], b[7]);
      mpz_addmul (tmp, a[6], b[6]);
      mpz_addmul (tmp, a[7], b[5]);
      mpz_addmul_ui (res_num, tmp, 1155);

      mpz_mul_ui (tmp, a[7], 3003);
      mpz_addmul (res_num, tmp, b[7]);
    }
    if (res_den != NULL) // We want the denominator
    {
       mpz_set_ui (res_den, 45045);
    }
  }
  else
  {
    fprintf (stderr, "Not implemented yet, deg =  %d\n", f->deg);
    abort();
  }
  mpz_clear (tmp);
}

/* Compute the square of the norm for the square L2 norm. */
void
mpz_poly_norm_sq_L2 (mpz_t res_num, mpz_t res_den, mpz_poly_srcptr f)
{
  mpz_poly_innerproduct_sq_L2 (res_num, res_den, f, f);
}

/* Compute inner product for the square L2 norm, with skewness. */
void
mpz_poly_skew_innerproduct_sq_L2 (mpz_t res_num, mpz_t res_den,
                                  mpz_poly_srcptr f, mpz_poly_srcptr g,
                                  mpz_t skew)
{
  FATAL_ERROR_CHECK (f->deg != g->deg, "f->deg should be equal to g->deg");

  mpz_t *a = f->coeff;
  mpz_t *b = g->coeff;

  mpz_t skew2, s, tmp;
  mpz_init (skew2);
  mpz_init_set_ui (s, 1);
  mpz_init (tmp);

  mpz_mul (skew2, skew, skew);

  if (f->deg == 0)
  {
    if (res_num != NULL) // We want the numerator
    {
      /* s = 1 = skew^0 */
      mpz_mul (res_num, a[0], b[0]);
    }
    if (res_den != NULL) // We want the denominator
    {
      mpz_pow_ui (res_den, skew, f->deg);
      mpz_mul_ui (res_den, res_den, 1);
    }
  }
  else if (f->deg == 1)
  {
    if (res_num != NULL) // We want the numerator
    {
      /* s = 1 = skew^0 */
      mpz_mul (res_num, a[0], b[0]);

      mpz_mul (s, s, skew2);        /* s = skew^2 */
      mpz_mul (tmp, a[1], b[1]);
      mpz_addmul (res_num, tmp, s);
    }
    if (res_den != NULL) // We want the denominator
    {
      mpz_pow_ui (res_den, skew, f->deg);
      mpz_mul_ui (res_den, res_den, 3);
    }
  }
  else if (f->deg == 2)
  {
    if (res_num != NULL) // We want the numerator
    {
      /* s = 1 = skew^0 */
      mpz_mul_ui (res_num, a[0], 9);
      mpz_mul (res_num, res_num, b[0]);

      mpz_mul (s, s, skew2);        /* s = skew^2 */
      mpz_mul (tmp, a[0], b[2]);
      mpz_addmul (tmp, a[1], b[1]);
      mpz_addmul (tmp, a[2], b[0]);
      mpz_mul_ui (tmp, tmp, 5);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^4 */
      mpz_mul (tmp, a[2], b[2]);
      mpz_mul_ui (tmp, tmp, 9);
      mpz_addmul (res_num, tmp, s);
    }
    if (res_den != NULL) // We want the denominator
    {
      mpz_pow_ui (res_den, skew, f->deg);
      mpz_mul_ui (res_den, res_den, 45);
    }
  }
  else if (f->deg == 3)
  {
    if (res_num != NULL) // We want the numerator
    {
      /* s = 1 = skew^0 */
      mpz_mul_ui (res_num, a[0], 15);
      mpz_mul (res_num, res_num, b[0]);

      mpz_mul (s, s, skew2);        /* s = skew^2 */
      mpz_mul (tmp, a[0], b[2]);
      mpz_addmul (tmp, a[1], b[1]);
      mpz_addmul (tmp, a[2], b[0]);
      mpz_mul_ui (tmp, tmp, 7);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^4 */
      mpz_mul (tmp, a[1], b[3]);
      mpz_addmul (tmp, a[2], b[2]);
      mpz_addmul (tmp, a[3], b[1]);
      mpz_mul_ui (tmp, tmp, 7);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^6 */
      mpz_mul (tmp, a[3], b[3]);
      mpz_mul_ui (tmp, tmp, 15);
      mpz_addmul (res_num, tmp, s);
    }
    if (res_den != NULL) // We want the denominator
    {
      mpz_pow_ui (res_den, skew, f->deg);
      mpz_mul_ui (res_den, res_den, 105);
    }
  }
  else if (f->deg == 4)
  {
    if (res_num != NULL) // We want the numerator
    {
      /* s = 1 = skew^0 */
      mpz_mul_ui (res_num, a[0], 175);
      mpz_mul (res_num, res_num, b[0]);

      mpz_mul (s, s, skew2);        /* s = skew^2 */
      mpz_mul (tmp, a[0], b[2]);
      mpz_addmul (tmp, a[1], b[1]);
      mpz_addmul (tmp, a[2], b[0]);
      mpz_mul_ui (tmp, tmp, 75);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^4 */
      mpz_mul (tmp, a[0], b[4]);
      mpz_addmul (tmp, a[1], b[3]);
      mpz_addmul (tmp, a[2], b[2]);
      mpz_addmul (tmp, a[3], b[1]);
      mpz_addmul (tmp, a[4], b[0]);
      mpz_mul_ui (tmp, tmp, 63);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^6 */
      mpz_mul (tmp, a[2], b[4]);
      mpz_addmul (tmp, a[3], b[3]);
      mpz_addmul (tmp, a[4], b[2]);
      mpz_mul_ui (tmp, tmp, 75);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^8 */
      mpz_mul (tmp, a[4], b[4]);
      mpz_mul_ui (tmp, tmp, 175);
      mpz_addmul (res_num, tmp, s);
    }
    if (res_den != NULL) // We want the denominator
    {
      mpz_pow_ui (res_den, skew, f->deg);
      mpz_mul_ui (res_den, res_den, 1575);
    }
  }
  else if (f->deg == 5)
  {
    if (res_num != NULL) // We want the numerator
    {
      /* s = 1 = skew^0 */
      mpz_mul_ui (res_num, a[0], 945);
      mpz_mul (res_num, res_num, b[0]);

      mpz_mul (s, s, skew2);        /* s = skew^2 */
      mpz_mul (tmp, a[0], b[2]);
      mpz_addmul (tmp, a[1], b[1]);
      mpz_addmul (tmp, a[2], b[0]);
      mpz_mul_ui (tmp, tmp, 385);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^4 */
      mpz_mul (tmp, a[0], b[4]);
      mpz_addmul (tmp, a[1], b[3]);
      mpz_addmul (tmp, a[2], b[2]);
      mpz_addmul (tmp, a[3], b[1]);
      mpz_addmul (tmp, a[4], b[0]);
      mpz_mul_ui (tmp, tmp, 297);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^6 */
      mpz_mul (tmp, a[1], b[5]);
      mpz_addmul (tmp, a[2], b[4]);
      mpz_addmul (tmp, a[3], b[3]);
      mpz_addmul (tmp, a[4], b[2]);
      mpz_addmul (tmp, a[5], b[1]);
      mpz_mul_ui (tmp, tmp, 297);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^8 */
      mpz_mul (tmp, a[3], b[5]);
      mpz_addmul (tmp, a[4], b[4]);
      mpz_addmul (tmp, a[5], b[3]);
      mpz_mul_ui (tmp, tmp, 385);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^10 */
      mpz_mul (tmp, a[5], b[5]);
      mpz_mul_ui (tmp, tmp, 945);
      mpz_addmul (res_num, tmp, s);
    }
    if (res_den != NULL) // We want the denominator
    {
      mpz_pow_ui (res_den, skew, f->deg);
      mpz_mul_ui (res_den, res_den, 10395);
    }
  }
  else if (f->deg == 6)
  {
    if (res_num != NULL) // We want the numerator
    {
      /* s = 1 = skew^0 */
      mpz_mul_ui (res_num, a[0], 24255);
      mpz_mul (res_num, res_num, b[0]);

      mpz_mul (s, s, skew2);        /* s = skew^2 */
      mpz_mul (tmp, a[0], b[2]);
      mpz_addmul (tmp, a[1], b[1]);
      mpz_addmul (tmp, a[2], b[0]);
      mpz_mul_ui (tmp, tmp, 9555);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^4 */
      mpz_mul (tmp, a[0], b[4]);
      mpz_addmul (tmp, a[1], b[3]);
      mpz_addmul (tmp, a[2], b[2]);
      mpz_addmul (tmp, a[3], b[1]);
      mpz_addmul (tmp, a[4], b[0]);
      mpz_mul_ui (tmp, tmp, 7007);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^6 */
      mpz_mul (tmp, a[0], b[6]);
      mpz_addmul (tmp, a[1], b[5]);
      mpz_addmul (tmp, a[2], b[4]);
      mpz_addmul (tmp, a[3], b[3]);
      mpz_addmul (tmp, a[4], b[2]);
      mpz_addmul (tmp, a[5], b[1]);
      mpz_addmul (tmp, a[6], b[0]);
      mpz_mul_ui (tmp, tmp, 6435);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^8 */
      mpz_mul (tmp, a[2], b[6]);
      mpz_addmul (tmp, a[3], b[5]);
      mpz_addmul (tmp, a[4], b[4]);
      mpz_addmul (tmp, a[5], b[3]);
      mpz_addmul (tmp, a[6], b[2]);
      mpz_mul_ui (tmp, tmp, 7007);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^10 */
      mpz_mul (tmp, a[4], b[6]);
      mpz_addmul (tmp, a[5], b[5]);
      mpz_addmul (tmp, a[6], b[4]);
      mpz_mul_ui (tmp, tmp, 9555);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^12 */
      mpz_mul (tmp, a[6], b[6]);
      mpz_mul_ui (tmp, tmp, 24255);
      mpz_addmul (res_num, tmp, s);
    }
    if (res_den != NULL) // We want the denominator
    {
      mpz_pow_ui (res_den, skew, f->deg);
      mpz_mul_ui (res_den, res_den, 315315);
    }
  }
  else if (f->deg == 7)
  {
    if (res_num != NULL) // We want the numerator
    {
      /* s = 1 = skew^0 */
      mpz_mul_ui (res_num, a[0], 3003);
      mpz_mul (res_num, res_num, b[0]);

      mpz_mul (s, s, skew2);        /* s = skew^2 */
      mpz_mul (tmp, a[0], b[2]);
      mpz_addmul (tmp, a[1], b[1]);
      mpz_addmul (tmp, a[2], b[0]);
      mpz_mul_ui (tmp, tmp, 1155);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^4 */
      mpz_mul (tmp, a[0], b[4]);
      mpz_addmul (tmp, a[1], b[3]);
      mpz_addmul (tmp, a[2], b[2]);
      mpz_addmul (tmp, a[3], b[1]);
      mpz_addmul (tmp, a[4], b[0]);
      mpz_mul_ui (tmp, tmp, 819);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^6 */
      mpz_mul (tmp, a[0], b[6]);
      mpz_addmul (tmp, a[1], b[5]);
      mpz_addmul (tmp, a[2], b[4]);
      mpz_addmul (tmp, a[3], b[3]);
      mpz_addmul (tmp, a[4], b[2]);
      mpz_addmul (tmp, a[5], b[1]);
      mpz_addmul (tmp, a[6], b[0]);
      mpz_mul_ui (tmp, tmp, 715);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^8 */
      mpz_mul (tmp, a[1], b[7]);
      mpz_addmul (tmp, a[2], b[6]);
      mpz_addmul (tmp, a[3], b[5]);
      mpz_addmul (tmp, a[4], b[4]);
      mpz_addmul (tmp, a[5], b[3]);
      mpz_addmul (tmp, a[6], b[2]);
      mpz_addmul (tmp, a[7], b[1]);
      mpz_mul_ui (tmp, tmp, 715);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^10 */
      mpz_mul (tmp, a[3], b[7]);
      mpz_addmul (tmp, a[4], b[6]);
      mpz_addmul (tmp, a[5], b[5]);
      mpz_addmul (tmp, a[6], b[4]);
      mpz_addmul (tmp, a[7], b[3]);
      mpz_mul_ui (tmp, tmp, 819);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^12 */
      mpz_mul (tmp, a[5], b[7]);
      mpz_addmul (tmp, a[6], b[6]);
      mpz_addmul (tmp, a[7], b[5]);
      mpz_mul_ui (tmp, tmp, 1155);
      mpz_addmul (res_num, tmp, s);

      mpz_mul (s, s, skew2);        /* s = skew^14 */
      mpz_mul (tmp, a[7], b[7]);
      mpz_mul_ui (tmp, tmp, 3003);
      mpz_addmul (res_num, tmp, s);
    }
    if (res_den != NULL) // We want the denominator
    {
      mpz_pow_ui (res_den, skew, f->deg);
      mpz_mul_ui (res_den, res_den, 45045);
    }
  }
  else
  {
    fprintf (stderr, "Not implemented yet, deg =  %d\n", f->deg);
    abort();
  }

  mpz_clear (tmp);
  mpz_clear (skew2);
  mpz_clear (s);
}

/* Compute the square of the norm for the square L2 norm, with skewness. */
void
mpz_poly_skew_norm_sq_L2 (mpz_t res_num, mpz_t res_den, mpz_poly_srcptr f, mpz_t skew)
{
  mpz_poly_skew_innerproduct_sq_L2 (res_num, res_den, f, f, skew);
}

/* Compute inner product for the square L2 norm. */
void
mpz_vector_innerproduct_sq_L2 (mpz_t res_num, mpz_t res_den, mpz_vector_srcptr u,
                               mpz_vector_srcptr v)
{
  mpz_poly_t f,g;
  f->coeff = u->c;
  f->deg = u->dim-1;
  g->coeff = v->c;
  g->deg = v->dim-1;
  mpz_poly_innerproduct_sq_L2 (res_num, res_den, f, g);
}

/* Compute the square of the norm for the square L2 norm. */
void
mpz_vector_norm_sq_L2 (mpz_t res_num, mpz_t res_den, mpz_vector_srcptr u)
{
  mpz_poly_t f;
  f->coeff = u->c;
  f->deg = u->dim-1;
  mpz_poly_norm_sq_L2 (res_num, res_den, f);
}

/* Compute inner product for the square L2 norm, with skewness. */
void
mpz_vector_skew_innerproduct_sq_L2 (mpz_t res_num, mpz_t res_den, mpz_vector_srcptr u,
                               mpz_vector_srcptr v ,mpz_t skew)
{
  mpz_poly_t f,g;
  f->coeff = u->c;
  f->deg = u->dim-1;
  g->coeff = v->c;
  g->deg = v->dim-1;
  mpz_poly_skew_innerproduct_sq_L2 (res_num, res_den, f, g, skew);
}

/* Compute the square of the norm for the square L2 norm, with skewness. */
void
mpz_vector_skew_norm_sq_L2 (mpz_t res_num, mpz_t res_den, mpz_vector_srcptr u,
                               mpz_t skew)
{
  mpz_poly_t f;
  f->coeff = u->c;
  f->deg = u->dim-1;
  mpz_poly_skew_norm_sq_L2 (res_num, res_den, f, skew);
}


/*******************************************************************/
/* Functions for double poly.                                      */
/*******************************************************************/


/* Compute inner product for the classical L2 norm. */
/* Assume deg(f) == deg(g) */
double
double_poly_innerproduct_L2 (double_poly_srcptr f, double_poly_srcptr g)
{
  FATAL_ERROR_CHECK (f->deg != g->deg, "f->deg should be equal to g->deg");
  double *a = f->coeff;
  double *b = g->coeff;
  if (f->deg == 0)
  {
    return ( a[0]*b[0] );
  }
  else if (f->deg == 1)
  {
    return ( a[1]*b[1] + a[0]*b[0] );
  }
  else if (f->deg == 2)
  {
    return ( a[2]*b[2] + a[1]*b[1] + a[0]*b[0] );
  }
  else if (f->deg == 3)
  {
    return ( a[3]*b[3] + a[2]*b[2] + a[1]*b[1] + a[0]*b[0] );
  }
  else if (f->deg == 4)
  {
    return ( a[4]*b[4] + a[3]*b[3] + a[2]*b[2] + a[1]*b[1] + a[0]*b[0] );
  }
  else if (f->deg == 5)
  {
    return ( a[5]*b[5] + a[4]*b[4] + a[3]*b[3] + a[2]*b[2] + a[1]*b[1] + a[0]*b[0] );
  }
  else if (f->deg == 6)
  {
    return ( a[6]*b[6] + a[5]*b[5] + a[4]*b[4] + a[3]*b[3] + a[2]*b[2] + a[1]*b[1] + a[0]*b[0] );
  }
  else if (f->deg == 7)
  {
    return ( a[7]*b[7] + a[6]*b[6] + a[5]*b[5] + a[4]*b[4] + a[3]*b[3] + a[2]*b[2] + a[1]*b[1] + a[0]*b[0] );
  }
  else
  {
    fprintf (stderr, "Not implemented yet, deg = %u\n", f->deg);
    abort();
  }
}

/* Compute the square of the norm for the classical L2 norm. */
/* Assume deg(f) == deg(g) */
double
double_poly_norm_L2 (double_poly_srcptr f)
{
  return double_poly_innerproduct_L2 (f, f);
}

/* Compute inner product for the classical L2 norm, with skewness */
/* Assume deg(f) == deg(g) */
double
double_poly_skew_innerproduct_L2 (double_poly_srcptr f, double_poly_srcptr g, double skew)
{
  FATAL_ERROR_CHECK (f->deg != g->deg, "f->deg should be equal to g->deg");
  double *a = f->coeff;
  double *b = g->coeff;
  double skew2 = skew*skew;
  double ret, s = 1;
  ret = a[0]*b[0];
  for (int i = 1; i <= f->deg; i++)
  {
    s = s*skew2;             /* s = skew^(2*i) */
    ret += s*a[i]*b[i];
  }
  if (f->deg == 1)
    ret = ret / (skew);
  if (f->deg == 2)
    ret = ret / (skew*skew);
  if (f->deg == 3)
    ret = ret / (skew*skew*skew);
  if (f->deg == 4)
    ret = ret / (skew*skew*skew*skew);
  if (f->deg == 5)
    ret = ret / (skew*skew*skew*skew*skew);
  if (f->deg == 6)
    ret = ret / (skew*skew*skew*skew*skew*skew);
  if (f->deg == 7)
    ret = ret / (skew*skew*skew*skew*skew*skew*skew);

  return ret;
}

/* Return the square of the skew-norm of the polynomial f. */
double
double_poly_skew_norm_L2 (double_poly_srcptr f, double skew)
{
  return double_poly_skew_innerproduct_L2 (f, f, skew);
}

/* Return the log of the skew-norm of the polynomial f. */
double
double_poly_lognorm_L2 (double_poly_srcptr f, double skew)
{
  double n = double_poly_skew_norm_L2 (f, skew);
  return 0.5 * log (n);
}


#define PI 3.14159265358979323846

/* Compute inner product for the circular L2 norm. */
/* Assume deg(f) == deg(g) */
double
double_poly_innerproduct_cir_L2 (double_poly_srcptr f, double_poly_srcptr g)
{
  FATAL_ERROR_CHECK (f->deg != g->deg, "f->deg should be equal to g->deg");
  double *a = f->coeff;
  double *b = g->coeff;
  if (f->deg == 0)
  {
    return ( a[0]*b[0] );
  }
  else if (f->deg == 1)
  {
    return ( a[1]*b[1] + a[0]*b[0] )/4.0;
  }
  else if (f->deg == 2)
  {
    return ( 3.0*a[2]*b[2] + a[1]*b[1] + 2.0*a[0]*b[2] + 3.0*a[0]*b[0] )/24.0;
  }
  else if (f->deg == 3)
  {
    return ( 5.0*a[3]*b[3] + a[2]*b[2] + 2.0*a[1]*b[3] + a[1]*b[1] + 2.0*a[0]*b[2] + 5.0*a[0]*b[0] )/64.0;
  }
  else if (f->deg == 4)
  {
    return ( 35.0*a[4]*b[4] + 5.0*a[3]*b[3] + 10.0*a[2]*b[4] + 3.0*a[2]*b[2] + 6.0*a[1]*b[3] + 5.0*a[1]*b[1] + 6.0*a[0]*b[4] + 10.0*a[0]*b[2] + 35.0*a[0]*b[0] )/640.0;
  }
  else if (f->deg == 5)
  {
    return ( 63.0*a[5]*b[5] + 7.0*a[4]*b[4] + 14.0*a[3]*b[5] + 3.0*a[3]*b[3] + 6.0*a[2]*b[4] + 3.0*a[2]*b[2] + 6.0*a[1]*b[5] + 6.0*a[1]*b[3] + 7.0*a[1]*b[1] + 6.0*a[0]*b[4] + 14.0*a[0]*b[2] + 63.0*a[0]*b[0] )/1536.0;
  }
  else if (f->deg == 6)
  {
    return ( 231.0*a[6]*b[6] + 21.0*a[5]*b[5] + 42.0*a[4]*b[6] + 7.0*a[4]*b[4] + 14.0*a[3]*b[5] + 5.0*a[3]*b[3] + 14.0*a[2]*b[6] + 10.0*a[2]*b[4] + 7.0*a[2]*b[2] + 10.0*a[1]*b[5] + 14.0*a[1]*b[3] + 21.0*a[1]*b[1] + 10.0*a[0]*b[6] + 14.0*a[0]*b[4] + 42.0*a[0]*b[2] + 231.0*a[0]*b[0] )/7168.0;
  }
  else if (f->deg == 7)
  {
    return ( 429.0*a[7]*b[7] + 33.0*a[6]*b[6] + 66.0*a[5]*b[7] + 9.0*a[5]*b[5] + 18.0*a[4]*b[6] + 5.0*a[4]*b[4] + 18.0*a[3]*b[7] + 10.0*a[3]*b[5] + 5.0*a[3]*b[3] + 10.0*a[2]*b[6] + 10.0*a[2]*b[4] + 9.0*a[2]*b[2] + 10.0*a[1]*b[7] + 10.0*a[1]*b[5] + 18.0*a[1]*b[3] + 33.0*a[1]*b[1] + 10.0*a[0]*b[6] + 18.0*a[0]*b[4] + 66.0*a[0]*b[2] + 429.0*a[0]*b[0] )/16384.0;
  }
  else
  {
    fprintf (stderr, "Not implemented yet, deg = %u\n", f->deg);
    abort();
  }
}

/* Compute the square of the norm for the circular L2 norm. */
/* Assume deg(f) == deg(g) */
double
double_poly_norm_cir_L2 (double_poly_srcptr f)
{
  return double_poly_innerproduct_cir_L2 (f, f);
}

/* Compute inner product for the circular L2 norm, with skewness */
/* Assume deg(f) == deg(g) */
double
double_poly_skew_innerproduct_cir_L2 (double_poly_srcptr f, double_poly_srcptr g, double skew)
{
  FATAL_ERROR_CHECK (f->deg != g->deg, "f->deg should be equal to g->deg");
  double *a = f->coeff;
  double *b = g->coeff;
  double skew2 = skew*skew;
  double ret, s = 1;
  if (f->deg == 0)
  {
    /* s = 1 = skew^0 */
    ret = 1.0*a[0]*b[0];

  }
  else if (f->deg == 1)
  {
    /* s = 1 = skew^0 */
    ret = 1.0*a[0]*b[0];

    s = s*skew2;        /* s = skew^2 */
    ret += 1.0*s*( a[1]*b[1] );

    ret = ret / (4.0*skew);
  }
  else if (f->deg == 2)
  {
    /* s = 1 = skew^0 */
    ret = 3.0*a[0]*b[0];

    s = s*skew2;        /* s = skew^2 */
    ret += 1.0*s*( a[0]*b[2] + a[1]*b[1] + a[2]*b[0] );

    s = s*skew2;        /* s = skew^4 */
    ret += 3.0*s*( a[2]*b[2] );

    ret = ret / (24.0*skew*skew);
  }
  else if (f->deg == 3)
  {
    /* s = 1 = skew^0 */
    ret = 5.0*a[0]*b[0];

    s = s*skew2;        /* s = skew^2 */
    ret += 1.0*s*( a[0]*b[2] + a[1]*b[1] + a[2]*b[0] );

    s = s*skew2;        /* s = skew^4 */
    ret += 1.0*s*( a[1]*b[3] + a[2]*b[2] + a[3]*b[1] );

    s = s*skew2;        /* s = skew^6 */
    ret += 5.0*s*( a[3]*b[3] );

    ret = ret / (64.0*skew*skew*skew);
  }
  else if (f->deg == 4)
  {
    /* s = 1 = skew^0 */
    ret = 35.0*a[0]*b[0];

    s = s*skew2;        /* s = skew^2 */
    ret += 5.0*s*( a[0]*b[2] + a[1]*b[1] + a[2]*b[0] );

    s = s*skew2;        /* s = skew^4 */
    ret += 3.0*s*( a[0]*b[4] + a[1]*b[3] + a[2]*b[2] + a[3]*b[1] + a[4]*b[0] );

    s = s*skew2;        /* s = skew^6 */
    ret += 5.0*s*( a[2]*b[4] + a[3]*b[3] + a[4]*b[2] );

    s = s*skew2;        /* s = skew^8 */
    ret += 35.0*s*( a[4]*b[4] );

    ret = ret / (640.0*skew*skew*skew*skew);
  }
  else if (f->deg == 5)
  {
    /* s = 1 = skew^0 */
    ret = 63.0*a[0]*b[0];

    s = s*skew2;        /* s = skew^2 */
    ret += 7.0*s*( a[0]*b[2] + a[1]*b[1] + a[2]*b[0] );

    s = s*skew2;        /* s = skew^4 */
    ret += 3.0*s*( a[0]*b[4] + a[1]*b[3] + a[2]*b[2] + a[3]*b[1] + a[4]*b[0] );

    s = s*skew2;        /* s = skew^6 */
    ret += 3.0*s*( a[1]*b[5] + a[2]*b[4] + a[3]*b[3] + a[4]*b[2] + a[5]*b[1] );

    s = s*skew2;        /* s = skew^8 */
    ret += 7.0*s*( a[3]*b[5] + a[4]*b[4] + a[5]*b[3] );

    s = s*skew2;        /* s = skew^10 */
    ret += 63.0*s*( a[5]*b[5] );

    ret = ret / (1536.0*skew*skew*skew*skew*skew);
  }
  else if (f->deg == 6)
  {
    /* s = 1 = skew^0 */
    ret = 231.0*a[0]*b[0];

    s = s*skew2;        /* s = skew^2 */
    ret += 21.0*s*( a[0]*b[2] + a[1]*b[1] + a[2]*b[0] );

    s = s*skew2;        /* s = skew^4 */
    ret += 7.0*s*( a[0]*b[4] + a[1]*b[3] + a[2]*b[2] + a[3]*b[1] + a[4]*b[0] );

    s = s*skew2;        /* s = skew^6 */
    ret += 5.0*s*( a[0]*b[6] + a[1]*b[5] + a[2]*b[4] + a[3]*b[3] + a[4]*b[2] + a[5]*b[1] + a[6]*b[0] );

    s = s*skew2;        /* s = skew^8 */
    ret += 7.0*s*( a[2]*b[6] + a[3]*b[5] + a[4]*b[4] + a[5]*b[3] + a[6]*b[2] );

    s = s*skew2;        /* s = skew^10 */
    ret += 21.0*s*( a[4]*b[6] + a[5]*b[5] + a[6]*b[4] );

    s = s*skew2;        /* s = skew^12 */
    ret += 231.0*s*( a[6]*b[6] );

    ret = ret / (7168.0*skew*skew*skew*skew*skew*skew);
  }
  else if (f->deg == 7)
  {
    /* s = 1 = skew^0 */
    ret = 429.0*a[0]*b[0];

    s = s*skew2;        /* s = skew^2 */
    ret += 33.0*s*( a[0]*b[2] + a[1]*b[1] + a[2]*b[0] );

    s = s*skew2;        /* s = skew^4 */
    ret += 9.0*s*( a[0]*b[4] + a[1]*b[3] + a[2]*b[2] + a[3]*b[1] + a[4]*b[0] );

    s = s*skew2;        /* s = skew^6 */
    ret += 5.0*s*( a[0]*b[6] + a[1]*b[5] + a[2]*b[4] + a[3]*b[3] + a[4]*b[2] + a[5]*b[1] + a[6]*b[0] );

    s = s*skew2;        /* s = skew^8 */
    ret += 5.0*s*( a[1]*b[7] + a[2]*b[6] + a[3]*b[5] + a[4]*b[4] + a[5]*b[3] + a[6]*b[2] + a[7]*b[1] );

    s = s*skew2;        /* s = skew^10 */
    ret += 9.0*s*( a[3]*b[7] + a[4]*b[6] + a[5]*b[5] + a[6]*b[4] + a[7]*b[3] );

    s = s*skew2;        /* s = skew^12 */
    ret += 33.0*s*( a[5]*b[7] + a[6]*b[6] + a[7]*b[5] );

    s = s*skew2;        /* s = skew^14 */
    ret += 429.0*s*( a[7]*b[7] );

    ret = ret / (16384.0*skew*skew*skew*skew*skew*skew*skew);
  }
  else
  {
    fprintf (stderr, "Not implemented yet, deg =  %d\n", f->deg);
    abort();
  }

  /* In CADO, there is no 1/pi factor in the definition of the
     circular norm*/
  return PI * ret;
}

/* Return the square of the skew-norm of the polynomial f. */
double
double_poly_skew_norm_cir_L2 (double_poly_srcptr f, double skew)
{
  return double_poly_skew_innerproduct_cir_L2 (f, f, skew);
}

/* Return the log of the skew-norm of the polynomial f. */
double
double_poly_lognorm_cir_L2 (double_poly_srcptr f, double skew)
{
  double n = double_poly_skew_norm_cir_L2 (f, skew);
  return 0.5 * log (n);
}



/* Compute inner product for the square L2 norm. */
/* Assume deg(f) == deg(g) */
double
double_poly_innerproduct_sq_L2 (double_poly_srcptr f, double_poly_srcptr g)
{
  FATAL_ERROR_CHECK (f->deg != g->deg, "f->deg should be equal to g->deg");
  double *a = f->coeff;
  double *b = g->coeff;
  if (f->deg == 0)
  {
    return ( a[0]*b[0] );
  }
  else if (f->deg == 1)
  {
    return ( a[1]*b[1] + a[0]*b[0] )/3.0;
  }
  else if (f->deg == 2)
  {
    return ( 9.0*a[2]*b[2] + 5.0*a[1]*b[1] + 10.0*a[0]*b[2] + 9.0*a[0]*b[0] )/45.0;
  }
  else if (f->deg == 3)
  {
    return ( 15.0*a[3]*b[3] + 7.0*a[2]*b[2] + 14.0*a[1]*b[3] + 7.0*a[1]*b[1] + 14.0*a[0]*b[2] + 15.0*a[0]*b[0] )/105.0;
  }
  else if (f->deg == 4)
  {
    return ( 175.0*a[4]*b[4] + 75.0*a[3]*b[3] + 150.0*a[2]*b[4] + 63.0*a[2]*b[2] + 126.0*a[1]*b[3] + 75.0*a[1]*b[1] + 126.0*a[0]*b[4] + 150.0*a[0]*b[2] + 175.0*a[0]*b[0] )/1575.0;
  }
  else if (f->deg == 5)
  {
    return ( 945.0*a[5]*b[5] + 385.0*a[4]*b[4] + 770.0*a[3]*b[5] + 297.0*a[3]*b[3] + 594.0*a[2]*b[4] + 297.0*a[2]*b[2] + 594.0*a[1]*b[5] + 594.0*a[1]*b[3] + 385.0*a[1]*b[1] + 594.0*a[0]*b[4] + 770.0*a[0]*b[2] + 945.0*a[0]*b[0] )/10395.0;
  }
  else if (f->deg == 6)
  {
    return ( 24255.0*a[6]*b[6] + 9555.0*a[5]*b[5] + 19110.0*a[4]*b[6] + 7007.0*a[4]*b[4] + 14014.0*a[3]*b[5] + 6435.0*a[3]*b[3] + 14014.0*a[2]*b[6] + 12870.0*a[2]*b[4] + 7007.0*a[2]*b[2] + 12870.0*a[1]*b[5] + 14014.0*a[1]*b[3] + 9555.0*a[1]*b[1] + 12870.0*a[0]*b[6] + 14014.0*a[0]*b[4] + 19110.0*a[0]*b[2] + 24255.0*a[0]*b[0] )/315315.0;
  }
  else if (f->deg == 7)
  {
    return ( 3003.0*a[7]*b[7] + 1155.0*a[6]*b[6] + 2310.0*a[5]*b[7] + 819.0*a[5]*b[5] + 1638.0*a[4]*b[6] + 715.0*a[4]*b[4] + 1638.0*a[3]*b[7] + 1430.0*a[3]*b[5] + 715.0*a[3]*b[3] + 1430.0*a[2]*b[6] + 1430.0*a[2]*b[4] + 819.0*a[2]*b[2] + 1430.0*a[1]*b[7] + 1430.0*a[1]*b[5] + 1638.0*a[1]*b[3] + 1155.0*a[1]*b[1] + 1430.0*a[0]*b[6] + 1638.0*a[0]*b[4] + 2310.0*a[0]*b[2] + 3003.0*a[0]*b[0] )/45045.0;
  }
  else
  {
    fprintf (stderr, "Not implemented yet, deg = %u\n", f->deg);
    abort();
  }
}

/* Compute the square of the norm for the square L2 norm. */
/* Assume deg(f) == deg(g) */
double
double_poly_norm_sq_L2 (double_poly_srcptr f)
{
  return double_poly_innerproduct_sq_L2 (f, f);
}

/* Compute inner product for the square L2 norm, with skewness */
/* Assume deg(f) == deg(g) */
double
double_poly_skew_innerproduct_sq_L2 (double_poly_srcptr f, double_poly_srcptr g, double skew)
{
  FATAL_ERROR_CHECK (f->deg != g->deg, "f->deg should be equal to g->deg");
  double *a = f->coeff;
  double *b = g->coeff;
  double skew2 = skew*skew;
  double ret, s = 1;
  if (f->deg == 0)
  {
    /* s = 1 = skew^0 */
    ret = 1.0*a[0]*b[0];

  }
  else if (f->deg == 1)
  {
    /* s = 1 = skew^0 */
    ret = 1.0*a[0]*b[0];

    s = s*skew2;        /* s = skew^2 */
    ret += 1.0*s*( a[1]*b[1] );

    ret = ret / (3.0*skew);
  }
  else if (f->deg == 2)
  {
    /* s = 1 = skew^0 */
    ret = 9.0*a[0]*b[0];

    s = s*skew2;        /* s = skew^2 */
    ret += 5.0*s*( a[0]*b[2] + a[1]*b[1] + a[2]*b[0] );

    s = s*skew2;        /* s = skew^4 */
    ret += 9.0*s*( a[2]*b[2] );

    ret = ret / (45.0*skew*skew);
  }
  else if (f->deg == 3)
  {
    /* s = 1 = skew^0 */
    ret = 15.0*a[0]*b[0];

    s = s*skew2;        /* s = skew^2 */
    ret += 7.0*s*( a[0]*b[2] + a[1]*b[1] + a[2]*b[0] );

    s = s*skew2;        /* s = skew^4 */
    ret += 7.0*s*( a[1]*b[3] + a[2]*b[2] + a[3]*b[1] );

    s = s*skew2;        /* s = skew^6 */
    ret += 15.0*s*( a[3]*b[3] );

    ret = ret / (105.0*skew*skew*skew);
  }
  else if (f->deg == 4)
  {
    /* s = 1 = skew^0 */
    ret = 175.0*a[0]*b[0];

    s = s*skew2;        /* s = skew^2 */
    ret += 75.0*s*( a[0]*b[2] + a[1]*b[1] + a[2]*b[0] );

    s = s*skew2;        /* s = skew^4 */
    ret += 63.0*s*( a[0]*b[4] + a[1]*b[3] + a[2]*b[2] + a[3]*b[1] + a[4]*b[0] );

    s = s*skew2;        /* s = skew^6 */
    ret += 75.0*s*( a[2]*b[4] + a[3]*b[3] + a[4]*b[2] );

    s = s*skew2;        /* s = skew^8 */
    ret += 175.0*s*( a[4]*b[4] );

    ret = ret / (1575.0*skew*skew*skew*skew);
  }
  else if (f->deg == 5)
  {
    /* s = 1 = skew^0 */
    ret = 945.0*a[0]*b[0];

    s = s*skew2;        /* s = skew^2 */
    ret += 385.0*s*( a[0]*b[2] + a[1]*b[1] + a[2]*b[0] );

    s = s*skew2;        /* s = skew^4 */
    ret += 297.0*s*( a[0]*b[4] + a[1]*b[3] + a[2]*b[2] + a[3]*b[1] + a[4]*b[0] );

    s = s*skew2;        /* s = skew^6 */
    ret += 297.0*s*( a[1]*b[5] + a[2]*b[4] + a[3]*b[3] + a[4]*b[2] + a[5]*b[1] );

    s = s*skew2;        /* s = skew^8 */
    ret += 385.0*s*( a[3]*b[5] + a[4]*b[4] + a[5]*b[3] );

    s = s*skew2;        /* s = skew^10 */
    ret += 945.0*s*( a[5]*b[5] );

    ret = ret / (10395.0*skew*skew*skew*skew*skew);
  }
  else if (f->deg == 6)
  {
    /* s = 1 = skew^0 */
    ret = 24255.0*a[0]*b[0];

    s = s*skew2;        /* s = skew^2 */
    ret += 9555.0*s*( a[0]*b[2] + a[1]*b[1] + a[2]*b[0] );

    s = s*skew2;        /* s = skew^4 */
    ret += 7007.0*s*( a[0]*b[4] + a[1]*b[3] + a[2]*b[2] + a[3]*b[1] + a[4]*b[0] );

    s = s*skew2;        /* s = skew^6 */
    ret += 6435.0*s*( a[0]*b[6] + a[1]*b[5] + a[2]*b[4] + a[3]*b[3] + a[4]*b[2] + a[5]*b[1] + a[6]*b[0] );

    s = s*skew2;        /* s = skew^8 */
    ret += 7007.0*s*( a[2]*b[6] + a[3]*b[5] + a[4]*b[4] + a[5]*b[3] + a[6]*b[2] );

    s = s*skew2;        /* s = skew^10 */
    ret += 9555.0*s*( a[4]*b[6] + a[5]*b[5] + a[6]*b[4] );

    s = s*skew2;        /* s = skew^12 */
    ret += 24255.0*s*( a[6]*b[6] );

    ret = ret / (315315.0*skew*skew*skew*skew*skew*skew);
  }
  else if (f->deg == 7)
  {
    /* s = 1 = skew^0 */
    ret = 3003.0*a[0]*b[0];

    s = s*skew2;        /* s = skew^2 */
    ret += 1155.0*s*( a[0]*b[2] + a[1]*b[1] + a[2]*b[0] );

    s = s*skew2;        /* s = skew^4 */
    ret += 819.0*s*( a[0]*b[4] + a[1]*b[3] + a[2]*b[2] + a[3]*b[1] + a[4]*b[0] );

    s = s*skew2;        /* s = skew^6 */
    ret += 715.0*s*( a[0]*b[6] + a[1]*b[5] + a[2]*b[4] + a[3]*b[3] + a[4]*b[2] + a[5]*b[1] + a[6]*b[0] );

    s = s*skew2;        /* s = skew^8 */
    ret += 715.0*s*( a[1]*b[7] + a[2]*b[6] + a[3]*b[5] + a[4]*b[4] + a[5]*b[3] + a[6]*b[2] + a[7]*b[1] );

    s = s*skew2;        /* s = skew^10 */
    ret += 819.0*s*( a[3]*b[7] + a[4]*b[6] + a[5]*b[5] + a[6]*b[4] + a[7]*b[3] );

    s = s*skew2;        /* s = skew^12 */
    ret += 1155.0*s*( a[5]*b[7] + a[6]*b[6] + a[7]*b[5] );

    s = s*skew2;        /* s = skew^14 */
    ret += 3003.0*s*( a[7]*b[7] );

    ret = ret / (45045.0*skew*skew*skew*skew*skew*skew*skew);
  }
  else
  {
    fprintf (stderr, "Not implemented yet, deg =  %d\n", f->deg);
    abort();
  }

  return ret;
}

/* Return the square of the skew-norm of the polynomial f. */
double
double_poly_skew_norm_sq_L2 (double_poly_srcptr f, double skew)
{
  return double_poly_skew_innerproduct_sq_L2 (f, f, skew);
}

/* Return the log of the skew-norm of the polynomial f. */
double
double_poly_lognorm_sq_L2 (double_poly_srcptr f, double skew)
{
  double n = double_poly_skew_norm_sq_L2 (f, skew);
  return 0.5 * log (n);
}

