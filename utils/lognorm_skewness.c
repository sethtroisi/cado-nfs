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
