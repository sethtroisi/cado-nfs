#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "makefb.h"

/* Come from test/utils/test_mpz_poly.c */
static void mpz_poly_setcoeffs_si_var(mpz_poly_t f, int d, ...)
{
    va_list ap;
    va_start(ap, d);
    mpz_poly_realloc(f, d + 1);
    for(int i = 0 ; i <= d ; i++) {
        mpz_set_si(f->coeff[i], va_arg(ap, int));
    }
    mpz_poly_cleandeg(f, d);
    va_end(ap);
}

int main()
{
  mpz_poly_factor_list l;
  mpz_poly_factor_list_init(l);
  mpz_poly_t f;
  mpz_poly_init(f, 0);
  
  // 0+1*x^1+2*x^2+2*x^3+2*x^4  
  mpz_poly_setcoeffs_si_var(f, 4, 0, 1, 2, 2, 2);
  mpz_poly_factor2(l, f);
  ASSERT_ALWAYS(l->size == 1);
  mpz_t coeff;
  mpz_init(coeff);
  mpz_poly_getcoeff(coeff, 0, l->factors[0]->f);
  ASSERT_ALWAYS(mpz_cmp_ui(coeff, 0) == 0);
  mpz_poly_getcoeff(coeff, 1, l->factors[0]->f);
  ASSERT_ALWAYS(mpz_cmp_ui(coeff, 1) == 0);

  // 1+0*x^1+2*x^2+2*x^3+2*x^4  
  mpz_poly_setcoeffs_si_var(f, 4, 1, 0, 2, 2, 2);
  mpz_poly_factor2(l, f);
  ASSERT_ALWAYS(l->size == 0);

  mpz_poly_setcoeffs_si_var(f, 4, 0, 1, -2, -2, -2);
  mpz_poly_factor2(l, f);
  ASSERT_ALWAYS(l->size == 1);
  mpz_poly_getcoeff(coeff, 0, l->factors[0]->f);
  ASSERT_ALWAYS(mpz_cmp_ui(coeff, 0) == 0);
  mpz_poly_getcoeff(coeff, 1, l->factors[0]->f);
  ASSERT_ALWAYS(mpz_cmp_ui(coeff, 1) == 0);
 
  //x^10 + x^9 + x^8 + 1
  mpz_poly_setcoeffs_si_var(f, 10, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1);
  mpz_poly_factor2(l, f);
  ASSERT_ALWAYS(l->size == 3);

  //x^10 + x^8 + x^7 + x^6 + x^2 + x + 1
  mpz_poly_setcoeffs_si_var(f, 10, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1);
  mpz_poly_factor2(l, f);
  ASSERT_ALWAYS(l->size == 1);

  //Very long, just to test if you can factor a large polynomial.
  //x^50 + x^49 + x^47 + x^46 + x^44 + x^40 + x^39 + x^38 + x^34 + x^33 + x^32
  //+ x^31 + x^28 + x^26 + x^25 + x^24 + x^23 + x^22 + x^20 + x^12 + x^11
  //+ x^7 + x^4 + x
  mpz_poly_setcoeffs_si_var(f, 50, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1,
      1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1);
  mpz_poly_factor2(l, f);
  ASSERT_ALWAYS(l->size == 6);

  mpz_clear(coeff);
  mpz_poly_clear(f);
  mpz_poly_factor_list_clear(l);
  return 0;
}
