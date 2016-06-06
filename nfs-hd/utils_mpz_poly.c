#include "utils_mpz_poly.h"

/*
 * Add 1 to f. If the constant term is equal to 1, set this term to 0 and
 *  propagate the addition of 1 to the next coefficient, and so on.
 *
 * f: the polynomial on which the addition is computed, the modifications are
 *  made on f.
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
 * Factorise naively a mpz_poly mod 2.
 */
void mpz_poly_factor2(mpz_poly_factor_list_ptr list, mpz_poly_srcptr f)
{
  //set p to 2, because we work in F2.
  mpz_t p;
  mpz_init(p);
  mpz_set_ui(p, 2);

  //make a copy of f.
  mpz_poly fcopy;
  mpz_poly_init(fcopy, f->deg);
  mpz_poly_set(fcopy, f);

  //reduce all the coefficient mod p.
  mpz_t coeff;
  mpz_init(coeff);
  for (int i = 0; i <= f->deg; i++) {
    mpz_poly_getcoeff(coeff, i, f);
    mpz_mod(coeff, coeff, p);
    mpz_poly_setcoeff(fcopy, i, coeff);
  }

  //Purge list.
  mpz_poly_factor_list_flush(list);

  //If deg(f) in F2 is less than 1, we have the factor.
  if (fcopy->deg < 1) {
    mpz_clear(coeff);
    mpz_poly_clear(fcopy);
    mpz_clear(p);
    ASSERT(list->size == 0);
    return;
  } else if (fcopy->deg == 1) {
    mpz_clear(coeff);
    mpz_poly_factor_list_push(list, fcopy, 1);
    mpz_poly_clear(fcopy);
    mpz_clear(p);
    ASSERT(list->size == 1);
    return;
  }

  if (mpz_poly_is_irreducible(f, p)) {
    //If f is irreducible mod 2, fcopy is the factorisation of f mod 2.
    mpz_poly_factor_list_push(list, fcopy, 1);
  } else {
    //Create the first possible factor.
    mpz_poly tmp;
    mpz_poly_init(tmp, 1);
    mpz_poly_setcoeff_int64(tmp, 1, 1);

    //Enumerate all the possible factor of f mod 2.
    while (tmp->deg <= fcopy->deg) {
      //tmp is a possible factor.
      if (mpz_poly_is_irreducible(tmp, p)) {
        mpz_poly q;
        mpz_poly_init(q, 0);
        mpz_poly r;
        mpz_poly_init(r, 0);
        //Euclidean division of fcopy
        mpz_poly_div_qr(q, r, fcopy, tmp, p);
        //Power of the possible factor.
        unsigned int m = 0;
        //While fcopy is divisible by tmp.
        while (r->deg == -1) {
          //Increase the power of tmp.
          m++;
          mpz_poly_set(fcopy, q);
          if (fcopy->deg == 0 || fcopy->deg == -1) {
            //No other possible factor.
            break;
          }
          mpz_poly_div_qr(q, r, fcopy, tmp, p);
        }
        if (m != 0) {
          //Push tmp^m as a factor of f mod 2.
          mpz_poly_factor_list_push(list, tmp, m);
        }
        mpz_poly_clear(q);
        mpz_poly_clear(r);
      }
      //Go to the next possible polynomial in F2.
      mpz_poly_add_one_in_F2(tmp);
    }
    mpz_poly_clear(tmp);
  }

#ifndef NDBEBUG
  //Verify if the factorisation is good.
  mpz_poly_cleandeg(fcopy, -1);
  for (int i = 0; i <= f->deg; i++) {
    mpz_poly_getcoeff(coeff, i, f);
    mpz_mod(coeff, coeff, p);
    mpz_poly_setcoeff(fcopy, i, coeff);
  }

  mpz_poly fmul;
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

  ASSERT(mpz_poly_cmp(fcopy, fmul) == 0);

  mpz_poly_clear(fmul);
#endif // NDBEBUG

  mpz_clear(coeff);
  mpz_poly_clear(fcopy);
  mpz_clear(p);
}

