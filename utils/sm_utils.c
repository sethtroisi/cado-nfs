#include "cado.h"
#include <stdio.h>
#include "portability.h"
#include "utils.h"

/* Init polynomial rel and set it to a - b*x */
void
mpz_poly_init_set_ab (mpz_poly_ptr rel, int64_t a, uint64_t b)
{
  mpz_poly_init(rel, 1);
  mpz_poly_setcoeff_int64(rel, 0, a);
  mpz_poly_setcoeff_int64(rel, 1, -b);
}

/* compute the SM */
void
compute_sm (mpz_poly_t SM, mpz_poly_t num, const mpz_poly_t F, const mpz_t ell,
            const mpz_t smexp, const mpz_t ell2, const mpz_t invl2)
{
  mpz_poly_power_mod_f_mod_mpz_barrett(SM, num, F, smexp, ell2, invl2);
  mpz_poly_sub_ui(SM, SM, 1);

  for(int j = 0; j <= SM->deg; j++)
  {
    ASSERT_ALWAYS(mpz_divisible_p(SM->coeff[j], ell));
    mpz_divexact(SM->coeff[j], SM->coeff[j], ell);
    ASSERT_ALWAYS(mpz_cmp(ell, SM->coeff[j])>0);
  }
}

/* Assume nSM > 0 */
void
print_sm (FILE *f, mpz_poly_t SM, int nSM, int d)
{
  if (nSM == 0)
    return;
  d--;
  if (d > SM->deg)
    fprintf(f, "0");
  else
    gmp_fprintf(f, "%Zu", SM->coeff[d]);
  for(int j = 1; j < nSM; j++)
  {
    d--;
    if (d > SM->deg)
      fprintf(f, " 0");
    else
      gmp_fprintf(f, " %Zu", SM->coeff[d]);
  }
  //fprintf(f, "\n");
}

/* Compute the Shirokauer maps of a single pair (a,b). 
   SM must be allocated and is viewed as a polynomial of degree F->deg. */
void
sm_single_rel(mpz_poly_t *SM, int64_t a, uint64_t b, mpz_poly_ptr *F,
              const mpz_ptr *eps, const mpz_t ell, const mpz_t ell2,
              const mpz_t invl2)
{
  mpz_poly_t rel;

  mpz_poly_init_set_ab(rel, a, b);
  compute_sm (SM[0], rel, F[0], ell, eps[0], ell2, invl2);
  compute_sm (SM[1], rel, F[1], ell, eps[1], ell2, invl2);
  mpz_poly_clear(rel);
}

void
sm_relset_init (sm_relset_t r, int *d)
{
  for (int side = 0; side < 2; side++) {
    mpz_poly_init (r->num[side], d[side]);
    mpz_poly_init (r->denom[side], d[side]);
  }
}

void
sm_relset_clear (sm_relset_t r)
{
  for (int side = 0; side < 2; side++) {
    mpz_poly_clear (r->num[side]);
    mpz_poly_clear (r->denom[side]);
  }
}

/* Given an array of index of rows and an array of abpolys,
 * construct the polynomial that corresponds to the relation-set, i.e. 
 *          rel = prod(abpoly[rk]^ek)
 * where rk is the index of a row, ek its exponent 
 * and abpoly[rk] = bk*x - ak.
 * rel is build as a fraction (sm_relset_t) and should be initialized.
 */
void
sm_build_one_relset (sm_relset_ptr rel, uint64_t *r, int64_t *e, int len,
                 mpz_poly_t * abpolys, mpz_poly_srcptr *F, const mpz_t ell2)
{
  mpz_t ee;
  mpz_init(ee);  
  mpz_poly_t tmp[2];
  mpz_poly_init(tmp[0], F[0]->deg);
  mpz_poly_init(tmp[1], F[1]->deg);

  /* Set the initial fraction to 1 / 1 */
  for (int side = 0; side < 2; side++) {
    mpz_poly_set_zero (rel->num[side]);
    mpz_poly_setcoeff_si(rel->num[side], 0, 1); 
    mpz_poly_set_zero (rel->denom[side]);
    mpz_poly_setcoeff_si(rel->denom[side], 0, 1); 
  }

  for(int k = 0 ; k < len ; k++)
  {
    /* Should never happen! */
    ASSERT_ALWAYS(e[k] != 0);

    if (e[k] > 0)
    {
      mpz_set_si(ee, e[k]);
      /* TODO: mpz_poly_long_power_mod_f_mod_mpz */
      for (int s = 0; s < 2; ++s) {
        mpz_poly_power_mod_f_mod_mpz (tmp[s], abpolys[r[k]], F[s], ee, ell2);
        mpz_poly_mul_mod_f_mod_mpz (rel->num[s], rel->num[s], tmp[s], F[s],
            ell2, NULL);
      }
    }
    else
    {
      mpz_set_si(ee, -e[k]);
      /* TODO: mpz_poly_long_power_mod_f_mod_mpz */
      for (int s = 0; s < 2; ++s) {
        mpz_poly_power_mod_f_mod_mpz(tmp[s], abpolys[r[k]], F[s], ee, ell2);
        mpz_poly_mul_mod_f_mod_mpz(rel->denom[s], rel->denom[s], tmp[s], F[s],
            ell2, NULL);
      }
    }
  }
  for (int s = 0; s < 2; ++s) {
    mpz_poly_cleandeg(rel->num[s], F[s]->deg);
    mpz_poly_cleandeg(rel->denom[s], F[s]->deg);
  }

  mpz_clear (ee);
  mpz_poly_clear(tmp[0]);
  mpz_poly_clear(tmp[1]);
}
