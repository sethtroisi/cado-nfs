#include "cado.h"
#include <stdio.h>
#include "portability.h"
#include "utils.h"

/* Init polynomial rel and set it to b*x - a */
void
mpz_poly_init_set_ab (mpz_poly_ptr rel, int64_t a, uint64_t b)
{
  if (b == 0)
  {
    /* freerel */
    mpz_poly_init(rel, 0);
    mpz_poly_setcoeff_int64(rel, 0, a);
  }
  else
  {
    /* an (a,b)-pair is a degree-1 poly */
    mpz_poly_init(rel, 1);
    mpz_poly_setcoeff_int64(rel, 0, a);
    mpz_poly_setcoeff_int64(rel, 1, -b);
  }
}

/* compute the SM */
void
compute_sm (mpz_poly_t SM, mpz_poly_t num, const mpz_poly_t F, const mpz_t ell,
            const mpz_t smexp, const mpz_t ell2, const mpz_t invl2)
{
  mpz_poly_power_mod_f_mod_mpz_Barrett(SM, num, F, smexp, ell2, invl2);
  mpz_poly_sub_ui(SM, 1);

  for(int j = 0; j <= SM->deg; j++)
  {
    ASSERT_ALWAYS(mpz_divisible_p(SM->coeff[j], ell));
    mpz_divexact(SM->coeff[j], SM->coeff[j], ell);
    ASSERT_ALWAYS(mpz_cmp(ell, SM->coeff[j])>0);
  }
}

void
print_sm (FILE *f, mpz_poly_t SM, int nSM)
{
  for(int j=0; j<nSM; j++)
  {
    if (j > SM->deg)
      fprintf(f, "0 ");
    else
      gmp_fprintf(f, "%Zu ", SM->coeff[j]);
  }
  fprintf(f, "\n");
}

/* Computed the Shirokauer maps of a single pair (a,b). 
   SM must be allocated and is viewed as a polynomial of degree F->deg. */
void sm_single_rel(mpz_poly_t SM, int64_t a, uint64_t b, mpz_poly_t F,
                   const mpz_t eps, const mpz_t ell, const mpz_t ell2,
                   const mpz_t invl2)
{
  mpz_poly_t rel;

  mpz_poly_init_set_ab(rel, a, b);
  compute_sm (SM, rel, F, ell, eps, ell2, invl2);
}
