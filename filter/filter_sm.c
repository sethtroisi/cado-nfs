#include "cado.h"
#include <stdio.h>
#include "portability.h"
#include "utils.h"
#include "filter_utils.h"
#include "filter_sm.h"

/* Q = P^a mod f, mod p. Note, p is mpz_t */
void
mpz_poly_power_mod_f_mod_mpz_Barrett (mpz_poly_t Q, const mpz_poly_t P, const mpz_poly_t f,
                    const mpz_t a, const mpz_t p, MAYBE_UNUSED const mpz_t invp)
{
  int k = mpz_sizeinbase(a, 2);

  if (mpz_cmp_ui(a, 0) == 0) {
    Q->deg = 0;
    mpz_set_ui(Q->coeff[0], 1);
    return;
  }
  
  // Initialize Q to P
  mpz_poly_copy(Q, P);

  // Horner
  for (k -= 2; k >= 0; k--)
  {
    mpz_poly_sqr_mod_f_mod_mpz(Q, Q, f, p, NULL);  // R <- R^2
    if (mpz_tstbit(a, k))
      mpz_poly_mul_mod_f_mod_mpz(Q, Q, P, f, p, NULL);  // R <- R*P
  }
}

inline void
mpz_poly_init_set_ab (mpz_poly_ptr rel, int64_t a, uint64_t b)
{
  if (b == 0)
  {
    /* freerel */
    mpz_poly_init(rel, 0);
    mpz_poly_setcoeff_int64(rel, 0, a);
    rel->deg=0;
  }
  else
  {
    /* an (a,b)-pair is a degree-1 poly */
    mpz_poly_init(rel, 1);
    mpz_poly_setcoeff_int64(rel, 0, a);
    mpz_poly_setcoeff_int64(rel, 1, -b);
    rel->deg = 1;
  }
}


/*  Reduce frac(=num/denom) mod F mod m (the return value is in frac->num) */
inline void
mpz_poly_reduce_frac_mod_f_mod_mpz (sm_relset_ptr frac, const mpz_poly_t F,
                        const mpz_t m, mpz_t tmp, mpz_poly_t g, mpz_poly_t U, mpz_poly_t V)
{
  if (frac->denom->deg == 0)
  {
    mpz_invert(tmp, frac->denom->coeff[0], m);
    mpz_poly_mul_mpz(frac->num, frac->num, tmp);
    mpz_poly_reduce_mod_mpz(frac->num, frac->num, m);
  }
  else
  {
    mpz_poly_xgcd_mpz (g, F, frac->denom, U, V, m);
    mpz_poly_mul (frac->num, frac->num, V);
    int d = mpz_poly_mod_f_mod_mpz (frac->num->coeff, frac->num->deg, F->coeff,
                                F->deg, m, NULL);
    cleandeg(frac->num, d);
  }
}

/* compute the SM */
inline void
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

inline void
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
void sm_single_rel(mpz_poly_t SM, int64_t a, uint64_t b, mpz_poly_t F, const mpz_t eps, 
                   const mpz_t ell, const mpz_t ell2, const mpz_t invl2)
{
  mpz_poly_t rel;

  SM->deg = 0;
  mpz_poly_setcoeff_si(SM, 0, 1);
  
  mpz_poly_init_set_ab(rel, a, b);
  compute_sm (SM, rel, F, ell, eps, ell2, invl2);
}


/* Construct mpz_poly_t F from cado_poly pol (algebraic side) */
/* XXX Temporary until utils/poly* mess is fixed */
void mpz_poly_t_from_cado_poly_alg (mpz_poly_t F, cado_poly pol)
{
  int deg = pol->pols[ALGEBRAIC_SIDE]->degree;
  mpz_t *f = pol->pols[ALGEBRAIC_SIDE]->f;
  ASSERT_ALWAYS(deg > 1);
  mpz_poly_init (F, deg);
  for (int i = deg; i >= 0; --i)
    mpz_poly_setcoeff (F, i, f[i]);
}
