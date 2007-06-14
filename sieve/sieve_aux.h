/* Auxiliary inlined functions used in the siever and related modules */

#include "config.h"
#include "cado.h"

static inline fbroot_t
first_sieve_loc (const fbprime_t p, const fbroot_t r, const fbroot_t amin_p, 
                 const unsigned long b, const int odd)
{
  modulus m;
  residue r1, r2;
  fbroot_t d;

  /* Find first index d in sievearray where p on this root divides.
     So we want a/b == r (mod p) <=> a == br (mod p). Then we want
     d so that amin + d * (1 + odd) == a (mod p) 
     <=> d == (a - amin) / (1 + odd) (mod p) 
     == (br - amin) / (1 + odd) (mod p) */

  ASSERT (odd == 0 || odd == 1);
  ASSERT (!odd || (b & 1) == 0);
  ASSERT (r < p);
  ASSERT (amin_p < p);
  ASSERT_EXPENSIVE (b % p != 0);

  mod_initmod_ul (m, p); /* Most of the mod_*() calls are no-ops */
  mod_init (r1, m);
  mod_init (r2, m);
  mod_set_ul (r1, b, m); /* Modular reduction */
  mod_set_ul_reduced (r2, r, m);
  mod_mul (r1, r1, r2, m); /* Multiply and mod reduction. If we keep 
    r_i in Montgomery representation, a single mul/REDC will compute 
    the product, reduce it mod p and return it as an integer. However,
    this requires storing p^{-1} (mod 2^32) in the factor base. */
  mod_sub_ul (r1, r1, amin_p, m);
  if (odd)
    mod_div2 (r1, r1, m);
  d = mod_get_ul (r1, m); 
  ASSERT (d < p);
  mod_clear (r1, m);
  mod_clear (r2, m);
  mod_clearmod (m);

#ifdef PARI
  printf ("(" FBROOT_FORMAT " + " FBROOT_FORMAT " * (1 + %d)) %% " FBPRIME_FORMAT
          " == (%lu * " FBROOT_FORMAT ") %% " FBPRIME_FORMAT " /* PARI */\n",
          amin_p, d, odd, p, b, r, p);
#endif

  return d;
}
