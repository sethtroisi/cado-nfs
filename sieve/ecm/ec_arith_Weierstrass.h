#include "ec_arith_common.h"

/*
  implemented functions:
    - weierstrass_add (R:weierstrass, P:weierstrass, Q:weierstrass)
        R <- P+Q

    - weierstrass_dbl (R:weierstrass, P:weierstrass)
        R <- 2*P
*/

static int
weierstrass_add (ec_point_t R, const ec_point_t P, const ec_point_t Q,
		 const residue_t a, const modulus_t m)
{
  return 0;
}

static int
weierstrass_dbl (ec_point_t R, const ec_point_t P,
		 const residue_t a, const modulus_t m)
{
  return 0;
}


/* Computes R = [e]P (mod m)  */
MAYBE_UNUSED 
static int
weierstrass_smul_ui (ec_point_t P, const unsigned long e,
		     const residue_t a, const modulus_t m)
{
  return 0;
}


