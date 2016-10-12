#include "getprime.h"

#define M_MAX 253

/* Bytecode an addition chain for k = prod (primes) < B1 */

double addchain_bytecode (const unsigned int B1, const double dbl_cost, const double dbladd_cost, const double dbl_precomp_cost, double add_precomp_cost, bc_state_t *state)
{
  double chaincost;
  mpz_t E;

  /* E = product of all primes between 3 and B1 */
  mpz_init (E);
  mpz_set_ui (E, 1UL);
  prime_info pi;
  prime_info_init (pi);
  p = (unsigned int) getprime_mt (pi);
  ASSERT (p == 3);
  
  for ( ; p <= B1; p = (unsigned int) getprime_mt (pi))
    {
      unsigned long q;
      /* Uses p^k s.t. (p-1)p^(k-1) <= B1, except for p=2 because our 
         base 2 is a QR for primes == 1 (mod 8) already */
      for (q = 1; q <= B1 / (p - 1); q *= p)
        mpz_mul_ui (E, E, p);
    }
  prime_info_clear (pi);
  
  if (verbose)
    gmp_printf ("E = %Zd;\n", E);

  
  /* Computes the bytecode for E */

  size_t nbits;
  nbits = mpz_sizeinbase (E, 2) + 1;   // CHECK!
  *bc = (char *) malloc (nbits * sizeof (char));
  
  double cost, mincost;
  unsigned int m, best_m = 1;

  mincost = f2 (NULL, E, 1, addcost, doublecost);

  /* loop on m */
  for (m = 3 ; m < M_MAX ; m += 2)
    {
      cost = f2 (NULL, E, m, addcost, doublecost);
      if (cost < mincost)
	{
	  mincost = cost;
	  best_m = m;
	}
    }
  
  
  chaincost = f2 (*bc, E, best_m, addcost, doublecost);
	

  mpz_clear (E);

  return chaincost;
}


