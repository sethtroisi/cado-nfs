/* auxiliary routines on GMP data-types */

#include <stdint.h>
#include "cado.h"

/* Set z to q. Warning: on 32-bit machines, we cannot use mpz_set_ui! */
void
mpz_set_uint64 (mpz_t z, uint64_t q)
{
  if (sizeof (unsigned long int) == 8)
    mpz_set_ui (z, (unsigned long int) q);
  else
    {
      ASSERT_ALWAYS (sizeof (unsigned long int) == 4);
      mpz_set_ui (z, (unsigned long int) (q >> 32));
      mpz_mul_2exp (z, z, 32);
      mpz_add_ui (z, z, (unsigned long int) (q & 4294967295UL));
    }
}


