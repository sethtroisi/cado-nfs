#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <gmp.h>
#include "ularith.h"
#include "modredc_ul.h"
#include "rootfinder.h"
#include "modredc_ul_default.h"

void omega (residue_t o, residue_t b, const unsigned long k, const modulus_t pp);

int main(int argc, char **argv) {
  unsigned long p, minp = 3, maxp = ~0UL;
  enumeratediv_t div;
  
  if (argc > 1)
    minp = strtoul (argv[1], NULL, 10);
  if (argc > 2)
    maxp = strtoul (argv[2], NULL, 10);

  printf ("minp = %lu, maxp = %lu\n", minp, maxp);


  for (p = minp; p <= maxp; p += 2UL) {
    modulus_t pp;
    residue_t o, b, pow;
    unsigned long k;

    if (p > 3 && p % 3 == 0)
      continue;
    
    mod_initmod_ul (pp, p);
    
    if (!mod_isprime(pp))
      continue;
    
    mod_init (b, pp);
    mod_init (o, pp);
    mod_init (pow, pp);

#if 1
    enumeratediv_init (&div, p-1);
    printf ("divisors(%lu) == vecsort([1", p-1);
    while ((k = enumeratediv(&div)) != 0)
      if (k > 1)
        printf (", %lu", k);
    printf ("])\n");
#endif

    enumeratediv_init (&div, p-1);
    while ((k = enumeratediv(&div)) != 0) {
      omega (o, b, k, pp);

#if 1
      /* Test using output to be piped into Pari/GP */
      printf ("znorder(Mod(%lu, %lu)) == %lu\n", 
              mod_get_ul(o, pp), mod_getmod_ul (pp), k);
#else
      /* Slow but reliable test */
      for (i = 1; i < k; i++) {
        mod_pow_ul (pow, o, i, pp);
        ASSERT_ALWAYS (!mod_is1(pow, pp));
        mod_pow_ul (pow, o, k, pp);
        ASSERT_ALWAYS (mod_is1(pow, pp));
      }
#endif
    }

    mod_clear (b, pp);
    mod_clear (o, pp);
    mod_clear (pow, pp);
    mod_clearmod (pp);
  }
  exit (EXIT_SUCCESS);
}
