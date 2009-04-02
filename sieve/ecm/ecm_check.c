/* Usage: ecm_check B1 B2 sigma < in > out

   Compile with:

   gcc -O2 -g -I. -I../../utils ecm_check.c -o ecm_check libfacul.a ../../utils/libutils.a -lgmp

   'in' is a file with one prime number per line
   
   Puts in 'out' the primes which are not found with the curve (B1,B2,sigma).
 */

#include "../../utils/mod_ul_default.h"
#include "facul.h"
#include "ecm.h"

/* return non-zero iff p is found */
static int
tryecm (const unsigned long p, const ecm_plan_t *plan)
{
  modulusredcul_t m;
  modint_t f;
  
  modredcul_initmod_uls (m, &p);
  ecm_ul (f, m, plan);
  modredcul_clearmod (m);
  return mod_intequal_ul (f, p);
}

int
main (int argc, char *argv[])
{
  unsigned long B1, B2, sigma, p;
  ecm_plan_t plan[1];

  B1 = atoi (argv[1]);
  B2 = atoi (argv[2]);
  sigma = atoi (argv[3]);

  ecm_make_plan (plan, B1, B2, BRENT12, sigma, 0);
  while (!feof (stdin))
    {
      scanf ("%lu\n", &p);
      if (p != 2 && p != 3 && tryecm (p, plan) == 0)
        printf ("%lu\n", p);
    }
  return 0;
}

