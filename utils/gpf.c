#include "cado.h"
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include "getprime.h"
#include "ularith.h"
#include "gpf.h"

unsigned int gpf_max = 0;
unsigned int *gpf = NULL;

void gpf_init(unsigned int m)
{
      if (m <= gpf_max)
          return;
      if (gpf == NULL) {
        gpf = malloc ((m + 1) * sizeof(unsigned int));
      } else {
        gpf = realloc (gpf, (m + 1) * sizeof(unsigned int));
      }
      gpf_max = m;
      
      /* When we increase gpf_max we could and probably should keep
         the old gpf data and sieve only over the newly added area.
         For now, Q&D: sieve everything again. */
      for (size_t i = 0; i <= m; i++) {
        gpf[i] = i;
      }
      prime_info pi;
      prime_info_init (pi);
      const unsigned int max_sieve = sqrt(m + 1);
      
      /* Do 2 separately so the compiler can use bit shift for division */
      if (2 <= max_sieve) {
          for (unsigned int i = 2; i <= m ; i += 2) {
              ASSERT(gpf[i] % 2 == 0);
              while (gpf[i] > 2) {
                gpf[i] /= 2;
                if (gpf[i] % 2 != 0)
                    break;
              }
          }
      }
      
      for (unsigned int p = getprime_mt (pi) ; p <= max_sieve; p = getprime_mt (pi)) {
          const unsigned int inv_p = ularith_invmod(p);
          const unsigned int lim_p = UINT_MAX / p;
          for (unsigned int i = 2*p; i <= m ; i += p) {
              ASSERT_EXPENSIVE(gpf[i] % p == 0);
              ASSERT_EXPENSIVE(gpf[i] / p == gpf[i] * inv_p);

              while (gpf[i] > p) {
                unsigned int candidate = gpf[i] * inv_p;
                if (candidate > lim_p)
                  break;
                gpf[i] = candidate;
              }
          }
      }
      prime_info_clear (pi);
}

void gpf_clear()
{
    free(gpf);
    gpf = NULL;
    gpf_max = 0;
}
