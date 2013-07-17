#include "cado.h"
#include <stdlib.h>

#include "las-unsieve.h"

void sieve_info_init_unsieve_data(sieve_info_ptr si)
{
  /* Store largest prime factor of k in si->us->lpf[k], 0 for k=0, 1 for k=1 */
  si->us->lpf = (unsigned int *) malloc (sizeof (unsigned int) << si->conf->logI);
  FATAL_ERROR_CHECK(si->us->lpf == NULL, "malloc failed");
  ASSERT_ALWAYS (si->us->lpf != NULL);
  si->us->lpf[0] = 0U;
  si->us->lpf[1] = 1U;
  for (unsigned int k = 2U; k < si->I; k++)
    {
      unsigned int p, c = k;
      for (p = 2U; p * p <= c; p += 1U + p % 2U)
        {
          while (c % p == 0U)
            c /= p;
          if (c == 1U)
            break;
        }
      si->us->lpf[k] = (c == 1U) ? p : c;
    }
}

void sieve_info_clear_unsieve_data(sieve_info_ptr si)
{
  free (si->us->lpf);
}

#ifdef UNSIEVE_NOT_COPRIME
static void
unsieve_one_prime (unsigned char *line_start, const unsigned int p, 
                   const unsigned int y, sieve_info_srcptr si)
{
  unsigned int x, np = p; /* if 2|y, np=2p, else np=p */

  x = (si->I / 2U) % p;
  if (y % 2U == 0U)
    {
      np += p;
      if (x % 2U == 0U)
        x += p;
    }
  for ( ; x < si->I; x += np)
    line_start[x] = 255;
}


/* Set locations where gcd(i,j) != 1 to 255*/
void 
unsieve_not_coprime (unsigned char *S, const int N, sieve_info_srcptr si)
{
  unsigned int y; /* Line coordinate within bucket region */
  for (y = 0U + (N == 0U ? 1U : 0U); 
       y < 1U << (LOG_BUCKET_REGION - si->logI); y++)
    {
      unsigned int c = y + (N << (LOG_BUCKET_REGION - si->logI));
      unsigned int p;
      unsigned char *line_start = S + (y << si->logI);

      p = si->us->lpf[c]; /* set p to largest prime factor of c */
      if (p > 3U)
        {
          ASSERT (c % p == 0U);
          unsieve_one_prime (line_start, p, y, si);
          do {c /= p;} while (c % p == 0U);
        }
      
      while (c % 2U == 0U) 
        c >>= 1;
      
      if (c % 3U == 0U)
        {
          const unsigned int I_ul  = si->I / sizeof (unsigned long);
          unsigned int i, x;
          unsigned long s[3] = {0UL, 0UL, 0UL};
          unsigned long * restrict ul_line_start = (unsigned long *) line_start;
          
          /* Sieve only the small pattern */
          for (x = (si->I / 2U) % 3U; x < 3U * sizeof(unsigned long); x += 3U)
            ((unsigned char *) s)[x] = 255;
          
          /* Then apply pattern to array */
          for (i = 0U; i < I_ul - 2U; i += 3U)
            {
              ul_line_start[i] |= s[0];
              ul_line_start[i + 1] |= s[1];
              ul_line_start[i + 2] |= s[2];
            }
          if (i < I_ul - 1U)
            line_start[i] |= s[0];
          if (i < I_ul)
            line_start[i + 1] |= s[1];
          
          do {c /= 3U;} while (c % 3U == 0U);
        }

      for (p = 5U; p * p <= c; p += 2U)
        if (c % p == 0U)
          {
            unsieve_one_prime (line_start, p, y, si);
            do {c /= p;} while (c % p == 0U);
          }

      /* Now c == 1 or c is a prime */
      if (c != 1U)
        {
          ASSERT(c > 3U && si->us->lpf[c] == c);
          unsieve_one_prime (line_start, c, y, si);
        }
    }
}
#endif /* ifdef UNSIEVE_NOT_COPRIME */

