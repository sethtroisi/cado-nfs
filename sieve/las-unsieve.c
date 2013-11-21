#include "cado.h"
#include <stdlib.h>

#include "las-unsieve.h"

void sieve_info_init_unsieve_data(sieve_info_ptr si)
{
  /* Store largest prime factor of k in si->us->lpf[k], 0 for k=0, 1 for k=1 */
  si->us->entries = (unsieve_entry_t *) malloc (sizeof (unsieve_entry_t) << si->conf->logI);
  FATAL_ERROR_CHECK(si->us->entries == NULL, "malloc failed");
  si->us->entries[0].lpf = 0U;
  si->us->entries[1].lpf = 1U;
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
      p = (c == 1U) ? p : c;
      c = k; do {c /= p;} while (c % p == 0);
      si->us->entries[k].lpf = p;
      si->us->entries[k].cof = c;
      si->us->entries[k].start = (si->I / 2U) % p;
    }

    /* Create pattern for sieving 3 */
    for (size_t x = 0; x < 3 * sizeof(unsigned long); x += 3)
      ((unsigned char *) si->us->pattern3)[x] = 255;
          
}

void sieve_info_clear_unsieve_data(sieve_info_ptr si)
{
  free (si->us->entries);
}

#ifdef UNSIEVE_NOT_COPRIME
static inline void
unsieve_one_prime (unsigned char *line_start, const unsigned int p, 
                   const unsigned int y, const unsigned int start_idx,
                   sieve_info_srcptr si)
{
  unsigned int x, np = p; /* if 2|y, np=2p, else np=p */

  x = start_idx;
  if (y % 2U == 0U)
    {
      np += p;
      if (x % 2U == 0U)
        x += p;
    }
  for ( ; x < si->I; x += np)
    line_start[x] = 255;
}


/* Set locations where gcd(i,j) != 1 to 255 */
void 
unsieve_not_coprime (unsigned char *S, const int N, sieve_info_srcptr si)
{
  unsigned int y; /* Line coordinate within bucket region */
  for (y = 0U + (N == 0U ? 1U : 0U); 
       y < 1U << (LOG_BUCKET_REGION - si->conf->logI); y++)
    {
      unsigned int c = y + (N << (LOG_BUCKET_REGION - si->conf->logI));
      unsigned int p, start_idx;
      unsigned char *line_start = S + (y << si->conf->logI);

      while (c % 2U == 0U) 
        c >>= 1;
      
      while (1)
        {
          p = si->us->entries[c].lpf; /* set p to largest prime factor of c */
          start_idx = si->us->entries[c].start;
          c = si->us->entries[c].cof;
          if (p <= 3)
            break;
          unsieve_one_prime (line_start, p, y, start_idx, si);
        }
      
      if (p == 3U)
        {
          const unsigned int I_ul  = si->I / sizeof (unsigned long);
          unsigned int i;
          unsigned long p0, p1, p2;
          unsigned long * restrict ul_line_start = (unsigned long *) line_start;

          /* If start_idx == 0, we want p0 to contain the pattern which has a
             hit at index 0, which is pattern3[0].
             If start_idx == 1, we want p0 to contain the pattern which has a
             hit at index 1, which is pattern3[1], because 2^k == 1 (mod 3)
             for k > 1.
             If start_idx == 2, we want the sole remaining case p0 = pattern3[2]. */
          p0 = si->us->pattern3[start_idx];
          p1 = si->us->pattern3[(start_idx + 1) % 3];
          p2 = si->us->pattern3[(start_idx + 2) % 3];

          /* Apply pattern to array */
          for (i = 0U; i < I_ul - 2U; i += 3U)
            {
              ul_line_start[i] |= p0;
              ul_line_start[i + 1] |= p1;
              ul_line_start[i + 2] |= p2;
            }
          if (i < I_ul - 1U)
            line_start[i] |= p0;
          if (i < I_ul)
            line_start[i + 1] |= p1;
        }
      ASSERT_ALWAYS(c <= 1);
    }
}
#endif /* ifdef UNSIEVE_NOT_COPRIME */
