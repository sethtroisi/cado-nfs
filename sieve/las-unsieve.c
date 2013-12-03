#include "cado.h"
#include <stdlib.h>
#include <string.h>

#include "las-unsieve.h"

/* Set every stride-th byte, starting at index 0, to 255 in an array of
   stride unsieve_pattern_t's */
static void
minisieve(unsieve_pattern_t *array, size_t stride)
{
    memset (array, 0, stride * sizeof(unsieve_pattern_t));
    for (size_t i = 0; i < stride * sizeof(unsieve_pattern_t); i += stride)
      ((unsigned char *) array)[i] = 255;
}

void sieve_info_init_unsieve_data(sieve_info_ptr si)
{
  /* Store largest prime factor of k in si->us->lpf[k], 0 for k=0, 1 for k=1 */
  si->us->entries = (unsieve_entry_t *) malloc (sizeof (unsieve_entry_t) << si->conf->logI);
  FATAL_ERROR_CHECK(si->us->entries == NULL, "malloc failed");
  si->us->entries[0].lpf = 0U;
  si->us->entries[0].cof = 0U;
  si->us->entries[0].start = 0U;
  si->us->entries[1].lpf = 1U;
  si->us->entries[1].cof = 1U;
  si->us->entries[1].start = 0U;
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
    minisieve(si->us->pattern3, 3);
    /* Create pattern for sieving 5 */
    minisieve(si->us->pattern5, 5);
    /* Create pattern for sieving 7 */
    minisieve(si->us->pattern7, 7);
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


static inline void
unsieve_3(unsigned char *line_start, const unsigned int start_idx,
          sieve_info_srcptr si)
{
  const unsigned int I_upt  = si->I / sizeof (unsieve_pattern_t);
  unsigned int i, pattern_idx;
  unsieve_pattern_t p0, p1, p2;
  unsieve_pattern_t * restrict ul_line_start = (unsieve_pattern_t *) line_start;

  if (sizeof(unsieve_pattern_t) == 4)
    /* -4^(-1) == 2 (mod 3) */
    pattern_idx = (2 * start_idx) % 3;
  else if (sizeof(unsieve_pattern_t) == 8)
    /* -8^(-1) == 1 (mod 3) */
    pattern_idx = start_idx;
  else if (sizeof(unsieve_pattern_t) == 16)
    /* -16^(-1) == 2 (mod 3) */
    pattern_idx = (2 * start_idx) % 3;
  else
    abort();
  
  p0 = si->us->pattern3[pattern_idx];
  p1 = si->us->pattern3[(pattern_idx + 1) % 3];
  p2 = si->us->pattern3[(pattern_idx + 2) % 3];

  ASSERT_ALWAYS(((unsigned char *)&p0)[start_idx] == 255);
  
  /* Apply pattern to array */
  for (i = 0U; i < I_upt - 2U; i += 3U)
    {
      UNSIEVE_OR(ul_line_start[i], p0);
      UNSIEVE_OR(ul_line_start[i + 1], p1);
      UNSIEVE_OR(ul_line_start[i + 2], p2);
    }
  if (i < I_upt)
    UNSIEVE_OR(ul_line_start[i], p0);
  if (i + 1 < I_upt)
    UNSIEVE_OR(ul_line_start[i + 1], p1);
}


static inline void
unsieve_5(unsigned char *line_start, const unsigned int start_idx,
          sieve_info_srcptr si)
{
  const unsigned int I_upt  = si->I / sizeof (unsieve_pattern_t);
  unsigned int i;
  unsieve_pattern_t p0, p1, p2, p3, p4;
  unsieve_pattern_t * restrict ul_line_start = (unsieve_pattern_t *) line_start;
  size_t pattern_idx;

  if (sizeof(unsieve_pattern_t) == 4) {
    /* -4^(-1) == 1 (mod 5) */
    pattern_idx = start_idx;
  } else if (sizeof(unsieve_pattern_t) == 8) {
    /* -8^(-1) == 3 (mod 5) */
    pattern_idx = (3 * start_idx) % 5;
  } else if (sizeof(unsieve_pattern_t) == 16) {
    /* -16^(-1) == -1 (mod 5) */
    pattern_idx = (5 - start_idx) % 5;
  } else
    abort();
  
  p0 = si->us->pattern5[pattern_idx];
  p1 = si->us->pattern5[(pattern_idx + 1) % 5];
  p2 = si->us->pattern5[(pattern_idx + 2) % 5];
  p3 = si->us->pattern5[(pattern_idx + 3) % 5];
  p4 = si->us->pattern5[(pattern_idx + 4) % 5];

  ASSERT_ALWAYS(((unsigned char *)&p0)[start_idx] == 255);
  
  /* Apply pattern to array */
  for (i = 0U; i < I_upt - 4U; i += 5U)
    {
      UNSIEVE_OR(ul_line_start[i], p0);
      UNSIEVE_OR(ul_line_start[i + 1], p1);
      UNSIEVE_OR(ul_line_start[i + 2], p2);
      UNSIEVE_OR(ul_line_start[i + 3], p3);
      UNSIEVE_OR(ul_line_start[i + 4], p4);
    }
  if (i < I_upt)
    UNSIEVE_OR(ul_line_start[i], p0);
  if (i + 1 < I_upt)
    UNSIEVE_OR(ul_line_start[i + 1], p1);
  if (i + 2 < I_upt)
    UNSIEVE_OR(ul_line_start[i + 2], p2);
  if (i + 3 < I_upt)
    UNSIEVE_OR(ul_line_start[i + 3], p3);
}

static inline void
unsieve_7(unsigned char *line_start, const unsigned int start_idx,
          sieve_info_srcptr si)
{
  const unsigned int I_upt  = si->I / sizeof (unsieve_pattern_t);
  unsigned int i;
  unsieve_pattern_t p0, p1, p2, p3, p4, p5, p6;
  unsieve_pattern_t * restrict ul_line_start = (unsieve_pattern_t *) line_start;
  size_t pattern_idx;

  if (sizeof(unsieve_pattern_t) == 4) {
    /* -4^(-1) == 5 (mod 7) */
    pattern_idx = (5 * start_idx) % 7;
  } else if (sizeof(unsieve_pattern_t) == 8) {
    /* -8^(-1) == -1 (mod 7) */
    pattern_idx = (7 - start_idx) % 7;
  } else if (sizeof(unsieve_pattern_t) == 16) {
    /* -16^(-1) == 3 (mod 7) */
    pattern_idx = (3 * start_idx) % 7;
  } else
    abort();
  
  p0 = si->us->pattern7[pattern_idx];
  p1 = si->us->pattern7[(pattern_idx + 1) % 7];
  p2 = si->us->pattern7[(pattern_idx + 2) % 7];
  p3 = si->us->pattern7[(pattern_idx + 3) % 7];
  p4 = si->us->pattern7[(pattern_idx + 4) % 7];
  p5 = si->us->pattern7[(pattern_idx + 5) % 7];
  p6 = si->us->pattern7[(pattern_idx + 6) % 7];

  ASSERT_ALWAYS(((unsigned char *)&p0)[start_idx] == 255);
  
  /* Apply pattern to array */
  for (i = 0U; i < I_upt - 6U; i += 7U)
    {
      UNSIEVE_OR(ul_line_start[i], p0);
      UNSIEVE_OR(ul_line_start[i + 1], p1);
      UNSIEVE_OR(ul_line_start[i + 2], p2);
      UNSIEVE_OR(ul_line_start[i + 3], p3);
      UNSIEVE_OR(ul_line_start[i + 4], p4);
      UNSIEVE_OR(ul_line_start[i + 5], p5);
      UNSIEVE_OR(ul_line_start[i + 6], p6);
    }
  if (i < I_upt)
    UNSIEVE_OR(ul_line_start[i], p0);
  if (i + 1 < I_upt)
    UNSIEVE_OR(ul_line_start[i + 1], p1);
  if (i + 2 < I_upt)
    UNSIEVE_OR(ul_line_start[i + 2], p2);
  if (i + 3 < I_upt)
    UNSIEVE_OR(ul_line_start[i + 3], p3);
  if (i + 4 < I_upt)
    UNSIEVE_OR(ul_line_start[i + 4], p4);
  if (i + 5 < I_upt)
    UNSIEVE_OR(ul_line_start[i + 5], p5);
}

/* Set locations where gcd(i,j) != 1 to 255 */
void 
unsieve_not_coprime (unsigned char *S, const int N, sieve_info_srcptr si)
{
  const unsigned int lines_per_region = 1U << (LOG_BUCKET_REGION - si->conf->logI);
  unsigned int y; /* Line coordinate within bucket region */
  for (y = (N == 0U ? 1U : 0U); y < lines_per_region; y++)
    {
      unsigned int c = y + N * lines_per_region;
      unsigned int p, start_idx;
      unsigned char *line_start = S + (y << si->conf->logI);

      while (c % 2U == 0U) 
        c >>= 1;

      while (1)
        {
          p = si->us->entries[c].lpf; /* set p to largest prime factor of c */
          start_idx = si->us->entries[c].start;
          c = si->us->entries[c].cof;
          if (p <= 7)
            break;
          unsieve_one_prime (line_start, p, y, start_idx, si);
        }
      
      if (p == 7U)
        {
          unsieve_7(line_start, start_idx, si);
          p = si->us->entries[c].lpf;
          start_idx = si->us->entries[c].start;
          c = si->us->entries[c].cof;
        }

      if (p == 5U)
        {
          unsieve_5(line_start, start_idx, si);
          p = si->us->entries[c].lpf;
          start_idx = si->us->entries[c].start;
          c = si->us->entries[c].cof;
        }

      if (p == 3U)
        {
          unsieve_3(line_start, start_idx, si);
        }
      ASSERT_ALWAYS(c <= 1);
    }
}
#endif /* ifdef UNSIEVE_NOT_COPRIME */
