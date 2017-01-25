#ifndef LAS_UNSIEVE_H_
#define LAS_UNSIEVE_H_

#include <stdint.h>
#include <vector>
#include "ularith.h"

typedef struct {
    unsigned int lpf, cof, start;
} unsieve_entry_t;

#ifdef HAVE_SSE2
#include <xmmintrin.h>
typedef __m128 unsieve_pattern_t;
#define UNSIEVE_OR(a,b) do{(a) = _mm_or_ps((a),(b));}while(0)
#else
typedef unsigned long unsieve_pattern_t;
#define UNSIEVE_OR(a,b) do{(a) |= (b);}while(0)
#endif

struct unsieve_aux_data_s {
    /* entry[i].lpf is largest prime factor of i, for i < I,
       cof is the cofactor i/lpf^k s.t. lfp^k || i,
       start is (I/2) % lpf */
    unsieve_entry_t *entries;
    unsieve_pattern_t pattern3[3];
    unsieve_pattern_t pattern5[5];
    unsieve_pattern_t pattern7[7];
};
typedef struct unsieve_aux_data_s unsieve_aux_data[1];
typedef struct unsieve_aux_data_s * unsieve_aux_data_ptr;
typedef const struct unsieve_aux_data_s * unsieve_aux_data_srcptr;

/* {{{ j_div
 * A structure for factoring the line-coordinate j into distinct odd primes
 * and performing fast divisibility tests of i by those primes.
 */
struct j_div_s {
  unsigned int p,     /* The odd prime that divides this entry */
               cof,   /* The cofactor after dividing out p as often as possible */
               inv,   /* p^(-1) (mod 2^32) */
               bound; /* (2^32-1) / p */
};
typedef struct j_div_s * j_div_ptr;
typedef const struct j_div_s * j_div_srcptr;


/* Helper function shared by las-unsieve.cpp and las-unsieve-sse2.cpp.
   Polluting global namespace is a little unfortunate here - #ifdef? */
static inline unsigned int 
extract_j_div(unsigned int (*div)[2], const unsigned int j, j_div_srcptr j_div, 
              const unsigned int pmin, const unsigned int pmax)
{
    unsigned int c, nr_div = 0;
    /* For each distict odd prime factor p of j, if pmin <= p <= pmax,
       store the inverse and bound in array */
    c = j;
    c >>= ularith_ctz(c);
    while (c > 1) {
      unsigned int p = j_div[c].p;
      if (p < pmin)
          break;
      if (p <= pmax) {
          div[nr_div][0] = j_div[c].inv;
          div[nr_div++][1] = j_div[c].bound;
      }
      c = j_div[c].cof;
    }
    return nr_div;
}


unsieve_aux_data_srcptr init_unsieve_data(uint32_t);
void clear_unsieve_data(unsieve_aux_data_srcptr);

j_div_srcptr init_j_div(uint32_t);
void clear_j_div(j_div_srcptr);
void search_survivors_in_line(unsigned char * const restrict[2],
        const unsigned char[2], unsigned int, unsigned int, int,
        j_div_srcptr, unsigned int, unsieve_aux_data_srcptr,
        std::vector<uint32_t> &, bool);
#ifdef HAVE_SSE2 
void search_survivors_in_line_sse2(unsigned char * const restrict[2],
        const unsigned char[2], unsigned int, unsigned int, int, j_div_srcptr,
        unsigned int, std::vector<uint32_t> &);
#endif

#endif	/* LAS_UNSIEVE_H_ */
