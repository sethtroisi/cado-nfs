#ifndef LAS_UNSIEVE_H_
#define LAS_UNSIEVE_H_

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

#include "las-types.h"

#ifdef __cplusplus
extern "C" {
#endif

void sieve_info_init_unsieve_data(sieve_info_ptr si);
void sieve_info_clear_unsieve_data(sieve_info_ptr si);
void unsieve_not_coprime_line(unsigned char *, unsigned int, unsigned int,
                              unsigned int, unsieve_aux_data_srcptr);
void unsieve_not_coprime (unsigned char *S, const int N, sieve_info_srcptr si);

#ifdef __cplusplus
}
#endif

#endif	/* LAS_UNSIEVE_H_ */
