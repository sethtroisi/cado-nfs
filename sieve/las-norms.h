#ifndef LAS_NORMS_H_
#define LAS_NORMS_H_

#include <stdint.h>
#include "las-types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*  initializing norms */
/* Knowing the norm on the rational side is bounded by 2^(2^k), compute
   lognorms approximations for k bits of exponent + NORM_BITS-k bits
   of mantissa */
void
init_norms (sieve_info_ptr si, int side);


/*  initialize norms for bucket regions */
/* Initialize lognorms on the rational side for the bucket_region
 * number N.
 * For the moment, nothing clever, wrt discarding (a,b) pairs that are
 * not coprime, except for the line j=0.
 */
void
init_rat_norms_bucket_region (unsigned char *S, unsigned int N, sieve_info_ptr si);

/* Initialize lognorms on the algebraic side for the bucket
 * number N.
 * Only the survivors of the rational sieve will be initialized, the
 * others are set to 255. Case GCD(i,j)!=1 also gets 255.
 * return nothing because the number of reports (= number of norm
 * initialisations) is algorithm dependent of ALG_RAT & ALG_LAZY.
 */
void
init_alg_norms_bucket_region (unsigned char *alg_S, 
                              unsigned char *rat_S,  unsigned int N, 
                              sieve_info_ptr si);

/* This prepares the auxiliary data which is used by
 * init_rat_norms_bucket_region and init_alg_norms_bucket_region
 */
void sieve_info_init_norm_data(FILE * output, sieve_info_ptr si, double q0d, int qside);

void sieve_info_clear_norm_data(sieve_info_ptr si);

void sieve_info_update_norm_data(sieve_info_ptr si, int nb_threads);
void sieve_info_init_norm_data_sq (sieve_info_ptr si, unsigned long q);

static inline int 
sieve_info_test_lognorm (const unsigned char C1, 
                         const unsigned char C2, 
                         const unsigned char S1,
                         const unsigned char S2)
{
  return S1 <= C1 && S2 <= C2;
}

#ifdef HAVE_SSE2
static inline int 
sieve_info_test_lognorm_sse2(__m128i *alg_S, const __m128i alg_pattern,
                             __m128i *rat_S, const __m128i rat_pattern)
{
    const __m128i zero = _mm_set1_epi8(0);
    const __m128i sign_conversion = _mm_set1_epi8(-128);
    __m128i a = *alg_S;
    __m128i r = *rat_S;
    __m128i m1, m2;
    /* _mm_cmpgt_epi8() performs signed comparison, but we have unsigned
    bytes. We can switch to signed in a way that preserves ordering by
    flipping the MSB, e.g., 0xFF (255 unsigned) becomes 0x7F (+127 signed), 
    and 0x00 (0 unsigned) becomes 0x80 (-128 signed).
    If a byte in the first operand is greater than the corresponding byte in
    the second operand, the corresponding byte in the result is set to all 1s
    (i.e., to 0xFF); otherwise, it is set to all 0s. */
    m1 = _mm_cmpgt_epi8 (_mm_xor_si128(a, sign_conversion), alg_pattern);
    m2 = _mm_cmpgt_epi8 (_mm_xor_si128(r, sign_conversion), rat_pattern);
    /* Logically OR the two masks: if at least one was 255, the sieve entry
    should be set to 255 */
    m1 = _mm_or_si128(m1, m2);
    *alg_S = _mm_or_si128(a, m1);
    *rat_S = _mm_or_si128(r, m1);
    /* Compute number of non-zero bytes. First sign flip: 0xFF -> 0x1 */
    m1 = _mm_sub_epi8(zero, m1);
    /* Using Sum of Absolute Differences with 0, which for us gives the
    number of non-zero bytes. This Sum of Absolute Differences uses
    unsigned arithmetic, thus we needed the sign flip first */
    m1 = _mm_sad_epu8(m1, zero);
    /* Sum is stored in two parts */
    int nr_set = _mm_extract_epi16(m1, 0) + _mm_extract_epi16(m1, 4);
    /* Return number of bytes that were not set to 255 */
    return 16 - nr_set;
}
#endif

#ifdef __cplusplus
}
#endif

#endif	/* LAS_NORMS_H_ */
