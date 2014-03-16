#ifndef MPFQ_GFP_COMMON_H_
#define MPFQ_GFP_COMMON_H_

#include "gmp.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpfq.h"

#ifdef __cplusplus
extern "C" {
#endif

/***  A general type for (all possible?) prime fields ***/

// Info for Montgomery representation
typedef struct {
  mp_limb_t *invP;   // -1/p mod R   (R = 2^(n*GMP_LIMB_BITS))
  mp_limb_t *invR;  // 1/R mod p
} mgy_info_struct;


// Info for Tonelli-Shanks. 
// Let q be the field cardinality (odd).
// Write q - 1 = h*2^e, with h odd.
// A value of e=0 means that the data is not initialized
typedef struct {
    int e;
    void * z; // a Generator of the 2-Sylow, castable into an elt.
    mp_limb_t * hh; // (h-1)/2, stored as an mpn of same length as q.
} ts_info_struct;


typedef struct {
  mpz_t p;
  mpz_t bigmul_p; // largest multiple of p that fits in an ur_elt
  long url_margin;  // number of adds of unreduced elts that are allowed (>=500)
  mgy_info_struct mgy_info;
  ts_info_struct ts_info;
  mpz_t factor;
  int io_base;
} mpfq_p_field_struct;

typedef mpfq_p_field_struct mpfq_p_field[1];
typedef const mpfq_p_field_struct * mpfq_p_src_field;
typedef mpfq_p_field_struct * mpfq_p_dst_field;


#ifdef __cplusplus
}
#endif

#endif	/* MPFQ_GFP_COMMON_H_ */
