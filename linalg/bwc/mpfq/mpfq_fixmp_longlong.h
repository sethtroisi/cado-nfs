#ifndef MPFQ_FIXMP_LONGLONG_H_
#define MPFQ_FIXMP_LONGLONG_H_

/* MPFQ generated file -- do not edit */

#include <gmp.h>
#include <limits.h>
#ifdef	MPFQ_LAST_GENERATED_TAG
#undef	MPFQ_LAST_GENERATED_TAG
#endif
#define MPFQ_LAST_GENERATED_TAG      fixmp

/* Options used:{ features={  }, tag=fixmp, w=0, } */


#ifdef  __cplusplus
extern "C" {
#endif
static inline
mp_limb_t mpfq_fixmp_1_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_1_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_1_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_1_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_1_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_1_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_1_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_1_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_1_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_1_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_1_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_1_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_1_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_1_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_1_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_1_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_1_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_1_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_1_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_1_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_1_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_1_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_1_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_1_lshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_1_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_1_long_lshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_1_long_rshift(mp_limb_t *, int, int);
static inline
mp_limb_t mpfq_fixmp_2_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_2_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_2_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_2_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_2_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_2_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_2_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_2_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_2_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_2_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_2_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_2_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_2_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_2_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_2_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_2_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_2_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_2_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_2_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_2_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_2_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_2_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_2_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_2_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_2_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_2_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_2_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_3_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_3_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_3_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_3_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_3_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_3_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_3_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_3_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_3_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_3_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_3_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_3_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_3_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_3_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_3_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_3_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_3_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_3_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_3_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_3_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_3_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_3_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_3_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_3_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_3_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_3_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_3_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_4_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_4_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_4_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_4_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_4_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_4_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_4_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_4_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_4_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_4_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_4_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_4_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_4_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_4_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_4_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_4_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_4_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_4_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_4_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_4_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_4_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_4_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_4_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_4_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_4_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_4_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_4_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_5_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_5_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_5_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_5_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_5_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_5_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_5_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_5_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_5_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_5_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_5_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_5_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_5_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_5_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_5_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_5_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_5_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_5_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_6_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_6_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_6_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_6_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_6_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_6_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_6_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_6_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_6_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_6_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_6_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_6_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_6_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_6_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_6_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_6_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_6_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_6_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_6_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_6_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_6_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_6_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_6_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_6_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_6_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_6_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_6_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_7_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_7_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_7_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_7_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_7_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_7_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_7_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_7_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_7_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_7_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_7_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_7_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_7_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_7_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_7_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_7_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_7_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_7_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_7_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_7_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_7_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_7_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_7_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_7_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_7_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_7_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_7_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_8_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_8_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_8_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_8_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_8_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_8_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_8_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_8_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_8_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_8_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_8_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_8_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_8_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_8_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_8_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_8_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_8_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_8_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_8_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_8_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_8_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_8_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_8_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_8_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_8_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_8_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_8_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_9_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_9_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_9_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_9_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_9_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_9_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_9_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_9_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_9_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_9_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_9_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_9_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_9_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_9_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_9_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_9_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_9_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_9_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_9_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_9_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_9_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_9_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_9_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_9_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_9_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_9_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_9_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_10_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_10_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_10_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_10_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_10_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_10_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_10_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_10_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_10_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_10_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_10_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_10_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_10_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_10_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_10_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_10_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_10_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_10_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_10_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_10_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_10_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_10_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_10_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_10_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_10_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_10_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_10_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_11_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_11_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_11_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_11_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_11_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_11_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_11_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_11_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_11_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_11_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_11_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_11_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_11_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_11_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_11_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_11_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_11_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_11_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_11_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_11_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_11_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_11_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_11_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_11_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_11_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_11_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_11_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_12_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_12_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_12_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_12_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_12_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_12_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_12_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_12_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_12_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_12_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_12_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_12_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_12_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_12_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_12_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_12_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_12_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_12_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_12_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_12_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_12_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_12_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_12_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_12_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_12_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_12_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_12_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_13_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_13_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_13_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_13_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_13_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_13_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_13_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_13_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_13_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_13_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_13_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_13_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_13_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_13_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_13_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_13_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_13_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_13_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_13_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_13_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_13_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_13_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_13_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_13_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_13_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_13_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_13_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_14_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_14_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_14_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_14_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_14_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_14_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_14_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_14_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_14_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_14_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_14_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_14_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_14_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_14_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_14_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_14_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_14_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_14_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_14_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_14_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_14_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_14_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_14_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_14_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_14_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_14_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_14_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_15_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_15_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_15_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_15_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_15_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_15_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_15_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_15_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_15_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_15_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_15_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_15_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_15_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_15_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_15_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_15_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_15_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_15_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_15_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_15_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_15_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_15_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_15_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_15_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_15_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_15_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_15_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_0_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_0_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_0_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_0_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_0_5_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_0_5_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_0_5_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_0_5_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_0_5_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_0_5_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_0_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_0_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_0_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_0_5_addmul05_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_0_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_0_5_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_0_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_0_5_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_0_5_addmul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_0_5_mul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_0_5_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_0_5_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_0_5_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_0_5_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_0_5_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_0_5_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_0_5_lshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_0_5_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_0_5_long_lshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_0_5_long_rshift(mp_limb_t *, int, int);
static inline
mp_limb_t mpfq_fixmp_1_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_1_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_1_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_1_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_1_5_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_1_5_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_1_5_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_1_5_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_1_5_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_1_5_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_1_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_1_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_1_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_1_5_addmul05_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_1_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_1_5_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_1_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_1_5_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_1_5_addmul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_1_5_mul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_1_5_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_1_5_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_1_5_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_1_5_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_1_5_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_1_5_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_1_5_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_1_5_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_1_5_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_1_5_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_2_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_2_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_2_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_2_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_2_5_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_2_5_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_2_5_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_2_5_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_2_5_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_2_5_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_2_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_2_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_2_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_2_5_addmul05_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_2_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_2_5_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_2_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_2_5_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_2_5_addmul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_2_5_mul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_2_5_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_2_5_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_2_5_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_2_5_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_2_5_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_2_5_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_2_5_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_2_5_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_2_5_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_2_5_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_3_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_3_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_3_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_3_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_3_5_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_3_5_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_3_5_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_3_5_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_3_5_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_3_5_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_3_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_3_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_3_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_3_5_addmul05_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_3_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_3_5_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_3_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_3_5_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_3_5_addmul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_3_5_mul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_3_5_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_3_5_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_3_5_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_3_5_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_3_5_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_3_5_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_3_5_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_3_5_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_3_5_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_3_5_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_4_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_4_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_4_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_4_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_4_5_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_4_5_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_4_5_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_4_5_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_4_5_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_4_5_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_4_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_4_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_4_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_4_5_addmul05_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_4_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_4_5_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_4_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_4_5_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_4_5_addmul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_4_5_mul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_4_5_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_4_5_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_4_5_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_4_5_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_4_5_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_4_5_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_4_5_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_4_5_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_4_5_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_4_5_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_5_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_5_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_5_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_5_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_5_5_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_5_5_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_5_5_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_5_5_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_5_5_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_5_5_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_5_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_5_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_5_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_5_5_addmul05_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_5_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_5_5_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_5_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_5_5_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_5_5_addmul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_5_5_mul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_5_5_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_5_5_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_5_5_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_5_5_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_5_5_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_5_5_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_5_5_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_5_5_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_5_5_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_5_5_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_6_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_6_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_6_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_6_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_6_5_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_6_5_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_6_5_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_6_5_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_6_5_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_6_5_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_6_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_6_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_6_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_6_5_addmul05_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_6_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_6_5_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_6_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_6_5_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_6_5_addmul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_6_5_mul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_6_5_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_6_5_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_6_5_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_6_5_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_6_5_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_6_5_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_6_5_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_6_5_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_6_5_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_6_5_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_7_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_7_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_7_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_7_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_7_5_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_7_5_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_7_5_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_7_5_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_7_5_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_7_5_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_7_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_7_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_7_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_7_5_addmul05_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_7_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_7_5_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_7_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_7_5_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_7_5_addmul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_7_5_mul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_7_5_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_7_5_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_7_5_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_7_5_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_7_5_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_7_5_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_7_5_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_7_5_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_7_5_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_7_5_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_8_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_8_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_8_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_8_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_8_5_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_8_5_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_8_5_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_8_5_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_8_5_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_8_5_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_8_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_8_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_8_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_8_5_addmul05_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_8_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_8_5_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_8_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_8_5_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_8_5_addmul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_8_5_mul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_8_5_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_8_5_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_8_5_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_8_5_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_8_5_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_8_5_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_8_5_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_8_5_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_8_5_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_8_5_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_9_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_9_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_9_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_9_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_9_5_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_9_5_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_9_5_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_9_5_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_9_5_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_9_5_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_9_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_9_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_9_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_9_5_addmul05_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_9_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_9_5_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_9_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_9_5_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_9_5_addmul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_9_5_mul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_9_5_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_9_5_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_9_5_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_9_5_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_9_5_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_9_5_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_9_5_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_9_5_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_9_5_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_9_5_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_10_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_10_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_10_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_10_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_10_5_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_10_5_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_10_5_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_10_5_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_10_5_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_10_5_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_10_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_10_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_10_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_10_5_addmul05_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_10_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_10_5_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_10_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_10_5_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_10_5_addmul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_10_5_mul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_10_5_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_10_5_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_10_5_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_10_5_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_10_5_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_10_5_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_10_5_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_10_5_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_10_5_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_10_5_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_11_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_11_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_11_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_11_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_11_5_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_11_5_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_11_5_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_11_5_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_11_5_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_11_5_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_11_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_11_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_11_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_11_5_addmul05_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_11_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_11_5_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_11_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_11_5_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_11_5_addmul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_11_5_mul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_11_5_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_11_5_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_11_5_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_11_5_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_11_5_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_11_5_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_11_5_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_11_5_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_11_5_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_11_5_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_12_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_12_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_12_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_12_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_12_5_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_12_5_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_12_5_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_12_5_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_12_5_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_12_5_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_12_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_12_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_12_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_12_5_addmul05_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_12_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_12_5_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_12_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_12_5_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_12_5_addmul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_12_5_mul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_12_5_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_12_5_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_12_5_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_12_5_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_12_5_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_12_5_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_12_5_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_12_5_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_12_5_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_12_5_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_13_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_13_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_13_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_13_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_13_5_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_13_5_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_13_5_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_13_5_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_13_5_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_13_5_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_13_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_13_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_13_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_13_5_addmul05_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_13_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_13_5_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_13_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_13_5_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_13_5_addmul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_13_5_mul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_13_5_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_13_5_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_13_5_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_13_5_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_13_5_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_13_5_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_13_5_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_13_5_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_13_5_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_13_5_lshift(mp_limb_t *, int);
static inline
mp_limb_t mpfq_fixmp_14_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_14_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_14_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_14_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_14_5_add_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_14_5_sub_ui(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_14_5_add_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_14_5_sub_ui_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
int mpfq_fixmp_14_5_cmp(const mp_limb_t *, const mp_limb_t *);
static inline
int mpfq_fixmp_14_5_cmp_ui(const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_14_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_14_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
mp_limb_t mpfq_fixmp_14_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_14_5_addmul05_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_14_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_14_5_sqr(mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_14_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_14_5_shortmul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
mp_limb_t mpfq_fixmp_14_5_addmul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_14_5_mul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
static inline
void mpfq_fixmp_14_5_mod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_14_5_rshift(mp_limb_t *, int);
static inline
void mpfq_fixmp_14_5_long_rshift(mp_limb_t *, int, int);
static inline
void mpfq_fixmp_14_5_long_lshift(mp_limb_t *, int, int);
static inline
int mpfq_fixmp_14_5_invmod(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_14_5_redc(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_14_5_redc_ur(mp_limb_t *, mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_14_5_mgy_encode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_14_5_mgy_decode(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
static inline
void mpfq_fixmp_14_5_lshift(mp_limb_t *, int);
#ifdef  __cplusplus
}
#endif

/* Implementations for inlines */
#if !defined(HAVE_native_mpfq_fixmp_1_add)
/* x, y, and z have 1 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 1_redc_ur */
static inline
mp_limb_t mpfq_fixmp_1_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_add) */

#if !defined(HAVE_native_mpfq_fixmp_1_sub)
/* x, y, and z have 1 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
static inline
mp_limb_t mpfq_fixmp_1_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_sub) */

#if !defined(HAVE_native_mpfq_fixmp_1_add_nc)
/* x, y, and z have 1 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_1_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_1_sub_nc)
/* x, y, and z have 1 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_1_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_1_add_ui)
/* x, y, and z have 1 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_1_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_1_sub_ui)
/* x, y, and z have 1 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_1_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_1_add_ui_nc)
/* x, y, and z have 1 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_1_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_1_sub_ui_nc)
/* x, y, and z have 1 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_1_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_1_cmp)
/* x and y have 1 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
static inline
int mpfq_fixmp_1_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 1-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_1_cmp_ui)
/* x has 1 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
static inline
int mpfq_fixmp_1_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_1_addmul1)
/* x has 1 words, z has 3.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
/* Triggered by: 1_redc, 1_redc_ur */
static inline
mp_limb_t mpfq_fixmp_1_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    z[1] += carry;
    return (z[1]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_1_addmul1_nc)
/* x has 1 words, z has 3.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 1_mul, 1_mgy_decode, 2_sqr, 2_shortmul, 3_sqr, 3_shortmul, 4_sqr, 4_shortmul, 5_sqr, 5_shortmul, 6_sqr, 6_shortmul, 7_sqr, 7_shortmul, 8_sqr, 8_shortmul, 9_sqr, 9_shortmul, 10_sqr, 10_shortmul, 11_sqr, 11_shortmul, 12_sqr, 12_shortmul, 13_sqr, 13_shortmul, 14_sqr, 14_shortmul, 15_sqr, 15_shortmul, 1_5_sqr, 1_5_shortmul, 2_5_sqr, 2_5_shortmul, 3_5_sqr, 3_5_shortmul, 4_5_sqr, 4_5_shortmul, 5_5_sqr, 5_5_shortmul, 6_5_sqr, 6_5_shortmul, 7_5_sqr, 7_5_shortmul, 8_5_sqr, 8_5_shortmul, 9_5_sqr, 9_5_shortmul, 10_5_sqr, 10_5_shortmul, 11_5_sqr, 11_5_shortmul, 12_5_sqr, 12_5_shortmul, 13_5_sqr, 13_5_shortmul, 14_5_sqr, 14_5_shortmul */
static inline
void mpfq_fixmp_1_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    z[1] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_1_addmul1_shortz)
/* x has 1 words, z has 2.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
static inline
mp_limb_t mpfq_fixmp_1_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_1_mul)
/* x and y have 1 words, z has 4. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 1_mgy_decode */
static inline
void mpfq_fixmp_1_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 2; z[i++] = 0) ;
    mpfq_fixmp_1_addmul1_nc (z + 0, x, y[0]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_mul) */

#if !defined(HAVE_native_mpfq_fixmp_1_sqr)
/* x has 1 words, z has 4. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_1_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[2] = {0,};
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpn_lshift(buf, buf, 2, 1);
    mpn_add_n(z, z, buf, 2);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_1_mul1)
/* x has 1 words, z has 3. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_1_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    z[1] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_1_shortmul)
/* x and y have 1 words, z has 2.
 * Put the low 2 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_1_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 1);
    z[1-1] += x[0]*y[1-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_1_mod)
/* x has 4 words. z and p have 1 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 1_mgy_decode */
static inline
void mpfq_fixmp_1_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[1+1], r[1];
    assert (p[1-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 2, p, 1);
    mpfq_copy(z, r, 1);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_mod) */

#if !defined(HAVE_native_mpfq_fixmp_1_invmod)
/* x, z, and p have 1 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_1_invmod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t a, b, u, v, fix;
      int t, lsh;
    
      a = *x;
      b = *p;
    
      if (a == 0 || a == b) {
        *z=0;
        return 0;
      }
      /* b must be odd and >a */
    
      fix = (b+1)>>1;
    
      assert (a < b);
    
      u = 1; v = 0; t = 0;
      
      /* compute u and t such that u*a_orig = 2^t mod b */
    
      /* we maintain:
       *    u*a_orig - (not kept)*b_orig = 2^t*a
       *    v*a_orig - (not kept)*b_orig = -2^t*b
       * a and b are both odd.
       * An update consists in reducing the largest by the smallest,
       * and then adjusting the valuation.  */
    
      lsh = mpfq_ctzl(a);
      a >>= lsh;
      t += lsh;
      v <<= lsh;
      do {
        do {
          b -= a; v += u;
          lsh = mpfq_ctzl(b);
          b >>= lsh;
          t += lsh;
          u <<= lsh;
        } while (a<b);
        if (a == b)
          break;
        do {
          a -= b; u += v;
          lsh = mpfq_ctzl(a);
          a >>= lsh;
          t += lsh;
          v <<= lsh;
        } while (b < a);
      } while (a != b);
      if (a != 1) {
        *z = a;
        return 0;
      }
      while (t>0) {
        mp_limb_t sig = u & 1UL;
        u >>= 1;
        if (sig)
          u += fix;
        --t;
      } 
      *z = u;
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_1_redc)
/* x has 4 words, z and p have 1 words.
 * only one word is read from invp.
 * Assuming R=W^2 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_1_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t t = x[0]*mip[0];
    mp_limb_t cy = mpfq_fixmp_1_addmul1(x, p, t);
    if (cy || (x[1]>=p[0])) {
        z[0] = x[1] - p[0];
    } else {
        z[0] = x[1];
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_redc) */

#if !defined(HAVE_native_mpfq_fixmp_1_redc_ur)
/* x has 5 words, z and p have 1 words.
 * only one word is read from invp.
 * Assuming R=W^2 is the redc modulus, we expect that x verifies:
 *  x < W*W^1*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_1_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[2];
    for (int i = 0; i < 1; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_1_addmul1(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy=mpfq_fixmp_1_add(x+1+1, x+1+1, x);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W+1)*p
    */
    if (cy) {
        /* x'/R-p < W*p, which fits in n+1 words */
        mpn_sub(x+1,x+1,1+1,p,1);
    }
    mpn_tdiv_qr(q, z, 0, x+1, 1+1, p, 1);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_1_mgy_encode)
/* x, z, and p have 1 words.
 * Assuming R=W^2 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_1_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[2] = { 0, x[0] };
    mpfq_fixmp_1_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_1_mgy_decode)
/* x, z, invR, and p have 1 words.
 * Assuming R=W^2 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_1_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[2];
    mpfq_fixmp_1_mul(t, x, invR);
    mpfq_fixmp_1_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_1_lshift)
/* a has 1 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_1_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_1_rshift)
/* a has 1 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_rshift */
static inline
void mpfq_fixmp_1_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    a[1-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_1_long_lshift)
/* a has 1 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_long_lshift */
static inline
void mpfq_fixmp_1_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        a[0] <<= cnt;
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_1_long_rshift)
/* a has 1 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_long_rshift */
static inline
void mpfq_fixmp_1_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        a[0] >>= cnt;
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_2_add)
/* x, y, and z have 2 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 2_invmod, 2_redc, 2_redc_ur */
static inline
mp_limb_t mpfq_fixmp_2_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_add) */

#if !defined(HAVE_native_mpfq_fixmp_2_sub)
/* x, y, and z have 2 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 2_invmod, 2_redc */
static inline
mp_limb_t mpfq_fixmp_2_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_sub) */

#if !defined(HAVE_native_mpfq_fixmp_2_add_nc)
/* x, y, and z have 2 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_2_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_2_sub_nc)
/* x, y, and z have 2 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_2_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_2_add_ui)
/* x, y, and z have 2 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_2_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_2_sub_ui)
/* x, y, and z have 2 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_2_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_2_add_ui_nc)
/* x, y, and z have 2 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_2_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_2_sub_ui_nc)
/* x, y, and z have 2 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_2_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_2_cmp)
/* x and y have 2 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 2_invmod, 2_redc */
static inline
int mpfq_fixmp_2_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 2-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_2_cmp_ui)
/* x has 2 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 2_invmod */
static inline
int mpfq_fixmp_2_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 2-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_2_addmul1)
/* x has 2 words, z has 4.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
/* Triggered by: 2_redc_ur */
static inline
mp_limb_t mpfq_fixmp_2_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    z[2] += carry;
    return (z[2]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_2_addmul1_nc)
/* x has 2 words, z has 4.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 2_mul, 2_mgy_decode, 3_sqr, 3_shortmul, 4_sqr, 4_shortmul, 5_sqr, 5_shortmul, 6_sqr, 6_shortmul, 7_sqr, 7_shortmul, 8_sqr, 8_shortmul, 9_sqr, 9_shortmul, 10_sqr, 10_shortmul, 11_sqr, 11_shortmul, 12_sqr, 12_shortmul, 13_sqr, 13_shortmul, 14_sqr, 14_shortmul, 15_sqr, 15_shortmul, 2_5_sqr, 2_5_shortmul, 3_5_sqr, 3_5_shortmul, 4_5_sqr, 4_5_shortmul, 5_5_sqr, 5_5_shortmul, 6_5_sqr, 6_5_shortmul, 7_5_sqr, 7_5_shortmul, 8_5_sqr, 8_5_shortmul, 9_5_sqr, 9_5_shortmul, 10_5_sqr, 10_5_shortmul, 11_5_sqr, 11_5_shortmul, 12_5_sqr, 12_5_shortmul, 13_5_sqr, 13_5_shortmul, 14_5_sqr, 14_5_shortmul */
static inline
void mpfq_fixmp_2_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    z[2] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_2_addmul1_shortz)
/* x has 2 words, z has 3.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 2_redc */
static inline
mp_limb_t mpfq_fixmp_2_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_2_mul)
/* x and y have 2 words, z has 6. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 2_mgy_decode */
static inline
void mpfq_fixmp_2_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 4; z[i++] = 0) ;
    mpfq_fixmp_2_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_2_addmul1_nc (z + 1, x, y[1]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_mul) */

#if !defined(HAVE_native_mpfq_fixmp_2_sqr)
/* x has 2 words, z has 6. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_2_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[4] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpn_lshift(buf, buf, 4, 1);
    mpn_add_n(z, z, buf, 4);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_2_mul1)
/* x has 2 words, z has 4. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_2_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    z[2] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_2_shortmul)
/* x and y have 2 words, z has 3.
 * Put the low 3 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_2_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 2);
    mpfq_fixmp_1_addmul1_nc (z+0, x, y[0]);
    z[2-1] += x[1]*y[0];
    z[2-1] += x[0]*y[2-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_2_mod)
/* x has 6 words. z and p have 2 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 2_mgy_decode */
static inline
void mpfq_fixmp_2_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[2+1], r[2];
    assert (p[2-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 4, p, 2);
    mpfq_copy(z, r, 2);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_mod) */

#if !defined(HAVE_native_mpfq_fixmp_2_rshift)
/* a has 2 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 2_invmod */
static inline
void mpfq_fixmp_2_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 2-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[2-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_2_long_rshift)
/* a has 2 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 2_invmod */
static inline
void mpfq_fixmp_2_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 2 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[2-off-1] = a[2-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 2 - off);
    }
    mpfq_zero(a + 2 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_2_long_lshift)
/* a has 2 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 2_invmod */
static inline
void mpfq_fixmp_2_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 2-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 2 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_2_invmod)
/* x, z, and p have 2 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_2_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[2], v[2], a[2], b[2], fix[2];
      int i, t, lsh;
    
      mpfq_zero(u, 2);
      mpfq_zero(v, 2);
      mpfq_copy(a, x, 2);
      mpfq_copy(b, p, 2);
      u[0] = 1UL;
      
      if (mpfq_fixmp_2_cmp(a, v) == 0 || mpfq_fixmp_2_cmp(a, b) == 0) {
        mpfq_zero(res, 2);
        return 0;
      }
    
      mpfq_fixmp_2_add(fix, b, u);
      mpfq_fixmp_2_rshift(fix, 1);
    
      assert (mpfq_fixmp_2_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 2);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_2_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_2_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_2_sub(b, b, a);
          mpfq_fixmp_2_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 2);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_2_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_2_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_2_cmp(a,b) < 0);
        if (mpfq_fixmp_2_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_2_sub(a, a, b);
          mpfq_fixmp_2_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 2);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_2_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_2_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_2_cmp(b,a)<0);
      } while (mpfq_fixmp_2_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_2_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 2);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_2_rshift(u, 1);
        if (sig)
          mpfq_fixmp_2_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 2);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_2_redc)
/* x has 6 words, z and p have 2 words.
 * only one word is read from invp.
 * Assuming R=W^3 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_2_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 2; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_2_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_2_add(x, x, x + 2);
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_2_cmp(x, p) >= 0) {
        mpfq_fixmp_2_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 2);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_redc) */

#if !defined(HAVE_native_mpfq_fixmp_2_redc_ur)
/* x has 7 words, z and p have 2 words.
 * only one word is read from invp.
 * Assuming R=W^3 is the redc modulus, we expect that x verifies:
 *  x < W*W^2*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_2_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[2];
    for (int i = 0; i < 2; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_2_addmul1(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy=mpfq_fixmp_2_add(x+2+1, x+2+1, x);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W+1)*p
    */
    if (cy) {
        /* x'/R-p < W*p, which fits in n+1 words */
        mpn_sub(x+2,x+2,2+1,p,2);
    }
    mpn_tdiv_qr(q, z, 0, x+2, 2+1, p, 2);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_2_mgy_encode)
/* x, z, and p have 2 words.
 * Assuming R=W^3 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_2_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[4] = { 0, 0, x[0], x[1] };
    mpfq_fixmp_2_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_2_mgy_decode)
/* x, z, invR, and p have 2 words.
 * Assuming R=W^3 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_2_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[4];
    mpfq_fixmp_2_mul(t, x, invR);
    mpfq_fixmp_2_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_2_lshift)
/* a has 2 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_2_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 2-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_3_add)
/* x, y, and z have 3 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 3_invmod, 3_redc, 3_redc_ur */
static inline
mp_limb_t mpfq_fixmp_3_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_add) */

#if !defined(HAVE_native_mpfq_fixmp_3_sub)
/* x, y, and z have 3 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 3_invmod, 3_redc */
static inline
mp_limb_t mpfq_fixmp_3_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_sub) */

#if !defined(HAVE_native_mpfq_fixmp_3_add_nc)
/* x, y, and z have 3 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_3_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_3_sub_nc)
/* x, y, and z have 3 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_3_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_3_add_ui)
/* x, y, and z have 3 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_3_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_3_sub_ui)
/* x, y, and z have 3 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_3_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_3_add_ui_nc)
/* x, y, and z have 3 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_3_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_3_sub_ui_nc)
/* x, y, and z have 3 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_3_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_3_cmp)
/* x and y have 3 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 3_invmod, 3_redc */
static inline
int mpfq_fixmp_3_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 3-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_3_cmp_ui)
/* x has 3 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 3_invmod */
static inline
int mpfq_fixmp_3_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 3-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_3_addmul1)
/* x has 3 words, z has 5.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
/* Triggered by: 3_redc_ur */
static inline
mp_limb_t mpfq_fixmp_3_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    z[3] += carry;
    return (z[3]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_3_addmul1_nc)
/* x has 3 words, z has 5.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 3_mul, 3_mgy_decode, 4_sqr, 4_shortmul, 5_sqr, 5_shortmul, 6_sqr, 6_shortmul, 7_sqr, 7_shortmul, 8_sqr, 8_shortmul, 9_sqr, 9_shortmul, 10_sqr, 10_shortmul, 11_sqr, 11_shortmul, 12_sqr, 12_shortmul, 13_sqr, 13_shortmul, 14_sqr, 14_shortmul, 15_sqr, 15_shortmul, 3_5_sqr, 3_5_shortmul, 4_5_sqr, 4_5_shortmul, 5_5_sqr, 5_5_shortmul, 6_5_sqr, 6_5_shortmul, 7_5_sqr, 7_5_shortmul, 8_5_sqr, 8_5_shortmul, 9_5_sqr, 9_5_shortmul, 10_5_sqr, 10_5_shortmul, 11_5_sqr, 11_5_shortmul, 12_5_sqr, 12_5_shortmul, 13_5_sqr, 13_5_shortmul, 14_5_sqr, 14_5_shortmul */
static inline
void mpfq_fixmp_3_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    z[3] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_3_addmul1_shortz)
/* x has 3 words, z has 4.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 3_redc */
static inline
mp_limb_t mpfq_fixmp_3_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_3_mul)
/* x and y have 3 words, z has 8. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 3_mgy_decode */
static inline
void mpfq_fixmp_3_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 6; z[i++] = 0) ;
    mpfq_fixmp_3_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_3_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_3_addmul1_nc (z + 2, x, y[2]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_mul) */

#if !defined(HAVE_native_mpfq_fixmp_3_sqr)
/* x has 3 words, z has 8. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_3_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[6] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpfq_umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
    mpn_lshift(buf, buf, 6, 1);
    mpn_add_n(z, z, buf, 6);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_3_mul1)
/* x has 3 words, z has 5. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_3_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    z[3] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_3_shortmul)
/* x and y have 3 words, z has 4.
 * Put the low 4 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_3_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 3);
    mpfq_fixmp_2_addmul1_nc (z+0, x, y[0]);
    z[3-1] += x[2]*y[0];
    mpfq_fixmp_1_addmul1_nc (z+1, x, y[1]);
    z[3-1] += x[1]*y[1];
    z[3-1] += x[0]*y[3-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_3_mod)
/* x has 8 words. z and p have 3 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 3_mgy_decode */
static inline
void mpfq_fixmp_3_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[3+1], r[3];
    assert (p[3-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 6, p, 3);
    mpfq_copy(z, r, 3);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_mod) */

#if !defined(HAVE_native_mpfq_fixmp_3_rshift)
/* a has 3 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 3_invmod */
static inline
void mpfq_fixmp_3_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 3-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[3-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_3_long_rshift)
/* a has 3 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 3_invmod */
static inline
void mpfq_fixmp_3_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 3 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[3-off-1] = a[3-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 3 - off);
    }
    mpfq_zero(a + 3 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_3_long_lshift)
/* a has 3 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 3_invmod */
static inline
void mpfq_fixmp_3_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 3-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 3 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_3_invmod)
/* x, z, and p have 3 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_3_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[3], v[3], a[3], b[3], fix[3];
      int i, t, lsh;
    
      mpfq_zero(u, 3);
      mpfq_zero(v, 3);
      mpfq_copy(a, x, 3);
      mpfq_copy(b, p, 3);
      u[0] = 1UL;
      
      if (mpfq_fixmp_3_cmp(a, v) == 0 || mpfq_fixmp_3_cmp(a, b) == 0) {
        mpfq_zero(res, 3);
        return 0;
      }
    
      mpfq_fixmp_3_add(fix, b, u);
      mpfq_fixmp_3_rshift(fix, 1);
    
      assert (mpfq_fixmp_3_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 3);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_3_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_3_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_3_sub(b, b, a);
          mpfq_fixmp_3_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 3);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_3_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_3_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_3_cmp(a,b) < 0);
        if (mpfq_fixmp_3_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_3_sub(a, a, b);
          mpfq_fixmp_3_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 3);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_3_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_3_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_3_cmp(b,a)<0);
      } while (mpfq_fixmp_3_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_3_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 3);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_3_rshift(u, 1);
        if (sig)
          mpfq_fixmp_3_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 3);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_3_redc)
/* x has 8 words, z and p have 3 words.
 * only one word is read from invp.
 * Assuming R=W^4 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_3_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 3; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_3_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_3_add(x, x, x + 3);
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_3_cmp(x, p) >= 0) {
        mpfq_fixmp_3_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 3);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_redc) */

#if !defined(HAVE_native_mpfq_fixmp_3_redc_ur)
/* x has 9 words, z and p have 3 words.
 * only one word is read from invp.
 * Assuming R=W^4 is the redc modulus, we expect that x verifies:
 *  x < W*W^3*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_3_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[2];
    for (int i = 0; i < 3; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_3_addmul1(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy=mpfq_fixmp_3_add(x+3+1, x+3+1, x);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W+1)*p
    */
    if (cy) {
        /* x'/R-p < W*p, which fits in n+1 words */
        mpn_sub(x+3,x+3,3+1,p,3);
    }
    mpn_tdiv_qr(q, z, 0, x+3, 3+1, p, 3);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_3_mgy_encode)
/* x, z, and p have 3 words.
 * Assuming R=W^4 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_3_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[6] = { 0, 0, 0, x[0], x[1], x[2] };
    mpfq_fixmp_3_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_3_mgy_decode)
/* x, z, invR, and p have 3 words.
 * Assuming R=W^4 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_3_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[6];
    mpfq_fixmp_3_mul(t, x, invR);
    mpfq_fixmp_3_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_3_lshift)
/* a has 3 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_3_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 3-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_4_add)
/* x, y, and z have 4 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 4_invmod, 4_redc, 4_redc_ur */
static inline
mp_limb_t mpfq_fixmp_4_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_add) */

#if !defined(HAVE_native_mpfq_fixmp_4_sub)
/* x, y, and z have 4 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 4_invmod, 4_redc */
static inline
mp_limb_t mpfq_fixmp_4_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_sub) */

#if !defined(HAVE_native_mpfq_fixmp_4_add_nc)
/* x, y, and z have 4 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_4_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_4_sub_nc)
/* x, y, and z have 4 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_4_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_4_add_ui)
/* x, y, and z have 4 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_4_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_4_sub_ui)
/* x, y, and z have 4 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_4_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_4_add_ui_nc)
/* x, y, and z have 4 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_4_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_4_sub_ui_nc)
/* x, y, and z have 4 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_4_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_4_cmp)
/* x and y have 4 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 4_invmod, 4_redc */
static inline
int mpfq_fixmp_4_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 4-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_4_cmp_ui)
/* x has 4 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 4_invmod */
static inline
int mpfq_fixmp_4_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 4-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_4_addmul1)
/* x has 4 words, z has 6.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
/* Triggered by: 4_redc_ur */
static inline
mp_limb_t mpfq_fixmp_4_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    z[4] += carry;
    return (z[4]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_4_addmul1_nc)
/* x has 4 words, z has 6.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 4_mul, 4_mgy_decode, 5_sqr, 5_shortmul, 6_sqr, 6_shortmul, 7_sqr, 7_shortmul, 8_sqr, 8_shortmul, 9_sqr, 9_shortmul, 10_sqr, 10_shortmul, 11_sqr, 11_shortmul, 12_sqr, 12_shortmul, 13_sqr, 13_shortmul, 14_sqr, 14_shortmul, 15_sqr, 15_shortmul, 4_5_sqr, 4_5_shortmul, 5_5_sqr, 5_5_shortmul, 6_5_sqr, 6_5_shortmul, 7_5_sqr, 7_5_shortmul, 8_5_sqr, 8_5_shortmul, 9_5_sqr, 9_5_shortmul, 10_5_sqr, 10_5_shortmul, 11_5_sqr, 11_5_shortmul, 12_5_sqr, 12_5_shortmul, 13_5_sqr, 13_5_shortmul, 14_5_sqr, 14_5_shortmul */
static inline
void mpfq_fixmp_4_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    z[4] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_4_addmul1_shortz)
/* x has 4 words, z has 5.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 4_redc */
static inline
mp_limb_t mpfq_fixmp_4_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_4_mul)
/* x and y have 4 words, z has 10. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 4_mgy_decode */
static inline
void mpfq_fixmp_4_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 8; z[i++] = 0) ;
    mpfq_fixmp_4_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_4_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_4_addmul1_nc (z + 2, x, y[2]);
    mpfq_fixmp_4_addmul1_nc (z + 3, x, y[3]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_mul) */

#if !defined(HAVE_native_mpfq_fixmp_4_sqr)
/* x has 4 words, z has 10. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_4_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[8] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_fixmp_3_addmul1_nc(buf + 3, x, x[3]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpfq_umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
    mpfq_umul_ppmm(z[2*3+1], z[2*3], x[3], x[3]);
    mpn_lshift(buf, buf, 8, 1);
    mpn_add_n(z, z, buf, 8);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_4_mul1)
/* x has 4 words, z has 6. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_4_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    z[4] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_4_shortmul)
/* x and y have 4 words, z has 5.
 * Put the low 5 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_4_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 4);
    mpfq_fixmp_3_addmul1_nc (z+0, x, y[0]);
    z[4-1] += x[3]*y[0];
    mpfq_fixmp_2_addmul1_nc (z+1, x, y[1]);
    z[4-1] += x[2]*y[1];
    mpfq_fixmp_1_addmul1_nc (z+2, x, y[2]);
    z[4-1] += x[1]*y[2];
    z[4-1] += x[0]*y[4-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_4_mod)
/* x has 10 words. z and p have 4 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 4_mgy_decode */
static inline
void mpfq_fixmp_4_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[4+1], r[4];
    assert (p[4-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 8, p, 4);
    mpfq_copy(z, r, 4);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_mod) */

#if !defined(HAVE_native_mpfq_fixmp_4_rshift)
/* a has 4 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 4_invmod */
static inline
void mpfq_fixmp_4_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 4-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[4-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_4_long_rshift)
/* a has 4 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 4_invmod */
static inline
void mpfq_fixmp_4_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 4 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[4-off-1] = a[4-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 4 - off);
    }
    mpfq_zero(a + 4 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_4_long_lshift)
/* a has 4 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 4_invmod */
static inline
void mpfq_fixmp_4_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 4-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 4 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_4_invmod)
/* x, z, and p have 4 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_4_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[4], v[4], a[4], b[4], fix[4];
      int i, t, lsh;
    
      mpfq_zero(u, 4);
      mpfq_zero(v, 4);
      mpfq_copy(a, x, 4);
      mpfq_copy(b, p, 4);
      u[0] = 1UL;
      
      if (mpfq_fixmp_4_cmp(a, v) == 0 || mpfq_fixmp_4_cmp(a, b) == 0) {
        mpfq_zero(res, 4);
        return 0;
      }
    
      mpfq_fixmp_4_add(fix, b, u);
      mpfq_fixmp_4_rshift(fix, 1);
    
      assert (mpfq_fixmp_4_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 4);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_4_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_4_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_4_sub(b, b, a);
          mpfq_fixmp_4_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 4);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_4_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_4_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_4_cmp(a,b) < 0);
        if (mpfq_fixmp_4_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_4_sub(a, a, b);
          mpfq_fixmp_4_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 4);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_4_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_4_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_4_cmp(b,a)<0);
      } while (mpfq_fixmp_4_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_4_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 4);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_4_rshift(u, 1);
        if (sig)
          mpfq_fixmp_4_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 4);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_4_redc)
/* x has 10 words, z and p have 4 words.
 * only one word is read from invp.
 * Assuming R=W^5 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_4_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 4; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_4_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_4_add(x, x, x + 4);
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_4_cmp(x, p) >= 0) {
        mpfq_fixmp_4_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 4);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_redc) */

#if !defined(HAVE_native_mpfq_fixmp_4_redc_ur)
/* x has 11 words, z and p have 4 words.
 * only one word is read from invp.
 * Assuming R=W^5 is the redc modulus, we expect that x verifies:
 *  x < W*W^4*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_4_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[2];
    for (int i = 0; i < 4; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_4_addmul1(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy=mpfq_fixmp_4_add(x+4+1, x+4+1, x);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W+1)*p
    */
    if (cy) {
        /* x'/R-p < W*p, which fits in n+1 words */
        mpn_sub(x+4,x+4,4+1,p,4);
    }
    mpn_tdiv_qr(q, z, 0, x+4, 4+1, p, 4);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_4_mgy_encode)
/* x, z, and p have 4 words.
 * Assuming R=W^5 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_4_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[8] = { 0, 0, 0, 0, x[0], x[1], x[2], x[3] };
    mpfq_fixmp_4_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_4_mgy_decode)
/* x, z, invR, and p have 4 words.
 * Assuming R=W^5 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_4_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[8];
    mpfq_fixmp_4_mul(t, x, invR);
    mpfq_fixmp_4_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_4_lshift)
/* a has 4 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_4_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 4-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_5_add)
/* x, y, and z have 5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 5_invmod, 5_redc, 5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_add) */

#if !defined(HAVE_native_mpfq_fixmp_5_sub)
/* x, y, and z have 5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 5_invmod, 5_redc */
static inline
mp_limb_t mpfq_fixmp_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_sub) */

#if !defined(HAVE_native_mpfq_fixmp_5_add_nc)
/* x, y, and z have 5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_5_sub_nc)
/* x, y, and z have 5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_5_add_ui)
/* x, y, and z have 5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_5_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_5_sub_ui)
/* x, y, and z have 5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_5_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_5_add_ui_nc)
/* x, y, and z have 5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_5_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_5_sub_ui_nc)
/* x, y, and z have 5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_5_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_5_cmp)
/* x and y have 5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 5_invmod, 5_redc */
static inline
int mpfq_fixmp_5_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 5-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_5_cmp_ui)
/* x has 5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 5_invmod */
static inline
int mpfq_fixmp_5_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 5-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_5_addmul1)
/* x has 5 words, z has 7.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
/* Triggered by: 5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    z[5] += carry;
    return (z[5]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_5_addmul1_nc)
/* x has 5 words, z has 7.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 5_mul, 5_mgy_decode, 6_sqr, 6_shortmul, 7_sqr, 7_shortmul, 8_sqr, 8_shortmul, 9_sqr, 9_shortmul, 10_sqr, 10_shortmul, 11_sqr, 11_shortmul, 12_sqr, 12_shortmul, 13_sqr, 13_shortmul, 14_sqr, 14_shortmul, 15_sqr, 15_shortmul, 5_5_sqr, 5_5_shortmul, 6_5_sqr, 6_5_shortmul, 7_5_sqr, 7_5_shortmul, 8_5_sqr, 8_5_shortmul, 9_5_sqr, 9_5_shortmul, 10_5_sqr, 10_5_shortmul, 11_5_sqr, 11_5_shortmul, 12_5_sqr, 12_5_shortmul, 13_5_sqr, 13_5_shortmul, 14_5_sqr, 14_5_shortmul */
static inline
void mpfq_fixmp_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    z[5] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_5_addmul1_shortz)
/* x has 5 words, z has 6.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 5_redc */
static inline
mp_limb_t mpfq_fixmp_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_5_mul)
/* x and y have 5 words, z has 12. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 5_mgy_decode */
static inline
void mpfq_fixmp_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 10; z[i++] = 0) ;
    mpfq_fixmp_5_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_5_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_5_addmul1_nc (z + 2, x, y[2]);
    mpfq_fixmp_5_addmul1_nc (z + 3, x, y[3]);
    mpfq_fixmp_5_addmul1_nc (z + 4, x, y[4]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_mul) */

#if !defined(HAVE_native_mpfq_fixmp_5_sqr)
/* x has 5 words, z has 12. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[10] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_fixmp_3_addmul1_nc(buf + 3, x, x[3]);
    mpfq_fixmp_4_addmul1_nc(buf + 4, x, x[4]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpfq_umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
    mpfq_umul_ppmm(z[2*3+1], z[2*3], x[3], x[3]);
    mpfq_umul_ppmm(z[2*4+1], z[2*4], x[4], x[4]);
    mpn_lshift(buf, buf, 10, 1);
    mpn_add_n(z, z, buf, 10);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_5_mul1)
/* x has 5 words, z has 7. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    z[5] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_5_shortmul)
/* x and y have 5 words, z has 6.
 * Put the low 6 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_5_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 5);
    mpfq_fixmp_4_addmul1_nc (z+0, x, y[0]);
    z[5-1] += x[4]*y[0];
    mpfq_fixmp_3_addmul1_nc (z+1, x, y[1]);
    z[5-1] += x[3]*y[1];
    mpfq_fixmp_2_addmul1_nc (z+2, x, y[2]);
    z[5-1] += x[2]*y[2];
    mpfq_fixmp_1_addmul1_nc (z+3, x, y[3]);
    z[5-1] += x[1]*y[3];
    z[5-1] += x[0]*y[5-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_5_mod)
/* x has 12 words. z and p have 5 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 5_mgy_decode */
static inline
void mpfq_fixmp_5_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[5+1], r[5];
    assert (p[5-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 10, p, 5);
    mpfq_copy(z, r, 5);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_mod) */

#if !defined(HAVE_native_mpfq_fixmp_5_rshift)
/* a has 5 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 5_invmod */
static inline
void mpfq_fixmp_5_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 5-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[5-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_5_long_rshift)
/* a has 5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 5_invmod */
static inline
void mpfq_fixmp_5_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 5 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[5-off-1] = a[5-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 5 - off);
    }
    mpfq_zero(a + 5 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_5_long_lshift)
/* a has 5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 5_invmod */
static inline
void mpfq_fixmp_5_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 5-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 5 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_5_invmod)
/* x, z, and p have 5 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_5_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[5], v[5], a[5], b[5], fix[5];
      int i, t, lsh;
    
      mpfq_zero(u, 5);
      mpfq_zero(v, 5);
      mpfq_copy(a, x, 5);
      mpfq_copy(b, p, 5);
      u[0] = 1UL;
      
      if (mpfq_fixmp_5_cmp(a, v) == 0 || mpfq_fixmp_5_cmp(a, b) == 0) {
        mpfq_zero(res, 5);
        return 0;
      }
    
      mpfq_fixmp_5_add(fix, b, u);
      mpfq_fixmp_5_rshift(fix, 1);
    
      assert (mpfq_fixmp_5_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 5);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_5_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_5_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_5_sub(b, b, a);
          mpfq_fixmp_5_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 5);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_5_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_5_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_5_cmp(a,b) < 0);
        if (mpfq_fixmp_5_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_5_sub(a, a, b);
          mpfq_fixmp_5_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 5);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_5_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_5_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_5_cmp(b,a)<0);
      } while (mpfq_fixmp_5_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_5_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 5);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_5_rshift(u, 1);
        if (sig)
          mpfq_fixmp_5_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 5);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_5_redc)
/* x has 12 words, z and p have 5 words.
 * only one word is read from invp.
 * Assuming R=W^6 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_5_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 5; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_5_add(x, x, x + 5);
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_5_cmp(x, p) >= 0) {
        mpfq_fixmp_5_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 5);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_redc) */

#if !defined(HAVE_native_mpfq_fixmp_5_redc_ur)
/* x has 13 words, z and p have 5 words.
 * only one word is read from invp.
 * Assuming R=W^6 is the redc modulus, we expect that x verifies:
 *  x < W*W^5*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_5_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[2];
    for (int i = 0; i < 5; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_5_addmul1(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy=mpfq_fixmp_5_add(x+5+1, x+5+1, x);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W+1)*p
    */
    if (cy) {
        /* x'/R-p < W*p, which fits in n+1 words */
        mpn_sub(x+5,x+5,5+1,p,5);
    }
    mpn_tdiv_qr(q, z, 0, x+5, 5+1, p, 5);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_5_mgy_encode)
/* x, z, and p have 5 words.
 * Assuming R=W^6 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_5_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[10] = { 0, 0, 0, 0, 0, x[0], x[1], x[2], x[3], x[4] };
    mpfq_fixmp_5_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_5_mgy_decode)
/* x, z, invR, and p have 5 words.
 * Assuming R=W^6 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_5_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[10];
    mpfq_fixmp_5_mul(t, x, invR);
    mpfq_fixmp_5_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_5_lshift)
/* a has 5 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_5_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 5-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_6_add)
/* x, y, and z have 6 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 6_invmod, 6_redc, 6_redc_ur */
static inline
mp_limb_t mpfq_fixmp_6_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_add) */

#if !defined(HAVE_native_mpfq_fixmp_6_sub)
/* x, y, and z have 6 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 6_invmod, 6_redc */
static inline
mp_limb_t mpfq_fixmp_6_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_sub) */

#if !defined(HAVE_native_mpfq_fixmp_6_add_nc)
/* x, y, and z have 6 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_6_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_6_sub_nc)
/* x, y, and z have 6 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_6_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_6_add_ui)
/* x, y, and z have 6 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_6_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_6_sub_ui)
/* x, y, and z have 6 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_6_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_6_add_ui_nc)
/* x, y, and z have 6 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_6_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_6_sub_ui_nc)
/* x, y, and z have 6 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_6_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_6_cmp)
/* x and y have 6 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 6_invmod, 6_redc */
static inline
int mpfq_fixmp_6_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 6-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_6_cmp_ui)
/* x has 6 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 6_invmod */
static inline
int mpfq_fixmp_6_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 6-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_6_addmul1)
/* x has 6 words, z has 8.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
/* Triggered by: 6_redc_ur */
static inline
mp_limb_t mpfq_fixmp_6_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    z[6] += carry;
    return (z[6]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_6_addmul1_nc)
/* x has 6 words, z has 8.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 6_mul, 6_mgy_decode, 7_sqr, 7_shortmul, 8_sqr, 8_shortmul, 9_sqr, 9_shortmul, 10_sqr, 10_shortmul, 11_sqr, 11_shortmul, 12_sqr, 12_shortmul, 13_sqr, 13_shortmul, 14_sqr, 14_shortmul, 15_sqr, 15_shortmul, 6_5_sqr, 6_5_shortmul, 7_5_sqr, 7_5_shortmul, 8_5_sqr, 8_5_shortmul, 9_5_sqr, 9_5_shortmul, 10_5_sqr, 10_5_shortmul, 11_5_sqr, 11_5_shortmul, 12_5_sqr, 12_5_shortmul, 13_5_sqr, 13_5_shortmul, 14_5_sqr, 14_5_shortmul */
static inline
void mpfq_fixmp_6_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    z[6] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_6_addmul1_shortz)
/* x has 6 words, z has 7.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 6_redc */
static inline
mp_limb_t mpfq_fixmp_6_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_6_mul)
/* x and y have 6 words, z has 14. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 6_mgy_decode */
static inline
void mpfq_fixmp_6_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 12; z[i++] = 0) ;
    mpfq_fixmp_6_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_6_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_6_addmul1_nc (z + 2, x, y[2]);
    mpfq_fixmp_6_addmul1_nc (z + 3, x, y[3]);
    mpfq_fixmp_6_addmul1_nc (z + 4, x, y[4]);
    mpfq_fixmp_6_addmul1_nc (z + 5, x, y[5]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_mul) */

#if !defined(HAVE_native_mpfq_fixmp_6_sqr)
/* x has 6 words, z has 14. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_6_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[12] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_fixmp_3_addmul1_nc(buf + 3, x, x[3]);
    mpfq_fixmp_4_addmul1_nc(buf + 4, x, x[4]);
    mpfq_fixmp_5_addmul1_nc(buf + 5, x, x[5]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpfq_umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
    mpfq_umul_ppmm(z[2*3+1], z[2*3], x[3], x[3]);
    mpfq_umul_ppmm(z[2*4+1], z[2*4], x[4], x[4]);
    mpfq_umul_ppmm(z[2*5+1], z[2*5], x[5], x[5]);
    mpn_lshift(buf, buf, 12, 1);
    mpn_add_n(z, z, buf, 12);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_6_mul1)
/* x has 6 words, z has 8. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_6_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    z[6] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_6_shortmul)
/* x and y have 6 words, z has 7.
 * Put the low 7 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_6_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 6);
    mpfq_fixmp_5_addmul1_nc (z+0, x, y[0]);
    z[6-1] += x[5]*y[0];
    mpfq_fixmp_4_addmul1_nc (z+1, x, y[1]);
    z[6-1] += x[4]*y[1];
    mpfq_fixmp_3_addmul1_nc (z+2, x, y[2]);
    z[6-1] += x[3]*y[2];
    mpfq_fixmp_2_addmul1_nc (z+3, x, y[3]);
    z[6-1] += x[2]*y[3];
    mpfq_fixmp_1_addmul1_nc (z+4, x, y[4]);
    z[6-1] += x[1]*y[4];
    z[6-1] += x[0]*y[6-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_6_mod)
/* x has 14 words. z and p have 6 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 6_mgy_decode */
static inline
void mpfq_fixmp_6_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[6+1], r[6];
    assert (p[6-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 12, p, 6);
    mpfq_copy(z, r, 6);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_mod) */

#if !defined(HAVE_native_mpfq_fixmp_6_rshift)
/* a has 6 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 6_invmod */
static inline
void mpfq_fixmp_6_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 6-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[6-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_6_long_rshift)
/* a has 6 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 6_invmod */
static inline
void mpfq_fixmp_6_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 6 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[6-off-1] = a[6-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 6 - off);
    }
    mpfq_zero(a + 6 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_6_long_lshift)
/* a has 6 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 6_invmod */
static inline
void mpfq_fixmp_6_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 6-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 6 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_6_invmod)
/* x, z, and p have 6 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_6_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[6], v[6], a[6], b[6], fix[6];
      int i, t, lsh;
    
      mpfq_zero(u, 6);
      mpfq_zero(v, 6);
      mpfq_copy(a, x, 6);
      mpfq_copy(b, p, 6);
      u[0] = 1UL;
      
      if (mpfq_fixmp_6_cmp(a, v) == 0 || mpfq_fixmp_6_cmp(a, b) == 0) {
        mpfq_zero(res, 6);
        return 0;
      }
    
      mpfq_fixmp_6_add(fix, b, u);
      mpfq_fixmp_6_rshift(fix, 1);
    
      assert (mpfq_fixmp_6_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 6);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_6_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_6_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_6_sub(b, b, a);
          mpfq_fixmp_6_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 6);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_6_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_6_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_6_cmp(a,b) < 0);
        if (mpfq_fixmp_6_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_6_sub(a, a, b);
          mpfq_fixmp_6_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 6);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_6_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_6_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_6_cmp(b,a)<0);
      } while (mpfq_fixmp_6_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_6_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 6);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_6_rshift(u, 1);
        if (sig)
          mpfq_fixmp_6_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 6);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_6_redc)
/* x has 14 words, z and p have 6 words.
 * only one word is read from invp.
 * Assuming R=W^7 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_6_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 6; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_6_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_6_add(x, x, x + 6);
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_6_cmp(x, p) >= 0) {
        mpfq_fixmp_6_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 6);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_redc) */

#if !defined(HAVE_native_mpfq_fixmp_6_redc_ur)
/* x has 15 words, z and p have 6 words.
 * only one word is read from invp.
 * Assuming R=W^7 is the redc modulus, we expect that x verifies:
 *  x < W*W^6*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_6_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[2];
    for (int i = 0; i < 6; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_6_addmul1(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy=mpfq_fixmp_6_add(x+6+1, x+6+1, x);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W+1)*p
    */
    if (cy) {
        /* x'/R-p < W*p, which fits in n+1 words */
        mpn_sub(x+6,x+6,6+1,p,6);
    }
    mpn_tdiv_qr(q, z, 0, x+6, 6+1, p, 6);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_6_mgy_encode)
/* x, z, and p have 6 words.
 * Assuming R=W^7 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_6_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[12] = { 0, 0, 0, 0, 0, 0, x[0], x[1], x[2], x[3], x[4], x[5] };
    mpfq_fixmp_6_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_6_mgy_decode)
/* x, z, invR, and p have 6 words.
 * Assuming R=W^7 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_6_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[12];
    mpfq_fixmp_6_mul(t, x, invR);
    mpfq_fixmp_6_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_6_lshift)
/* a has 6 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_6_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 6-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_7_add)
/* x, y, and z have 7 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 7_invmod, 7_redc, 7_redc_ur */
static inline
mp_limb_t mpfq_fixmp_7_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_add) */

#if !defined(HAVE_native_mpfq_fixmp_7_sub)
/* x, y, and z have 7 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 7_invmod, 7_redc */
static inline
mp_limb_t mpfq_fixmp_7_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_sub) */

#if !defined(HAVE_native_mpfq_fixmp_7_add_nc)
/* x, y, and z have 7 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_7_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_7_sub_nc)
/* x, y, and z have 7 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_7_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_7_add_ui)
/* x, y, and z have 7 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_7_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_7_sub_ui)
/* x, y, and z have 7 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_7_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_7_add_ui_nc)
/* x, y, and z have 7 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_7_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_7_sub_ui_nc)
/* x, y, and z have 7 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_7_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_7_cmp)
/* x and y have 7 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 7_invmod, 7_redc */
static inline
int mpfq_fixmp_7_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 7-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_7_cmp_ui)
/* x has 7 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 7_invmod */
static inline
int mpfq_fixmp_7_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 7-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_7_addmul1)
/* x has 7 words, z has 9.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
/* Triggered by: 7_redc_ur */
static inline
mp_limb_t mpfq_fixmp_7_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    z[7] += carry;
    return (z[7]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_7_addmul1_nc)
/* x has 7 words, z has 9.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 7_mul, 7_mgy_decode, 8_sqr, 8_shortmul, 9_sqr, 9_shortmul, 10_sqr, 10_shortmul, 11_sqr, 11_shortmul, 12_sqr, 12_shortmul, 13_sqr, 13_shortmul, 14_sqr, 14_shortmul, 15_sqr, 15_shortmul, 7_5_sqr, 7_5_shortmul, 8_5_sqr, 8_5_shortmul, 9_5_sqr, 9_5_shortmul, 10_5_sqr, 10_5_shortmul, 11_5_sqr, 11_5_shortmul, 12_5_sqr, 12_5_shortmul, 13_5_sqr, 13_5_shortmul, 14_5_sqr, 14_5_shortmul */
static inline
void mpfq_fixmp_7_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    z[7] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_7_addmul1_shortz)
/* x has 7 words, z has 8.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 7_redc */
static inline
mp_limb_t mpfq_fixmp_7_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_7_mul)
/* x and y have 7 words, z has 16. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 7_mgy_decode */
static inline
void mpfq_fixmp_7_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 14; z[i++] = 0) ;
    mpfq_fixmp_7_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_7_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_7_addmul1_nc (z + 2, x, y[2]);
    mpfq_fixmp_7_addmul1_nc (z + 3, x, y[3]);
    mpfq_fixmp_7_addmul1_nc (z + 4, x, y[4]);
    mpfq_fixmp_7_addmul1_nc (z + 5, x, y[5]);
    mpfq_fixmp_7_addmul1_nc (z + 6, x, y[6]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_mul) */

#if !defined(HAVE_native_mpfq_fixmp_7_sqr)
/* x has 7 words, z has 16. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_7_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[14] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_fixmp_3_addmul1_nc(buf + 3, x, x[3]);
    mpfq_fixmp_4_addmul1_nc(buf + 4, x, x[4]);
    mpfq_fixmp_5_addmul1_nc(buf + 5, x, x[5]);
    mpfq_fixmp_6_addmul1_nc(buf + 6, x, x[6]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpfq_umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
    mpfq_umul_ppmm(z[2*3+1], z[2*3], x[3], x[3]);
    mpfq_umul_ppmm(z[2*4+1], z[2*4], x[4], x[4]);
    mpfq_umul_ppmm(z[2*5+1], z[2*5], x[5], x[5]);
    mpfq_umul_ppmm(z[2*6+1], z[2*6], x[6], x[6]);
    mpn_lshift(buf, buf, 14, 1);
    mpn_add_n(z, z, buf, 14);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_7_mul1)
/* x has 7 words, z has 9. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_7_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    z[7] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_7_shortmul)
/* x and y have 7 words, z has 8.
 * Put the low 8 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_7_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 7);
    mpfq_fixmp_6_addmul1_nc (z+0, x, y[0]);
    z[7-1] += x[6]*y[0];
    mpfq_fixmp_5_addmul1_nc (z+1, x, y[1]);
    z[7-1] += x[5]*y[1];
    mpfq_fixmp_4_addmul1_nc (z+2, x, y[2]);
    z[7-1] += x[4]*y[2];
    mpfq_fixmp_3_addmul1_nc (z+3, x, y[3]);
    z[7-1] += x[3]*y[3];
    mpfq_fixmp_2_addmul1_nc (z+4, x, y[4]);
    z[7-1] += x[2]*y[4];
    mpfq_fixmp_1_addmul1_nc (z+5, x, y[5]);
    z[7-1] += x[1]*y[5];
    z[7-1] += x[0]*y[7-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_7_mod)
/* x has 16 words. z and p have 7 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 7_mgy_decode */
static inline
void mpfq_fixmp_7_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[7+1], r[7];
    assert (p[7-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 14, p, 7);
    mpfq_copy(z, r, 7);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_mod) */

#if !defined(HAVE_native_mpfq_fixmp_7_rshift)
/* a has 7 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 7_invmod */
static inline
void mpfq_fixmp_7_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 7-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[7-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_7_long_rshift)
/* a has 7 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 7_invmod */
static inline
void mpfq_fixmp_7_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 7 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[7-off-1] = a[7-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 7 - off);
    }
    mpfq_zero(a + 7 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_7_long_lshift)
/* a has 7 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 7_invmod */
static inline
void mpfq_fixmp_7_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 7-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 7 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_7_invmod)
/* x, z, and p have 7 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_7_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[7], v[7], a[7], b[7], fix[7];
      int i, t, lsh;
    
      mpfq_zero(u, 7);
      mpfq_zero(v, 7);
      mpfq_copy(a, x, 7);
      mpfq_copy(b, p, 7);
      u[0] = 1UL;
      
      if (mpfq_fixmp_7_cmp(a, v) == 0 || mpfq_fixmp_7_cmp(a, b) == 0) {
        mpfq_zero(res, 7);
        return 0;
      }
    
      mpfq_fixmp_7_add(fix, b, u);
      mpfq_fixmp_7_rshift(fix, 1);
    
      assert (mpfq_fixmp_7_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 7);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_7_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_7_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_7_sub(b, b, a);
          mpfq_fixmp_7_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 7);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_7_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_7_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_7_cmp(a,b) < 0);
        if (mpfq_fixmp_7_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_7_sub(a, a, b);
          mpfq_fixmp_7_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 7);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_7_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_7_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_7_cmp(b,a)<0);
      } while (mpfq_fixmp_7_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_7_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 7);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_7_rshift(u, 1);
        if (sig)
          mpfq_fixmp_7_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 7);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_7_redc)
/* x has 16 words, z and p have 7 words.
 * only one word is read from invp.
 * Assuming R=W^8 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_7_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 7; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_7_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_7_add(x, x, x + 7);
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_7_cmp(x, p) >= 0) {
        mpfq_fixmp_7_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 7);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_redc) */

#if !defined(HAVE_native_mpfq_fixmp_7_redc_ur)
/* x has 17 words, z and p have 7 words.
 * only one word is read from invp.
 * Assuming R=W^8 is the redc modulus, we expect that x verifies:
 *  x < W*W^7*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_7_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[2];
    for (int i = 0; i < 7; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_7_addmul1(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy=mpfq_fixmp_7_add(x+7+1, x+7+1, x);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W+1)*p
    */
    if (cy) {
        /* x'/R-p < W*p, which fits in n+1 words */
        mpn_sub(x+7,x+7,7+1,p,7);
    }
    mpn_tdiv_qr(q, z, 0, x+7, 7+1, p, 7);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_7_mgy_encode)
/* x, z, and p have 7 words.
 * Assuming R=W^8 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_7_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[14] = { 0, 0, 0, 0, 0, 0, 0, x[0], x[1], x[2], x[3], x[4], x[5], x[6] };
    mpfq_fixmp_7_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_7_mgy_decode)
/* x, z, invR, and p have 7 words.
 * Assuming R=W^8 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_7_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[14];
    mpfq_fixmp_7_mul(t, x, invR);
    mpfq_fixmp_7_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_7_lshift)
/* a has 7 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_7_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 7-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_8_add)
/* x, y, and z have 8 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 8_invmod, 8_redc, 8_redc_ur */
static inline
mp_limb_t mpfq_fixmp_8_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_add) */

#if !defined(HAVE_native_mpfq_fixmp_8_sub)
/* x, y, and z have 8 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 8_invmod, 8_redc */
static inline
mp_limb_t mpfq_fixmp_8_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_sub) */

#if !defined(HAVE_native_mpfq_fixmp_8_add_nc)
/* x, y, and z have 8 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_8_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_8_sub_nc)
/* x, y, and z have 8 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_8_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_8_add_ui)
/* x, y, and z have 8 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_8_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_8_sub_ui)
/* x, y, and z have 8 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_8_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_8_add_ui_nc)
/* x, y, and z have 8 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_8_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_8_sub_ui_nc)
/* x, y, and z have 8 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_8_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_8_cmp)
/* x and y have 8 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 8_invmod, 8_redc */
static inline
int mpfq_fixmp_8_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 8-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_8_cmp_ui)
/* x has 8 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 8_invmod */
static inline
int mpfq_fixmp_8_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 8-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_8_addmul1)
/* x has 8 words, z has 10.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
/* Triggered by: 8_redc_ur */
static inline
mp_limb_t mpfq_fixmp_8_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    z[8] += carry;
    return (z[8]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_8_addmul1_nc)
/* x has 8 words, z has 10.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 8_mul, 8_mgy_decode, 9_sqr, 9_shortmul, 10_sqr, 10_shortmul, 11_sqr, 11_shortmul, 12_sqr, 12_shortmul, 13_sqr, 13_shortmul, 14_sqr, 14_shortmul, 15_sqr, 15_shortmul, 8_5_sqr, 8_5_shortmul, 9_5_sqr, 9_5_shortmul, 10_5_sqr, 10_5_shortmul, 11_5_sqr, 11_5_shortmul, 12_5_sqr, 12_5_shortmul, 13_5_sqr, 13_5_shortmul, 14_5_sqr, 14_5_shortmul */
static inline
void mpfq_fixmp_8_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    z[8] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_8_addmul1_shortz)
/* x has 8 words, z has 9.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 8_redc */
static inline
mp_limb_t mpfq_fixmp_8_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_8_mul)
/* x and y have 8 words, z has 18. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 8_mgy_decode */
static inline
void mpfq_fixmp_8_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 16; z[i++] = 0) ;
    mpfq_fixmp_8_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_8_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_8_addmul1_nc (z + 2, x, y[2]);
    mpfq_fixmp_8_addmul1_nc (z + 3, x, y[3]);
    mpfq_fixmp_8_addmul1_nc (z + 4, x, y[4]);
    mpfq_fixmp_8_addmul1_nc (z + 5, x, y[5]);
    mpfq_fixmp_8_addmul1_nc (z + 6, x, y[6]);
    mpfq_fixmp_8_addmul1_nc (z + 7, x, y[7]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_mul) */

#if !defined(HAVE_native_mpfq_fixmp_8_sqr)
/* x has 8 words, z has 18. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_8_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[16] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_fixmp_3_addmul1_nc(buf + 3, x, x[3]);
    mpfq_fixmp_4_addmul1_nc(buf + 4, x, x[4]);
    mpfq_fixmp_5_addmul1_nc(buf + 5, x, x[5]);
    mpfq_fixmp_6_addmul1_nc(buf + 6, x, x[6]);
    mpfq_fixmp_7_addmul1_nc(buf + 7, x, x[7]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpfq_umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
    mpfq_umul_ppmm(z[2*3+1], z[2*3], x[3], x[3]);
    mpfq_umul_ppmm(z[2*4+1], z[2*4], x[4], x[4]);
    mpfq_umul_ppmm(z[2*5+1], z[2*5], x[5], x[5]);
    mpfq_umul_ppmm(z[2*6+1], z[2*6], x[6], x[6]);
    mpfq_umul_ppmm(z[2*7+1], z[2*7], x[7], x[7]);
    mpn_lshift(buf, buf, 16, 1);
    mpn_add_n(z, z, buf, 16);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_8_mul1)
/* x has 8 words, z has 10. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_8_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[7] = lo;
    z[8] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_8_shortmul)
/* x and y have 8 words, z has 9.
 * Put the low 9 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_8_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 8);
    mpfq_fixmp_7_addmul1_nc (z+0, x, y[0]);
    z[8-1] += x[7]*y[0];
    mpfq_fixmp_6_addmul1_nc (z+1, x, y[1]);
    z[8-1] += x[6]*y[1];
    mpfq_fixmp_5_addmul1_nc (z+2, x, y[2]);
    z[8-1] += x[5]*y[2];
    mpfq_fixmp_4_addmul1_nc (z+3, x, y[3]);
    z[8-1] += x[4]*y[3];
    mpfq_fixmp_3_addmul1_nc (z+4, x, y[4]);
    z[8-1] += x[3]*y[4];
    mpfq_fixmp_2_addmul1_nc (z+5, x, y[5]);
    z[8-1] += x[2]*y[5];
    mpfq_fixmp_1_addmul1_nc (z+6, x, y[6]);
    z[8-1] += x[1]*y[6];
    z[8-1] += x[0]*y[8-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_8_mod)
/* x has 18 words. z and p have 8 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 8_mgy_decode */
static inline
void mpfq_fixmp_8_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[8+1], r[8];
    assert (p[8-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 16, p, 8);
    mpfq_copy(z, r, 8);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_mod) */

#if !defined(HAVE_native_mpfq_fixmp_8_rshift)
/* a has 8 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 8_invmod */
static inline
void mpfq_fixmp_8_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 8-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[8-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_8_long_rshift)
/* a has 8 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 8_invmod */
static inline
void mpfq_fixmp_8_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 8 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[8-off-1] = a[8-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 8 - off);
    }
    mpfq_zero(a + 8 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_8_long_lshift)
/* a has 8 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 8_invmod */
static inline
void mpfq_fixmp_8_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 8-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 8 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_8_invmod)
/* x, z, and p have 8 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_8_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[8], v[8], a[8], b[8], fix[8];
      int i, t, lsh;
    
      mpfq_zero(u, 8);
      mpfq_zero(v, 8);
      mpfq_copy(a, x, 8);
      mpfq_copy(b, p, 8);
      u[0] = 1UL;
      
      if (mpfq_fixmp_8_cmp(a, v) == 0 || mpfq_fixmp_8_cmp(a, b) == 0) {
        mpfq_zero(res, 8);
        return 0;
      }
    
      mpfq_fixmp_8_add(fix, b, u);
      mpfq_fixmp_8_rshift(fix, 1);
    
      assert (mpfq_fixmp_8_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 8);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_8_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_8_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_8_sub(b, b, a);
          mpfq_fixmp_8_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 8);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_8_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_8_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_8_cmp(a,b) < 0);
        if (mpfq_fixmp_8_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_8_sub(a, a, b);
          mpfq_fixmp_8_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 8);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_8_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_8_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_8_cmp(b,a)<0);
      } while (mpfq_fixmp_8_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_8_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 8);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_8_rshift(u, 1);
        if (sig)
          mpfq_fixmp_8_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 8);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_8_redc)
/* x has 18 words, z and p have 8 words.
 * only one word is read from invp.
 * Assuming R=W^9 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_8_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 8; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_8_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_8_add(x, x, x + 8);
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_8_cmp(x, p) >= 0) {
        mpfq_fixmp_8_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 8);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_redc) */

#if !defined(HAVE_native_mpfq_fixmp_8_redc_ur)
/* x has 19 words, z and p have 8 words.
 * only one word is read from invp.
 * Assuming R=W^9 is the redc modulus, we expect that x verifies:
 *  x < W*W^8*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_8_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[2];
    for (int i = 0; i < 8; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_8_addmul1(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy=mpfq_fixmp_8_add(x+8+1, x+8+1, x);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W+1)*p
    */
    if (cy) {
        /* x'/R-p < W*p, which fits in n+1 words */
        mpn_sub(x+8,x+8,8+1,p,8);
    }
    mpn_tdiv_qr(q, z, 0, x+8, 8+1, p, 8);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_8_mgy_encode)
/* x, z, and p have 8 words.
 * Assuming R=W^9 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_8_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[16] = { 0, 0, 0, 0, 0, 0, 0, 0, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7] };
    mpfq_fixmp_8_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_8_mgy_decode)
/* x, z, invR, and p have 8 words.
 * Assuming R=W^9 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_8_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[16];
    mpfq_fixmp_8_mul(t, x, invR);
    mpfq_fixmp_8_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_8_lshift)
/* a has 8 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_8_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 8-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_9_add)
/* x, y, and z have 9 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 9_invmod, 9_redc, 9_redc_ur */
static inline
mp_limb_t mpfq_fixmp_9_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_add) */

#if !defined(HAVE_native_mpfq_fixmp_9_sub)
/* x, y, and z have 9 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 9_invmod, 9_redc */
static inline
mp_limb_t mpfq_fixmp_9_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_sub) */

#if !defined(HAVE_native_mpfq_fixmp_9_add_nc)
/* x, y, and z have 9 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_9_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_9_sub_nc)
/* x, y, and z have 9 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_9_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_9_add_ui)
/* x, y, and z have 9 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_9_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_9_sub_ui)
/* x, y, and z have 9 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_9_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_9_add_ui_nc)
/* x, y, and z have 9 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_9_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_9_sub_ui_nc)
/* x, y, and z have 9 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_9_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_9_cmp)
/* x and y have 9 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 9_invmod, 9_redc */
static inline
int mpfq_fixmp_9_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 9-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_9_cmp_ui)
/* x has 9 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 9_invmod */
static inline
int mpfq_fixmp_9_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 9-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_9_addmul1)
/* x has 9 words, z has 11.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
/* Triggered by: 9_redc_ur */
static inline
mp_limb_t mpfq_fixmp_9_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    z[9] += carry;
    return (z[9]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_9_addmul1_nc)
/* x has 9 words, z has 11.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 9_mul, 9_mgy_decode, 10_sqr, 10_shortmul, 11_sqr, 11_shortmul, 12_sqr, 12_shortmul, 13_sqr, 13_shortmul, 14_sqr, 14_shortmul, 15_sqr, 15_shortmul, 9_5_sqr, 9_5_shortmul, 10_5_sqr, 10_5_shortmul, 11_5_sqr, 11_5_shortmul, 12_5_sqr, 12_5_shortmul, 13_5_sqr, 13_5_shortmul, 14_5_sqr, 14_5_shortmul */
static inline
void mpfq_fixmp_9_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    z[9] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_9_addmul1_shortz)
/* x has 9 words, z has 10.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 9_redc */
static inline
mp_limb_t mpfq_fixmp_9_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_9_mul)
/* x and y have 9 words, z has 20. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 9_mgy_decode */
static inline
void mpfq_fixmp_9_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 18; z[i++] = 0) ;
    mpfq_fixmp_9_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_9_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_9_addmul1_nc (z + 2, x, y[2]);
    mpfq_fixmp_9_addmul1_nc (z + 3, x, y[3]);
    mpfq_fixmp_9_addmul1_nc (z + 4, x, y[4]);
    mpfq_fixmp_9_addmul1_nc (z + 5, x, y[5]);
    mpfq_fixmp_9_addmul1_nc (z + 6, x, y[6]);
    mpfq_fixmp_9_addmul1_nc (z + 7, x, y[7]);
    mpfq_fixmp_9_addmul1_nc (z + 8, x, y[8]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_mul) */

#if !defined(HAVE_native_mpfq_fixmp_9_sqr)
/* x has 9 words, z has 20. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_9_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[18] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_fixmp_3_addmul1_nc(buf + 3, x, x[3]);
    mpfq_fixmp_4_addmul1_nc(buf + 4, x, x[4]);
    mpfq_fixmp_5_addmul1_nc(buf + 5, x, x[5]);
    mpfq_fixmp_6_addmul1_nc(buf + 6, x, x[6]);
    mpfq_fixmp_7_addmul1_nc(buf + 7, x, x[7]);
    mpfq_fixmp_8_addmul1_nc(buf + 8, x, x[8]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpfq_umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
    mpfq_umul_ppmm(z[2*3+1], z[2*3], x[3], x[3]);
    mpfq_umul_ppmm(z[2*4+1], z[2*4], x[4], x[4]);
    mpfq_umul_ppmm(z[2*5+1], z[2*5], x[5], x[5]);
    mpfq_umul_ppmm(z[2*6+1], z[2*6], x[6], x[6]);
    mpfq_umul_ppmm(z[2*7+1], z[2*7], x[7], x[7]);
    mpfq_umul_ppmm(z[2*8+1], z[2*8], x[8], x[8]);
    mpn_lshift(buf, buf, 18, 1);
    mpn_add_n(z, z, buf, 18);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_9_mul1)
/* x has 9 words, z has 11. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_9_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[8] = lo;
    z[9] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_9_shortmul)
/* x and y have 9 words, z has 10.
 * Put the low 10 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_9_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 9);
    mpfq_fixmp_8_addmul1_nc (z+0, x, y[0]);
    z[9-1] += x[8]*y[0];
    mpfq_fixmp_7_addmul1_nc (z+1, x, y[1]);
    z[9-1] += x[7]*y[1];
    mpfq_fixmp_6_addmul1_nc (z+2, x, y[2]);
    z[9-1] += x[6]*y[2];
    mpfq_fixmp_5_addmul1_nc (z+3, x, y[3]);
    z[9-1] += x[5]*y[3];
    mpfq_fixmp_4_addmul1_nc (z+4, x, y[4]);
    z[9-1] += x[4]*y[4];
    mpfq_fixmp_3_addmul1_nc (z+5, x, y[5]);
    z[9-1] += x[3]*y[5];
    mpfq_fixmp_2_addmul1_nc (z+6, x, y[6]);
    z[9-1] += x[2]*y[6];
    mpfq_fixmp_1_addmul1_nc (z+7, x, y[7]);
    z[9-1] += x[1]*y[7];
    z[9-1] += x[0]*y[9-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_9_mod)
/* x has 20 words. z and p have 9 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 9_mgy_decode */
static inline
void mpfq_fixmp_9_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[9+1], r[9];
    assert (p[9-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 18, p, 9);
    mpfq_copy(z, r, 9);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_mod) */

#if !defined(HAVE_native_mpfq_fixmp_9_rshift)
/* a has 9 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 9_invmod */
static inline
void mpfq_fixmp_9_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 9-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[9-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_9_long_rshift)
/* a has 9 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 9_invmod */
static inline
void mpfq_fixmp_9_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 9 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[9-off-1] = a[9-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 9 - off);
    }
    mpfq_zero(a + 9 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_9_long_lshift)
/* a has 9 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 9_invmod */
static inline
void mpfq_fixmp_9_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 9-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 9 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_9_invmod)
/* x, z, and p have 9 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_9_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[9], v[9], a[9], b[9], fix[9];
      int i, t, lsh;
    
      mpfq_zero(u, 9);
      mpfq_zero(v, 9);
      mpfq_copy(a, x, 9);
      mpfq_copy(b, p, 9);
      u[0] = 1UL;
      
      if (mpfq_fixmp_9_cmp(a, v) == 0 || mpfq_fixmp_9_cmp(a, b) == 0) {
        mpfq_zero(res, 9);
        return 0;
      }
    
      mpfq_fixmp_9_add(fix, b, u);
      mpfq_fixmp_9_rshift(fix, 1);
    
      assert (mpfq_fixmp_9_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 9);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_9_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_9_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_9_sub(b, b, a);
          mpfq_fixmp_9_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 9);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_9_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_9_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_9_cmp(a,b) < 0);
        if (mpfq_fixmp_9_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_9_sub(a, a, b);
          mpfq_fixmp_9_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 9);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_9_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_9_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_9_cmp(b,a)<0);
      } while (mpfq_fixmp_9_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_9_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 9);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_9_rshift(u, 1);
        if (sig)
          mpfq_fixmp_9_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 9);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_9_redc)
/* x has 20 words, z and p have 9 words.
 * only one word is read from invp.
 * Assuming R=W^10 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_9_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 9; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_9_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_9_add(x, x, x + 9);
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_9_cmp(x, p) >= 0) {
        mpfq_fixmp_9_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 9);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_redc) */

#if !defined(HAVE_native_mpfq_fixmp_9_redc_ur)
/* x has 21 words, z and p have 9 words.
 * only one word is read from invp.
 * Assuming R=W^10 is the redc modulus, we expect that x verifies:
 *  x < W*W^9*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_9_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[2];
    for (int i = 0; i < 9; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_9_addmul1(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy=mpfq_fixmp_9_add(x+9+1, x+9+1, x);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W+1)*p
    */
    if (cy) {
        /* x'/R-p < W*p, which fits in n+1 words */
        mpn_sub(x+9,x+9,9+1,p,9);
    }
    mpn_tdiv_qr(q, z, 0, x+9, 9+1, p, 9);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_9_mgy_encode)
/* x, z, and p have 9 words.
 * Assuming R=W^10 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_9_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[18] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8] };
    mpfq_fixmp_9_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_9_mgy_decode)
/* x, z, invR, and p have 9 words.
 * Assuming R=W^10 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_9_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[18];
    mpfq_fixmp_9_mul(t, x, invR);
    mpfq_fixmp_9_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_9_lshift)
/* a has 9 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_9_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 9-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_10_add)
/* x, y, and z have 10 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 10_invmod, 10_redc, 10_redc_ur */
static inline
mp_limb_t mpfq_fixmp_10_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r + y[9];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[9] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_add) */

#if !defined(HAVE_native_mpfq_fixmp_10_sub)
/* x, y, and z have 10 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 10_invmod, 10_redc */
static inline
mp_limb_t mpfq_fixmp_10_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r - y[9];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[9] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_sub) */

#if !defined(HAVE_native_mpfq_fixmp_10_add_nc)
/* x, y, and z have 10 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_10_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r + y[9];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[9] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_10_sub_nc)
/* x, y, and z have 10 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_10_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r - y[9];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[9] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_10_add_ui)
/* x, y, and z have 10 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_10_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    s = x[9];
    t = s + cy;
    cy = t < s;
    z[9] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_10_sub_ui)
/* x, y, and z have 10 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_10_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    s = x[9];
    t = s - cy;
    cy = t > s;
    z[9] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_10_add_ui_nc)
/* x, y, and z have 10 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_10_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    s = x[9];
    t = s + cy;
    cy = t < s;
    z[9] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_10_sub_ui_nc)
/* x, y, and z have 10 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_10_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    s = x[9];
    t = s - cy;
    cy = t > s;
    z[9] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_10_cmp)
/* x and y have 10 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 10_invmod, 10_redc */
static inline
int mpfq_fixmp_10_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 10-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_10_cmp_ui)
/* x has 10 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 10_invmod */
static inline
int mpfq_fixmp_10_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 10-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_10_addmul1)
/* x has 10 words, z has 12.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
/* Triggered by: 10_redc_ur */
static inline
mp_limb_t mpfq_fixmp_10_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    z[10] += carry;
    return (z[10]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_10_addmul1_nc)
/* x has 10 words, z has 12.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 10_mul, 10_mgy_decode, 11_sqr, 11_shortmul, 12_sqr, 12_shortmul, 13_sqr, 13_shortmul, 14_sqr, 14_shortmul, 15_sqr, 15_shortmul, 10_5_sqr, 10_5_shortmul, 11_5_sqr, 11_5_shortmul, 12_5_sqr, 12_5_shortmul, 13_5_sqr, 13_5_shortmul, 14_5_sqr, 14_5_shortmul */
static inline
void mpfq_fixmp_10_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    z[10] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_10_addmul1_shortz)
/* x has 10 words, z has 11.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 10_redc */
static inline
mp_limb_t mpfq_fixmp_10_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_10_mul)
/* x and y have 10 words, z has 22. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 10_mgy_decode */
static inline
void mpfq_fixmp_10_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 20; z[i++] = 0) ;
    mpfq_fixmp_10_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_10_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_10_addmul1_nc (z + 2, x, y[2]);
    mpfq_fixmp_10_addmul1_nc (z + 3, x, y[3]);
    mpfq_fixmp_10_addmul1_nc (z + 4, x, y[4]);
    mpfq_fixmp_10_addmul1_nc (z + 5, x, y[5]);
    mpfq_fixmp_10_addmul1_nc (z + 6, x, y[6]);
    mpfq_fixmp_10_addmul1_nc (z + 7, x, y[7]);
    mpfq_fixmp_10_addmul1_nc (z + 8, x, y[8]);
    mpfq_fixmp_10_addmul1_nc (z + 9, x, y[9]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_mul) */

#if !defined(HAVE_native_mpfq_fixmp_10_sqr)
/* x has 10 words, z has 22. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_10_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[20] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_fixmp_3_addmul1_nc(buf + 3, x, x[3]);
    mpfq_fixmp_4_addmul1_nc(buf + 4, x, x[4]);
    mpfq_fixmp_5_addmul1_nc(buf + 5, x, x[5]);
    mpfq_fixmp_6_addmul1_nc(buf + 6, x, x[6]);
    mpfq_fixmp_7_addmul1_nc(buf + 7, x, x[7]);
    mpfq_fixmp_8_addmul1_nc(buf + 8, x, x[8]);
    mpfq_fixmp_9_addmul1_nc(buf + 9, x, x[9]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpfq_umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
    mpfq_umul_ppmm(z[2*3+1], z[2*3], x[3], x[3]);
    mpfq_umul_ppmm(z[2*4+1], z[2*4], x[4], x[4]);
    mpfq_umul_ppmm(z[2*5+1], z[2*5], x[5], x[5]);
    mpfq_umul_ppmm(z[2*6+1], z[2*6], x[6], x[6]);
    mpfq_umul_ppmm(z[2*7+1], z[2*7], x[7], x[7]);
    mpfq_umul_ppmm(z[2*8+1], z[2*8], x[8], x[8]);
    mpfq_umul_ppmm(z[2*9+1], z[2*9], x[9], x[9]);
    mpn_lshift(buf, buf, 20, 1);
    mpn_add_n(z, z, buf, 20);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_10_mul1)
/* x has 10 words, z has 12. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_10_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[9] = lo;
    z[10] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_10_shortmul)
/* x and y have 10 words, z has 11.
 * Put the low 11 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_10_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 10);
    mpfq_fixmp_9_addmul1_nc (z+0, x, y[0]);
    z[10-1] += x[9]*y[0];
    mpfq_fixmp_8_addmul1_nc (z+1, x, y[1]);
    z[10-1] += x[8]*y[1];
    mpfq_fixmp_7_addmul1_nc (z+2, x, y[2]);
    z[10-1] += x[7]*y[2];
    mpfq_fixmp_6_addmul1_nc (z+3, x, y[3]);
    z[10-1] += x[6]*y[3];
    mpfq_fixmp_5_addmul1_nc (z+4, x, y[4]);
    z[10-1] += x[5]*y[4];
    mpfq_fixmp_4_addmul1_nc (z+5, x, y[5]);
    z[10-1] += x[4]*y[5];
    mpfq_fixmp_3_addmul1_nc (z+6, x, y[6]);
    z[10-1] += x[3]*y[6];
    mpfq_fixmp_2_addmul1_nc (z+7, x, y[7]);
    z[10-1] += x[2]*y[7];
    mpfq_fixmp_1_addmul1_nc (z+8, x, y[8]);
    z[10-1] += x[1]*y[8];
    z[10-1] += x[0]*y[10-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_10_mod)
/* x has 22 words. z and p have 10 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 10_mgy_decode */
static inline
void mpfq_fixmp_10_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[10+1], r[10];
    assert (p[10-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 20, p, 10);
    mpfq_copy(z, r, 10);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_mod) */

#if !defined(HAVE_native_mpfq_fixmp_10_rshift)
/* a has 10 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 10_invmod */
static inline
void mpfq_fixmp_10_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 10-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[10-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_10_long_rshift)
/* a has 10 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 10_invmod */
static inline
void mpfq_fixmp_10_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 10 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[10-off-1] = a[10-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 10 - off);
    }
    mpfq_zero(a + 10 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_10_long_lshift)
/* a has 10 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 10_invmod */
static inline
void mpfq_fixmp_10_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 10-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 10 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_10_invmod)
/* x, z, and p have 10 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_10_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[10], v[10], a[10], b[10], fix[10];
      int i, t, lsh;
    
      mpfq_zero(u, 10);
      mpfq_zero(v, 10);
      mpfq_copy(a, x, 10);
      mpfq_copy(b, p, 10);
      u[0] = 1UL;
      
      if (mpfq_fixmp_10_cmp(a, v) == 0 || mpfq_fixmp_10_cmp(a, b) == 0) {
        mpfq_zero(res, 10);
        return 0;
      }
    
      mpfq_fixmp_10_add(fix, b, u);
      mpfq_fixmp_10_rshift(fix, 1);
    
      assert (mpfq_fixmp_10_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 10);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_10_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_10_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_10_sub(b, b, a);
          mpfq_fixmp_10_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 10);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_10_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_10_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_10_cmp(a,b) < 0);
        if (mpfq_fixmp_10_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_10_sub(a, a, b);
          mpfq_fixmp_10_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 10);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_10_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_10_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_10_cmp(b,a)<0);
      } while (mpfq_fixmp_10_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_10_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 10);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_10_rshift(u, 1);
        if (sig)
          mpfq_fixmp_10_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 10);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_10_redc)
/* x has 22 words, z and p have 10 words.
 * only one word is read from invp.
 * Assuming R=W^11 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_10_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 10; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_10_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_10_add(x, x, x + 10);
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_10_cmp(x, p) >= 0) {
        mpfq_fixmp_10_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 10);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_redc) */

#if !defined(HAVE_native_mpfq_fixmp_10_redc_ur)
/* x has 23 words, z and p have 10 words.
 * only one word is read from invp.
 * Assuming R=W^11 is the redc modulus, we expect that x verifies:
 *  x < W*W^10*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_10_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[2];
    for (int i = 0; i < 10; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_10_addmul1(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy=mpfq_fixmp_10_add(x+10+1, x+10+1, x);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W+1)*p
    */
    if (cy) {
        /* x'/R-p < W*p, which fits in n+1 words */
        mpn_sub(x+10,x+10,10+1,p,10);
    }
    mpn_tdiv_qr(q, z, 0, x+10, 10+1, p, 10);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_10_mgy_encode)
/* x, z, and p have 10 words.
 * Assuming R=W^11 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_10_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[20] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9] };
    mpfq_fixmp_10_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_10_mgy_decode)
/* x, z, invR, and p have 10 words.
 * Assuming R=W^11 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_10_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[20];
    mpfq_fixmp_10_mul(t, x, invR);
    mpfq_fixmp_10_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_10_lshift)
/* a has 10 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_10_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 10-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_11_add)
/* x, y, and z have 11 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 11_invmod, 11_redc, 11_redc_ur */
static inline
mp_limb_t mpfq_fixmp_11_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r + y[9];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r + y[10];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[10] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_add) */

#if !defined(HAVE_native_mpfq_fixmp_11_sub)
/* x, y, and z have 11 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 11_invmod, 11_redc */
static inline
mp_limb_t mpfq_fixmp_11_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r - y[9];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r - y[10];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[10] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_sub) */

#if !defined(HAVE_native_mpfq_fixmp_11_add_nc)
/* x, y, and z have 11 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_11_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r + y[9];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r + y[10];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[10] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_11_sub_nc)
/* x, y, and z have 11 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_11_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r - y[9];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r - y[10];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[10] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_11_add_ui)
/* x, y, and z have 11 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_11_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    s = x[9];
    t = s + cy;
    cy = t < s;
    z[9] = t;
    s = x[10];
    t = s + cy;
    cy = t < s;
    z[10] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_11_sub_ui)
/* x, y, and z have 11 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_11_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    s = x[9];
    t = s - cy;
    cy = t > s;
    z[9] = t;
    s = x[10];
    t = s - cy;
    cy = t > s;
    z[10] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_11_add_ui_nc)
/* x, y, and z have 11 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_11_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    s = x[9];
    t = s + cy;
    cy = t < s;
    z[9] = t;
    s = x[10];
    t = s + cy;
    cy = t < s;
    z[10] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_11_sub_ui_nc)
/* x, y, and z have 11 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_11_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    s = x[9];
    t = s - cy;
    cy = t > s;
    z[9] = t;
    s = x[10];
    t = s - cy;
    cy = t > s;
    z[10] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_11_cmp)
/* x and y have 11 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 11_invmod, 11_redc */
static inline
int mpfq_fixmp_11_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 11-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_11_cmp_ui)
/* x has 11 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 11_invmod */
static inline
int mpfq_fixmp_11_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 11-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_11_addmul1)
/* x has 11 words, z has 13.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
/* Triggered by: 11_redc_ur */
static inline
mp_limb_t mpfq_fixmp_11_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    z[11] += carry;
    return (z[11]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_11_addmul1_nc)
/* x has 11 words, z has 13.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 11_mul, 11_mgy_decode, 12_sqr, 12_shortmul, 13_sqr, 13_shortmul, 14_sqr, 14_shortmul, 15_sqr, 15_shortmul, 11_5_sqr, 11_5_shortmul, 12_5_sqr, 12_5_shortmul, 13_5_sqr, 13_5_shortmul, 14_5_sqr, 14_5_shortmul */
static inline
void mpfq_fixmp_11_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    z[11] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_11_addmul1_shortz)
/* x has 11 words, z has 12.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 11_redc */
static inline
mp_limb_t mpfq_fixmp_11_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_11_mul)
/* x and y have 11 words, z has 24. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 11_mgy_decode */
static inline
void mpfq_fixmp_11_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 22; z[i++] = 0) ;
    mpfq_fixmp_11_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_11_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_11_addmul1_nc (z + 2, x, y[2]);
    mpfq_fixmp_11_addmul1_nc (z + 3, x, y[3]);
    mpfq_fixmp_11_addmul1_nc (z + 4, x, y[4]);
    mpfq_fixmp_11_addmul1_nc (z + 5, x, y[5]);
    mpfq_fixmp_11_addmul1_nc (z + 6, x, y[6]);
    mpfq_fixmp_11_addmul1_nc (z + 7, x, y[7]);
    mpfq_fixmp_11_addmul1_nc (z + 8, x, y[8]);
    mpfq_fixmp_11_addmul1_nc (z + 9, x, y[9]);
    mpfq_fixmp_11_addmul1_nc (z + 10, x, y[10]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_mul) */

#if !defined(HAVE_native_mpfq_fixmp_11_sqr)
/* x has 11 words, z has 24. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_11_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[22] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_fixmp_3_addmul1_nc(buf + 3, x, x[3]);
    mpfq_fixmp_4_addmul1_nc(buf + 4, x, x[4]);
    mpfq_fixmp_5_addmul1_nc(buf + 5, x, x[5]);
    mpfq_fixmp_6_addmul1_nc(buf + 6, x, x[6]);
    mpfq_fixmp_7_addmul1_nc(buf + 7, x, x[7]);
    mpfq_fixmp_8_addmul1_nc(buf + 8, x, x[8]);
    mpfq_fixmp_9_addmul1_nc(buf + 9, x, x[9]);
    mpfq_fixmp_10_addmul1_nc(buf + 10, x, x[10]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpfq_umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
    mpfq_umul_ppmm(z[2*3+1], z[2*3], x[3], x[3]);
    mpfq_umul_ppmm(z[2*4+1], z[2*4], x[4], x[4]);
    mpfq_umul_ppmm(z[2*5+1], z[2*5], x[5], x[5]);
    mpfq_umul_ppmm(z[2*6+1], z[2*6], x[6], x[6]);
    mpfq_umul_ppmm(z[2*7+1], z[2*7], x[7], x[7]);
    mpfq_umul_ppmm(z[2*8+1], z[2*8], x[8], x[8]);
    mpfq_umul_ppmm(z[2*9+1], z[2*9], x[9], x[9]);
    mpfq_umul_ppmm(z[2*10+1], z[2*10], x[10], x[10]);
    mpn_lshift(buf, buf, 22, 1);
    mpn_add_n(z, z, buf, 22);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_11_mul1)
/* x has 11 words, z has 13. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_11_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[10] = lo;
    z[11] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_11_shortmul)
/* x and y have 11 words, z has 12.
 * Put the low 12 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_11_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 11);
    mpfq_fixmp_10_addmul1_nc (z+0, x, y[0]);
    z[11-1] += x[10]*y[0];
    mpfq_fixmp_9_addmul1_nc (z+1, x, y[1]);
    z[11-1] += x[9]*y[1];
    mpfq_fixmp_8_addmul1_nc (z+2, x, y[2]);
    z[11-1] += x[8]*y[2];
    mpfq_fixmp_7_addmul1_nc (z+3, x, y[3]);
    z[11-1] += x[7]*y[3];
    mpfq_fixmp_6_addmul1_nc (z+4, x, y[4]);
    z[11-1] += x[6]*y[4];
    mpfq_fixmp_5_addmul1_nc (z+5, x, y[5]);
    z[11-1] += x[5]*y[5];
    mpfq_fixmp_4_addmul1_nc (z+6, x, y[6]);
    z[11-1] += x[4]*y[6];
    mpfq_fixmp_3_addmul1_nc (z+7, x, y[7]);
    z[11-1] += x[3]*y[7];
    mpfq_fixmp_2_addmul1_nc (z+8, x, y[8]);
    z[11-1] += x[2]*y[8];
    mpfq_fixmp_1_addmul1_nc (z+9, x, y[9]);
    z[11-1] += x[1]*y[9];
    z[11-1] += x[0]*y[11-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_11_mod)
/* x has 24 words. z and p have 11 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 11_mgy_decode */
static inline
void mpfq_fixmp_11_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[11+1], r[11];
    assert (p[11-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 22, p, 11);
    mpfq_copy(z, r, 11);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_mod) */

#if !defined(HAVE_native_mpfq_fixmp_11_rshift)
/* a has 11 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 11_invmod */
static inline
void mpfq_fixmp_11_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 11-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[11-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_11_long_rshift)
/* a has 11 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 11_invmod */
static inline
void mpfq_fixmp_11_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 11 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[11-off-1] = a[11-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 11 - off);
    }
    mpfq_zero(a + 11 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_11_long_lshift)
/* a has 11 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 11_invmod */
static inline
void mpfq_fixmp_11_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 11-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 11 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_11_invmod)
/* x, z, and p have 11 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_11_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[11], v[11], a[11], b[11], fix[11];
      int i, t, lsh;
    
      mpfq_zero(u, 11);
      mpfq_zero(v, 11);
      mpfq_copy(a, x, 11);
      mpfq_copy(b, p, 11);
      u[0] = 1UL;
      
      if (mpfq_fixmp_11_cmp(a, v) == 0 || mpfq_fixmp_11_cmp(a, b) == 0) {
        mpfq_zero(res, 11);
        return 0;
      }
    
      mpfq_fixmp_11_add(fix, b, u);
      mpfq_fixmp_11_rshift(fix, 1);
    
      assert (mpfq_fixmp_11_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 11);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_11_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_11_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_11_sub(b, b, a);
          mpfq_fixmp_11_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 11);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_11_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_11_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_11_cmp(a,b) < 0);
        if (mpfq_fixmp_11_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_11_sub(a, a, b);
          mpfq_fixmp_11_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 11);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_11_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_11_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_11_cmp(b,a)<0);
      } while (mpfq_fixmp_11_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_11_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 11);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_11_rshift(u, 1);
        if (sig)
          mpfq_fixmp_11_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 11);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_11_redc)
/* x has 24 words, z and p have 11 words.
 * only one word is read from invp.
 * Assuming R=W^12 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_11_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 11; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_11_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_11_add(x, x, x + 11);
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_11_cmp(x, p) >= 0) {
        mpfq_fixmp_11_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 11);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_redc) */

#if !defined(HAVE_native_mpfq_fixmp_11_redc_ur)
/* x has 25 words, z and p have 11 words.
 * only one word is read from invp.
 * Assuming R=W^12 is the redc modulus, we expect that x verifies:
 *  x < W*W^11*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_11_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[2];
    for (int i = 0; i < 11; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_11_addmul1(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy=mpfq_fixmp_11_add(x+11+1, x+11+1, x);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W+1)*p
    */
    if (cy) {
        /* x'/R-p < W*p, which fits in n+1 words */
        mpn_sub(x+11,x+11,11+1,p,11);
    }
    mpn_tdiv_qr(q, z, 0, x+11, 11+1, p, 11);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_11_mgy_encode)
/* x, z, and p have 11 words.
 * Assuming R=W^12 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_11_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[22] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10] };
    mpfq_fixmp_11_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_11_mgy_decode)
/* x, z, invR, and p have 11 words.
 * Assuming R=W^12 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_11_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[22];
    mpfq_fixmp_11_mul(t, x, invR);
    mpfq_fixmp_11_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_11_lshift)
/* a has 11 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_11_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 11-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_12_add)
/* x, y, and z have 12 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 12_invmod, 12_redc, 12_redc_ur */
static inline
mp_limb_t mpfq_fixmp_12_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r + y[9];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r + y[10];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r + y[11];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[11] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_add) */

#if !defined(HAVE_native_mpfq_fixmp_12_sub)
/* x, y, and z have 12 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 12_invmod, 12_redc */
static inline
mp_limb_t mpfq_fixmp_12_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r - y[9];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r - y[10];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r - y[11];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[11] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_sub) */

#if !defined(HAVE_native_mpfq_fixmp_12_add_nc)
/* x, y, and z have 12 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_12_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r + y[9];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r + y[10];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r + y[11];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[11] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_12_sub_nc)
/* x, y, and z have 12 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_12_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r - y[9];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r - y[10];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r - y[11];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[11] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_12_add_ui)
/* x, y, and z have 12 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_12_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    s = x[9];
    t = s + cy;
    cy = t < s;
    z[9] = t;
    s = x[10];
    t = s + cy;
    cy = t < s;
    z[10] = t;
    s = x[11];
    t = s + cy;
    cy = t < s;
    z[11] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_12_sub_ui)
/* x, y, and z have 12 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_12_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    s = x[9];
    t = s - cy;
    cy = t > s;
    z[9] = t;
    s = x[10];
    t = s - cy;
    cy = t > s;
    z[10] = t;
    s = x[11];
    t = s - cy;
    cy = t > s;
    z[11] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_12_add_ui_nc)
/* x, y, and z have 12 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_12_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    s = x[9];
    t = s + cy;
    cy = t < s;
    z[9] = t;
    s = x[10];
    t = s + cy;
    cy = t < s;
    z[10] = t;
    s = x[11];
    t = s + cy;
    cy = t < s;
    z[11] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_12_sub_ui_nc)
/* x, y, and z have 12 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_12_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    s = x[9];
    t = s - cy;
    cy = t > s;
    z[9] = t;
    s = x[10];
    t = s - cy;
    cy = t > s;
    z[10] = t;
    s = x[11];
    t = s - cy;
    cy = t > s;
    z[11] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_12_cmp)
/* x and y have 12 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 12_invmod, 12_redc */
static inline
int mpfq_fixmp_12_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 12-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_12_cmp_ui)
/* x has 12 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 12_invmod */
static inline
int mpfq_fixmp_12_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 12-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_12_addmul1)
/* x has 12 words, z has 14.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
/* Triggered by: 12_redc_ur */
static inline
mp_limb_t mpfq_fixmp_12_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    z[12] += carry;
    return (z[12]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_12_addmul1_nc)
/* x has 12 words, z has 14.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 12_mul, 12_mgy_decode, 13_sqr, 13_shortmul, 14_sqr, 14_shortmul, 15_sqr, 15_shortmul, 12_5_sqr, 12_5_shortmul, 13_5_sqr, 13_5_shortmul, 14_5_sqr, 14_5_shortmul */
static inline
void mpfq_fixmp_12_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    z[12] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_12_addmul1_shortz)
/* x has 12 words, z has 13.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 12_redc */
static inline
mp_limb_t mpfq_fixmp_12_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_12_mul)
/* x and y have 12 words, z has 26. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 12_mgy_decode */
static inline
void mpfq_fixmp_12_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 24; z[i++] = 0) ;
    mpfq_fixmp_12_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_12_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_12_addmul1_nc (z + 2, x, y[2]);
    mpfq_fixmp_12_addmul1_nc (z + 3, x, y[3]);
    mpfq_fixmp_12_addmul1_nc (z + 4, x, y[4]);
    mpfq_fixmp_12_addmul1_nc (z + 5, x, y[5]);
    mpfq_fixmp_12_addmul1_nc (z + 6, x, y[6]);
    mpfq_fixmp_12_addmul1_nc (z + 7, x, y[7]);
    mpfq_fixmp_12_addmul1_nc (z + 8, x, y[8]);
    mpfq_fixmp_12_addmul1_nc (z + 9, x, y[9]);
    mpfq_fixmp_12_addmul1_nc (z + 10, x, y[10]);
    mpfq_fixmp_12_addmul1_nc (z + 11, x, y[11]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_mul) */

#if !defined(HAVE_native_mpfq_fixmp_12_sqr)
/* x has 12 words, z has 26. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_12_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[24] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_fixmp_3_addmul1_nc(buf + 3, x, x[3]);
    mpfq_fixmp_4_addmul1_nc(buf + 4, x, x[4]);
    mpfq_fixmp_5_addmul1_nc(buf + 5, x, x[5]);
    mpfq_fixmp_6_addmul1_nc(buf + 6, x, x[6]);
    mpfq_fixmp_7_addmul1_nc(buf + 7, x, x[7]);
    mpfq_fixmp_8_addmul1_nc(buf + 8, x, x[8]);
    mpfq_fixmp_9_addmul1_nc(buf + 9, x, x[9]);
    mpfq_fixmp_10_addmul1_nc(buf + 10, x, x[10]);
    mpfq_fixmp_11_addmul1_nc(buf + 11, x, x[11]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpfq_umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
    mpfq_umul_ppmm(z[2*3+1], z[2*3], x[3], x[3]);
    mpfq_umul_ppmm(z[2*4+1], z[2*4], x[4], x[4]);
    mpfq_umul_ppmm(z[2*5+1], z[2*5], x[5], x[5]);
    mpfq_umul_ppmm(z[2*6+1], z[2*6], x[6], x[6]);
    mpfq_umul_ppmm(z[2*7+1], z[2*7], x[7], x[7]);
    mpfq_umul_ppmm(z[2*8+1], z[2*8], x[8], x[8]);
    mpfq_umul_ppmm(z[2*9+1], z[2*9], x[9], x[9]);
    mpfq_umul_ppmm(z[2*10+1], z[2*10], x[10], x[10]);
    mpfq_umul_ppmm(z[2*11+1], z[2*11], x[11], x[11]);
    mpn_lshift(buf, buf, 24, 1);
    mpn_add_n(z, z, buf, 24);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_12_mul1)
/* x has 12 words, z has 14. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_12_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[11] = lo;
    z[12] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_12_shortmul)
/* x and y have 12 words, z has 13.
 * Put the low 13 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_12_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 12);
    mpfq_fixmp_11_addmul1_nc (z+0, x, y[0]);
    z[12-1] += x[11]*y[0];
    mpfq_fixmp_10_addmul1_nc (z+1, x, y[1]);
    z[12-1] += x[10]*y[1];
    mpfq_fixmp_9_addmul1_nc (z+2, x, y[2]);
    z[12-1] += x[9]*y[2];
    mpfq_fixmp_8_addmul1_nc (z+3, x, y[3]);
    z[12-1] += x[8]*y[3];
    mpfq_fixmp_7_addmul1_nc (z+4, x, y[4]);
    z[12-1] += x[7]*y[4];
    mpfq_fixmp_6_addmul1_nc (z+5, x, y[5]);
    z[12-1] += x[6]*y[5];
    mpfq_fixmp_5_addmul1_nc (z+6, x, y[6]);
    z[12-1] += x[5]*y[6];
    mpfq_fixmp_4_addmul1_nc (z+7, x, y[7]);
    z[12-1] += x[4]*y[7];
    mpfq_fixmp_3_addmul1_nc (z+8, x, y[8]);
    z[12-1] += x[3]*y[8];
    mpfq_fixmp_2_addmul1_nc (z+9, x, y[9]);
    z[12-1] += x[2]*y[9];
    mpfq_fixmp_1_addmul1_nc (z+10, x, y[10]);
    z[12-1] += x[1]*y[10];
    z[12-1] += x[0]*y[12-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_12_mod)
/* x has 26 words. z and p have 12 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 12_mgy_decode */
static inline
void mpfq_fixmp_12_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[12+1], r[12];
    assert (p[12-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 24, p, 12);
    mpfq_copy(z, r, 12);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_mod) */

#if !defined(HAVE_native_mpfq_fixmp_12_rshift)
/* a has 12 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 12_invmod */
static inline
void mpfq_fixmp_12_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 12-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[12-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_12_long_rshift)
/* a has 12 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 12_invmod */
static inline
void mpfq_fixmp_12_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 12 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[12-off-1] = a[12-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 12 - off);
    }
    mpfq_zero(a + 12 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_12_long_lshift)
/* a has 12 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 12_invmod */
static inline
void mpfq_fixmp_12_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 12-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 12 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_12_invmod)
/* x, z, and p have 12 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_12_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[12], v[12], a[12], b[12], fix[12];
      int i, t, lsh;
    
      mpfq_zero(u, 12);
      mpfq_zero(v, 12);
      mpfq_copy(a, x, 12);
      mpfq_copy(b, p, 12);
      u[0] = 1UL;
      
      if (mpfq_fixmp_12_cmp(a, v) == 0 || mpfq_fixmp_12_cmp(a, b) == 0) {
        mpfq_zero(res, 12);
        return 0;
      }
    
      mpfq_fixmp_12_add(fix, b, u);
      mpfq_fixmp_12_rshift(fix, 1);
    
      assert (mpfq_fixmp_12_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 12);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_12_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_12_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_12_sub(b, b, a);
          mpfq_fixmp_12_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 12);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_12_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_12_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_12_cmp(a,b) < 0);
        if (mpfq_fixmp_12_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_12_sub(a, a, b);
          mpfq_fixmp_12_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 12);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_12_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_12_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_12_cmp(b,a)<0);
      } while (mpfq_fixmp_12_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_12_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 12);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_12_rshift(u, 1);
        if (sig)
          mpfq_fixmp_12_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 12);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_12_redc)
/* x has 26 words, z and p have 12 words.
 * only one word is read from invp.
 * Assuming R=W^13 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_12_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 12; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_12_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_12_add(x, x, x + 12);
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_12_cmp(x, p) >= 0) {
        mpfq_fixmp_12_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 12);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_redc) */

#if !defined(HAVE_native_mpfq_fixmp_12_redc_ur)
/* x has 27 words, z and p have 12 words.
 * only one word is read from invp.
 * Assuming R=W^13 is the redc modulus, we expect that x verifies:
 *  x < W*W^12*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_12_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[2];
    for (int i = 0; i < 12; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_12_addmul1(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy=mpfq_fixmp_12_add(x+12+1, x+12+1, x);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W+1)*p
    */
    if (cy) {
        /* x'/R-p < W*p, which fits in n+1 words */
        mpn_sub(x+12,x+12,12+1,p,12);
    }
    mpn_tdiv_qr(q, z, 0, x+12, 12+1, p, 12);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_12_mgy_encode)
/* x, z, and p have 12 words.
 * Assuming R=W^13 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_12_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[24] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11] };
    mpfq_fixmp_12_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_12_mgy_decode)
/* x, z, invR, and p have 12 words.
 * Assuming R=W^13 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_12_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[24];
    mpfq_fixmp_12_mul(t, x, invR);
    mpfq_fixmp_12_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_12_lshift)
/* a has 12 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_12_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 12-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_13_add)
/* x, y, and z have 13 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 13_invmod, 13_redc, 13_redc_ur */
static inline
mp_limb_t mpfq_fixmp_13_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r + y[9];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r + y[10];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r + y[11];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[11] = t;
    r = x[12];
    s = r + y[12];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[12] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_add) */

#if !defined(HAVE_native_mpfq_fixmp_13_sub)
/* x, y, and z have 13 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 13_invmod, 13_redc */
static inline
mp_limb_t mpfq_fixmp_13_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r - y[9];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r - y[10];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r - y[11];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[11] = t;
    r = x[12];
    s = r - y[12];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[12] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_sub) */

#if !defined(HAVE_native_mpfq_fixmp_13_add_nc)
/* x, y, and z have 13 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_13_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r + y[9];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r + y[10];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r + y[11];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[11] = t;
    r = x[12];
    s = r + y[12];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[12] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_13_sub_nc)
/* x, y, and z have 13 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_13_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r - y[9];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r - y[10];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r - y[11];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[11] = t;
    r = x[12];
    s = r - y[12];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[12] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_13_add_ui)
/* x, y, and z have 13 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_13_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    s = x[9];
    t = s + cy;
    cy = t < s;
    z[9] = t;
    s = x[10];
    t = s + cy;
    cy = t < s;
    z[10] = t;
    s = x[11];
    t = s + cy;
    cy = t < s;
    z[11] = t;
    s = x[12];
    t = s + cy;
    cy = t < s;
    z[12] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_13_sub_ui)
/* x, y, and z have 13 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_13_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    s = x[9];
    t = s - cy;
    cy = t > s;
    z[9] = t;
    s = x[10];
    t = s - cy;
    cy = t > s;
    z[10] = t;
    s = x[11];
    t = s - cy;
    cy = t > s;
    z[11] = t;
    s = x[12];
    t = s - cy;
    cy = t > s;
    z[12] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_13_add_ui_nc)
/* x, y, and z have 13 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_13_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    s = x[9];
    t = s + cy;
    cy = t < s;
    z[9] = t;
    s = x[10];
    t = s + cy;
    cy = t < s;
    z[10] = t;
    s = x[11];
    t = s + cy;
    cy = t < s;
    z[11] = t;
    s = x[12];
    t = s + cy;
    cy = t < s;
    z[12] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_13_sub_ui_nc)
/* x, y, and z have 13 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_13_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    s = x[9];
    t = s - cy;
    cy = t > s;
    z[9] = t;
    s = x[10];
    t = s - cy;
    cy = t > s;
    z[10] = t;
    s = x[11];
    t = s - cy;
    cy = t > s;
    z[11] = t;
    s = x[12];
    t = s - cy;
    cy = t > s;
    z[12] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_13_cmp)
/* x and y have 13 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 13_invmod, 13_redc */
static inline
int mpfq_fixmp_13_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 13-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_13_cmp_ui)
/* x has 13 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 13_invmod */
static inline
int mpfq_fixmp_13_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 13-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_13_addmul1)
/* x has 13 words, z has 15.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
/* Triggered by: 13_redc_ur */
static inline
mp_limb_t mpfq_fixmp_13_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[12];
    lo += buf;
    carry += (lo<buf);
    z[12] = lo;
    z[13] += carry;
    return (z[13]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_13_addmul1_nc)
/* x has 13 words, z has 15.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 13_mul, 13_mgy_decode, 14_sqr, 14_shortmul, 15_sqr, 15_shortmul, 13_5_sqr, 13_5_shortmul, 14_5_sqr, 14_5_shortmul */
static inline
void mpfq_fixmp_13_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[12];
    lo += buf;
    carry += (lo<buf);
    z[12] = lo;
    z[13] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_13_addmul1_shortz)
/* x has 13 words, z has 14.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 13_redc */
static inline
mp_limb_t mpfq_fixmp_13_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[12];
    lo += buf;
    carry += (lo<buf);
    z[12] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_13_mul)
/* x and y have 13 words, z has 28. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 13_mgy_decode */
static inline
void mpfq_fixmp_13_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 26; z[i++] = 0) ;
    mpfq_fixmp_13_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_13_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_13_addmul1_nc (z + 2, x, y[2]);
    mpfq_fixmp_13_addmul1_nc (z + 3, x, y[3]);
    mpfq_fixmp_13_addmul1_nc (z + 4, x, y[4]);
    mpfq_fixmp_13_addmul1_nc (z + 5, x, y[5]);
    mpfq_fixmp_13_addmul1_nc (z + 6, x, y[6]);
    mpfq_fixmp_13_addmul1_nc (z + 7, x, y[7]);
    mpfq_fixmp_13_addmul1_nc (z + 8, x, y[8]);
    mpfq_fixmp_13_addmul1_nc (z + 9, x, y[9]);
    mpfq_fixmp_13_addmul1_nc (z + 10, x, y[10]);
    mpfq_fixmp_13_addmul1_nc (z + 11, x, y[11]);
    mpfq_fixmp_13_addmul1_nc (z + 12, x, y[12]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_mul) */

#if !defined(HAVE_native_mpfq_fixmp_13_sqr)
/* x has 13 words, z has 28. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_13_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[26] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_fixmp_3_addmul1_nc(buf + 3, x, x[3]);
    mpfq_fixmp_4_addmul1_nc(buf + 4, x, x[4]);
    mpfq_fixmp_5_addmul1_nc(buf + 5, x, x[5]);
    mpfq_fixmp_6_addmul1_nc(buf + 6, x, x[6]);
    mpfq_fixmp_7_addmul1_nc(buf + 7, x, x[7]);
    mpfq_fixmp_8_addmul1_nc(buf + 8, x, x[8]);
    mpfq_fixmp_9_addmul1_nc(buf + 9, x, x[9]);
    mpfq_fixmp_10_addmul1_nc(buf + 10, x, x[10]);
    mpfq_fixmp_11_addmul1_nc(buf + 11, x, x[11]);
    mpfq_fixmp_12_addmul1_nc(buf + 12, x, x[12]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpfq_umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
    mpfq_umul_ppmm(z[2*3+1], z[2*3], x[3], x[3]);
    mpfq_umul_ppmm(z[2*4+1], z[2*4], x[4], x[4]);
    mpfq_umul_ppmm(z[2*5+1], z[2*5], x[5], x[5]);
    mpfq_umul_ppmm(z[2*6+1], z[2*6], x[6], x[6]);
    mpfq_umul_ppmm(z[2*7+1], z[2*7], x[7], x[7]);
    mpfq_umul_ppmm(z[2*8+1], z[2*8], x[8], x[8]);
    mpfq_umul_ppmm(z[2*9+1], z[2*9], x[9], x[9]);
    mpfq_umul_ppmm(z[2*10+1], z[2*10], x[10], x[10]);
    mpfq_umul_ppmm(z[2*11+1], z[2*11], x[11], x[11]);
    mpfq_umul_ppmm(z[2*12+1], z[2*12], x[12], x[12]);
    mpn_lshift(buf, buf, 26, 1);
    mpn_add_n(z, z, buf, 26);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_13_mul1)
/* x has 13 words, z has 15. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_13_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[12] = lo;
    z[13] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_13_shortmul)
/* x and y have 13 words, z has 14.
 * Put the low 14 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_13_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 13);
    mpfq_fixmp_12_addmul1_nc (z+0, x, y[0]);
    z[13-1] += x[12]*y[0];
    mpfq_fixmp_11_addmul1_nc (z+1, x, y[1]);
    z[13-1] += x[11]*y[1];
    mpfq_fixmp_10_addmul1_nc (z+2, x, y[2]);
    z[13-1] += x[10]*y[2];
    mpfq_fixmp_9_addmul1_nc (z+3, x, y[3]);
    z[13-1] += x[9]*y[3];
    mpfq_fixmp_8_addmul1_nc (z+4, x, y[4]);
    z[13-1] += x[8]*y[4];
    mpfq_fixmp_7_addmul1_nc (z+5, x, y[5]);
    z[13-1] += x[7]*y[5];
    mpfq_fixmp_6_addmul1_nc (z+6, x, y[6]);
    z[13-1] += x[6]*y[6];
    mpfq_fixmp_5_addmul1_nc (z+7, x, y[7]);
    z[13-1] += x[5]*y[7];
    mpfq_fixmp_4_addmul1_nc (z+8, x, y[8]);
    z[13-1] += x[4]*y[8];
    mpfq_fixmp_3_addmul1_nc (z+9, x, y[9]);
    z[13-1] += x[3]*y[9];
    mpfq_fixmp_2_addmul1_nc (z+10, x, y[10]);
    z[13-1] += x[2]*y[10];
    mpfq_fixmp_1_addmul1_nc (z+11, x, y[11]);
    z[13-1] += x[1]*y[11];
    z[13-1] += x[0]*y[13-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_13_mod)
/* x has 28 words. z and p have 13 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 13_mgy_decode */
static inline
void mpfq_fixmp_13_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[13+1], r[13];
    assert (p[13-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 26, p, 13);
    mpfq_copy(z, r, 13);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_mod) */

#if !defined(HAVE_native_mpfq_fixmp_13_rshift)
/* a has 13 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 13_invmod */
static inline
void mpfq_fixmp_13_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 13-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[13-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_13_long_rshift)
/* a has 13 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 13_invmod */
static inline
void mpfq_fixmp_13_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 13 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[13-off-1] = a[13-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 13 - off);
    }
    mpfq_zero(a + 13 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_13_long_lshift)
/* a has 13 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 13_invmod */
static inline
void mpfq_fixmp_13_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 13-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 13 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_13_invmod)
/* x, z, and p have 13 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_13_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[13], v[13], a[13], b[13], fix[13];
      int i, t, lsh;
    
      mpfq_zero(u, 13);
      mpfq_zero(v, 13);
      mpfq_copy(a, x, 13);
      mpfq_copy(b, p, 13);
      u[0] = 1UL;
      
      if (mpfq_fixmp_13_cmp(a, v) == 0 || mpfq_fixmp_13_cmp(a, b) == 0) {
        mpfq_zero(res, 13);
        return 0;
      }
    
      mpfq_fixmp_13_add(fix, b, u);
      mpfq_fixmp_13_rshift(fix, 1);
    
      assert (mpfq_fixmp_13_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 13);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_13_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_13_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_13_sub(b, b, a);
          mpfq_fixmp_13_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 13);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_13_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_13_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_13_cmp(a,b) < 0);
        if (mpfq_fixmp_13_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_13_sub(a, a, b);
          mpfq_fixmp_13_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 13);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_13_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_13_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_13_cmp(b,a)<0);
      } while (mpfq_fixmp_13_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_13_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 13);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_13_rshift(u, 1);
        if (sig)
          mpfq_fixmp_13_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 13);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_13_redc)
/* x has 28 words, z and p have 13 words.
 * only one word is read from invp.
 * Assuming R=W^14 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_13_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 13; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_13_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_13_add(x, x, x + 13);
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_13_cmp(x, p) >= 0) {
        mpfq_fixmp_13_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 13);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_redc) */

#if !defined(HAVE_native_mpfq_fixmp_13_redc_ur)
/* x has 29 words, z and p have 13 words.
 * only one word is read from invp.
 * Assuming R=W^14 is the redc modulus, we expect that x verifies:
 *  x < W*W^13*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_13_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[2];
    for (int i = 0; i < 13; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_13_addmul1(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy=mpfq_fixmp_13_add(x+13+1, x+13+1, x);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W+1)*p
    */
    if (cy) {
        /* x'/R-p < W*p, which fits in n+1 words */
        mpn_sub(x+13,x+13,13+1,p,13);
    }
    mpn_tdiv_qr(q, z, 0, x+13, 13+1, p, 13);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_13_mgy_encode)
/* x, z, and p have 13 words.
 * Assuming R=W^14 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_13_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[26] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12] };
    mpfq_fixmp_13_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_13_mgy_decode)
/* x, z, invR, and p have 13 words.
 * Assuming R=W^14 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_13_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[26];
    mpfq_fixmp_13_mul(t, x, invR);
    mpfq_fixmp_13_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_13_lshift)
/* a has 13 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_13_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 13-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_14_add)
/* x, y, and z have 14 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 14_invmod, 14_redc, 14_redc_ur */
static inline
mp_limb_t mpfq_fixmp_14_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r + y[9];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r + y[10];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r + y[11];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[11] = t;
    r = x[12];
    s = r + y[12];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[12] = t;
    r = x[13];
    s = r + y[13];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[13] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_add) */

#if !defined(HAVE_native_mpfq_fixmp_14_sub)
/* x, y, and z have 14 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 14_invmod, 14_redc */
static inline
mp_limb_t mpfq_fixmp_14_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r - y[9];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r - y[10];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r - y[11];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[11] = t;
    r = x[12];
    s = r - y[12];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[12] = t;
    r = x[13];
    s = r - y[13];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[13] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_sub) */

#if !defined(HAVE_native_mpfq_fixmp_14_add_nc)
/* x, y, and z have 14 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_14_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r + y[9];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r + y[10];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r + y[11];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[11] = t;
    r = x[12];
    s = r + y[12];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[12] = t;
    r = x[13];
    s = r + y[13];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[13] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_14_sub_nc)
/* x, y, and z have 14 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_14_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r - y[9];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r - y[10];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r - y[11];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[11] = t;
    r = x[12];
    s = r - y[12];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[12] = t;
    r = x[13];
    s = r - y[13];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[13] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_14_add_ui)
/* x, y, and z have 14 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_14_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    s = x[9];
    t = s + cy;
    cy = t < s;
    z[9] = t;
    s = x[10];
    t = s + cy;
    cy = t < s;
    z[10] = t;
    s = x[11];
    t = s + cy;
    cy = t < s;
    z[11] = t;
    s = x[12];
    t = s + cy;
    cy = t < s;
    z[12] = t;
    s = x[13];
    t = s + cy;
    cy = t < s;
    z[13] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_14_sub_ui)
/* x, y, and z have 14 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_14_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    s = x[9];
    t = s - cy;
    cy = t > s;
    z[9] = t;
    s = x[10];
    t = s - cy;
    cy = t > s;
    z[10] = t;
    s = x[11];
    t = s - cy;
    cy = t > s;
    z[11] = t;
    s = x[12];
    t = s - cy;
    cy = t > s;
    z[12] = t;
    s = x[13];
    t = s - cy;
    cy = t > s;
    z[13] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_14_add_ui_nc)
/* x, y, and z have 14 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_14_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    s = x[9];
    t = s + cy;
    cy = t < s;
    z[9] = t;
    s = x[10];
    t = s + cy;
    cy = t < s;
    z[10] = t;
    s = x[11];
    t = s + cy;
    cy = t < s;
    z[11] = t;
    s = x[12];
    t = s + cy;
    cy = t < s;
    z[12] = t;
    s = x[13];
    t = s + cy;
    cy = t < s;
    z[13] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_14_sub_ui_nc)
/* x, y, and z have 14 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_14_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    s = x[9];
    t = s - cy;
    cy = t > s;
    z[9] = t;
    s = x[10];
    t = s - cy;
    cy = t > s;
    z[10] = t;
    s = x[11];
    t = s - cy;
    cy = t > s;
    z[11] = t;
    s = x[12];
    t = s - cy;
    cy = t > s;
    z[12] = t;
    s = x[13];
    t = s - cy;
    cy = t > s;
    z[13] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_14_cmp)
/* x and y have 14 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 14_invmod, 14_redc */
static inline
int mpfq_fixmp_14_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 14-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_14_cmp_ui)
/* x has 14 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 14_invmod */
static inline
int mpfq_fixmp_14_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 14-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_14_addmul1)
/* x has 14 words, z has 16.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
/* Triggered by: 14_redc_ur */
static inline
mp_limb_t mpfq_fixmp_14_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[12];
    lo += buf;
    carry += (lo<buf);
    z[12] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[13]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[13];
    lo += buf;
    carry += (lo<buf);
    z[13] = lo;
    z[14] += carry;
    return (z[14]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_14_addmul1_nc)
/* x has 14 words, z has 16.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 14_mul, 14_mgy_decode, 15_sqr, 15_shortmul, 14_5_sqr, 14_5_shortmul */
static inline
void mpfq_fixmp_14_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[12];
    lo += buf;
    carry += (lo<buf);
    z[12] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[13]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[13];
    lo += buf;
    carry += (lo<buf);
    z[13] = lo;
    z[14] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_14_addmul1_shortz)
/* x has 14 words, z has 15.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 14_redc */
static inline
mp_limb_t mpfq_fixmp_14_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[12];
    lo += buf;
    carry += (lo<buf);
    z[12] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[13]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[13];
    lo += buf;
    carry += (lo<buf);
    z[13] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_14_mul)
/* x and y have 14 words, z has 30. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 14_mgy_decode */
static inline
void mpfq_fixmp_14_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 28; z[i++] = 0) ;
    mpfq_fixmp_14_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_14_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_14_addmul1_nc (z + 2, x, y[2]);
    mpfq_fixmp_14_addmul1_nc (z + 3, x, y[3]);
    mpfq_fixmp_14_addmul1_nc (z + 4, x, y[4]);
    mpfq_fixmp_14_addmul1_nc (z + 5, x, y[5]);
    mpfq_fixmp_14_addmul1_nc (z + 6, x, y[6]);
    mpfq_fixmp_14_addmul1_nc (z + 7, x, y[7]);
    mpfq_fixmp_14_addmul1_nc (z + 8, x, y[8]);
    mpfq_fixmp_14_addmul1_nc (z + 9, x, y[9]);
    mpfq_fixmp_14_addmul1_nc (z + 10, x, y[10]);
    mpfq_fixmp_14_addmul1_nc (z + 11, x, y[11]);
    mpfq_fixmp_14_addmul1_nc (z + 12, x, y[12]);
    mpfq_fixmp_14_addmul1_nc (z + 13, x, y[13]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_mul) */

#if !defined(HAVE_native_mpfq_fixmp_14_sqr)
/* x has 14 words, z has 30. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_14_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[28] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_fixmp_3_addmul1_nc(buf + 3, x, x[3]);
    mpfq_fixmp_4_addmul1_nc(buf + 4, x, x[4]);
    mpfq_fixmp_5_addmul1_nc(buf + 5, x, x[5]);
    mpfq_fixmp_6_addmul1_nc(buf + 6, x, x[6]);
    mpfq_fixmp_7_addmul1_nc(buf + 7, x, x[7]);
    mpfq_fixmp_8_addmul1_nc(buf + 8, x, x[8]);
    mpfq_fixmp_9_addmul1_nc(buf + 9, x, x[9]);
    mpfq_fixmp_10_addmul1_nc(buf + 10, x, x[10]);
    mpfq_fixmp_11_addmul1_nc(buf + 11, x, x[11]);
    mpfq_fixmp_12_addmul1_nc(buf + 12, x, x[12]);
    mpfq_fixmp_13_addmul1_nc(buf + 13, x, x[13]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpfq_umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
    mpfq_umul_ppmm(z[2*3+1], z[2*3], x[3], x[3]);
    mpfq_umul_ppmm(z[2*4+1], z[2*4], x[4], x[4]);
    mpfq_umul_ppmm(z[2*5+1], z[2*5], x[5], x[5]);
    mpfq_umul_ppmm(z[2*6+1], z[2*6], x[6], x[6]);
    mpfq_umul_ppmm(z[2*7+1], z[2*7], x[7], x[7]);
    mpfq_umul_ppmm(z[2*8+1], z[2*8], x[8], x[8]);
    mpfq_umul_ppmm(z[2*9+1], z[2*9], x[9], x[9]);
    mpfq_umul_ppmm(z[2*10+1], z[2*10], x[10], x[10]);
    mpfq_umul_ppmm(z[2*11+1], z[2*11], x[11], x[11]);
    mpfq_umul_ppmm(z[2*12+1], z[2*12], x[12], x[12]);
    mpfq_umul_ppmm(z[2*13+1], z[2*13], x[13], x[13]);
    mpn_lshift(buf, buf, 28, 1);
    mpn_add_n(z, z, buf, 28);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_14_mul1)
/* x has 14 words, z has 16. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_14_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[12] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[13]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[13] = lo;
    z[14] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_14_shortmul)
/* x and y have 14 words, z has 15.
 * Put the low 15 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_14_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 14);
    mpfq_fixmp_13_addmul1_nc (z+0, x, y[0]);
    z[14-1] += x[13]*y[0];
    mpfq_fixmp_12_addmul1_nc (z+1, x, y[1]);
    z[14-1] += x[12]*y[1];
    mpfq_fixmp_11_addmul1_nc (z+2, x, y[2]);
    z[14-1] += x[11]*y[2];
    mpfq_fixmp_10_addmul1_nc (z+3, x, y[3]);
    z[14-1] += x[10]*y[3];
    mpfq_fixmp_9_addmul1_nc (z+4, x, y[4]);
    z[14-1] += x[9]*y[4];
    mpfq_fixmp_8_addmul1_nc (z+5, x, y[5]);
    z[14-1] += x[8]*y[5];
    mpfq_fixmp_7_addmul1_nc (z+6, x, y[6]);
    z[14-1] += x[7]*y[6];
    mpfq_fixmp_6_addmul1_nc (z+7, x, y[7]);
    z[14-1] += x[6]*y[7];
    mpfq_fixmp_5_addmul1_nc (z+8, x, y[8]);
    z[14-1] += x[5]*y[8];
    mpfq_fixmp_4_addmul1_nc (z+9, x, y[9]);
    z[14-1] += x[4]*y[9];
    mpfq_fixmp_3_addmul1_nc (z+10, x, y[10]);
    z[14-1] += x[3]*y[10];
    mpfq_fixmp_2_addmul1_nc (z+11, x, y[11]);
    z[14-1] += x[2]*y[11];
    mpfq_fixmp_1_addmul1_nc (z+12, x, y[12]);
    z[14-1] += x[1]*y[12];
    z[14-1] += x[0]*y[14-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_14_mod)
/* x has 30 words. z and p have 14 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 14_mgy_decode */
static inline
void mpfq_fixmp_14_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[14+1], r[14];
    assert (p[14-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 28, p, 14);
    mpfq_copy(z, r, 14);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_mod) */

#if !defined(HAVE_native_mpfq_fixmp_14_rshift)
/* a has 14 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 14_invmod */
static inline
void mpfq_fixmp_14_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 14-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[14-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_14_long_rshift)
/* a has 14 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 14_invmod */
static inline
void mpfq_fixmp_14_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 14 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[14-off-1] = a[14-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 14 - off);
    }
    mpfq_zero(a + 14 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_14_long_lshift)
/* a has 14 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 14_invmod */
static inline
void mpfq_fixmp_14_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 14-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 14 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_14_invmod)
/* x, z, and p have 14 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_14_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[14], v[14], a[14], b[14], fix[14];
      int i, t, lsh;
    
      mpfq_zero(u, 14);
      mpfq_zero(v, 14);
      mpfq_copy(a, x, 14);
      mpfq_copy(b, p, 14);
      u[0] = 1UL;
      
      if (mpfq_fixmp_14_cmp(a, v) == 0 || mpfq_fixmp_14_cmp(a, b) == 0) {
        mpfq_zero(res, 14);
        return 0;
      }
    
      mpfq_fixmp_14_add(fix, b, u);
      mpfq_fixmp_14_rshift(fix, 1);
    
      assert (mpfq_fixmp_14_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 14);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_14_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_14_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_14_sub(b, b, a);
          mpfq_fixmp_14_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 14);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_14_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_14_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_14_cmp(a,b) < 0);
        if (mpfq_fixmp_14_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_14_sub(a, a, b);
          mpfq_fixmp_14_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 14);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_14_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_14_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_14_cmp(b,a)<0);
      } while (mpfq_fixmp_14_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_14_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 14);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_14_rshift(u, 1);
        if (sig)
          mpfq_fixmp_14_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 14);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_14_redc)
/* x has 30 words, z and p have 14 words.
 * only one word is read from invp.
 * Assuming R=W^15 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_14_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 14; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_14_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_14_add(x, x, x + 14);
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_14_cmp(x, p) >= 0) {
        mpfq_fixmp_14_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 14);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_redc) */

#if !defined(HAVE_native_mpfq_fixmp_14_redc_ur)
/* x has 31 words, z and p have 14 words.
 * only one word is read from invp.
 * Assuming R=W^15 is the redc modulus, we expect that x verifies:
 *  x < W*W^14*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_14_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[2];
    for (int i = 0; i < 14; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_14_addmul1(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy=mpfq_fixmp_14_add(x+14+1, x+14+1, x);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W+1)*p
    */
    if (cy) {
        /* x'/R-p < W*p, which fits in n+1 words */
        mpn_sub(x+14,x+14,14+1,p,14);
    }
    mpn_tdiv_qr(q, z, 0, x+14, 14+1, p, 14);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_14_mgy_encode)
/* x, z, and p have 14 words.
 * Assuming R=W^15 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_14_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[28] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13] };
    mpfq_fixmp_14_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_14_mgy_decode)
/* x, z, invR, and p have 14 words.
 * Assuming R=W^15 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_14_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[28];
    mpfq_fixmp_14_mul(t, x, invR);
    mpfq_fixmp_14_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_14_lshift)
/* a has 14 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_14_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 14-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_15_add)
/* x, y, and z have 15 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 15_invmod, 15_redc, 15_redc_ur */
static inline
mp_limb_t mpfq_fixmp_15_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r + y[9];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r + y[10];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r + y[11];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[11] = t;
    r = x[12];
    s = r + y[12];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[12] = t;
    r = x[13];
    s = r + y[13];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[13] = t;
    r = x[14];
    s = r + y[14];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[14] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_add) */

#if !defined(HAVE_native_mpfq_fixmp_15_sub)
/* x, y, and z have 15 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 15_invmod, 15_redc */
static inline
mp_limb_t mpfq_fixmp_15_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r - y[9];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r - y[10];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r - y[11];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[11] = t;
    r = x[12];
    s = r - y[12];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[12] = t;
    r = x[13];
    s = r - y[13];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[13] = t;
    r = x[14];
    s = r - y[14];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[14] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_sub) */

#if !defined(HAVE_native_mpfq_fixmp_15_add_nc)
/* x, y, and z have 15 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_15_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r + y[9];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r + y[10];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r + y[11];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[11] = t;
    r = x[12];
    s = r + y[12];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[12] = t;
    r = x[13];
    s = r + y[13];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[13] = t;
    r = x[14];
    s = r + y[14];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[14] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_15_sub_nc)
/* x, y, and z have 15 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_15_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r - y[9];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r - y[10];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r - y[11];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[11] = t;
    r = x[12];
    s = r - y[12];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[12] = t;
    r = x[13];
    s = r - y[13];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[13] = t;
    r = x[14];
    s = r - y[14];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[14] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_15_add_ui)
/* x, y, and z have 15 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_15_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    s = x[9];
    t = s + cy;
    cy = t < s;
    z[9] = t;
    s = x[10];
    t = s + cy;
    cy = t < s;
    z[10] = t;
    s = x[11];
    t = s + cy;
    cy = t < s;
    z[11] = t;
    s = x[12];
    t = s + cy;
    cy = t < s;
    z[12] = t;
    s = x[13];
    t = s + cy;
    cy = t < s;
    z[13] = t;
    s = x[14];
    t = s + cy;
    cy = t < s;
    z[14] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_15_sub_ui)
/* x, y, and z have 15 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_15_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    s = x[9];
    t = s - cy;
    cy = t > s;
    z[9] = t;
    s = x[10];
    t = s - cy;
    cy = t > s;
    z[10] = t;
    s = x[11];
    t = s - cy;
    cy = t > s;
    z[11] = t;
    s = x[12];
    t = s - cy;
    cy = t > s;
    z[12] = t;
    s = x[13];
    t = s - cy;
    cy = t > s;
    z[13] = t;
    s = x[14];
    t = s - cy;
    cy = t > s;
    z[14] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_15_add_ui_nc)
/* x, y, and z have 15 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_15_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    s = x[9];
    t = s + cy;
    cy = t < s;
    z[9] = t;
    s = x[10];
    t = s + cy;
    cy = t < s;
    z[10] = t;
    s = x[11];
    t = s + cy;
    cy = t < s;
    z[11] = t;
    s = x[12];
    t = s + cy;
    cy = t < s;
    z[12] = t;
    s = x[13];
    t = s + cy;
    cy = t < s;
    z[13] = t;
    s = x[14];
    t = s + cy;
    cy = t < s;
    z[14] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_15_sub_ui_nc)
/* x, y, and z have 15 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_15_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    s = x[9];
    t = s - cy;
    cy = t > s;
    z[9] = t;
    s = x[10];
    t = s - cy;
    cy = t > s;
    z[10] = t;
    s = x[11];
    t = s - cy;
    cy = t > s;
    z[11] = t;
    s = x[12];
    t = s - cy;
    cy = t > s;
    z[12] = t;
    s = x[13];
    t = s - cy;
    cy = t > s;
    z[13] = t;
    s = x[14];
    t = s - cy;
    cy = t > s;
    z[14] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_15_cmp)
/* x and y have 15 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 15_invmod, 15_redc */
static inline
int mpfq_fixmp_15_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 15-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_15_cmp_ui)
/* x has 15 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 15_invmod */
static inline
int mpfq_fixmp_15_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 15-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_15_addmul1)
/* x has 15 words, z has 17.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
/* Triggered by: 15_redc_ur */
static inline
mp_limb_t mpfq_fixmp_15_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[12];
    lo += buf;
    carry += (lo<buf);
    z[12] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[13]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[13];
    lo += buf;
    carry += (lo<buf);
    z[13] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[14]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[14];
    lo += buf;
    carry += (lo<buf);
    z[14] = lo;
    z[15] += carry;
    return (z[15]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_15_addmul1_nc)
/* x has 15 words, z has 17.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 15_mul, 15_mgy_decode */
static inline
void mpfq_fixmp_15_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[12];
    lo += buf;
    carry += (lo<buf);
    z[12] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[13]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[13];
    lo += buf;
    carry += (lo<buf);
    z[13] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[14]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[14];
    lo += buf;
    carry += (lo<buf);
    z[14] = lo;
    z[15] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_15_addmul1_shortz)
/* x has 15 words, z has 16.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 15_redc */
static inline
mp_limb_t mpfq_fixmp_15_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[12];
    lo += buf;
    carry += (lo<buf);
    z[12] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[13]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[13];
    lo += buf;
    carry += (lo<buf);
    z[13] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[14]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[14];
    lo += buf;
    carry += (lo<buf);
    z[14] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_15_mul)
/* x and y have 15 words, z has 32. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 15_mgy_decode */
static inline
void mpfq_fixmp_15_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 30; z[i++] = 0) ;
    mpfq_fixmp_15_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_15_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_15_addmul1_nc (z + 2, x, y[2]);
    mpfq_fixmp_15_addmul1_nc (z + 3, x, y[3]);
    mpfq_fixmp_15_addmul1_nc (z + 4, x, y[4]);
    mpfq_fixmp_15_addmul1_nc (z + 5, x, y[5]);
    mpfq_fixmp_15_addmul1_nc (z + 6, x, y[6]);
    mpfq_fixmp_15_addmul1_nc (z + 7, x, y[7]);
    mpfq_fixmp_15_addmul1_nc (z + 8, x, y[8]);
    mpfq_fixmp_15_addmul1_nc (z + 9, x, y[9]);
    mpfq_fixmp_15_addmul1_nc (z + 10, x, y[10]);
    mpfq_fixmp_15_addmul1_nc (z + 11, x, y[11]);
    mpfq_fixmp_15_addmul1_nc (z + 12, x, y[12]);
    mpfq_fixmp_15_addmul1_nc (z + 13, x, y[13]);
    mpfq_fixmp_15_addmul1_nc (z + 14, x, y[14]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_mul) */

#if !defined(HAVE_native_mpfq_fixmp_15_sqr)
/* x has 15 words, z has 32. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_15_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[30] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_fixmp_3_addmul1_nc(buf + 3, x, x[3]);
    mpfq_fixmp_4_addmul1_nc(buf + 4, x, x[4]);
    mpfq_fixmp_5_addmul1_nc(buf + 5, x, x[5]);
    mpfq_fixmp_6_addmul1_nc(buf + 6, x, x[6]);
    mpfq_fixmp_7_addmul1_nc(buf + 7, x, x[7]);
    mpfq_fixmp_8_addmul1_nc(buf + 8, x, x[8]);
    mpfq_fixmp_9_addmul1_nc(buf + 9, x, x[9]);
    mpfq_fixmp_10_addmul1_nc(buf + 10, x, x[10]);
    mpfq_fixmp_11_addmul1_nc(buf + 11, x, x[11]);
    mpfq_fixmp_12_addmul1_nc(buf + 12, x, x[12]);
    mpfq_fixmp_13_addmul1_nc(buf + 13, x, x[13]);
    mpfq_fixmp_14_addmul1_nc(buf + 14, x, x[14]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpfq_umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
    mpfq_umul_ppmm(z[2*3+1], z[2*3], x[3], x[3]);
    mpfq_umul_ppmm(z[2*4+1], z[2*4], x[4], x[4]);
    mpfq_umul_ppmm(z[2*5+1], z[2*5], x[5], x[5]);
    mpfq_umul_ppmm(z[2*6+1], z[2*6], x[6], x[6]);
    mpfq_umul_ppmm(z[2*7+1], z[2*7], x[7], x[7]);
    mpfq_umul_ppmm(z[2*8+1], z[2*8], x[8], x[8]);
    mpfq_umul_ppmm(z[2*9+1], z[2*9], x[9], x[9]);
    mpfq_umul_ppmm(z[2*10+1], z[2*10], x[10], x[10]);
    mpfq_umul_ppmm(z[2*11+1], z[2*11], x[11], x[11]);
    mpfq_umul_ppmm(z[2*12+1], z[2*12], x[12], x[12]);
    mpfq_umul_ppmm(z[2*13+1], z[2*13], x[13], x[13]);
    mpfq_umul_ppmm(z[2*14+1], z[2*14], x[14], x[14]);
    mpn_lshift(buf, buf, 30, 1);
    mpn_add_n(z, z, buf, 30);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_15_mul1)
/* x has 15 words, z has 17. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_15_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[12] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[13]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[13] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[14]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[14] = lo;
    z[15] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_15_shortmul)
/* x and y have 15 words, z has 16.
 * Put the low 16 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_15_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 15);
    mpfq_fixmp_14_addmul1_nc (z+0, x, y[0]);
    z[15-1] += x[14]*y[0];
    mpfq_fixmp_13_addmul1_nc (z+1, x, y[1]);
    z[15-1] += x[13]*y[1];
    mpfq_fixmp_12_addmul1_nc (z+2, x, y[2]);
    z[15-1] += x[12]*y[2];
    mpfq_fixmp_11_addmul1_nc (z+3, x, y[3]);
    z[15-1] += x[11]*y[3];
    mpfq_fixmp_10_addmul1_nc (z+4, x, y[4]);
    z[15-1] += x[10]*y[4];
    mpfq_fixmp_9_addmul1_nc (z+5, x, y[5]);
    z[15-1] += x[9]*y[5];
    mpfq_fixmp_8_addmul1_nc (z+6, x, y[6]);
    z[15-1] += x[8]*y[6];
    mpfq_fixmp_7_addmul1_nc (z+7, x, y[7]);
    z[15-1] += x[7]*y[7];
    mpfq_fixmp_6_addmul1_nc (z+8, x, y[8]);
    z[15-1] += x[6]*y[8];
    mpfq_fixmp_5_addmul1_nc (z+9, x, y[9]);
    z[15-1] += x[5]*y[9];
    mpfq_fixmp_4_addmul1_nc (z+10, x, y[10]);
    z[15-1] += x[4]*y[10];
    mpfq_fixmp_3_addmul1_nc (z+11, x, y[11]);
    z[15-1] += x[3]*y[11];
    mpfq_fixmp_2_addmul1_nc (z+12, x, y[12]);
    z[15-1] += x[2]*y[12];
    mpfq_fixmp_1_addmul1_nc (z+13, x, y[13]);
    z[15-1] += x[1]*y[13];
    z[15-1] += x[0]*y[15-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_15_mod)
/* x has 32 words. z and p have 15 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 15_mgy_decode */
static inline
void mpfq_fixmp_15_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[15+1], r[15];
    assert (p[15-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 30, p, 15);
    mpfq_copy(z, r, 15);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_mod) */

#if !defined(HAVE_native_mpfq_fixmp_15_rshift)
/* a has 15 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 15_invmod */
static inline
void mpfq_fixmp_15_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 15-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[15-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_15_long_rshift)
/* a has 15 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 15_invmod */
static inline
void mpfq_fixmp_15_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 15 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[15-off-1] = a[15-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 15 - off);
    }
    mpfq_zero(a + 15 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_15_long_lshift)
/* a has 15 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 15_invmod */
static inline
void mpfq_fixmp_15_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 15-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 15 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_15_invmod)
/* x, z, and p have 15 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_15_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[15], v[15], a[15], b[15], fix[15];
      int i, t, lsh;
    
      mpfq_zero(u, 15);
      mpfq_zero(v, 15);
      mpfq_copy(a, x, 15);
      mpfq_copy(b, p, 15);
      u[0] = 1UL;
      
      if (mpfq_fixmp_15_cmp(a, v) == 0 || mpfq_fixmp_15_cmp(a, b) == 0) {
        mpfq_zero(res, 15);
        return 0;
      }
    
      mpfq_fixmp_15_add(fix, b, u);
      mpfq_fixmp_15_rshift(fix, 1);
    
      assert (mpfq_fixmp_15_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 15);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_15_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_15_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_15_sub(b, b, a);
          mpfq_fixmp_15_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 15);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_15_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_15_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_15_cmp(a,b) < 0);
        if (mpfq_fixmp_15_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_15_sub(a, a, b);
          mpfq_fixmp_15_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 15);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_15_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_15_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_15_cmp(b,a)<0);
      } while (mpfq_fixmp_15_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_15_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 15);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_15_rshift(u, 1);
        if (sig)
          mpfq_fixmp_15_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 15);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_15_redc)
/* x has 32 words, z and p have 15 words.
 * only one word is read from invp.
 * Assuming R=W^16 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_15_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 15; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_15_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_15_add(x, x, x + 15);
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_15_cmp(x, p) >= 0) {
        mpfq_fixmp_15_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 15);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_redc) */

#if !defined(HAVE_native_mpfq_fixmp_15_redc_ur)
/* x has 33 words, z and p have 15 words.
 * only one word is read from invp.
 * Assuming R=W^16 is the redc modulus, we expect that x verifies:
 *  x < W*W^15*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_15_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[2];
    for (int i = 0; i < 15; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_15_addmul1(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy=mpfq_fixmp_15_add(x+15+1, x+15+1, x);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W+1)*p
    */
    if (cy) {
        /* x'/R-p < W*p, which fits in n+1 words */
        mpn_sub(x+15,x+15,15+1,p,15);
    }
    mpn_tdiv_qr(q, z, 0, x+15, 15+1, p, 15);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_15_mgy_encode)
/* x, z, and p have 15 words.
 * Assuming R=W^16 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_15_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[30] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14] };
    mpfq_fixmp_15_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_15_mgy_decode)
/* x, z, invR, and p have 15 words.
 * Assuming R=W^16 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_15_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[30];
    mpfq_fixmp_15_mul(t, x, invR);
    mpfq_fixmp_15_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_15_lshift)
/* a has 15 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_15_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 15-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_15_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_add)
/* x, y, and z have 0.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 0_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_0_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_add) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_sub)
/* x, y, and z have 0.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 0_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_0_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_sub) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_add_nc)
/* x, y, and z have 0.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_0_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_sub_nc)
/* x, y, and z have 0.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_0_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_add_ui)
/* x, y, and z have 0.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_0_5_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_sub_ui)
/* x, y, and z have 0.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_0_5_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_add_ui_nc)
/* x, y, and z have 0.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_0_5_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_sub_ui_nc)
/* x, y, and z have 0.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_0_5_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_cmp)
/* x and y have 0.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 0_5_redc_ur */
static inline
int mpfq_fixmp_0_5_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 1-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_cmp_ui)
/* x has 0.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
static inline
int mpfq_fixmp_0_5_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_addmul1)
/* x has 0.5 words, z has 2.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_0_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    z[1] += carry;
    return (z[1]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_addmul1_nc)
/* x has 0.5 words, z has 2.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
static inline
void mpfq_fixmp_0_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    z[1] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_addmul1_shortz)
/* x has 0.5 words, z has 1.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 0_5_redc, 0_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_0_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_addmul05_nc)
/* x has 0.5 words, z has 1. c is 0.5 word.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 0_5_mul, 0_5_mgy_decode */
static inline
void mpfq_fixmp_0_5_addmul05_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t lo, carry;
    carry = 0;
    lo = c*x[0] + carry;
    assert(lo >= carry);
    z[0] += lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_addmul05_nc) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_mul)
/* x and y have 0.5 words, z has 1. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 0_5_mgy_decode */
static inline
void mpfq_fixmp_0_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 1; z[i++] = 0) ;
    mpfq_fixmp_0_5_addmul05_nc (z + 0, x, y[0]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_mul) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_sqr)
/* x has 0.5 words, z has 1. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_0_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[1] = {0,};
    z[2*0] = x[0] * x[0];
    mpn_lshift(buf, buf, 1, 1);
    mpn_add_n(z, z, buf, 1);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_mul1)
/* x has 0.5 words, z has 2. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_0_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    z[1] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_shortmul)
/* x and y have 0.5 words, z has 1.
 * Put the low 1 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_0_5_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 1);
    z[1-1] += x[0]*y[1-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_addmul05)
/* x has 0.5 words, z has 1. c is 0.5 word.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul05 */
static inline
mp_limb_t mpfq_fixmp_0_5_addmul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t lo, carry;
    carry = 0;
    lo = c*x[0] + carry;
    assert(lo >= carry);
    z[0] += lo;
    return z[0] < lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_addmul05) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_mul05)
/* x has 0.5 words, z has 1. c is 0.5 word.
 * Put (x*c) in z. No carry. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul05 */
static inline
void mpfq_fixmp_0_5_mul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t lo, carry;
    carry = 0;
    lo = c*x[0] + carry;
    assert(lo >= carry);
    z[0] = lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_mul05) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_mod)
/* x has 1 words. z and p have 0.5 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 0_5_mgy_decode */
static inline
void mpfq_fixmp_0_5_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[0+1], r[1];
    assert (p[1-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 1, p, 1);
    mpfq_copy(z, r, 1);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_mod) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_invmod)
/* x, z, and p have 0.5 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_0_5_invmod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t a, b, u, v, fix;
      int t, lsh;
    
      a = *x;
      b = *p;
    
      if (a == 0 || a == b) {
        *z=0;
        return 0;
      }
      /* b must be odd and >a */
    
      fix = (b+1)>>1;
    
      assert (a < b);
    
      u = 1; v = 0; t = 0;
      
      /* compute u and t such that u*a_orig = 2^t mod b */
    
      /* we maintain:
       *    u*a_orig - (not kept)*b_orig = 2^t*a
       *    v*a_orig - (not kept)*b_orig = -2^t*b
       * a and b are both odd.
       * An update consists in reducing the largest by the smallest,
       * and then adjusting the valuation.  */
    
      lsh = mpfq_ctzl(a);
      a >>= lsh;
      t += lsh;
      v <<= lsh;
      do {
        do {
          b -= a; v += u;
          lsh = mpfq_ctzl(b);
          b >>= lsh;
          t += lsh;
          u <<= lsh;
        } while (a<b);
        if (a == b)
          break;
        do {
          a -= b; u += v;
          lsh = mpfq_ctzl(a);
          a >>= lsh;
          t += lsh;
          v <<= lsh;
        } while (b < a);
      } while (a != b);
      if (a != 1) {
        *z = a;
        return 0;
      }
      while (t>0) {
        mp_limb_t sig = u & 1UL;
        u >>= 1;
        if (sig)
          u += fix;
        --t;
      } 
      *z = u;
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_redc)
/* x has 1 words, z and p have 0.5 words.
 * only one word is read from invp.
 * Assuming R=W^1 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_0_5_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t t = x[0]*mip[0];
    mp_limb_t cy = mpfq_fixmp_0_5_addmul1_shortz(x, p, t);
    if (cy >= p[0]) {
        z[0] = cy - p[0];
    } else {
        z[0] = cy;
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_redc) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_redc_ur)
/* x has 2 words, z and p have 0.5 words.
 * only one word is read from invp.
 * Assuming R=W^1 is the redc modulus, we expect that x verifies:
 *  x < W*W^0.5*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_0_5_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[1];
    for(int i = 0; i < 1; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_0_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_0_5_add(x + 1, x, x + 1);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W^0.5*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W^0.5+1)*p
    */
    if (cy) {
        /* x'/R-p < W^0.5*p, which fits in n words. */
        mpfq_fixmp_0_5_sub(x + 1, x + 1, p);
    }
    mpn_tdiv_qr(q, z, 0, x + 1, 1, p, 1);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_mgy_encode)
/* x, z, and p have 0.5 words.
 * Assuming R=W^1 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_0_5_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[2] = { 0, x[0] };
    mp_limb_t qq[1+1];
    mpn_tdiv_qr(qq, z, 0, t, 2, p, 1);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_mgy_decode)
/* x, z, invR, and p have 0.5 words.
 * Assuming R=W^1 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_0_5_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[2];
    mpfq_fixmp_0_5_mul(t, x, invR);
    mpfq_fixmp_0_5_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_lshift)
/* a has 0.5 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_0_5_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_rshift)
/* a has 0.5 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_rshift */
static inline
void mpfq_fixmp_0_5_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    a[1-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_long_lshift)
/* a has 0.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_long_lshift */
static inline
void mpfq_fixmp_0_5_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        a[0] <<= cnt;
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_0_5_long_rshift)
/* a has 0.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_long_rshift */
static inline
void mpfq_fixmp_0_5_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        a[0] >>= cnt;
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_0_5_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_add)
/* x, y, and z have 1.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 1_5_invmod, 1_5_redc, 1_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_1_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_add) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_sub)
/* x, y, and z have 1.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 1_5_invmod, 1_5_redc, 1_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_1_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_sub) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_add_nc)
/* x, y, and z have 1.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_1_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_sub_nc)
/* x, y, and z have 1.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_1_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_add_ui)
/* x, y, and z have 1.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_1_5_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_sub_ui)
/* x, y, and z have 1.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_1_5_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_add_ui_nc)
/* x, y, and z have 1.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_1_5_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_sub_ui_nc)
/* x, y, and z have 1.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_1_5_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_cmp)
/* x and y have 1.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 1_5_invmod, 1_5_redc, 1_5_redc_ur */
static inline
int mpfq_fixmp_1_5_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 2-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_cmp_ui)
/* x has 1.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 1_5_invmod */
static inline
int mpfq_fixmp_1_5_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 2-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_addmul1)
/* x has 1.5 words, z has 3.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_1_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    z[2] += carry;
    return (z[2]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_addmul1_nc)
/* x has 1.5 words, z has 3.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 1_5_mul, 1_5_mgy_decode */
static inline
void mpfq_fixmp_1_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    z[2] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_addmul1_shortz)
/* x has 1.5 words, z has 2.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 1_5_redc, 1_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_1_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_addmul05_nc)
/* x has 1.5 words, z has 2. c is 0.5 word.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 1_5_mul, 1_5_mgy_decode */
static inline
void mpfq_fixmp_1_5_addmul05_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    lo = c*x[1] + carry;
    assert(lo >= carry);
    z[1] += lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_addmul05_nc) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_mul)
/* x and y have 1.5 words, z has 3. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 1_5_mgy_decode */
static inline
void mpfq_fixmp_1_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 3; z[i++] = 0) ;
    mpfq_fixmp_1_5_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_1_5_addmul05_nc (z + 1, x, y[1]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_mul) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_sqr)
/* x has 1.5 words, z has 3. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_1_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[3] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    z[2*1] = x[1] * x[1];
    mpn_lshift(buf, buf, 3, 1);
    mpn_add_n(z, z, buf, 3);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_mul1)
/* x has 1.5 words, z has 3. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_1_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    z[2] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_shortmul)
/* x and y have 1.5 words, z has 2.
 * Put the low 2 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_1_5_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 2);
    mpfq_fixmp_1_addmul1_nc (z+0, x, y[0]);
    z[2-1] += x[1]*y[0];
    z[2-1] += x[0]*y[2-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_addmul05)
/* x has 1.5 words, z has 2. c is 0.5 word.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul05 */
static inline
mp_limb_t mpfq_fixmp_1_5_addmul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    lo = c*x[1] + carry;
    assert(lo >= carry);
    z[1] += lo;
    return z[1] < lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_addmul05) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_mul05)
/* x has 1.5 words, z has 2. c is 0.5 word.
 * Put (x*c) in z. No carry. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul05 */
static inline
void mpfq_fixmp_1_5_mul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    lo = c*x[1] + carry;
    assert(lo >= carry);
    z[1] = lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_mul05) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_mod)
/* x has 3 words. z and p have 1.5 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 1_5_mgy_decode */
static inline
void mpfq_fixmp_1_5_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[1+1], r[2];
    assert (p[2-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 3, p, 2);
    mpfq_copy(z, r, 2);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_mod) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_rshift)
/* a has 1.5 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 1_5_invmod */
static inline
void mpfq_fixmp_1_5_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 2-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[2-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_long_rshift)
/* a has 1.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 1_5_invmod */
static inline
void mpfq_fixmp_1_5_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 2 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[2-off-1] = a[2-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 2 - off);
    }
    mpfq_zero(a + 2 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_long_lshift)
/* a has 1.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 1_5_invmod */
static inline
void mpfq_fixmp_1_5_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 2-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 2 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_invmod)
/* x, z, and p have 1.5 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_1_5_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[2], v[2], a[2], b[2], fix[2];
      int i, t, lsh;
    
      mpfq_zero(u, 2);
      mpfq_zero(v, 2);
      mpfq_copy(a, x, 2);
      mpfq_copy(b, p, 2);
      u[0] = 1UL;
      
      if (mpfq_fixmp_1_5_cmp(a, v) == 0 || mpfq_fixmp_1_5_cmp(a, b) == 0) {
        mpfq_zero(res, 2);
        return 0;
      }
    
      mpfq_fixmp_1_5_add(fix, b, u);
      mpfq_fixmp_1_5_rshift(fix, 1);
    
      assert (mpfq_fixmp_1_5_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 2);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_1_5_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_1_5_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_1_5_sub(b, b, a);
          mpfq_fixmp_1_5_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 2);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_1_5_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_1_5_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_1_5_cmp(a,b) < 0);
        if (mpfq_fixmp_1_5_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_1_5_sub(a, a, b);
          mpfq_fixmp_1_5_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 2);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_1_5_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_1_5_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_1_5_cmp(b,a)<0);
      } while (mpfq_fixmp_1_5_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_1_5_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 2);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_1_5_rshift(u, 1);
        if (sig)
          mpfq_fixmp_1_5_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 2);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_redc)
/* x has 3 words, z and p have 1.5 words.
 * only one word is read from invp.
 * Assuming R=W^2 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_1_5_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 2; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_1_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    {
        mp_limb_t ret[2] = { x[2], 0 };
        cy = mpfq_fixmp_1_5_add(x, x, ret);
    }
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_1_5_cmp(x, p) >= 0) {
        mpfq_fixmp_1_5_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 2);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_redc) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_redc_ur)
/* x has 4 words, z and p have 1.5 words.
 * only one word is read from invp.
 * Assuming R=W^2 is the redc modulus, we expect that x verifies:
 *  x < W*W^1.5*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_1_5_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[1];
    for(int i = 0; i < 2; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_1_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_1_5_add(x + 2, x, x + 2);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W^0.5*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W^0.5+1)*p
    */
    if (cy) {
        /* x'/R-p < W^0.5*p, which fits in n words. */
        mpfq_fixmp_1_5_sub(x + 2, x + 2, p);
    }
    mpn_tdiv_qr(q, z, 0, x + 2, 2, p, 2);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_mgy_encode)
/* x, z, and p have 1.5 words.
 * Assuming R=W^2 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_1_5_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[4] = { 0, 0, x[0], x[1] };
    mp_limb_t qq[2+1];
    mpn_tdiv_qr(qq, z, 0, t, 4, p, 2);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_mgy_decode)
/* x, z, invR, and p have 1.5 words.
 * Assuming R=W^2 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_1_5_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[4];
    mpfq_fixmp_1_5_mul(t, x, invR);
    mpfq_fixmp_1_5_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_1_5_lshift)
/* a has 1.5 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_1_5_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 2-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_1_5_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_add)
/* x, y, and z have 2.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 2_5_invmod, 2_5_redc, 2_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_2_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_add) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_sub)
/* x, y, and z have 2.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 2_5_invmod, 2_5_redc, 2_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_2_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_sub) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_add_nc)
/* x, y, and z have 2.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_2_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_sub_nc)
/* x, y, and z have 2.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_2_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_add_ui)
/* x, y, and z have 2.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_2_5_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_sub_ui)
/* x, y, and z have 2.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_2_5_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_add_ui_nc)
/* x, y, and z have 2.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_2_5_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_sub_ui_nc)
/* x, y, and z have 2.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_2_5_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_cmp)
/* x and y have 2.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 2_5_invmod, 2_5_redc, 2_5_redc_ur */
static inline
int mpfq_fixmp_2_5_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 3-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_cmp_ui)
/* x has 2.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 2_5_invmod */
static inline
int mpfq_fixmp_2_5_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 3-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_addmul1)
/* x has 2.5 words, z has 4.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_2_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    z[3] += carry;
    return (z[3]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_addmul1_nc)
/* x has 2.5 words, z has 4.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 2_5_mul, 2_5_mgy_decode */
static inline
void mpfq_fixmp_2_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    z[3] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_addmul1_shortz)
/* x has 2.5 words, z has 3.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 2_5_redc, 2_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_2_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_addmul05_nc)
/* x has 2.5 words, z has 3. c is 0.5 word.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 2_5_mul, 2_5_mgy_decode */
static inline
void mpfq_fixmp_2_5_addmul05_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    lo = c*x[2] + carry;
    assert(lo >= carry);
    z[2] += lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_addmul05_nc) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_mul)
/* x and y have 2.5 words, z has 5. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 2_5_mgy_decode */
static inline
void mpfq_fixmp_2_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 5; z[i++] = 0) ;
    mpfq_fixmp_2_5_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_2_5_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_2_5_addmul05_nc (z + 2, x, y[2]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_mul) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_sqr)
/* x has 2.5 words, z has 5. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_2_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[5] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    z[2*2] = x[2] * x[2];
    mpn_lshift(buf, buf, 5, 1);
    mpn_add_n(z, z, buf, 5);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_mul1)
/* x has 2.5 words, z has 4. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_2_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    z[3] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_shortmul)
/* x and y have 2.5 words, z has 3.
 * Put the low 3 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_2_5_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 3);
    mpfq_fixmp_2_addmul1_nc (z+0, x, y[0]);
    z[3-1] += x[2]*y[0];
    mpfq_fixmp_1_addmul1_nc (z+1, x, y[1]);
    z[3-1] += x[1]*y[1];
    z[3-1] += x[0]*y[3-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_addmul05)
/* x has 2.5 words, z has 3. c is 0.5 word.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul05 */
static inline
mp_limb_t mpfq_fixmp_2_5_addmul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    lo = c*x[2] + carry;
    assert(lo >= carry);
    z[2] += lo;
    return z[2] < lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_addmul05) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_mul05)
/* x has 2.5 words, z has 3. c is 0.5 word.
 * Put (x*c) in z. No carry. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul05 */
static inline
void mpfq_fixmp_2_5_mul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    lo = c*x[2] + carry;
    assert(lo >= carry);
    z[2] = lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_mul05) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_mod)
/* x has 5 words. z and p have 2.5 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 2_5_mgy_decode */
static inline
void mpfq_fixmp_2_5_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[2+1], r[3];
    assert (p[3-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 5, p, 3);
    mpfq_copy(z, r, 3);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_mod) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_rshift)
/* a has 2.5 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 2_5_invmod */
static inline
void mpfq_fixmp_2_5_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 3-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[3-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_long_rshift)
/* a has 2.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 2_5_invmod */
static inline
void mpfq_fixmp_2_5_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 3 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[3-off-1] = a[3-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 3 - off);
    }
    mpfq_zero(a + 3 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_long_lshift)
/* a has 2.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 2_5_invmod */
static inline
void mpfq_fixmp_2_5_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 3-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 3 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_invmod)
/* x, z, and p have 2.5 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_2_5_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[3], v[3], a[3], b[3], fix[3];
      int i, t, lsh;
    
      mpfq_zero(u, 3);
      mpfq_zero(v, 3);
      mpfq_copy(a, x, 3);
      mpfq_copy(b, p, 3);
      u[0] = 1UL;
      
      if (mpfq_fixmp_2_5_cmp(a, v) == 0 || mpfq_fixmp_2_5_cmp(a, b) == 0) {
        mpfq_zero(res, 3);
        return 0;
      }
    
      mpfq_fixmp_2_5_add(fix, b, u);
      mpfq_fixmp_2_5_rshift(fix, 1);
    
      assert (mpfq_fixmp_2_5_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 3);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_2_5_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_2_5_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_2_5_sub(b, b, a);
          mpfq_fixmp_2_5_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 3);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_2_5_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_2_5_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_2_5_cmp(a,b) < 0);
        if (mpfq_fixmp_2_5_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_2_5_sub(a, a, b);
          mpfq_fixmp_2_5_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 3);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_2_5_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_2_5_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_2_5_cmp(b,a)<0);
      } while (mpfq_fixmp_2_5_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_2_5_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 3);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_2_5_rshift(u, 1);
        if (sig)
          mpfq_fixmp_2_5_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 3);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_redc)
/* x has 5 words, z and p have 2.5 words.
 * only one word is read from invp.
 * Assuming R=W^3 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_2_5_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 3; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_2_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    {
        mp_limb_t ret[3] = { x[3], x[4], 0 };
        cy = mpfq_fixmp_2_5_add(x, x, ret);
    }
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_2_5_cmp(x, p) >= 0) {
        mpfq_fixmp_2_5_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 3);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_redc) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_redc_ur)
/* x has 6 words, z and p have 2.5 words.
 * only one word is read from invp.
 * Assuming R=W^3 is the redc modulus, we expect that x verifies:
 *  x < W*W^2.5*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_2_5_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[1];
    for(int i = 0; i < 3; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_2_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_2_5_add(x + 3, x, x + 3);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W^0.5*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W^0.5+1)*p
    */
    if (cy) {
        /* x'/R-p < W^0.5*p, which fits in n words. */
        mpfq_fixmp_2_5_sub(x + 3, x + 3, p);
    }
    mpn_tdiv_qr(q, z, 0, x + 3, 3, p, 3);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_mgy_encode)
/* x, z, and p have 2.5 words.
 * Assuming R=W^3 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_2_5_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[6] = { 0, 0, 0, x[0], x[1], x[2] };
    mp_limb_t qq[3+1];
    mpn_tdiv_qr(qq, z, 0, t, 6, p, 3);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_mgy_decode)
/* x, z, invR, and p have 2.5 words.
 * Assuming R=W^3 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_2_5_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[6];
    mpfq_fixmp_2_5_mul(t, x, invR);
    mpfq_fixmp_2_5_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_2_5_lshift)
/* a has 2.5 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_2_5_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 3-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_2_5_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_add)
/* x, y, and z have 3.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 3_5_invmod, 3_5_redc, 3_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_3_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_add) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_sub)
/* x, y, and z have 3.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 3_5_invmod, 3_5_redc, 3_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_3_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_sub) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_add_nc)
/* x, y, and z have 3.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_3_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_sub_nc)
/* x, y, and z have 3.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_3_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_add_ui)
/* x, y, and z have 3.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_3_5_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_sub_ui)
/* x, y, and z have 3.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_3_5_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_add_ui_nc)
/* x, y, and z have 3.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_3_5_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_sub_ui_nc)
/* x, y, and z have 3.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_3_5_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_cmp)
/* x and y have 3.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 3_5_invmod, 3_5_redc, 3_5_redc_ur */
static inline
int mpfq_fixmp_3_5_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 4-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_cmp_ui)
/* x has 3.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 3_5_invmod */
static inline
int mpfq_fixmp_3_5_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 4-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_addmul1)
/* x has 3.5 words, z has 5.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_3_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    z[4] += carry;
    return (z[4]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_addmul1_nc)
/* x has 3.5 words, z has 5.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 3_5_mul, 3_5_mgy_decode */
static inline
void mpfq_fixmp_3_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    z[4] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_addmul1_shortz)
/* x has 3.5 words, z has 4.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 3_5_redc, 3_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_3_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_addmul05_nc)
/* x has 3.5 words, z has 4. c is 0.5 word.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 3_5_mul, 3_5_mgy_decode */
static inline
void mpfq_fixmp_3_5_addmul05_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    lo = c*x[3] + carry;
    assert(lo >= carry);
    z[3] += lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_addmul05_nc) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_mul)
/* x and y have 3.5 words, z has 7. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 3_5_mgy_decode */
static inline
void mpfq_fixmp_3_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 7; z[i++] = 0) ;
    mpfq_fixmp_3_5_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_3_5_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_3_5_addmul1_nc (z + 2, x, y[2]);
    mpfq_fixmp_3_5_addmul05_nc (z + 3, x, y[3]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_mul) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_sqr)
/* x has 3.5 words, z has 7. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_3_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[7] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_fixmp_3_addmul1_nc(buf + 3, x, x[3]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpfq_umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
    z[2*3] = x[3] * x[3];
    mpn_lshift(buf, buf, 7, 1);
    mpn_add_n(z, z, buf, 7);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_mul1)
/* x has 3.5 words, z has 5. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_3_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    z[4] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_shortmul)
/* x and y have 3.5 words, z has 4.
 * Put the low 4 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_3_5_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 4);
    mpfq_fixmp_3_addmul1_nc (z+0, x, y[0]);
    z[4-1] += x[3]*y[0];
    mpfq_fixmp_2_addmul1_nc (z+1, x, y[1]);
    z[4-1] += x[2]*y[1];
    mpfq_fixmp_1_addmul1_nc (z+2, x, y[2]);
    z[4-1] += x[1]*y[2];
    z[4-1] += x[0]*y[4-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_addmul05)
/* x has 3.5 words, z has 4. c is 0.5 word.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul05 */
static inline
mp_limb_t mpfq_fixmp_3_5_addmul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    lo = c*x[3] + carry;
    assert(lo >= carry);
    z[3] += lo;
    return z[3] < lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_addmul05) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_mul05)
/* x has 3.5 words, z has 4. c is 0.5 word.
 * Put (x*c) in z. No carry. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul05 */
static inline
void mpfq_fixmp_3_5_mul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    lo = c*x[3] + carry;
    assert(lo >= carry);
    z[3] = lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_mul05) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_mod)
/* x has 7 words. z and p have 3.5 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 3_5_mgy_decode */
static inline
void mpfq_fixmp_3_5_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[3+1], r[4];
    assert (p[4-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 7, p, 4);
    mpfq_copy(z, r, 4);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_mod) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_rshift)
/* a has 3.5 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 3_5_invmod */
static inline
void mpfq_fixmp_3_5_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 4-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[4-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_long_rshift)
/* a has 3.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 3_5_invmod */
static inline
void mpfq_fixmp_3_5_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 4 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[4-off-1] = a[4-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 4 - off);
    }
    mpfq_zero(a + 4 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_long_lshift)
/* a has 3.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 3_5_invmod */
static inline
void mpfq_fixmp_3_5_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 4-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 4 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_invmod)
/* x, z, and p have 3.5 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_3_5_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[4], v[4], a[4], b[4], fix[4];
      int i, t, lsh;
    
      mpfq_zero(u, 4);
      mpfq_zero(v, 4);
      mpfq_copy(a, x, 4);
      mpfq_copy(b, p, 4);
      u[0] = 1UL;
      
      if (mpfq_fixmp_3_5_cmp(a, v) == 0 || mpfq_fixmp_3_5_cmp(a, b) == 0) {
        mpfq_zero(res, 4);
        return 0;
      }
    
      mpfq_fixmp_3_5_add(fix, b, u);
      mpfq_fixmp_3_5_rshift(fix, 1);
    
      assert (mpfq_fixmp_3_5_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 4);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_3_5_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_3_5_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_3_5_sub(b, b, a);
          mpfq_fixmp_3_5_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 4);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_3_5_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_3_5_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_3_5_cmp(a,b) < 0);
        if (mpfq_fixmp_3_5_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_3_5_sub(a, a, b);
          mpfq_fixmp_3_5_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 4);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_3_5_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_3_5_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_3_5_cmp(b,a)<0);
      } while (mpfq_fixmp_3_5_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_3_5_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 4);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_3_5_rshift(u, 1);
        if (sig)
          mpfq_fixmp_3_5_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 4);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_redc)
/* x has 7 words, z and p have 3.5 words.
 * only one word is read from invp.
 * Assuming R=W^4 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_3_5_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 4; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_3_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    {
        mp_limb_t ret[4] = { x[4], x[5], x[6], 0 };
        cy = mpfq_fixmp_3_5_add(x, x, ret);
    }
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_3_5_cmp(x, p) >= 0) {
        mpfq_fixmp_3_5_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 4);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_redc) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_redc_ur)
/* x has 8 words, z and p have 3.5 words.
 * only one word is read from invp.
 * Assuming R=W^4 is the redc modulus, we expect that x verifies:
 *  x < W*W^3.5*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_3_5_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[1];
    for(int i = 0; i < 4; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_3_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_3_5_add(x + 4, x, x + 4);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W^0.5*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W^0.5+1)*p
    */
    if (cy) {
        /* x'/R-p < W^0.5*p, which fits in n words. */
        mpfq_fixmp_3_5_sub(x + 4, x + 4, p);
    }
    mpn_tdiv_qr(q, z, 0, x + 4, 4, p, 4);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_mgy_encode)
/* x, z, and p have 3.5 words.
 * Assuming R=W^4 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_3_5_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[8] = { 0, 0, 0, 0, x[0], x[1], x[2], x[3] };
    mp_limb_t qq[4+1];
    mpn_tdiv_qr(qq, z, 0, t, 8, p, 4);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_mgy_decode)
/* x, z, invR, and p have 3.5 words.
 * Assuming R=W^4 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_3_5_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[8];
    mpfq_fixmp_3_5_mul(t, x, invR);
    mpfq_fixmp_3_5_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_3_5_lshift)
/* a has 3.5 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_3_5_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 4-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_3_5_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_add)
/* x, y, and z have 4.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 4_5_invmod, 4_5_redc, 4_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_4_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_add) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_sub)
/* x, y, and z have 4.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 4_5_invmod, 4_5_redc, 4_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_4_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_sub) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_add_nc)
/* x, y, and z have 4.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_4_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_sub_nc)
/* x, y, and z have 4.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_4_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_add_ui)
/* x, y, and z have 4.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_4_5_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_sub_ui)
/* x, y, and z have 4.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_4_5_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_add_ui_nc)
/* x, y, and z have 4.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_4_5_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_sub_ui_nc)
/* x, y, and z have 4.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_4_5_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_cmp)
/* x and y have 4.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 4_5_invmod, 4_5_redc, 4_5_redc_ur */
static inline
int mpfq_fixmp_4_5_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 5-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_cmp_ui)
/* x has 4.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 4_5_invmod */
static inline
int mpfq_fixmp_4_5_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 5-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_addmul1)
/* x has 4.5 words, z has 6.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_4_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    z[5] += carry;
    return (z[5]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_addmul1_nc)
/* x has 4.5 words, z has 6.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 4_5_mul, 4_5_mgy_decode */
static inline
void mpfq_fixmp_4_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    z[5] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_addmul1_shortz)
/* x has 4.5 words, z has 5.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 4_5_redc, 4_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_4_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_addmul05_nc)
/* x has 4.5 words, z has 5. c is 0.5 word.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 4_5_mul, 4_5_mgy_decode */
static inline
void mpfq_fixmp_4_5_addmul05_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    lo = c*x[4] + carry;
    assert(lo >= carry);
    z[4] += lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_addmul05_nc) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_mul)
/* x and y have 4.5 words, z has 9. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 4_5_mgy_decode */
static inline
void mpfq_fixmp_4_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 9; z[i++] = 0) ;
    mpfq_fixmp_4_5_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_4_5_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_4_5_addmul1_nc (z + 2, x, y[2]);
    mpfq_fixmp_4_5_addmul1_nc (z + 3, x, y[3]);
    mpfq_fixmp_4_5_addmul05_nc (z + 4, x, y[4]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_mul) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_sqr)
/* x has 4.5 words, z has 9. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_4_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[9] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_fixmp_3_addmul1_nc(buf + 3, x, x[3]);
    mpfq_fixmp_4_addmul1_nc(buf + 4, x, x[4]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpfq_umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
    mpfq_umul_ppmm(z[2*3+1], z[2*3], x[3], x[3]);
    z[2*4] = x[4] * x[4];
    mpn_lshift(buf, buf, 9, 1);
    mpn_add_n(z, z, buf, 9);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_mul1)
/* x has 4.5 words, z has 6. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_4_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    z[5] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_shortmul)
/* x and y have 4.5 words, z has 5.
 * Put the low 5 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_4_5_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 5);
    mpfq_fixmp_4_addmul1_nc (z+0, x, y[0]);
    z[5-1] += x[4]*y[0];
    mpfq_fixmp_3_addmul1_nc (z+1, x, y[1]);
    z[5-1] += x[3]*y[1];
    mpfq_fixmp_2_addmul1_nc (z+2, x, y[2]);
    z[5-1] += x[2]*y[2];
    mpfq_fixmp_1_addmul1_nc (z+3, x, y[3]);
    z[5-1] += x[1]*y[3];
    z[5-1] += x[0]*y[5-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_addmul05)
/* x has 4.5 words, z has 5. c is 0.5 word.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul05 */
static inline
mp_limb_t mpfq_fixmp_4_5_addmul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    lo = c*x[4] + carry;
    assert(lo >= carry);
    z[4] += lo;
    return z[4] < lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_addmul05) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_mul05)
/* x has 4.5 words, z has 5. c is 0.5 word.
 * Put (x*c) in z. No carry. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul05 */
static inline
void mpfq_fixmp_4_5_mul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    lo = c*x[4] + carry;
    assert(lo >= carry);
    z[4] = lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_mul05) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_mod)
/* x has 9 words. z and p have 4.5 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 4_5_mgy_decode */
static inline
void mpfq_fixmp_4_5_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[4+1], r[5];
    assert (p[5-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 9, p, 5);
    mpfq_copy(z, r, 5);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_mod) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_rshift)
/* a has 4.5 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 4_5_invmod */
static inline
void mpfq_fixmp_4_5_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 5-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[5-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_long_rshift)
/* a has 4.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 4_5_invmod */
static inline
void mpfq_fixmp_4_5_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 5 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[5-off-1] = a[5-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 5 - off);
    }
    mpfq_zero(a + 5 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_long_lshift)
/* a has 4.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 4_5_invmod */
static inline
void mpfq_fixmp_4_5_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 5-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 5 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_invmod)
/* x, z, and p have 4.5 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_4_5_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[5], v[5], a[5], b[5], fix[5];
      int i, t, lsh;
    
      mpfq_zero(u, 5);
      mpfq_zero(v, 5);
      mpfq_copy(a, x, 5);
      mpfq_copy(b, p, 5);
      u[0] = 1UL;
      
      if (mpfq_fixmp_4_5_cmp(a, v) == 0 || mpfq_fixmp_4_5_cmp(a, b) == 0) {
        mpfq_zero(res, 5);
        return 0;
      }
    
      mpfq_fixmp_4_5_add(fix, b, u);
      mpfq_fixmp_4_5_rshift(fix, 1);
    
      assert (mpfq_fixmp_4_5_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 5);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_4_5_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_4_5_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_4_5_sub(b, b, a);
          mpfq_fixmp_4_5_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 5);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_4_5_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_4_5_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_4_5_cmp(a,b) < 0);
        if (mpfq_fixmp_4_5_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_4_5_sub(a, a, b);
          mpfq_fixmp_4_5_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 5);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_4_5_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_4_5_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_4_5_cmp(b,a)<0);
      } while (mpfq_fixmp_4_5_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_4_5_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 5);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_4_5_rshift(u, 1);
        if (sig)
          mpfq_fixmp_4_5_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 5);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_redc)
/* x has 9 words, z and p have 4.5 words.
 * only one word is read from invp.
 * Assuming R=W^5 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_4_5_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 5; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_4_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    {
        mp_limb_t ret[5] = { x[5], x[6], x[7], x[8], 0 };
        cy = mpfq_fixmp_4_5_add(x, x, ret);
    }
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_4_5_cmp(x, p) >= 0) {
        mpfq_fixmp_4_5_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 5);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_redc) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_redc_ur)
/* x has 10 words, z and p have 4.5 words.
 * only one word is read from invp.
 * Assuming R=W^5 is the redc modulus, we expect that x verifies:
 *  x < W*W^4.5*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_4_5_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[1];
    for(int i = 0; i < 5; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_4_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_4_5_add(x + 5, x, x + 5);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W^0.5*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W^0.5+1)*p
    */
    if (cy) {
        /* x'/R-p < W^0.5*p, which fits in n words. */
        mpfq_fixmp_4_5_sub(x + 5, x + 5, p);
    }
    mpn_tdiv_qr(q, z, 0, x + 5, 5, p, 5);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_mgy_encode)
/* x, z, and p have 4.5 words.
 * Assuming R=W^5 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_4_5_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[10] = { 0, 0, 0, 0, 0, x[0], x[1], x[2], x[3], x[4] };
    mp_limb_t qq[5+1];
    mpn_tdiv_qr(qq, z, 0, t, 10, p, 5);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_mgy_decode)
/* x, z, invR, and p have 4.5 words.
 * Assuming R=W^5 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_4_5_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[10];
    mpfq_fixmp_4_5_mul(t, x, invR);
    mpfq_fixmp_4_5_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_4_5_lshift)
/* a has 4.5 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_4_5_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 5-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_4_5_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_add)
/* x, y, and z have 5.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 5_5_invmod, 5_5_redc, 5_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_5_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_add) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_sub)
/* x, y, and z have 5.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 5_5_invmod, 5_5_redc, 5_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_5_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_sub) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_add_nc)
/* x, y, and z have 5.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_5_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_sub_nc)
/* x, y, and z have 5.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_5_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_add_ui)
/* x, y, and z have 5.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_5_5_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_sub_ui)
/* x, y, and z have 5.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_5_5_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_add_ui_nc)
/* x, y, and z have 5.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_5_5_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_sub_ui_nc)
/* x, y, and z have 5.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_5_5_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_cmp)
/* x and y have 5.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 5_5_invmod, 5_5_redc, 5_5_redc_ur */
static inline
int mpfq_fixmp_5_5_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 6-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_cmp_ui)
/* x has 5.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 5_5_invmod */
static inline
int mpfq_fixmp_5_5_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 6-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_addmul1)
/* x has 5.5 words, z has 7.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_5_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    z[6] += carry;
    return (z[6]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_addmul1_nc)
/* x has 5.5 words, z has 7.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 5_5_mul, 5_5_mgy_decode */
static inline
void mpfq_fixmp_5_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    z[6] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_addmul1_shortz)
/* x has 5.5 words, z has 6.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 5_5_redc, 5_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_5_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_addmul05_nc)
/* x has 5.5 words, z has 6. c is 0.5 word.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 5_5_mul, 5_5_mgy_decode */
static inline
void mpfq_fixmp_5_5_addmul05_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    lo = c*x[5] + carry;
    assert(lo >= carry);
    z[5] += lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_addmul05_nc) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_mul)
/* x and y have 5.5 words, z has 11. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 5_5_mgy_decode */
static inline
void mpfq_fixmp_5_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 11; z[i++] = 0) ;
    mpfq_fixmp_5_5_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_5_5_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_5_5_addmul1_nc (z + 2, x, y[2]);
    mpfq_fixmp_5_5_addmul1_nc (z + 3, x, y[3]);
    mpfq_fixmp_5_5_addmul1_nc (z + 4, x, y[4]);
    mpfq_fixmp_5_5_addmul05_nc (z + 5, x, y[5]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_mul) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_sqr)
/* x has 5.5 words, z has 11. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_5_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[11] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_fixmp_3_addmul1_nc(buf + 3, x, x[3]);
    mpfq_fixmp_4_addmul1_nc(buf + 4, x, x[4]);
    mpfq_fixmp_5_addmul1_nc(buf + 5, x, x[5]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpfq_umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
    mpfq_umul_ppmm(z[2*3+1], z[2*3], x[3], x[3]);
    mpfq_umul_ppmm(z[2*4+1], z[2*4], x[4], x[4]);
    z[2*5] = x[5] * x[5];
    mpn_lshift(buf, buf, 11, 1);
    mpn_add_n(z, z, buf, 11);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_mul1)
/* x has 5.5 words, z has 7. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_5_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    z[6] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_shortmul)
/* x and y have 5.5 words, z has 6.
 * Put the low 6 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_5_5_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 6);
    mpfq_fixmp_5_addmul1_nc (z+0, x, y[0]);
    z[6-1] += x[5]*y[0];
    mpfq_fixmp_4_addmul1_nc (z+1, x, y[1]);
    z[6-1] += x[4]*y[1];
    mpfq_fixmp_3_addmul1_nc (z+2, x, y[2]);
    z[6-1] += x[3]*y[2];
    mpfq_fixmp_2_addmul1_nc (z+3, x, y[3]);
    z[6-1] += x[2]*y[3];
    mpfq_fixmp_1_addmul1_nc (z+4, x, y[4]);
    z[6-1] += x[1]*y[4];
    z[6-1] += x[0]*y[6-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_addmul05)
/* x has 5.5 words, z has 6. c is 0.5 word.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul05 */
static inline
mp_limb_t mpfq_fixmp_5_5_addmul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    lo = c*x[5] + carry;
    assert(lo >= carry);
    z[5] += lo;
    return z[5] < lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_addmul05) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_mul05)
/* x has 5.5 words, z has 6. c is 0.5 word.
 * Put (x*c) in z. No carry. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul05 */
static inline
void mpfq_fixmp_5_5_mul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    lo = c*x[5] + carry;
    assert(lo >= carry);
    z[5] = lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_mul05) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_mod)
/* x has 11 words. z and p have 5.5 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 5_5_mgy_decode */
static inline
void mpfq_fixmp_5_5_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[5+1], r[6];
    assert (p[6-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 11, p, 6);
    mpfq_copy(z, r, 6);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_mod) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_rshift)
/* a has 5.5 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 5_5_invmod */
static inline
void mpfq_fixmp_5_5_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 6-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[6-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_long_rshift)
/* a has 5.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 5_5_invmod */
static inline
void mpfq_fixmp_5_5_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 6 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[6-off-1] = a[6-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 6 - off);
    }
    mpfq_zero(a + 6 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_long_lshift)
/* a has 5.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 5_5_invmod */
static inline
void mpfq_fixmp_5_5_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 6-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 6 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_invmod)
/* x, z, and p have 5.5 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_5_5_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[6], v[6], a[6], b[6], fix[6];
      int i, t, lsh;
    
      mpfq_zero(u, 6);
      mpfq_zero(v, 6);
      mpfq_copy(a, x, 6);
      mpfq_copy(b, p, 6);
      u[0] = 1UL;
      
      if (mpfq_fixmp_5_5_cmp(a, v) == 0 || mpfq_fixmp_5_5_cmp(a, b) == 0) {
        mpfq_zero(res, 6);
        return 0;
      }
    
      mpfq_fixmp_5_5_add(fix, b, u);
      mpfq_fixmp_5_5_rshift(fix, 1);
    
      assert (mpfq_fixmp_5_5_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 6);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_5_5_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_5_5_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_5_5_sub(b, b, a);
          mpfq_fixmp_5_5_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 6);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_5_5_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_5_5_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_5_5_cmp(a,b) < 0);
        if (mpfq_fixmp_5_5_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_5_5_sub(a, a, b);
          mpfq_fixmp_5_5_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 6);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_5_5_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_5_5_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_5_5_cmp(b,a)<0);
      } while (mpfq_fixmp_5_5_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_5_5_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 6);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_5_5_rshift(u, 1);
        if (sig)
          mpfq_fixmp_5_5_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 6);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_redc)
/* x has 11 words, z and p have 5.5 words.
 * only one word is read from invp.
 * Assuming R=W^6 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_5_5_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 6; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_5_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    {
        mp_limb_t ret[6] = { x[6], x[7], x[8], x[9], x[10], 0 };
        cy = mpfq_fixmp_5_5_add(x, x, ret);
    }
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_5_5_cmp(x, p) >= 0) {
        mpfq_fixmp_5_5_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 6);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_redc) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_redc_ur)
/* x has 12 words, z and p have 5.5 words.
 * only one word is read from invp.
 * Assuming R=W^6 is the redc modulus, we expect that x verifies:
 *  x < W*W^5.5*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_5_5_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[1];
    for(int i = 0; i < 6; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_5_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_5_5_add(x + 6, x, x + 6);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W^0.5*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W^0.5+1)*p
    */
    if (cy) {
        /* x'/R-p < W^0.5*p, which fits in n words. */
        mpfq_fixmp_5_5_sub(x + 6, x + 6, p);
    }
    mpn_tdiv_qr(q, z, 0, x + 6, 6, p, 6);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_mgy_encode)
/* x, z, and p have 5.5 words.
 * Assuming R=W^6 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_5_5_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[12] = { 0, 0, 0, 0, 0, 0, x[0], x[1], x[2], x[3], x[4], x[5] };
    mp_limb_t qq[6+1];
    mpn_tdiv_qr(qq, z, 0, t, 12, p, 6);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_mgy_decode)
/* x, z, invR, and p have 5.5 words.
 * Assuming R=W^6 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_5_5_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[12];
    mpfq_fixmp_5_5_mul(t, x, invR);
    mpfq_fixmp_5_5_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_5_5_lshift)
/* a has 5.5 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_5_5_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 6-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_5_5_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_add)
/* x, y, and z have 6.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 6_5_invmod, 6_5_redc, 6_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_6_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_add) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_sub)
/* x, y, and z have 6.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 6_5_invmod, 6_5_redc, 6_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_6_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_sub) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_add_nc)
/* x, y, and z have 6.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_6_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_sub_nc)
/* x, y, and z have 6.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_6_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_add_ui)
/* x, y, and z have 6.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_6_5_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_sub_ui)
/* x, y, and z have 6.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_6_5_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_add_ui_nc)
/* x, y, and z have 6.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_6_5_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_sub_ui_nc)
/* x, y, and z have 6.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_6_5_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_cmp)
/* x and y have 6.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 6_5_invmod, 6_5_redc, 6_5_redc_ur */
static inline
int mpfq_fixmp_6_5_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 7-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_cmp_ui)
/* x has 6.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 6_5_invmod */
static inline
int mpfq_fixmp_6_5_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 7-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_addmul1)
/* x has 6.5 words, z has 8.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_6_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    z[7] += carry;
    return (z[7]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_addmul1_nc)
/* x has 6.5 words, z has 8.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 6_5_mul, 6_5_mgy_decode */
static inline
void mpfq_fixmp_6_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    z[7] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_addmul1_shortz)
/* x has 6.5 words, z has 7.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 6_5_redc, 6_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_6_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_addmul05_nc)
/* x has 6.5 words, z has 7. c is 0.5 word.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 6_5_mul, 6_5_mgy_decode */
static inline
void mpfq_fixmp_6_5_addmul05_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    lo = c*x[6] + carry;
    assert(lo >= carry);
    z[6] += lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_addmul05_nc) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_mul)
/* x and y have 6.5 words, z has 13. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 6_5_mgy_decode */
static inline
void mpfq_fixmp_6_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 13; z[i++] = 0) ;
    mpfq_fixmp_6_5_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_6_5_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_6_5_addmul1_nc (z + 2, x, y[2]);
    mpfq_fixmp_6_5_addmul1_nc (z + 3, x, y[3]);
    mpfq_fixmp_6_5_addmul1_nc (z + 4, x, y[4]);
    mpfq_fixmp_6_5_addmul1_nc (z + 5, x, y[5]);
    mpfq_fixmp_6_5_addmul05_nc (z + 6, x, y[6]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_mul) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_sqr)
/* x has 6.5 words, z has 13. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_6_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[13] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_fixmp_3_addmul1_nc(buf + 3, x, x[3]);
    mpfq_fixmp_4_addmul1_nc(buf + 4, x, x[4]);
    mpfq_fixmp_5_addmul1_nc(buf + 5, x, x[5]);
    mpfq_fixmp_6_addmul1_nc(buf + 6, x, x[6]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpfq_umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
    mpfq_umul_ppmm(z[2*3+1], z[2*3], x[3], x[3]);
    mpfq_umul_ppmm(z[2*4+1], z[2*4], x[4], x[4]);
    mpfq_umul_ppmm(z[2*5+1], z[2*5], x[5], x[5]);
    z[2*6] = x[6] * x[6];
    mpn_lshift(buf, buf, 13, 1);
    mpn_add_n(z, z, buf, 13);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_mul1)
/* x has 6.5 words, z has 8. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_6_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    z[7] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_shortmul)
/* x and y have 6.5 words, z has 7.
 * Put the low 7 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_6_5_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 7);
    mpfq_fixmp_6_addmul1_nc (z+0, x, y[0]);
    z[7-1] += x[6]*y[0];
    mpfq_fixmp_5_addmul1_nc (z+1, x, y[1]);
    z[7-1] += x[5]*y[1];
    mpfq_fixmp_4_addmul1_nc (z+2, x, y[2]);
    z[7-1] += x[4]*y[2];
    mpfq_fixmp_3_addmul1_nc (z+3, x, y[3]);
    z[7-1] += x[3]*y[3];
    mpfq_fixmp_2_addmul1_nc (z+4, x, y[4]);
    z[7-1] += x[2]*y[4];
    mpfq_fixmp_1_addmul1_nc (z+5, x, y[5]);
    z[7-1] += x[1]*y[5];
    z[7-1] += x[0]*y[7-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_addmul05)
/* x has 6.5 words, z has 7. c is 0.5 word.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul05 */
static inline
mp_limb_t mpfq_fixmp_6_5_addmul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    lo = c*x[6] + carry;
    assert(lo >= carry);
    z[6] += lo;
    return z[6] < lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_addmul05) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_mul05)
/* x has 6.5 words, z has 7. c is 0.5 word.
 * Put (x*c) in z. No carry. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul05 */
static inline
void mpfq_fixmp_6_5_mul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    lo = c*x[6] + carry;
    assert(lo >= carry);
    z[6] = lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_mul05) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_mod)
/* x has 13 words. z and p have 6.5 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 6_5_mgy_decode */
static inline
void mpfq_fixmp_6_5_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[6+1], r[7];
    assert (p[7-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 13, p, 7);
    mpfq_copy(z, r, 7);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_mod) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_rshift)
/* a has 6.5 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 6_5_invmod */
static inline
void mpfq_fixmp_6_5_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 7-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[7-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_long_rshift)
/* a has 6.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 6_5_invmod */
static inline
void mpfq_fixmp_6_5_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 7 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[7-off-1] = a[7-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 7 - off);
    }
    mpfq_zero(a + 7 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_long_lshift)
/* a has 6.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 6_5_invmod */
static inline
void mpfq_fixmp_6_5_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 7-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 7 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_invmod)
/* x, z, and p have 6.5 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_6_5_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[7], v[7], a[7], b[7], fix[7];
      int i, t, lsh;
    
      mpfq_zero(u, 7);
      mpfq_zero(v, 7);
      mpfq_copy(a, x, 7);
      mpfq_copy(b, p, 7);
      u[0] = 1UL;
      
      if (mpfq_fixmp_6_5_cmp(a, v) == 0 || mpfq_fixmp_6_5_cmp(a, b) == 0) {
        mpfq_zero(res, 7);
        return 0;
      }
    
      mpfq_fixmp_6_5_add(fix, b, u);
      mpfq_fixmp_6_5_rshift(fix, 1);
    
      assert (mpfq_fixmp_6_5_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 7);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_6_5_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_6_5_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_6_5_sub(b, b, a);
          mpfq_fixmp_6_5_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 7);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_6_5_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_6_5_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_6_5_cmp(a,b) < 0);
        if (mpfq_fixmp_6_5_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_6_5_sub(a, a, b);
          mpfq_fixmp_6_5_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 7);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_6_5_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_6_5_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_6_5_cmp(b,a)<0);
      } while (mpfq_fixmp_6_5_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_6_5_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 7);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_6_5_rshift(u, 1);
        if (sig)
          mpfq_fixmp_6_5_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 7);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_redc)
/* x has 13 words, z and p have 6.5 words.
 * only one word is read from invp.
 * Assuming R=W^7 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_6_5_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 7; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_6_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    {
        mp_limb_t ret[7] = { x[7], x[8], x[9], x[10], x[11], x[12], 0 };
        cy = mpfq_fixmp_6_5_add(x, x, ret);
    }
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_6_5_cmp(x, p) >= 0) {
        mpfq_fixmp_6_5_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 7);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_redc) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_redc_ur)
/* x has 14 words, z and p have 6.5 words.
 * only one word is read from invp.
 * Assuming R=W^7 is the redc modulus, we expect that x verifies:
 *  x < W*W^6.5*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_6_5_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[1];
    for(int i = 0; i < 7; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_6_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_6_5_add(x + 7, x, x + 7);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W^0.5*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W^0.5+1)*p
    */
    if (cy) {
        /* x'/R-p < W^0.5*p, which fits in n words. */
        mpfq_fixmp_6_5_sub(x + 7, x + 7, p);
    }
    mpn_tdiv_qr(q, z, 0, x + 7, 7, p, 7);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_mgy_encode)
/* x, z, and p have 6.5 words.
 * Assuming R=W^7 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_6_5_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[14] = { 0, 0, 0, 0, 0, 0, 0, x[0], x[1], x[2], x[3], x[4], x[5], x[6] };
    mp_limb_t qq[7+1];
    mpn_tdiv_qr(qq, z, 0, t, 14, p, 7);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_mgy_decode)
/* x, z, invR, and p have 6.5 words.
 * Assuming R=W^7 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_6_5_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[14];
    mpfq_fixmp_6_5_mul(t, x, invR);
    mpfq_fixmp_6_5_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_6_5_lshift)
/* a has 6.5 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_6_5_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 7-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_6_5_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_add)
/* x, y, and z have 7.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 7_5_invmod, 7_5_redc, 7_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_7_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_add) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_sub)
/* x, y, and z have 7.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 7_5_invmod, 7_5_redc, 7_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_7_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_sub) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_add_nc)
/* x, y, and z have 7.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_7_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_sub_nc)
/* x, y, and z have 7.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_7_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_add_ui)
/* x, y, and z have 7.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_7_5_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_sub_ui)
/* x, y, and z have 7.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_7_5_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_add_ui_nc)
/* x, y, and z have 7.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_7_5_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_sub_ui_nc)
/* x, y, and z have 7.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_7_5_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_cmp)
/* x and y have 7.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 7_5_invmod, 7_5_redc, 7_5_redc_ur */
static inline
int mpfq_fixmp_7_5_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 8-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_cmp_ui)
/* x has 7.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 7_5_invmod */
static inline
int mpfq_fixmp_7_5_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 8-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_addmul1)
/* x has 7.5 words, z has 9.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_7_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    z[8] += carry;
    return (z[8]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_addmul1_nc)
/* x has 7.5 words, z has 9.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 7_5_mul, 7_5_mgy_decode */
static inline
void mpfq_fixmp_7_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    z[8] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_addmul1_shortz)
/* x has 7.5 words, z has 8.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 7_5_redc, 7_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_7_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_addmul05_nc)
/* x has 7.5 words, z has 8. c is 0.5 word.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 7_5_mul, 7_5_mgy_decode */
static inline
void mpfq_fixmp_7_5_addmul05_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    lo = c*x[7] + carry;
    assert(lo >= carry);
    z[7] += lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_addmul05_nc) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_mul)
/* x and y have 7.5 words, z has 15. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 7_5_mgy_decode */
static inline
void mpfq_fixmp_7_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 15; z[i++] = 0) ;
    mpfq_fixmp_7_5_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_7_5_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_7_5_addmul1_nc (z + 2, x, y[2]);
    mpfq_fixmp_7_5_addmul1_nc (z + 3, x, y[3]);
    mpfq_fixmp_7_5_addmul1_nc (z + 4, x, y[4]);
    mpfq_fixmp_7_5_addmul1_nc (z + 5, x, y[5]);
    mpfq_fixmp_7_5_addmul1_nc (z + 6, x, y[6]);
    mpfq_fixmp_7_5_addmul05_nc (z + 7, x, y[7]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_mul) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_sqr)
/* x has 7.5 words, z has 15. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_7_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[15] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_fixmp_3_addmul1_nc(buf + 3, x, x[3]);
    mpfq_fixmp_4_addmul1_nc(buf + 4, x, x[4]);
    mpfq_fixmp_5_addmul1_nc(buf + 5, x, x[5]);
    mpfq_fixmp_6_addmul1_nc(buf + 6, x, x[6]);
    mpfq_fixmp_7_addmul1_nc(buf + 7, x, x[7]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpfq_umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
    mpfq_umul_ppmm(z[2*3+1], z[2*3], x[3], x[3]);
    mpfq_umul_ppmm(z[2*4+1], z[2*4], x[4], x[4]);
    mpfq_umul_ppmm(z[2*5+1], z[2*5], x[5], x[5]);
    mpfq_umul_ppmm(z[2*6+1], z[2*6], x[6], x[6]);
    z[2*7] = x[7] * x[7];
    mpn_lshift(buf, buf, 15, 1);
    mpn_add_n(z, z, buf, 15);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_mul1)
/* x has 7.5 words, z has 9. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_7_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[7] = lo;
    z[8] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_shortmul)
/* x and y have 7.5 words, z has 8.
 * Put the low 8 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_7_5_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 8);
    mpfq_fixmp_7_addmul1_nc (z+0, x, y[0]);
    z[8-1] += x[7]*y[0];
    mpfq_fixmp_6_addmul1_nc (z+1, x, y[1]);
    z[8-1] += x[6]*y[1];
    mpfq_fixmp_5_addmul1_nc (z+2, x, y[2]);
    z[8-1] += x[5]*y[2];
    mpfq_fixmp_4_addmul1_nc (z+3, x, y[3]);
    z[8-1] += x[4]*y[3];
    mpfq_fixmp_3_addmul1_nc (z+4, x, y[4]);
    z[8-1] += x[3]*y[4];
    mpfq_fixmp_2_addmul1_nc (z+5, x, y[5]);
    z[8-1] += x[2]*y[5];
    mpfq_fixmp_1_addmul1_nc (z+6, x, y[6]);
    z[8-1] += x[1]*y[6];
    z[8-1] += x[0]*y[8-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_addmul05)
/* x has 7.5 words, z has 8. c is 0.5 word.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul05 */
static inline
mp_limb_t mpfq_fixmp_7_5_addmul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    lo = c*x[7] + carry;
    assert(lo >= carry);
    z[7] += lo;
    return z[7] < lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_addmul05) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_mul05)
/* x has 7.5 words, z has 8. c is 0.5 word.
 * Put (x*c) in z. No carry. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul05 */
static inline
void mpfq_fixmp_7_5_mul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    lo = c*x[7] + carry;
    assert(lo >= carry);
    z[7] = lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_mul05) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_mod)
/* x has 15 words. z and p have 7.5 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 7_5_mgy_decode */
static inline
void mpfq_fixmp_7_5_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[7+1], r[8];
    assert (p[8-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 15, p, 8);
    mpfq_copy(z, r, 8);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_mod) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_rshift)
/* a has 7.5 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 7_5_invmod */
static inline
void mpfq_fixmp_7_5_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 8-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[8-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_long_rshift)
/* a has 7.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 7_5_invmod */
static inline
void mpfq_fixmp_7_5_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 8 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[8-off-1] = a[8-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 8 - off);
    }
    mpfq_zero(a + 8 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_long_lshift)
/* a has 7.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 7_5_invmod */
static inline
void mpfq_fixmp_7_5_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 8-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 8 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_invmod)
/* x, z, and p have 7.5 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_7_5_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[8], v[8], a[8], b[8], fix[8];
      int i, t, lsh;
    
      mpfq_zero(u, 8);
      mpfq_zero(v, 8);
      mpfq_copy(a, x, 8);
      mpfq_copy(b, p, 8);
      u[0] = 1UL;
      
      if (mpfq_fixmp_7_5_cmp(a, v) == 0 || mpfq_fixmp_7_5_cmp(a, b) == 0) {
        mpfq_zero(res, 8);
        return 0;
      }
    
      mpfq_fixmp_7_5_add(fix, b, u);
      mpfq_fixmp_7_5_rshift(fix, 1);
    
      assert (mpfq_fixmp_7_5_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 8);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_7_5_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_7_5_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_7_5_sub(b, b, a);
          mpfq_fixmp_7_5_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 8);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_7_5_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_7_5_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_7_5_cmp(a,b) < 0);
        if (mpfq_fixmp_7_5_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_7_5_sub(a, a, b);
          mpfq_fixmp_7_5_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 8);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_7_5_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_7_5_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_7_5_cmp(b,a)<0);
      } while (mpfq_fixmp_7_5_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_7_5_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 8);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_7_5_rshift(u, 1);
        if (sig)
          mpfq_fixmp_7_5_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 8);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_redc)
/* x has 15 words, z and p have 7.5 words.
 * only one word is read from invp.
 * Assuming R=W^8 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_7_5_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 8; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_7_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    {
        mp_limb_t ret[8] = { x[8], x[9], x[10], x[11], x[12], x[13], x[14], 0 };
        cy = mpfq_fixmp_7_5_add(x, x, ret);
    }
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_7_5_cmp(x, p) >= 0) {
        mpfq_fixmp_7_5_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 8);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_redc) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_redc_ur)
/* x has 16 words, z and p have 7.5 words.
 * only one word is read from invp.
 * Assuming R=W^8 is the redc modulus, we expect that x verifies:
 *  x < W*W^7.5*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_7_5_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[1];
    for(int i = 0; i < 8; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_7_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_7_5_add(x + 8, x, x + 8);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W^0.5*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W^0.5+1)*p
    */
    if (cy) {
        /* x'/R-p < W^0.5*p, which fits in n words. */
        mpfq_fixmp_7_5_sub(x + 8, x + 8, p);
    }
    mpn_tdiv_qr(q, z, 0, x + 8, 8, p, 8);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_mgy_encode)
/* x, z, and p have 7.5 words.
 * Assuming R=W^8 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_7_5_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[16] = { 0, 0, 0, 0, 0, 0, 0, 0, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7] };
    mp_limb_t qq[8+1];
    mpn_tdiv_qr(qq, z, 0, t, 16, p, 8);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_mgy_decode)
/* x, z, invR, and p have 7.5 words.
 * Assuming R=W^8 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_7_5_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[16];
    mpfq_fixmp_7_5_mul(t, x, invR);
    mpfq_fixmp_7_5_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_7_5_lshift)
/* a has 7.5 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_7_5_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 8-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_7_5_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_add)
/* x, y, and z have 8.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 8_5_invmod, 8_5_redc, 8_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_8_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_add) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_sub)
/* x, y, and z have 8.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 8_5_invmod, 8_5_redc, 8_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_8_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_sub) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_add_nc)
/* x, y, and z have 8.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_8_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_sub_nc)
/* x, y, and z have 8.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_8_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_add_ui)
/* x, y, and z have 8.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_8_5_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_sub_ui)
/* x, y, and z have 8.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_8_5_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_add_ui_nc)
/* x, y, and z have 8.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_8_5_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_sub_ui_nc)
/* x, y, and z have 8.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_8_5_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_cmp)
/* x and y have 8.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 8_5_invmod, 8_5_redc, 8_5_redc_ur */
static inline
int mpfq_fixmp_8_5_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 9-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_cmp_ui)
/* x has 8.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 8_5_invmod */
static inline
int mpfq_fixmp_8_5_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 9-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_addmul1)
/* x has 8.5 words, z has 10.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_8_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    z[9] += carry;
    return (z[9]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_addmul1_nc)
/* x has 8.5 words, z has 10.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 8_5_mul, 8_5_mgy_decode */
static inline
void mpfq_fixmp_8_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    z[9] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_addmul1_shortz)
/* x has 8.5 words, z has 9.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 8_5_redc, 8_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_8_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_addmul05_nc)
/* x has 8.5 words, z has 9. c is 0.5 word.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 8_5_mul, 8_5_mgy_decode */
static inline
void mpfq_fixmp_8_5_addmul05_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    lo = c*x[8] + carry;
    assert(lo >= carry);
    z[8] += lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_addmul05_nc) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_mul)
/* x and y have 8.5 words, z has 17. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 8_5_mgy_decode */
static inline
void mpfq_fixmp_8_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 17; z[i++] = 0) ;
    mpfq_fixmp_8_5_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_8_5_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_8_5_addmul1_nc (z + 2, x, y[2]);
    mpfq_fixmp_8_5_addmul1_nc (z + 3, x, y[3]);
    mpfq_fixmp_8_5_addmul1_nc (z + 4, x, y[4]);
    mpfq_fixmp_8_5_addmul1_nc (z + 5, x, y[5]);
    mpfq_fixmp_8_5_addmul1_nc (z + 6, x, y[6]);
    mpfq_fixmp_8_5_addmul1_nc (z + 7, x, y[7]);
    mpfq_fixmp_8_5_addmul05_nc (z + 8, x, y[8]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_mul) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_sqr)
/* x has 8.5 words, z has 17. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_8_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[17] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_fixmp_3_addmul1_nc(buf + 3, x, x[3]);
    mpfq_fixmp_4_addmul1_nc(buf + 4, x, x[4]);
    mpfq_fixmp_5_addmul1_nc(buf + 5, x, x[5]);
    mpfq_fixmp_6_addmul1_nc(buf + 6, x, x[6]);
    mpfq_fixmp_7_addmul1_nc(buf + 7, x, x[7]);
    mpfq_fixmp_8_addmul1_nc(buf + 8, x, x[8]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpfq_umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
    mpfq_umul_ppmm(z[2*3+1], z[2*3], x[3], x[3]);
    mpfq_umul_ppmm(z[2*4+1], z[2*4], x[4], x[4]);
    mpfq_umul_ppmm(z[2*5+1], z[2*5], x[5], x[5]);
    mpfq_umul_ppmm(z[2*6+1], z[2*6], x[6], x[6]);
    mpfq_umul_ppmm(z[2*7+1], z[2*7], x[7], x[7]);
    z[2*8] = x[8] * x[8];
    mpn_lshift(buf, buf, 17, 1);
    mpn_add_n(z, z, buf, 17);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_mul1)
/* x has 8.5 words, z has 10. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_8_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[8] = lo;
    z[9] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_shortmul)
/* x and y have 8.5 words, z has 9.
 * Put the low 9 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_8_5_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 9);
    mpfq_fixmp_8_addmul1_nc (z+0, x, y[0]);
    z[9-1] += x[8]*y[0];
    mpfq_fixmp_7_addmul1_nc (z+1, x, y[1]);
    z[9-1] += x[7]*y[1];
    mpfq_fixmp_6_addmul1_nc (z+2, x, y[2]);
    z[9-1] += x[6]*y[2];
    mpfq_fixmp_5_addmul1_nc (z+3, x, y[3]);
    z[9-1] += x[5]*y[3];
    mpfq_fixmp_4_addmul1_nc (z+4, x, y[4]);
    z[9-1] += x[4]*y[4];
    mpfq_fixmp_3_addmul1_nc (z+5, x, y[5]);
    z[9-1] += x[3]*y[5];
    mpfq_fixmp_2_addmul1_nc (z+6, x, y[6]);
    z[9-1] += x[2]*y[6];
    mpfq_fixmp_1_addmul1_nc (z+7, x, y[7]);
    z[9-1] += x[1]*y[7];
    z[9-1] += x[0]*y[9-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_addmul05)
/* x has 8.5 words, z has 9. c is 0.5 word.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul05 */
static inline
mp_limb_t mpfq_fixmp_8_5_addmul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    lo = c*x[8] + carry;
    assert(lo >= carry);
    z[8] += lo;
    return z[8] < lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_addmul05) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_mul05)
/* x has 8.5 words, z has 9. c is 0.5 word.
 * Put (x*c) in z. No carry. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul05 */
static inline
void mpfq_fixmp_8_5_mul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[7] = lo;
    lo = c*x[8] + carry;
    assert(lo >= carry);
    z[8] = lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_mul05) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_mod)
/* x has 17 words. z and p have 8.5 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 8_5_mgy_decode */
static inline
void mpfq_fixmp_8_5_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[8+1], r[9];
    assert (p[9-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 17, p, 9);
    mpfq_copy(z, r, 9);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_mod) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_rshift)
/* a has 8.5 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 8_5_invmod */
static inline
void mpfq_fixmp_8_5_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 9-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[9-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_long_rshift)
/* a has 8.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 8_5_invmod */
static inline
void mpfq_fixmp_8_5_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 9 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[9-off-1] = a[9-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 9 - off);
    }
    mpfq_zero(a + 9 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_long_lshift)
/* a has 8.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 8_5_invmod */
static inline
void mpfq_fixmp_8_5_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 9-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 9 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_invmod)
/* x, z, and p have 8.5 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_8_5_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[9], v[9], a[9], b[9], fix[9];
      int i, t, lsh;
    
      mpfq_zero(u, 9);
      mpfq_zero(v, 9);
      mpfq_copy(a, x, 9);
      mpfq_copy(b, p, 9);
      u[0] = 1UL;
      
      if (mpfq_fixmp_8_5_cmp(a, v) == 0 || mpfq_fixmp_8_5_cmp(a, b) == 0) {
        mpfq_zero(res, 9);
        return 0;
      }
    
      mpfq_fixmp_8_5_add(fix, b, u);
      mpfq_fixmp_8_5_rshift(fix, 1);
    
      assert (mpfq_fixmp_8_5_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 9);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_8_5_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_8_5_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_8_5_sub(b, b, a);
          mpfq_fixmp_8_5_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 9);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_8_5_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_8_5_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_8_5_cmp(a,b) < 0);
        if (mpfq_fixmp_8_5_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_8_5_sub(a, a, b);
          mpfq_fixmp_8_5_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 9);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_8_5_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_8_5_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_8_5_cmp(b,a)<0);
      } while (mpfq_fixmp_8_5_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_8_5_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 9);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_8_5_rshift(u, 1);
        if (sig)
          mpfq_fixmp_8_5_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 9);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_redc)
/* x has 17 words, z and p have 8.5 words.
 * only one word is read from invp.
 * Assuming R=W^9 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_8_5_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 9; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_8_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    {
        mp_limb_t ret[9] = { x[9], x[10], x[11], x[12], x[13], x[14], x[15], x[16], 0 };
        cy = mpfq_fixmp_8_5_add(x, x, ret);
    }
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_8_5_cmp(x, p) >= 0) {
        mpfq_fixmp_8_5_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 9);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_redc) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_redc_ur)
/* x has 18 words, z and p have 8.5 words.
 * only one word is read from invp.
 * Assuming R=W^9 is the redc modulus, we expect that x verifies:
 *  x < W*W^8.5*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_8_5_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[1];
    for(int i = 0; i < 9; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_8_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_8_5_add(x + 9, x, x + 9);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W^0.5*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W^0.5+1)*p
    */
    if (cy) {
        /* x'/R-p < W^0.5*p, which fits in n words. */
        mpfq_fixmp_8_5_sub(x + 9, x + 9, p);
    }
    mpn_tdiv_qr(q, z, 0, x + 9, 9, p, 9);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_mgy_encode)
/* x, z, and p have 8.5 words.
 * Assuming R=W^9 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_8_5_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[18] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8] };
    mp_limb_t qq[9+1];
    mpn_tdiv_qr(qq, z, 0, t, 18, p, 9);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_mgy_decode)
/* x, z, invR, and p have 8.5 words.
 * Assuming R=W^9 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_8_5_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[18];
    mpfq_fixmp_8_5_mul(t, x, invR);
    mpfq_fixmp_8_5_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_8_5_lshift)
/* a has 8.5 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_8_5_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 9-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_8_5_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_add)
/* x, y, and z have 9.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 9_5_invmod, 9_5_redc, 9_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_9_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r + y[9];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[9] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_add) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_sub)
/* x, y, and z have 9.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 9_5_invmod, 9_5_redc, 9_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_9_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r - y[9];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[9] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_sub) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_add_nc)
/* x, y, and z have 9.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_9_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r + y[9];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[9] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_sub_nc)
/* x, y, and z have 9.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_9_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r - y[9];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[9] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_add_ui)
/* x, y, and z have 9.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_9_5_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    s = x[9];
    t = s + cy;
    cy = t < s;
    z[9] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_sub_ui)
/* x, y, and z have 9.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_9_5_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    s = x[9];
    t = s - cy;
    cy = t > s;
    z[9] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_add_ui_nc)
/* x, y, and z have 9.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_9_5_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    s = x[9];
    t = s + cy;
    cy = t < s;
    z[9] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_sub_ui_nc)
/* x, y, and z have 9.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_9_5_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    s = x[9];
    t = s - cy;
    cy = t > s;
    z[9] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_cmp)
/* x and y have 9.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 9_5_invmod, 9_5_redc, 9_5_redc_ur */
static inline
int mpfq_fixmp_9_5_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 10-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_cmp_ui)
/* x has 9.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 9_5_invmod */
static inline
int mpfq_fixmp_9_5_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 10-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_addmul1)
/* x has 9.5 words, z has 11.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_9_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    z[10] += carry;
    return (z[10]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_addmul1_nc)
/* x has 9.5 words, z has 11.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 9_5_mul, 9_5_mgy_decode */
static inline
void mpfq_fixmp_9_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    z[10] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_addmul1_shortz)
/* x has 9.5 words, z has 10.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 9_5_redc, 9_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_9_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_addmul05_nc)
/* x has 9.5 words, z has 10. c is 0.5 word.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 9_5_mul, 9_5_mgy_decode */
static inline
void mpfq_fixmp_9_5_addmul05_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    lo = c*x[9] + carry;
    assert(lo >= carry);
    z[9] += lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_addmul05_nc) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_mul)
/* x and y have 9.5 words, z has 19. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 9_5_mgy_decode */
static inline
void mpfq_fixmp_9_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 19; z[i++] = 0) ;
    mpfq_fixmp_9_5_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_9_5_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_9_5_addmul1_nc (z + 2, x, y[2]);
    mpfq_fixmp_9_5_addmul1_nc (z + 3, x, y[3]);
    mpfq_fixmp_9_5_addmul1_nc (z + 4, x, y[4]);
    mpfq_fixmp_9_5_addmul1_nc (z + 5, x, y[5]);
    mpfq_fixmp_9_5_addmul1_nc (z + 6, x, y[6]);
    mpfq_fixmp_9_5_addmul1_nc (z + 7, x, y[7]);
    mpfq_fixmp_9_5_addmul1_nc (z + 8, x, y[8]);
    mpfq_fixmp_9_5_addmul05_nc (z + 9, x, y[9]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_mul) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_sqr)
/* x has 9.5 words, z has 19. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_9_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[19] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_fixmp_3_addmul1_nc(buf + 3, x, x[3]);
    mpfq_fixmp_4_addmul1_nc(buf + 4, x, x[4]);
    mpfq_fixmp_5_addmul1_nc(buf + 5, x, x[5]);
    mpfq_fixmp_6_addmul1_nc(buf + 6, x, x[6]);
    mpfq_fixmp_7_addmul1_nc(buf + 7, x, x[7]);
    mpfq_fixmp_8_addmul1_nc(buf + 8, x, x[8]);
    mpfq_fixmp_9_addmul1_nc(buf + 9, x, x[9]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpfq_umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
    mpfq_umul_ppmm(z[2*3+1], z[2*3], x[3], x[3]);
    mpfq_umul_ppmm(z[2*4+1], z[2*4], x[4], x[4]);
    mpfq_umul_ppmm(z[2*5+1], z[2*5], x[5], x[5]);
    mpfq_umul_ppmm(z[2*6+1], z[2*6], x[6], x[6]);
    mpfq_umul_ppmm(z[2*7+1], z[2*7], x[7], x[7]);
    mpfq_umul_ppmm(z[2*8+1], z[2*8], x[8], x[8]);
    z[2*9] = x[9] * x[9];
    mpn_lshift(buf, buf, 19, 1);
    mpn_add_n(z, z, buf, 19);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_mul1)
/* x has 9.5 words, z has 11. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_9_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[9] = lo;
    z[10] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_shortmul)
/* x and y have 9.5 words, z has 10.
 * Put the low 10 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_9_5_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 10);
    mpfq_fixmp_9_addmul1_nc (z+0, x, y[0]);
    z[10-1] += x[9]*y[0];
    mpfq_fixmp_8_addmul1_nc (z+1, x, y[1]);
    z[10-1] += x[8]*y[1];
    mpfq_fixmp_7_addmul1_nc (z+2, x, y[2]);
    z[10-1] += x[7]*y[2];
    mpfq_fixmp_6_addmul1_nc (z+3, x, y[3]);
    z[10-1] += x[6]*y[3];
    mpfq_fixmp_5_addmul1_nc (z+4, x, y[4]);
    z[10-1] += x[5]*y[4];
    mpfq_fixmp_4_addmul1_nc (z+5, x, y[5]);
    z[10-1] += x[4]*y[5];
    mpfq_fixmp_3_addmul1_nc (z+6, x, y[6]);
    z[10-1] += x[3]*y[6];
    mpfq_fixmp_2_addmul1_nc (z+7, x, y[7]);
    z[10-1] += x[2]*y[7];
    mpfq_fixmp_1_addmul1_nc (z+8, x, y[8]);
    z[10-1] += x[1]*y[8];
    z[10-1] += x[0]*y[10-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_addmul05)
/* x has 9.5 words, z has 10. c is 0.5 word.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul05 */
static inline
mp_limb_t mpfq_fixmp_9_5_addmul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    lo = c*x[9] + carry;
    assert(lo >= carry);
    z[9] += lo;
    return z[9] < lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_addmul05) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_mul05)
/* x has 9.5 words, z has 10. c is 0.5 word.
 * Put (x*c) in z. No carry. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul05 */
static inline
void mpfq_fixmp_9_5_mul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[8] = lo;
    lo = c*x[9] + carry;
    assert(lo >= carry);
    z[9] = lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_mul05) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_mod)
/* x has 19 words. z and p have 9.5 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 9_5_mgy_decode */
static inline
void mpfq_fixmp_9_5_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[9+1], r[10];
    assert (p[10-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 19, p, 10);
    mpfq_copy(z, r, 10);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_mod) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_rshift)
/* a has 9.5 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 9_5_invmod */
static inline
void mpfq_fixmp_9_5_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 10-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[10-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_long_rshift)
/* a has 9.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 9_5_invmod */
static inline
void mpfq_fixmp_9_5_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 10 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[10-off-1] = a[10-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 10 - off);
    }
    mpfq_zero(a + 10 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_long_lshift)
/* a has 9.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 9_5_invmod */
static inline
void mpfq_fixmp_9_5_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 10-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 10 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_invmod)
/* x, z, and p have 9.5 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_9_5_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[10], v[10], a[10], b[10], fix[10];
      int i, t, lsh;
    
      mpfq_zero(u, 10);
      mpfq_zero(v, 10);
      mpfq_copy(a, x, 10);
      mpfq_copy(b, p, 10);
      u[0] = 1UL;
      
      if (mpfq_fixmp_9_5_cmp(a, v) == 0 || mpfq_fixmp_9_5_cmp(a, b) == 0) {
        mpfq_zero(res, 10);
        return 0;
      }
    
      mpfq_fixmp_9_5_add(fix, b, u);
      mpfq_fixmp_9_5_rshift(fix, 1);
    
      assert (mpfq_fixmp_9_5_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 10);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_9_5_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_9_5_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_9_5_sub(b, b, a);
          mpfq_fixmp_9_5_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 10);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_9_5_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_9_5_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_9_5_cmp(a,b) < 0);
        if (mpfq_fixmp_9_5_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_9_5_sub(a, a, b);
          mpfq_fixmp_9_5_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 10);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_9_5_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_9_5_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_9_5_cmp(b,a)<0);
      } while (mpfq_fixmp_9_5_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_9_5_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 10);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_9_5_rshift(u, 1);
        if (sig)
          mpfq_fixmp_9_5_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 10);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_redc)
/* x has 19 words, z and p have 9.5 words.
 * only one word is read from invp.
 * Assuming R=W^10 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_9_5_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 10; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_9_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    {
        mp_limb_t ret[10] = { x[10], x[11], x[12], x[13], x[14], x[15], x[16], x[17], x[18], 0 };
        cy = mpfq_fixmp_9_5_add(x, x, ret);
    }
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_9_5_cmp(x, p) >= 0) {
        mpfq_fixmp_9_5_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 10);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_redc) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_redc_ur)
/* x has 20 words, z and p have 9.5 words.
 * only one word is read from invp.
 * Assuming R=W^10 is the redc modulus, we expect that x verifies:
 *  x < W*W^9.5*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_9_5_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[1];
    for(int i = 0; i < 10; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_9_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_9_5_add(x + 10, x, x + 10);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W^0.5*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W^0.5+1)*p
    */
    if (cy) {
        /* x'/R-p < W^0.5*p, which fits in n words. */
        mpfq_fixmp_9_5_sub(x + 10, x + 10, p);
    }
    mpn_tdiv_qr(q, z, 0, x + 10, 10, p, 10);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_mgy_encode)
/* x, z, and p have 9.5 words.
 * Assuming R=W^10 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_9_5_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[20] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9] };
    mp_limb_t qq[10+1];
    mpn_tdiv_qr(qq, z, 0, t, 20, p, 10);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_mgy_decode)
/* x, z, invR, and p have 9.5 words.
 * Assuming R=W^10 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_9_5_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[20];
    mpfq_fixmp_9_5_mul(t, x, invR);
    mpfq_fixmp_9_5_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_9_5_lshift)
/* a has 9.5 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_9_5_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 10-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_9_5_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_add)
/* x, y, and z have 10.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 10_5_invmod, 10_5_redc, 10_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_10_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r + y[9];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r + y[10];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[10] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_add) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_sub)
/* x, y, and z have 10.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 10_5_invmod, 10_5_redc, 10_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_10_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r - y[9];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r - y[10];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[10] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_sub) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_add_nc)
/* x, y, and z have 10.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_10_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r + y[9];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r + y[10];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[10] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_sub_nc)
/* x, y, and z have 10.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_10_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r - y[9];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r - y[10];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[10] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_add_ui)
/* x, y, and z have 10.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_10_5_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    s = x[9];
    t = s + cy;
    cy = t < s;
    z[9] = t;
    s = x[10];
    t = s + cy;
    cy = t < s;
    z[10] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_sub_ui)
/* x, y, and z have 10.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_10_5_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    s = x[9];
    t = s - cy;
    cy = t > s;
    z[9] = t;
    s = x[10];
    t = s - cy;
    cy = t > s;
    z[10] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_add_ui_nc)
/* x, y, and z have 10.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_10_5_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    s = x[9];
    t = s + cy;
    cy = t < s;
    z[9] = t;
    s = x[10];
    t = s + cy;
    cy = t < s;
    z[10] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_sub_ui_nc)
/* x, y, and z have 10.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_10_5_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    s = x[9];
    t = s - cy;
    cy = t > s;
    z[9] = t;
    s = x[10];
    t = s - cy;
    cy = t > s;
    z[10] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_cmp)
/* x and y have 10.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 10_5_invmod, 10_5_redc, 10_5_redc_ur */
static inline
int mpfq_fixmp_10_5_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 11-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_cmp_ui)
/* x has 10.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 10_5_invmod */
static inline
int mpfq_fixmp_10_5_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 11-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_addmul1)
/* x has 10.5 words, z has 12.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_10_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    z[11] += carry;
    return (z[11]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_addmul1_nc)
/* x has 10.5 words, z has 12.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 10_5_mul, 10_5_mgy_decode */
static inline
void mpfq_fixmp_10_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    z[11] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_addmul1_shortz)
/* x has 10.5 words, z has 11.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 10_5_redc, 10_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_10_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_addmul05_nc)
/* x has 10.5 words, z has 11. c is 0.5 word.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 10_5_mul, 10_5_mgy_decode */
static inline
void mpfq_fixmp_10_5_addmul05_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    lo = c*x[10] + carry;
    assert(lo >= carry);
    z[10] += lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_addmul05_nc) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_mul)
/* x and y have 10.5 words, z has 21. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 10_5_mgy_decode */
static inline
void mpfq_fixmp_10_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 21; z[i++] = 0) ;
    mpfq_fixmp_10_5_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_10_5_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_10_5_addmul1_nc (z + 2, x, y[2]);
    mpfq_fixmp_10_5_addmul1_nc (z + 3, x, y[3]);
    mpfq_fixmp_10_5_addmul1_nc (z + 4, x, y[4]);
    mpfq_fixmp_10_5_addmul1_nc (z + 5, x, y[5]);
    mpfq_fixmp_10_5_addmul1_nc (z + 6, x, y[6]);
    mpfq_fixmp_10_5_addmul1_nc (z + 7, x, y[7]);
    mpfq_fixmp_10_5_addmul1_nc (z + 8, x, y[8]);
    mpfq_fixmp_10_5_addmul1_nc (z + 9, x, y[9]);
    mpfq_fixmp_10_5_addmul05_nc (z + 10, x, y[10]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_mul) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_sqr)
/* x has 10.5 words, z has 21. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_10_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[21] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_fixmp_3_addmul1_nc(buf + 3, x, x[3]);
    mpfq_fixmp_4_addmul1_nc(buf + 4, x, x[4]);
    mpfq_fixmp_5_addmul1_nc(buf + 5, x, x[5]);
    mpfq_fixmp_6_addmul1_nc(buf + 6, x, x[6]);
    mpfq_fixmp_7_addmul1_nc(buf + 7, x, x[7]);
    mpfq_fixmp_8_addmul1_nc(buf + 8, x, x[8]);
    mpfq_fixmp_9_addmul1_nc(buf + 9, x, x[9]);
    mpfq_fixmp_10_addmul1_nc(buf + 10, x, x[10]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpfq_umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
    mpfq_umul_ppmm(z[2*3+1], z[2*3], x[3], x[3]);
    mpfq_umul_ppmm(z[2*4+1], z[2*4], x[4], x[4]);
    mpfq_umul_ppmm(z[2*5+1], z[2*5], x[5], x[5]);
    mpfq_umul_ppmm(z[2*6+1], z[2*6], x[6], x[6]);
    mpfq_umul_ppmm(z[2*7+1], z[2*7], x[7], x[7]);
    mpfq_umul_ppmm(z[2*8+1], z[2*8], x[8], x[8]);
    mpfq_umul_ppmm(z[2*9+1], z[2*9], x[9], x[9]);
    z[2*10] = x[10] * x[10];
    mpn_lshift(buf, buf, 21, 1);
    mpn_add_n(z, z, buf, 21);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_mul1)
/* x has 10.5 words, z has 12. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_10_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[10] = lo;
    z[11] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_shortmul)
/* x and y have 10.5 words, z has 11.
 * Put the low 11 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_10_5_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 11);
    mpfq_fixmp_10_addmul1_nc (z+0, x, y[0]);
    z[11-1] += x[10]*y[0];
    mpfq_fixmp_9_addmul1_nc (z+1, x, y[1]);
    z[11-1] += x[9]*y[1];
    mpfq_fixmp_8_addmul1_nc (z+2, x, y[2]);
    z[11-1] += x[8]*y[2];
    mpfq_fixmp_7_addmul1_nc (z+3, x, y[3]);
    z[11-1] += x[7]*y[3];
    mpfq_fixmp_6_addmul1_nc (z+4, x, y[4]);
    z[11-1] += x[6]*y[4];
    mpfq_fixmp_5_addmul1_nc (z+5, x, y[5]);
    z[11-1] += x[5]*y[5];
    mpfq_fixmp_4_addmul1_nc (z+6, x, y[6]);
    z[11-1] += x[4]*y[6];
    mpfq_fixmp_3_addmul1_nc (z+7, x, y[7]);
    z[11-1] += x[3]*y[7];
    mpfq_fixmp_2_addmul1_nc (z+8, x, y[8]);
    z[11-1] += x[2]*y[8];
    mpfq_fixmp_1_addmul1_nc (z+9, x, y[9]);
    z[11-1] += x[1]*y[9];
    z[11-1] += x[0]*y[11-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_addmul05)
/* x has 10.5 words, z has 11. c is 0.5 word.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul05 */
static inline
mp_limb_t mpfq_fixmp_10_5_addmul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    lo = c*x[10] + carry;
    assert(lo >= carry);
    z[10] += lo;
    return z[10] < lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_addmul05) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_mul05)
/* x has 10.5 words, z has 11. c is 0.5 word.
 * Put (x*c) in z. No carry. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul05 */
static inline
void mpfq_fixmp_10_5_mul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[9] = lo;
    lo = c*x[10] + carry;
    assert(lo >= carry);
    z[10] = lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_mul05) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_mod)
/* x has 21 words. z and p have 10.5 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 10_5_mgy_decode */
static inline
void mpfq_fixmp_10_5_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[10+1], r[11];
    assert (p[11-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 21, p, 11);
    mpfq_copy(z, r, 11);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_mod) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_rshift)
/* a has 10.5 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 10_5_invmod */
static inline
void mpfq_fixmp_10_5_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 11-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[11-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_long_rshift)
/* a has 10.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 10_5_invmod */
static inline
void mpfq_fixmp_10_5_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 11 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[11-off-1] = a[11-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 11 - off);
    }
    mpfq_zero(a + 11 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_long_lshift)
/* a has 10.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 10_5_invmod */
static inline
void mpfq_fixmp_10_5_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 11-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 11 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_invmod)
/* x, z, and p have 10.5 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_10_5_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[11], v[11], a[11], b[11], fix[11];
      int i, t, lsh;
    
      mpfq_zero(u, 11);
      mpfq_zero(v, 11);
      mpfq_copy(a, x, 11);
      mpfq_copy(b, p, 11);
      u[0] = 1UL;
      
      if (mpfq_fixmp_10_5_cmp(a, v) == 0 || mpfq_fixmp_10_5_cmp(a, b) == 0) {
        mpfq_zero(res, 11);
        return 0;
      }
    
      mpfq_fixmp_10_5_add(fix, b, u);
      mpfq_fixmp_10_5_rshift(fix, 1);
    
      assert (mpfq_fixmp_10_5_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 11);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_10_5_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_10_5_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_10_5_sub(b, b, a);
          mpfq_fixmp_10_5_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 11);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_10_5_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_10_5_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_10_5_cmp(a,b) < 0);
        if (mpfq_fixmp_10_5_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_10_5_sub(a, a, b);
          mpfq_fixmp_10_5_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 11);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_10_5_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_10_5_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_10_5_cmp(b,a)<0);
      } while (mpfq_fixmp_10_5_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_10_5_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 11);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_10_5_rshift(u, 1);
        if (sig)
          mpfq_fixmp_10_5_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 11);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_redc)
/* x has 21 words, z and p have 10.5 words.
 * only one word is read from invp.
 * Assuming R=W^11 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_10_5_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 11; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_10_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    {
        mp_limb_t ret[11] = { x[11], x[12], x[13], x[14], x[15], x[16], x[17], x[18], x[19], x[20], 0 };
        cy = mpfq_fixmp_10_5_add(x, x, ret);
    }
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_10_5_cmp(x, p) >= 0) {
        mpfq_fixmp_10_5_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 11);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_redc) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_redc_ur)
/* x has 22 words, z and p have 10.5 words.
 * only one word is read from invp.
 * Assuming R=W^11 is the redc modulus, we expect that x verifies:
 *  x < W*W^10.5*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_10_5_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[1];
    for(int i = 0; i < 11; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_10_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_10_5_add(x + 11, x, x + 11);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W^0.5*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W^0.5+1)*p
    */
    if (cy) {
        /* x'/R-p < W^0.5*p, which fits in n words. */
        mpfq_fixmp_10_5_sub(x + 11, x + 11, p);
    }
    mpn_tdiv_qr(q, z, 0, x + 11, 11, p, 11);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_mgy_encode)
/* x, z, and p have 10.5 words.
 * Assuming R=W^11 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_10_5_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[22] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10] };
    mp_limb_t qq[11+1];
    mpn_tdiv_qr(qq, z, 0, t, 22, p, 11);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_mgy_decode)
/* x, z, invR, and p have 10.5 words.
 * Assuming R=W^11 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_10_5_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[22];
    mpfq_fixmp_10_5_mul(t, x, invR);
    mpfq_fixmp_10_5_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_10_5_lshift)
/* a has 10.5 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_10_5_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 11-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_10_5_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_add)
/* x, y, and z have 11.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 11_5_invmod, 11_5_redc, 11_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_11_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r + y[9];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r + y[10];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r + y[11];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[11] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_add) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_sub)
/* x, y, and z have 11.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 11_5_invmod, 11_5_redc, 11_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_11_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r - y[9];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r - y[10];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r - y[11];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[11] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_sub) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_add_nc)
/* x, y, and z have 11.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_11_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r + y[9];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r + y[10];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r + y[11];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[11] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_sub_nc)
/* x, y, and z have 11.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_11_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r - y[9];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r - y[10];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r - y[11];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[11] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_add_ui)
/* x, y, and z have 11.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_11_5_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    s = x[9];
    t = s + cy;
    cy = t < s;
    z[9] = t;
    s = x[10];
    t = s + cy;
    cy = t < s;
    z[10] = t;
    s = x[11];
    t = s + cy;
    cy = t < s;
    z[11] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_sub_ui)
/* x, y, and z have 11.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_11_5_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    s = x[9];
    t = s - cy;
    cy = t > s;
    z[9] = t;
    s = x[10];
    t = s - cy;
    cy = t > s;
    z[10] = t;
    s = x[11];
    t = s - cy;
    cy = t > s;
    z[11] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_add_ui_nc)
/* x, y, and z have 11.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_11_5_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    s = x[9];
    t = s + cy;
    cy = t < s;
    z[9] = t;
    s = x[10];
    t = s + cy;
    cy = t < s;
    z[10] = t;
    s = x[11];
    t = s + cy;
    cy = t < s;
    z[11] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_sub_ui_nc)
/* x, y, and z have 11.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_11_5_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    s = x[9];
    t = s - cy;
    cy = t > s;
    z[9] = t;
    s = x[10];
    t = s - cy;
    cy = t > s;
    z[10] = t;
    s = x[11];
    t = s - cy;
    cy = t > s;
    z[11] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_cmp)
/* x and y have 11.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 11_5_invmod, 11_5_redc, 11_5_redc_ur */
static inline
int mpfq_fixmp_11_5_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 12-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_cmp_ui)
/* x has 11.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 11_5_invmod */
static inline
int mpfq_fixmp_11_5_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 12-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_addmul1)
/* x has 11.5 words, z has 13.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_11_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    z[12] += carry;
    return (z[12]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_addmul1_nc)
/* x has 11.5 words, z has 13.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 11_5_mul, 11_5_mgy_decode */
static inline
void mpfq_fixmp_11_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    z[12] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_addmul1_shortz)
/* x has 11.5 words, z has 12.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 11_5_redc, 11_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_11_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_addmul05_nc)
/* x has 11.5 words, z has 12. c is 0.5 word.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 11_5_mul, 11_5_mgy_decode */
static inline
void mpfq_fixmp_11_5_addmul05_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    lo = c*x[11] + carry;
    assert(lo >= carry);
    z[11] += lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_addmul05_nc) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_mul)
/* x and y have 11.5 words, z has 23. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 11_5_mgy_decode */
static inline
void mpfq_fixmp_11_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 23; z[i++] = 0) ;
    mpfq_fixmp_11_5_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_11_5_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_11_5_addmul1_nc (z + 2, x, y[2]);
    mpfq_fixmp_11_5_addmul1_nc (z + 3, x, y[3]);
    mpfq_fixmp_11_5_addmul1_nc (z + 4, x, y[4]);
    mpfq_fixmp_11_5_addmul1_nc (z + 5, x, y[5]);
    mpfq_fixmp_11_5_addmul1_nc (z + 6, x, y[6]);
    mpfq_fixmp_11_5_addmul1_nc (z + 7, x, y[7]);
    mpfq_fixmp_11_5_addmul1_nc (z + 8, x, y[8]);
    mpfq_fixmp_11_5_addmul1_nc (z + 9, x, y[9]);
    mpfq_fixmp_11_5_addmul1_nc (z + 10, x, y[10]);
    mpfq_fixmp_11_5_addmul05_nc (z + 11, x, y[11]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_mul) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_sqr)
/* x has 11.5 words, z has 23. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_11_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[23] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_fixmp_3_addmul1_nc(buf + 3, x, x[3]);
    mpfq_fixmp_4_addmul1_nc(buf + 4, x, x[4]);
    mpfq_fixmp_5_addmul1_nc(buf + 5, x, x[5]);
    mpfq_fixmp_6_addmul1_nc(buf + 6, x, x[6]);
    mpfq_fixmp_7_addmul1_nc(buf + 7, x, x[7]);
    mpfq_fixmp_8_addmul1_nc(buf + 8, x, x[8]);
    mpfq_fixmp_9_addmul1_nc(buf + 9, x, x[9]);
    mpfq_fixmp_10_addmul1_nc(buf + 10, x, x[10]);
    mpfq_fixmp_11_addmul1_nc(buf + 11, x, x[11]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpfq_umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
    mpfq_umul_ppmm(z[2*3+1], z[2*3], x[3], x[3]);
    mpfq_umul_ppmm(z[2*4+1], z[2*4], x[4], x[4]);
    mpfq_umul_ppmm(z[2*5+1], z[2*5], x[5], x[5]);
    mpfq_umul_ppmm(z[2*6+1], z[2*6], x[6], x[6]);
    mpfq_umul_ppmm(z[2*7+1], z[2*7], x[7], x[7]);
    mpfq_umul_ppmm(z[2*8+1], z[2*8], x[8], x[8]);
    mpfq_umul_ppmm(z[2*9+1], z[2*9], x[9], x[9]);
    mpfq_umul_ppmm(z[2*10+1], z[2*10], x[10], x[10]);
    z[2*11] = x[11] * x[11];
    mpn_lshift(buf, buf, 23, 1);
    mpn_add_n(z, z, buf, 23);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_mul1)
/* x has 11.5 words, z has 13. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_11_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[11] = lo;
    z[12] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_shortmul)
/* x and y have 11.5 words, z has 12.
 * Put the low 12 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_11_5_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 12);
    mpfq_fixmp_11_addmul1_nc (z+0, x, y[0]);
    z[12-1] += x[11]*y[0];
    mpfq_fixmp_10_addmul1_nc (z+1, x, y[1]);
    z[12-1] += x[10]*y[1];
    mpfq_fixmp_9_addmul1_nc (z+2, x, y[2]);
    z[12-1] += x[9]*y[2];
    mpfq_fixmp_8_addmul1_nc (z+3, x, y[3]);
    z[12-1] += x[8]*y[3];
    mpfq_fixmp_7_addmul1_nc (z+4, x, y[4]);
    z[12-1] += x[7]*y[4];
    mpfq_fixmp_6_addmul1_nc (z+5, x, y[5]);
    z[12-1] += x[6]*y[5];
    mpfq_fixmp_5_addmul1_nc (z+6, x, y[6]);
    z[12-1] += x[5]*y[6];
    mpfq_fixmp_4_addmul1_nc (z+7, x, y[7]);
    z[12-1] += x[4]*y[7];
    mpfq_fixmp_3_addmul1_nc (z+8, x, y[8]);
    z[12-1] += x[3]*y[8];
    mpfq_fixmp_2_addmul1_nc (z+9, x, y[9]);
    z[12-1] += x[2]*y[9];
    mpfq_fixmp_1_addmul1_nc (z+10, x, y[10]);
    z[12-1] += x[1]*y[10];
    z[12-1] += x[0]*y[12-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_addmul05)
/* x has 11.5 words, z has 12. c is 0.5 word.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul05 */
static inline
mp_limb_t mpfq_fixmp_11_5_addmul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    lo = c*x[11] + carry;
    assert(lo >= carry);
    z[11] += lo;
    return z[11] < lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_addmul05) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_mul05)
/* x has 11.5 words, z has 12. c is 0.5 word.
 * Put (x*c) in z. No carry. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul05 */
static inline
void mpfq_fixmp_11_5_mul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[10] = lo;
    lo = c*x[11] + carry;
    assert(lo >= carry);
    z[11] = lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_mul05) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_mod)
/* x has 23 words. z and p have 11.5 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 11_5_mgy_decode */
static inline
void mpfq_fixmp_11_5_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[11+1], r[12];
    assert (p[12-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 23, p, 12);
    mpfq_copy(z, r, 12);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_mod) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_rshift)
/* a has 11.5 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 11_5_invmod */
static inline
void mpfq_fixmp_11_5_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 12-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[12-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_long_rshift)
/* a has 11.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 11_5_invmod */
static inline
void mpfq_fixmp_11_5_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 12 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[12-off-1] = a[12-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 12 - off);
    }
    mpfq_zero(a + 12 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_long_lshift)
/* a has 11.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 11_5_invmod */
static inline
void mpfq_fixmp_11_5_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 12-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 12 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_invmod)
/* x, z, and p have 11.5 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_11_5_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[12], v[12], a[12], b[12], fix[12];
      int i, t, lsh;
    
      mpfq_zero(u, 12);
      mpfq_zero(v, 12);
      mpfq_copy(a, x, 12);
      mpfq_copy(b, p, 12);
      u[0] = 1UL;
      
      if (mpfq_fixmp_11_5_cmp(a, v) == 0 || mpfq_fixmp_11_5_cmp(a, b) == 0) {
        mpfq_zero(res, 12);
        return 0;
      }
    
      mpfq_fixmp_11_5_add(fix, b, u);
      mpfq_fixmp_11_5_rshift(fix, 1);
    
      assert (mpfq_fixmp_11_5_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 12);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_11_5_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_11_5_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_11_5_sub(b, b, a);
          mpfq_fixmp_11_5_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 12);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_11_5_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_11_5_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_11_5_cmp(a,b) < 0);
        if (mpfq_fixmp_11_5_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_11_5_sub(a, a, b);
          mpfq_fixmp_11_5_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 12);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_11_5_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_11_5_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_11_5_cmp(b,a)<0);
      } while (mpfq_fixmp_11_5_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_11_5_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 12);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_11_5_rshift(u, 1);
        if (sig)
          mpfq_fixmp_11_5_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 12);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_redc)
/* x has 23 words, z and p have 11.5 words.
 * only one word is read from invp.
 * Assuming R=W^12 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_11_5_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 12; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_11_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    {
        mp_limb_t ret[12] = { x[12], x[13], x[14], x[15], x[16], x[17], x[18], x[19], x[20], x[21], x[22], 0 };
        cy = mpfq_fixmp_11_5_add(x, x, ret);
    }
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_11_5_cmp(x, p) >= 0) {
        mpfq_fixmp_11_5_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 12);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_redc) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_redc_ur)
/* x has 24 words, z and p have 11.5 words.
 * only one word is read from invp.
 * Assuming R=W^12 is the redc modulus, we expect that x verifies:
 *  x < W*W^11.5*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_11_5_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[1];
    for(int i = 0; i < 12; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_11_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_11_5_add(x + 12, x, x + 12);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W^0.5*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W^0.5+1)*p
    */
    if (cy) {
        /* x'/R-p < W^0.5*p, which fits in n words. */
        mpfq_fixmp_11_5_sub(x + 12, x + 12, p);
    }
    mpn_tdiv_qr(q, z, 0, x + 12, 12, p, 12);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_mgy_encode)
/* x, z, and p have 11.5 words.
 * Assuming R=W^12 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_11_5_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[24] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11] };
    mp_limb_t qq[12+1];
    mpn_tdiv_qr(qq, z, 0, t, 24, p, 12);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_mgy_decode)
/* x, z, invR, and p have 11.5 words.
 * Assuming R=W^12 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_11_5_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[24];
    mpfq_fixmp_11_5_mul(t, x, invR);
    mpfq_fixmp_11_5_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_11_5_lshift)
/* a has 11.5 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_11_5_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 12-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_11_5_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_add)
/* x, y, and z have 12.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 12_5_invmod, 12_5_redc, 12_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_12_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r + y[9];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r + y[10];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r + y[11];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[11] = t;
    r = x[12];
    s = r + y[12];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[12] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_add) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_sub)
/* x, y, and z have 12.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 12_5_invmod, 12_5_redc, 12_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_12_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r - y[9];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r - y[10];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r - y[11];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[11] = t;
    r = x[12];
    s = r - y[12];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[12] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_sub) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_add_nc)
/* x, y, and z have 12.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_12_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r + y[9];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r + y[10];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r + y[11];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[11] = t;
    r = x[12];
    s = r + y[12];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[12] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_sub_nc)
/* x, y, and z have 12.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_12_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r - y[9];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r - y[10];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r - y[11];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[11] = t;
    r = x[12];
    s = r - y[12];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[12] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_add_ui)
/* x, y, and z have 12.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_12_5_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    s = x[9];
    t = s + cy;
    cy = t < s;
    z[9] = t;
    s = x[10];
    t = s + cy;
    cy = t < s;
    z[10] = t;
    s = x[11];
    t = s + cy;
    cy = t < s;
    z[11] = t;
    s = x[12];
    t = s + cy;
    cy = t < s;
    z[12] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_sub_ui)
/* x, y, and z have 12.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_12_5_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    s = x[9];
    t = s - cy;
    cy = t > s;
    z[9] = t;
    s = x[10];
    t = s - cy;
    cy = t > s;
    z[10] = t;
    s = x[11];
    t = s - cy;
    cy = t > s;
    z[11] = t;
    s = x[12];
    t = s - cy;
    cy = t > s;
    z[12] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_add_ui_nc)
/* x, y, and z have 12.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_12_5_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    s = x[9];
    t = s + cy;
    cy = t < s;
    z[9] = t;
    s = x[10];
    t = s + cy;
    cy = t < s;
    z[10] = t;
    s = x[11];
    t = s + cy;
    cy = t < s;
    z[11] = t;
    s = x[12];
    t = s + cy;
    cy = t < s;
    z[12] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_sub_ui_nc)
/* x, y, and z have 12.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_12_5_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    s = x[9];
    t = s - cy;
    cy = t > s;
    z[9] = t;
    s = x[10];
    t = s - cy;
    cy = t > s;
    z[10] = t;
    s = x[11];
    t = s - cy;
    cy = t > s;
    z[11] = t;
    s = x[12];
    t = s - cy;
    cy = t > s;
    z[12] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_cmp)
/* x and y have 12.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 12_5_invmod, 12_5_redc, 12_5_redc_ur */
static inline
int mpfq_fixmp_12_5_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 13-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_cmp_ui)
/* x has 12.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 12_5_invmod */
static inline
int mpfq_fixmp_12_5_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 13-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_addmul1)
/* x has 12.5 words, z has 14.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_12_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[12];
    lo += buf;
    carry += (lo<buf);
    z[12] = lo;
    z[13] += carry;
    return (z[13]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_addmul1_nc)
/* x has 12.5 words, z has 14.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 12_5_mul, 12_5_mgy_decode */
static inline
void mpfq_fixmp_12_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[12];
    lo += buf;
    carry += (lo<buf);
    z[12] = lo;
    z[13] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_addmul1_shortz)
/* x has 12.5 words, z has 13.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 12_5_redc, 12_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_12_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[12];
    lo += buf;
    carry += (lo<buf);
    z[12] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_addmul05_nc)
/* x has 12.5 words, z has 13. c is 0.5 word.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 12_5_mul, 12_5_mgy_decode */
static inline
void mpfq_fixmp_12_5_addmul05_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    lo = c*x[12] + carry;
    assert(lo >= carry);
    z[12] += lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_addmul05_nc) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_mul)
/* x and y have 12.5 words, z has 25. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 12_5_mgy_decode */
static inline
void mpfq_fixmp_12_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 25; z[i++] = 0) ;
    mpfq_fixmp_12_5_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_12_5_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_12_5_addmul1_nc (z + 2, x, y[2]);
    mpfq_fixmp_12_5_addmul1_nc (z + 3, x, y[3]);
    mpfq_fixmp_12_5_addmul1_nc (z + 4, x, y[4]);
    mpfq_fixmp_12_5_addmul1_nc (z + 5, x, y[5]);
    mpfq_fixmp_12_5_addmul1_nc (z + 6, x, y[6]);
    mpfq_fixmp_12_5_addmul1_nc (z + 7, x, y[7]);
    mpfq_fixmp_12_5_addmul1_nc (z + 8, x, y[8]);
    mpfq_fixmp_12_5_addmul1_nc (z + 9, x, y[9]);
    mpfq_fixmp_12_5_addmul1_nc (z + 10, x, y[10]);
    mpfq_fixmp_12_5_addmul1_nc (z + 11, x, y[11]);
    mpfq_fixmp_12_5_addmul05_nc (z + 12, x, y[12]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_mul) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_sqr)
/* x has 12.5 words, z has 25. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_12_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[25] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_fixmp_3_addmul1_nc(buf + 3, x, x[3]);
    mpfq_fixmp_4_addmul1_nc(buf + 4, x, x[4]);
    mpfq_fixmp_5_addmul1_nc(buf + 5, x, x[5]);
    mpfq_fixmp_6_addmul1_nc(buf + 6, x, x[6]);
    mpfq_fixmp_7_addmul1_nc(buf + 7, x, x[7]);
    mpfq_fixmp_8_addmul1_nc(buf + 8, x, x[8]);
    mpfq_fixmp_9_addmul1_nc(buf + 9, x, x[9]);
    mpfq_fixmp_10_addmul1_nc(buf + 10, x, x[10]);
    mpfq_fixmp_11_addmul1_nc(buf + 11, x, x[11]);
    mpfq_fixmp_12_addmul1_nc(buf + 12, x, x[12]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpfq_umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
    mpfq_umul_ppmm(z[2*3+1], z[2*3], x[3], x[3]);
    mpfq_umul_ppmm(z[2*4+1], z[2*4], x[4], x[4]);
    mpfq_umul_ppmm(z[2*5+1], z[2*5], x[5], x[5]);
    mpfq_umul_ppmm(z[2*6+1], z[2*6], x[6], x[6]);
    mpfq_umul_ppmm(z[2*7+1], z[2*7], x[7], x[7]);
    mpfq_umul_ppmm(z[2*8+1], z[2*8], x[8], x[8]);
    mpfq_umul_ppmm(z[2*9+1], z[2*9], x[9], x[9]);
    mpfq_umul_ppmm(z[2*10+1], z[2*10], x[10], x[10]);
    mpfq_umul_ppmm(z[2*11+1], z[2*11], x[11], x[11]);
    z[2*12] = x[12] * x[12];
    mpn_lshift(buf, buf, 25, 1);
    mpn_add_n(z, z, buf, 25);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_mul1)
/* x has 12.5 words, z has 14. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_12_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[12] = lo;
    z[13] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_shortmul)
/* x and y have 12.5 words, z has 13.
 * Put the low 13 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_12_5_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 13);
    mpfq_fixmp_12_addmul1_nc (z+0, x, y[0]);
    z[13-1] += x[12]*y[0];
    mpfq_fixmp_11_addmul1_nc (z+1, x, y[1]);
    z[13-1] += x[11]*y[1];
    mpfq_fixmp_10_addmul1_nc (z+2, x, y[2]);
    z[13-1] += x[10]*y[2];
    mpfq_fixmp_9_addmul1_nc (z+3, x, y[3]);
    z[13-1] += x[9]*y[3];
    mpfq_fixmp_8_addmul1_nc (z+4, x, y[4]);
    z[13-1] += x[8]*y[4];
    mpfq_fixmp_7_addmul1_nc (z+5, x, y[5]);
    z[13-1] += x[7]*y[5];
    mpfq_fixmp_6_addmul1_nc (z+6, x, y[6]);
    z[13-1] += x[6]*y[6];
    mpfq_fixmp_5_addmul1_nc (z+7, x, y[7]);
    z[13-1] += x[5]*y[7];
    mpfq_fixmp_4_addmul1_nc (z+8, x, y[8]);
    z[13-1] += x[4]*y[8];
    mpfq_fixmp_3_addmul1_nc (z+9, x, y[9]);
    z[13-1] += x[3]*y[9];
    mpfq_fixmp_2_addmul1_nc (z+10, x, y[10]);
    z[13-1] += x[2]*y[10];
    mpfq_fixmp_1_addmul1_nc (z+11, x, y[11]);
    z[13-1] += x[1]*y[11];
    z[13-1] += x[0]*y[13-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_addmul05)
/* x has 12.5 words, z has 13. c is 0.5 word.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul05 */
static inline
mp_limb_t mpfq_fixmp_12_5_addmul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    lo = c*x[12] + carry;
    assert(lo >= carry);
    z[12] += lo;
    return z[12] < lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_addmul05) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_mul05)
/* x has 12.5 words, z has 13. c is 0.5 word.
 * Put (x*c) in z. No carry. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul05 */
static inline
void mpfq_fixmp_12_5_mul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[11] = lo;
    lo = c*x[12] + carry;
    assert(lo >= carry);
    z[12] = lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_mul05) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_mod)
/* x has 25 words. z and p have 12.5 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 12_5_mgy_decode */
static inline
void mpfq_fixmp_12_5_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[12+1], r[13];
    assert (p[13-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 25, p, 13);
    mpfq_copy(z, r, 13);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_mod) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_rshift)
/* a has 12.5 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 12_5_invmod */
static inline
void mpfq_fixmp_12_5_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 13-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[13-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_long_rshift)
/* a has 12.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 12_5_invmod */
static inline
void mpfq_fixmp_12_5_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 13 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[13-off-1] = a[13-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 13 - off);
    }
    mpfq_zero(a + 13 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_long_lshift)
/* a has 12.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 12_5_invmod */
static inline
void mpfq_fixmp_12_5_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 13-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 13 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_invmod)
/* x, z, and p have 12.5 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_12_5_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[13], v[13], a[13], b[13], fix[13];
      int i, t, lsh;
    
      mpfq_zero(u, 13);
      mpfq_zero(v, 13);
      mpfq_copy(a, x, 13);
      mpfq_copy(b, p, 13);
      u[0] = 1UL;
      
      if (mpfq_fixmp_12_5_cmp(a, v) == 0 || mpfq_fixmp_12_5_cmp(a, b) == 0) {
        mpfq_zero(res, 13);
        return 0;
      }
    
      mpfq_fixmp_12_5_add(fix, b, u);
      mpfq_fixmp_12_5_rshift(fix, 1);
    
      assert (mpfq_fixmp_12_5_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 13);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_12_5_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_12_5_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_12_5_sub(b, b, a);
          mpfq_fixmp_12_5_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 13);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_12_5_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_12_5_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_12_5_cmp(a,b) < 0);
        if (mpfq_fixmp_12_5_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_12_5_sub(a, a, b);
          mpfq_fixmp_12_5_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 13);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_12_5_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_12_5_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_12_5_cmp(b,a)<0);
      } while (mpfq_fixmp_12_5_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_12_5_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 13);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_12_5_rshift(u, 1);
        if (sig)
          mpfq_fixmp_12_5_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 13);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_redc)
/* x has 25 words, z and p have 12.5 words.
 * only one word is read from invp.
 * Assuming R=W^13 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_12_5_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 13; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_12_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    {
        mp_limb_t ret[13] = { x[13], x[14], x[15], x[16], x[17], x[18], x[19], x[20], x[21], x[22], x[23], x[24], 0 };
        cy = mpfq_fixmp_12_5_add(x, x, ret);
    }
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_12_5_cmp(x, p) >= 0) {
        mpfq_fixmp_12_5_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 13);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_redc) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_redc_ur)
/* x has 26 words, z and p have 12.5 words.
 * only one word is read from invp.
 * Assuming R=W^13 is the redc modulus, we expect that x verifies:
 *  x < W*W^12.5*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_12_5_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[1];
    for(int i = 0; i < 13; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_12_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_12_5_add(x + 13, x, x + 13);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W^0.5*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W^0.5+1)*p
    */
    if (cy) {
        /* x'/R-p < W^0.5*p, which fits in n words. */
        mpfq_fixmp_12_5_sub(x + 13, x + 13, p);
    }
    mpn_tdiv_qr(q, z, 0, x + 13, 13, p, 13);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_mgy_encode)
/* x, z, and p have 12.5 words.
 * Assuming R=W^13 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_12_5_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[26] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12] };
    mp_limb_t qq[13+1];
    mpn_tdiv_qr(qq, z, 0, t, 26, p, 13);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_mgy_decode)
/* x, z, invR, and p have 12.5 words.
 * Assuming R=W^13 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_12_5_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[26];
    mpfq_fixmp_12_5_mul(t, x, invR);
    mpfq_fixmp_12_5_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_12_5_lshift)
/* a has 12.5 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_12_5_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 13-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_12_5_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_add)
/* x, y, and z have 13.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 13_5_invmod, 13_5_redc, 13_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_13_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r + y[9];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r + y[10];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r + y[11];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[11] = t;
    r = x[12];
    s = r + y[12];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[12] = t;
    r = x[13];
    s = r + y[13];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[13] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_add) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_sub)
/* x, y, and z have 13.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 13_5_invmod, 13_5_redc, 13_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_13_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r - y[9];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r - y[10];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r - y[11];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[11] = t;
    r = x[12];
    s = r - y[12];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[12] = t;
    r = x[13];
    s = r - y[13];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[13] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_sub) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_add_nc)
/* x, y, and z have 13.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_13_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r + y[9];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r + y[10];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r + y[11];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[11] = t;
    r = x[12];
    s = r + y[12];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[12] = t;
    r = x[13];
    s = r + y[13];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[13] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_sub_nc)
/* x, y, and z have 13.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_13_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r - y[9];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r - y[10];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r - y[11];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[11] = t;
    r = x[12];
    s = r - y[12];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[12] = t;
    r = x[13];
    s = r - y[13];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[13] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_add_ui)
/* x, y, and z have 13.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_13_5_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    s = x[9];
    t = s + cy;
    cy = t < s;
    z[9] = t;
    s = x[10];
    t = s + cy;
    cy = t < s;
    z[10] = t;
    s = x[11];
    t = s + cy;
    cy = t < s;
    z[11] = t;
    s = x[12];
    t = s + cy;
    cy = t < s;
    z[12] = t;
    s = x[13];
    t = s + cy;
    cy = t < s;
    z[13] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_sub_ui)
/* x, y, and z have 13.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_13_5_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    s = x[9];
    t = s - cy;
    cy = t > s;
    z[9] = t;
    s = x[10];
    t = s - cy;
    cy = t > s;
    z[10] = t;
    s = x[11];
    t = s - cy;
    cy = t > s;
    z[11] = t;
    s = x[12];
    t = s - cy;
    cy = t > s;
    z[12] = t;
    s = x[13];
    t = s - cy;
    cy = t > s;
    z[13] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_add_ui_nc)
/* x, y, and z have 13.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_13_5_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    s = x[9];
    t = s + cy;
    cy = t < s;
    z[9] = t;
    s = x[10];
    t = s + cy;
    cy = t < s;
    z[10] = t;
    s = x[11];
    t = s + cy;
    cy = t < s;
    z[11] = t;
    s = x[12];
    t = s + cy;
    cy = t < s;
    z[12] = t;
    s = x[13];
    t = s + cy;
    cy = t < s;
    z[13] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_sub_ui_nc)
/* x, y, and z have 13.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_13_5_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    s = x[9];
    t = s - cy;
    cy = t > s;
    z[9] = t;
    s = x[10];
    t = s - cy;
    cy = t > s;
    z[10] = t;
    s = x[11];
    t = s - cy;
    cy = t > s;
    z[11] = t;
    s = x[12];
    t = s - cy;
    cy = t > s;
    z[12] = t;
    s = x[13];
    t = s - cy;
    cy = t > s;
    z[13] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_cmp)
/* x and y have 13.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 13_5_invmod, 13_5_redc, 13_5_redc_ur */
static inline
int mpfq_fixmp_13_5_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 14-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_cmp_ui)
/* x has 13.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 13_5_invmod */
static inline
int mpfq_fixmp_13_5_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 14-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_addmul1)
/* x has 13.5 words, z has 15.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_13_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[12];
    lo += buf;
    carry += (lo<buf);
    z[12] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[13]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[13];
    lo += buf;
    carry += (lo<buf);
    z[13] = lo;
    z[14] += carry;
    return (z[14]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_addmul1_nc)
/* x has 13.5 words, z has 15.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 13_5_mul, 13_5_mgy_decode */
static inline
void mpfq_fixmp_13_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[12];
    lo += buf;
    carry += (lo<buf);
    z[12] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[13]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[13];
    lo += buf;
    carry += (lo<buf);
    z[13] = lo;
    z[14] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_addmul1_shortz)
/* x has 13.5 words, z has 14.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 13_5_redc, 13_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_13_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[12];
    lo += buf;
    carry += (lo<buf);
    z[12] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[13]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[13];
    lo += buf;
    carry += (lo<buf);
    z[13] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_addmul05_nc)
/* x has 13.5 words, z has 14. c is 0.5 word.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 13_5_mul, 13_5_mgy_decode */
static inline
void mpfq_fixmp_13_5_addmul05_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[12];
    lo += buf;
    carry += (lo<buf);
    z[12] = lo;
    lo = c*x[13] + carry;
    assert(lo >= carry);
    z[13] += lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_addmul05_nc) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_mul)
/* x and y have 13.5 words, z has 27. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 13_5_mgy_decode */
static inline
void mpfq_fixmp_13_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 27; z[i++] = 0) ;
    mpfq_fixmp_13_5_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_13_5_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_13_5_addmul1_nc (z + 2, x, y[2]);
    mpfq_fixmp_13_5_addmul1_nc (z + 3, x, y[3]);
    mpfq_fixmp_13_5_addmul1_nc (z + 4, x, y[4]);
    mpfq_fixmp_13_5_addmul1_nc (z + 5, x, y[5]);
    mpfq_fixmp_13_5_addmul1_nc (z + 6, x, y[6]);
    mpfq_fixmp_13_5_addmul1_nc (z + 7, x, y[7]);
    mpfq_fixmp_13_5_addmul1_nc (z + 8, x, y[8]);
    mpfq_fixmp_13_5_addmul1_nc (z + 9, x, y[9]);
    mpfq_fixmp_13_5_addmul1_nc (z + 10, x, y[10]);
    mpfq_fixmp_13_5_addmul1_nc (z + 11, x, y[11]);
    mpfq_fixmp_13_5_addmul1_nc (z + 12, x, y[12]);
    mpfq_fixmp_13_5_addmul05_nc (z + 13, x, y[13]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_mul) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_sqr)
/* x has 13.5 words, z has 27. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_13_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[27] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_fixmp_3_addmul1_nc(buf + 3, x, x[3]);
    mpfq_fixmp_4_addmul1_nc(buf + 4, x, x[4]);
    mpfq_fixmp_5_addmul1_nc(buf + 5, x, x[5]);
    mpfq_fixmp_6_addmul1_nc(buf + 6, x, x[6]);
    mpfq_fixmp_7_addmul1_nc(buf + 7, x, x[7]);
    mpfq_fixmp_8_addmul1_nc(buf + 8, x, x[8]);
    mpfq_fixmp_9_addmul1_nc(buf + 9, x, x[9]);
    mpfq_fixmp_10_addmul1_nc(buf + 10, x, x[10]);
    mpfq_fixmp_11_addmul1_nc(buf + 11, x, x[11]);
    mpfq_fixmp_12_addmul1_nc(buf + 12, x, x[12]);
    mpfq_fixmp_13_addmul1_nc(buf + 13, x, x[13]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpfq_umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
    mpfq_umul_ppmm(z[2*3+1], z[2*3], x[3], x[3]);
    mpfq_umul_ppmm(z[2*4+1], z[2*4], x[4], x[4]);
    mpfq_umul_ppmm(z[2*5+1], z[2*5], x[5], x[5]);
    mpfq_umul_ppmm(z[2*6+1], z[2*6], x[6], x[6]);
    mpfq_umul_ppmm(z[2*7+1], z[2*7], x[7], x[7]);
    mpfq_umul_ppmm(z[2*8+1], z[2*8], x[8], x[8]);
    mpfq_umul_ppmm(z[2*9+1], z[2*9], x[9], x[9]);
    mpfq_umul_ppmm(z[2*10+1], z[2*10], x[10], x[10]);
    mpfq_umul_ppmm(z[2*11+1], z[2*11], x[11], x[11]);
    mpfq_umul_ppmm(z[2*12+1], z[2*12], x[12], x[12]);
    z[2*13] = x[13] * x[13];
    mpn_lshift(buf, buf, 27, 1);
    mpn_add_n(z, z, buf, 27);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_mul1)
/* x has 13.5 words, z has 15. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_13_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[12] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[13]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[13] = lo;
    z[14] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_shortmul)
/* x and y have 13.5 words, z has 14.
 * Put the low 14 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_13_5_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 14);
    mpfq_fixmp_13_addmul1_nc (z+0, x, y[0]);
    z[14-1] += x[13]*y[0];
    mpfq_fixmp_12_addmul1_nc (z+1, x, y[1]);
    z[14-1] += x[12]*y[1];
    mpfq_fixmp_11_addmul1_nc (z+2, x, y[2]);
    z[14-1] += x[11]*y[2];
    mpfq_fixmp_10_addmul1_nc (z+3, x, y[3]);
    z[14-1] += x[10]*y[3];
    mpfq_fixmp_9_addmul1_nc (z+4, x, y[4]);
    z[14-1] += x[9]*y[4];
    mpfq_fixmp_8_addmul1_nc (z+5, x, y[5]);
    z[14-1] += x[8]*y[5];
    mpfq_fixmp_7_addmul1_nc (z+6, x, y[6]);
    z[14-1] += x[7]*y[6];
    mpfq_fixmp_6_addmul1_nc (z+7, x, y[7]);
    z[14-1] += x[6]*y[7];
    mpfq_fixmp_5_addmul1_nc (z+8, x, y[8]);
    z[14-1] += x[5]*y[8];
    mpfq_fixmp_4_addmul1_nc (z+9, x, y[9]);
    z[14-1] += x[4]*y[9];
    mpfq_fixmp_3_addmul1_nc (z+10, x, y[10]);
    z[14-1] += x[3]*y[10];
    mpfq_fixmp_2_addmul1_nc (z+11, x, y[11]);
    z[14-1] += x[2]*y[11];
    mpfq_fixmp_1_addmul1_nc (z+12, x, y[12]);
    z[14-1] += x[1]*y[12];
    z[14-1] += x[0]*y[14-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_addmul05)
/* x has 13.5 words, z has 14. c is 0.5 word.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul05 */
static inline
mp_limb_t mpfq_fixmp_13_5_addmul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[12];
    lo += buf;
    carry += (lo<buf);
    z[12] = lo;
    lo = c*x[13] + carry;
    assert(lo >= carry);
    z[13] += lo;
    return z[13] < lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_addmul05) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_mul05)
/* x has 13.5 words, z has 14. c is 0.5 word.
 * Put (x*c) in z. No carry. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul05 */
static inline
void mpfq_fixmp_13_5_mul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[12] = lo;
    lo = c*x[13] + carry;
    assert(lo >= carry);
    z[13] = lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_mul05) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_mod)
/* x has 27 words. z and p have 13.5 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 13_5_mgy_decode */
static inline
void mpfq_fixmp_13_5_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[13+1], r[14];
    assert (p[14-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 27, p, 14);
    mpfq_copy(z, r, 14);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_mod) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_rshift)
/* a has 13.5 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 13_5_invmod */
static inline
void mpfq_fixmp_13_5_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 14-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[14-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_long_rshift)
/* a has 13.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 13_5_invmod */
static inline
void mpfq_fixmp_13_5_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 14 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[14-off-1] = a[14-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 14 - off);
    }
    mpfq_zero(a + 14 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_long_lshift)
/* a has 13.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 13_5_invmod */
static inline
void mpfq_fixmp_13_5_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 14-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 14 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_invmod)
/* x, z, and p have 13.5 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_13_5_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[14], v[14], a[14], b[14], fix[14];
      int i, t, lsh;
    
      mpfq_zero(u, 14);
      mpfq_zero(v, 14);
      mpfq_copy(a, x, 14);
      mpfq_copy(b, p, 14);
      u[0] = 1UL;
      
      if (mpfq_fixmp_13_5_cmp(a, v) == 0 || mpfq_fixmp_13_5_cmp(a, b) == 0) {
        mpfq_zero(res, 14);
        return 0;
      }
    
      mpfq_fixmp_13_5_add(fix, b, u);
      mpfq_fixmp_13_5_rshift(fix, 1);
    
      assert (mpfq_fixmp_13_5_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 14);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_13_5_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_13_5_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_13_5_sub(b, b, a);
          mpfq_fixmp_13_5_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 14);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_13_5_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_13_5_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_13_5_cmp(a,b) < 0);
        if (mpfq_fixmp_13_5_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_13_5_sub(a, a, b);
          mpfq_fixmp_13_5_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 14);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_13_5_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_13_5_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_13_5_cmp(b,a)<0);
      } while (mpfq_fixmp_13_5_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_13_5_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 14);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_13_5_rshift(u, 1);
        if (sig)
          mpfq_fixmp_13_5_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 14);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_redc)
/* x has 27 words, z and p have 13.5 words.
 * only one word is read from invp.
 * Assuming R=W^14 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_13_5_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 14; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_13_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    {
        mp_limb_t ret[14] = { x[14], x[15], x[16], x[17], x[18], x[19], x[20], x[21], x[22], x[23], x[24], x[25], x[26], 0 };
        cy = mpfq_fixmp_13_5_add(x, x, ret);
    }
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_13_5_cmp(x, p) >= 0) {
        mpfq_fixmp_13_5_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 14);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_redc) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_redc_ur)
/* x has 28 words, z and p have 13.5 words.
 * only one word is read from invp.
 * Assuming R=W^14 is the redc modulus, we expect that x verifies:
 *  x < W*W^13.5*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_13_5_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[1];
    for(int i = 0; i < 14; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_13_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_13_5_add(x + 14, x, x + 14);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W^0.5*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W^0.5+1)*p
    */
    if (cy) {
        /* x'/R-p < W^0.5*p, which fits in n words. */
        mpfq_fixmp_13_5_sub(x + 14, x + 14, p);
    }
    mpn_tdiv_qr(q, z, 0, x + 14, 14, p, 14);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_mgy_encode)
/* x, z, and p have 13.5 words.
 * Assuming R=W^14 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_13_5_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[28] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13] };
    mp_limb_t qq[14+1];
    mpn_tdiv_qr(qq, z, 0, t, 28, p, 14);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_mgy_decode)
/* x, z, invR, and p have 13.5 words.
 * Assuming R=W^14 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_13_5_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[28];
    mpfq_fixmp_13_5_mul(t, x, invR);
    mpfq_fixmp_13_5_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_13_5_lshift)
/* a has 13.5 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_13_5_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 14-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_13_5_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_add)
/* x, y, and z have 14.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add */
/* Triggered by: 14_5_invmod, 14_5_redc, 14_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_14_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r + y[9];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r + y[10];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r + y[11];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[11] = t;
    r = x[12];
    s = r + y[12];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[12] = t;
    r = x[13];
    s = r + y[13];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[13] = t;
    r = x[14];
    s = r + y[14];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[14] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_add) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_sub)
/* x, y, and z have 14.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub */
/* Triggered by: 14_5_invmod, 14_5_redc, 14_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_14_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r - y[9];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r - y[10];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r - y[11];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[11] = t;
    r = x[12];
    s = r - y[12];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[12] = t;
    r = x[13];
    s = r - y[13];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[13] = t;
    r = x[14];
    s = r - y[14];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[14] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_sub) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_add_nc)
/* x, y, and z have 14.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_14_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y[0];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r + y[1];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r + y[2];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r + y[3];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r + y[4];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r + y[5];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r + y[6];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r + y[7];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r + y[8];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r + y[9];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r + y[10];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r + y[11];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[11] = t;
    r = x[12];
    s = r + y[12];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[12] = t;
    r = x[13];
    s = r + y[13];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[13] = t;
    r = x[14];
    s = r + y[14];
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[14] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_add_nc) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_sub_nc)
/* x, y, and z have 14.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_14_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y[0];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    r = x[1];
    s = r - y[1];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[1] = t;
    r = x[2];
    s = r - y[2];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[2] = t;
    r = x[3];
    s = r - y[3];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[3] = t;
    r = x[4];
    s = r - y[4];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[4] = t;
    r = x[5];
    s = r - y[5];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[5] = t;
    r = x[6];
    s = r - y[6];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[6] = t;
    r = x[7];
    s = r - y[7];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[7] = t;
    r = x[8];
    s = r - y[8];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[8] = t;
    r = x[9];
    s = r - y[9];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[9] = t;
    r = x[10];
    s = r - y[10];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[10] = t;
    r = x[11];
    s = r - y[11];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[11] = t;
    r = x[12];
    s = r - y[12];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[12] = t;
    r = x[13];
    s = r - y[13];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[13] = t;
    r = x[14];
    s = r - y[14];
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[14] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_sub_nc) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_add_ui)
/* x, y, and z have 14.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui */
static inline
mp_limb_t mpfq_fixmp_14_5_add_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    s = x[9];
    t = s + cy;
    cy = t < s;
    z[9] = t;
    s = x[10];
    t = s + cy;
    cy = t < s;
    z[10] = t;
    s = x[11];
    t = s + cy;
    cy = t < s;
    z[11] = t;
    s = x[12];
    t = s + cy;
    cy = t < s;
    z[12] = t;
    s = x[13];
    t = s + cy;
    cy = t < s;
    z[13] = t;
    s = x[14];
    t = s + cy;
    cy = t < s;
    z[14] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_add_ui) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_sub_ui)
/* x, y, and z have 14.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui */
static inline
mp_limb_t mpfq_fixmp_14_5_sub_ui(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    s = x[9];
    t = s - cy;
    cy = t > s;
    z[9] = t;
    s = x[10];
    t = s - cy;
    cy = t > s;
    z[10] = t;
    s = x[11];
    t = s - cy;
    cy = t > s;
    z[11] = t;
    s = x[12];
    t = s - cy;
    cy = t > s;
    z[12] = t;
    s = x[13];
    t = s - cy;
    cy = t > s;
    z[13] = t;
    s = x[14];
    t = s - cy;
    cy = t > s;
    z[14] = t;
    return cy;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_sub_ui) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_add_ui_nc)
/* x, y, and z have 14.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_add_ui_nc */
static inline
void mpfq_fixmp_14_5_add_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r + y;
    cy1 = s < r;
    t = s + cy;
    cy2 = t < s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s + cy;
    cy = t < s;
    z[1] = t;
    s = x[2];
    t = s + cy;
    cy = t < s;
    z[2] = t;
    s = x[3];
    t = s + cy;
    cy = t < s;
    z[3] = t;
    s = x[4];
    t = s + cy;
    cy = t < s;
    z[4] = t;
    s = x[5];
    t = s + cy;
    cy = t < s;
    z[5] = t;
    s = x[6];
    t = s + cy;
    cy = t < s;
    z[6] = t;
    s = x[7];
    t = s + cy;
    cy = t < s;
    z[7] = t;
    s = x[8];
    t = s + cy;
    cy = t < s;
    z[8] = t;
    s = x[9];
    t = s + cy;
    cy = t < s;
    z[9] = t;
    s = x[10];
    t = s + cy;
    cy = t < s;
    z[10] = t;
    s = x[11];
    t = s + cy;
    cy = t < s;
    z[11] = t;
    s = x[12];
    t = s + cy;
    cy = t < s;
    z[12] = t;
    s = x[13];
    t = s + cy;
    cy = t < s;
    z[13] = t;
    s = x[14];
    t = s + cy;
    cy = t < s;
    z[14] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_add_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_sub_ui_nc)
/* x, y, and z have 14.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sub_ui_nc */
static inline
void mpfq_fixmp_14_5_sub_ui_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t y)
{
    mp_limb_t r, s, t, cy, cy1, cy2;
    cy = 0;
    r = x[0];
    s = r - y;
    cy1 = s > r;
    t = s - cy;
    cy2 = t > s;
    cy = cy1 | cy2;
    z[0] = t;
    s = x[1];
    t = s - cy;
    cy = t > s;
    z[1] = t;
    s = x[2];
    t = s - cy;
    cy = t > s;
    z[2] = t;
    s = x[3];
    t = s - cy;
    cy = t > s;
    z[3] = t;
    s = x[4];
    t = s - cy;
    cy = t > s;
    z[4] = t;
    s = x[5];
    t = s - cy;
    cy = t > s;
    z[5] = t;
    s = x[6];
    t = s - cy;
    cy = t > s;
    z[6] = t;
    s = x[7];
    t = s - cy;
    cy = t > s;
    z[7] = t;
    s = x[8];
    t = s - cy;
    cy = t > s;
    z[8] = t;
    s = x[9];
    t = s - cy;
    cy = t > s;
    z[9] = t;
    s = x[10];
    t = s - cy;
    cy = t > s;
    z[10] = t;
    s = x[11];
    t = s - cy;
    cy = t > s;
    z[11] = t;
    s = x[12];
    t = s - cy;
    cy = t > s;
    z[12] = t;
    s = x[13];
    t = s - cy;
    cy = t > s;
    z[13] = t;
    s = x[14];
    t = s - cy;
    cy = t > s;
    z[14] = t;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_sub_ui_nc) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_cmp)
/* x and y have 14.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp */
/* Triggered by: 14_5_invmod, 14_5_redc, 14_5_redc_ur */
static inline
int mpfq_fixmp_14_5_cmp(const mp_limb_t * x, const mp_limb_t * y)
{
    for (int i = 15-1; i >= 0; --i) {
        if (x[i] > y[i]) return 1;
        if (x[i] < y[i]) return -1;
    }
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_cmp) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_cmp_ui)
/* x has 14.5 words. Return sign of difference x-y. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_cmp_ui */
/* Triggered by: 14_5_invmod */
static inline
int mpfq_fixmp_14_5_cmp_ui(const mp_limb_t * x, mp_limb_t y)
{
    for (int i = 15-1; i > 0; --i) {
        if (x[i]) return 1;
    }
    if (x[0]>y) return 1;
    if (x[0]<y) return -1;
    return 0;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_cmp_ui) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_addmul1)
/* x has 14.5 words, z has 16.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_14_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[12];
    lo += buf;
    carry += (lo<buf);
    z[12] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[13]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[13];
    lo += buf;
    carry += (lo<buf);
    z[13] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[14]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[14];
    lo += buf;
    carry += (lo<buf);
    z[14] = lo;
    z[15] += carry;
    return (z[15]<carry);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_addmul1) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_addmul1_nc)
/* x has 14.5 words, z has 16.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_nc */
/* Triggered by: 14_5_mul, 14_5_mgy_decode */
static inline
void mpfq_fixmp_14_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[12];
    lo += buf;
    carry += (lo<buf);
    z[12] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[13]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[13];
    lo += buf;
    carry += (lo<buf);
    z[13] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[14]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[14];
    lo += buf;
    carry += (lo<buf);
    z[14] = lo;
    z[15] += carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_addmul1_nc) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_addmul1_shortz)
/* x has 14.5 words, z has 15.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul1_shortz */
/* Triggered by: 14_5_redc, 14_5_redc_ur */
static inline
mp_limb_t mpfq_fixmp_14_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[12];
    lo += buf;
    carry += (lo<buf);
    z[12] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[13]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[13];
    lo += buf;
    carry += (lo<buf);
    z[13] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[14]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[14];
    lo += buf;
    carry += (lo<buf);
    z[14] = lo;
    return carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_addmul1_shortz) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_addmul05_nc)
/* x has 14.5 words, z has 15. c is 0.5 word.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 14_5_mul, 14_5_mgy_decode */
static inline
void mpfq_fixmp_14_5_addmul05_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[12];
    lo += buf;
    carry += (lo<buf);
    z[12] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[13]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[13];
    lo += buf;
    carry += (lo<buf);
    z[13] = lo;
    lo = c*x[14] + carry;
    assert(lo >= carry);
    z[14] += lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_addmul05_nc) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_mul)
/* x and y have 14.5 words, z has 29. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul */
/* Triggered by: 14_5_mgy_decode */
static inline
void mpfq_fixmp_14_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    assert(z != x && z != y);
    for (int i = 0; i < 29; z[i++] = 0) ;
    mpfq_fixmp_14_5_addmul1_nc (z + 0, x, y[0]);
    mpfq_fixmp_14_5_addmul1_nc (z + 1, x, y[1]);
    mpfq_fixmp_14_5_addmul1_nc (z + 2, x, y[2]);
    mpfq_fixmp_14_5_addmul1_nc (z + 3, x, y[3]);
    mpfq_fixmp_14_5_addmul1_nc (z + 4, x, y[4]);
    mpfq_fixmp_14_5_addmul1_nc (z + 5, x, y[5]);
    mpfq_fixmp_14_5_addmul1_nc (z + 6, x, y[6]);
    mpfq_fixmp_14_5_addmul1_nc (z + 7, x, y[7]);
    mpfq_fixmp_14_5_addmul1_nc (z + 8, x, y[8]);
    mpfq_fixmp_14_5_addmul1_nc (z + 9, x, y[9]);
    mpfq_fixmp_14_5_addmul1_nc (z + 10, x, y[10]);
    mpfq_fixmp_14_5_addmul1_nc (z + 11, x, y[11]);
    mpfq_fixmp_14_5_addmul1_nc (z + 12, x, y[12]);
    mpfq_fixmp_14_5_addmul1_nc (z + 13, x, y[13]);
    mpfq_fixmp_14_5_addmul05_nc (z + 14, x, y[14]);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_mul) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_sqr)
/* x has 14.5 words, z has 29. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_14_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mp_limb_t buf[29] = {0,};
    mpfq_fixmp_1_addmul1_nc(buf + 1, x, x[1]);
    mpfq_fixmp_2_addmul1_nc(buf + 2, x, x[2]);
    mpfq_fixmp_3_addmul1_nc(buf + 3, x, x[3]);
    mpfq_fixmp_4_addmul1_nc(buf + 4, x, x[4]);
    mpfq_fixmp_5_addmul1_nc(buf + 5, x, x[5]);
    mpfq_fixmp_6_addmul1_nc(buf + 6, x, x[6]);
    mpfq_fixmp_7_addmul1_nc(buf + 7, x, x[7]);
    mpfq_fixmp_8_addmul1_nc(buf + 8, x, x[8]);
    mpfq_fixmp_9_addmul1_nc(buf + 9, x, x[9]);
    mpfq_fixmp_10_addmul1_nc(buf + 10, x, x[10]);
    mpfq_fixmp_11_addmul1_nc(buf + 11, x, x[11]);
    mpfq_fixmp_12_addmul1_nc(buf + 12, x, x[12]);
    mpfq_fixmp_13_addmul1_nc(buf + 13, x, x[13]);
    mpfq_fixmp_14_addmul1_nc(buf + 14, x, x[14]);
    mpfq_umul_ppmm(z[2*0+1], z[2*0], x[0], x[0]);
    mpfq_umul_ppmm(z[2*1+1], z[2*1], x[1], x[1]);
    mpfq_umul_ppmm(z[2*2+1], z[2*2], x[2], x[2]);
    mpfq_umul_ppmm(z[2*3+1], z[2*3], x[3], x[3]);
    mpfq_umul_ppmm(z[2*4+1], z[2*4], x[4], x[4]);
    mpfq_umul_ppmm(z[2*5+1], z[2*5], x[5], x[5]);
    mpfq_umul_ppmm(z[2*6+1], z[2*6], x[6], x[6]);
    mpfq_umul_ppmm(z[2*7+1], z[2*7], x[7], x[7]);
    mpfq_umul_ppmm(z[2*8+1], z[2*8], x[8], x[8]);
    mpfq_umul_ppmm(z[2*9+1], z[2*9], x[9], x[9]);
    mpfq_umul_ppmm(z[2*10+1], z[2*10], x[10], x[10]);
    mpfq_umul_ppmm(z[2*11+1], z[2*11], x[11], x[11]);
    mpfq_umul_ppmm(z[2*12+1], z[2*12], x[12], x[12]);
    mpfq_umul_ppmm(z[2*13+1], z[2*13], x[13], x[13]);
    z[2*14] = x[14] * x[14];
    mpn_lshift(buf, buf, 29, 1);
    mpn_add_n(z, z, buf, 29);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_sqr) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_mul1)
/* x has 14.5 words, z has 16. Put x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_14_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[12] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[13]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[13] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[14]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[14] = lo;
    z[15] = carry;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_mul1) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_shortmul)
/* x and y have 14.5 words, z has 15.
 * Put the low 15 words of x*y in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_shortmul */
static inline
void mpfq_fixmp_14_5_shortmul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mpfq_zero(z, 15);
    mpfq_fixmp_14_addmul1_nc (z+0, x, y[0]);
    z[15-1] += x[14]*y[0];
    mpfq_fixmp_13_addmul1_nc (z+1, x, y[1]);
    z[15-1] += x[13]*y[1];
    mpfq_fixmp_12_addmul1_nc (z+2, x, y[2]);
    z[15-1] += x[12]*y[2];
    mpfq_fixmp_11_addmul1_nc (z+3, x, y[3]);
    z[15-1] += x[11]*y[3];
    mpfq_fixmp_10_addmul1_nc (z+4, x, y[4]);
    z[15-1] += x[10]*y[4];
    mpfq_fixmp_9_addmul1_nc (z+5, x, y[5]);
    z[15-1] += x[9]*y[5];
    mpfq_fixmp_8_addmul1_nc (z+6, x, y[6]);
    z[15-1] += x[8]*y[6];
    mpfq_fixmp_7_addmul1_nc (z+7, x, y[7]);
    z[15-1] += x[7]*y[7];
    mpfq_fixmp_6_addmul1_nc (z+8, x, y[8]);
    z[15-1] += x[6]*y[8];
    mpfq_fixmp_5_addmul1_nc (z+9, x, y[9]);
    z[15-1] += x[5]*y[9];
    mpfq_fixmp_4_addmul1_nc (z+10, x, y[10]);
    z[15-1] += x[4]*y[10];
    mpfq_fixmp_3_addmul1_nc (z+11, x, y[11]);
    z[15-1] += x[3]*y[11];
    mpfq_fixmp_2_addmul1_nc (z+12, x, y[12]);
    z[15-1] += x[2]*y[12];
    mpfq_fixmp_1_addmul1_nc (z+13, x, y[13]);
    z[15-1] += x[1]*y[13];
    z[15-1] += x[0]*y[15-1];
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_shortmul) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_addmul05)
/* x has 14.5 words, z has 15. c is 0.5 word.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_addmul05 */
static inline
mp_limb_t mpfq_fixmp_14_5_addmul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry, buf;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[0];
    lo += buf;
    carry += (lo<buf);
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[1];
    lo += buf;
    carry += (lo<buf);
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[2];
    lo += buf;
    carry += (lo<buf);
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[3];
    lo += buf;
    carry += (lo<buf);
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[4];
    lo += buf;
    carry += (lo<buf);
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[5];
    lo += buf;
    carry += (lo<buf);
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[6];
    lo += buf;
    carry += (lo<buf);
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[7];
    lo += buf;
    carry += (lo<buf);
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[8];
    lo += buf;
    carry += (lo<buf);
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[9];
    lo += buf;
    carry += (lo<buf);
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[10];
    lo += buf;
    carry += (lo<buf);
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[11];
    lo += buf;
    carry += (lo<buf);
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[12];
    lo += buf;
    carry += (lo<buf);
    z[12] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[13]);
    lo += carry;
    carry = (lo<carry) + hi;
    buf = z[13];
    lo += buf;
    carry += (lo<buf);
    z[13] = lo;
    lo = c*x[14] + carry;
    assert(lo >= carry);
    z[14] += lo;
    return z[14] < lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_addmul05) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_mul05)
/* x has 14.5 words, z has 15. c is 0.5 word.
 * Put (x*c) in z. No carry. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mul05 */
static inline
void mpfq_fixmp_14_5_mul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t hi, lo, carry;
    carry = 0;
    mpfq_umul_ppmm(hi,lo,c,x[0]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[0] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[1]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[1] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[2]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[2] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[3]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[3] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[4]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[4] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[5]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[5] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[6]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[6] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[7]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[7] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[8]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[8] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[9]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[9] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[10]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[10] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[11]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[11] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[12]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[12] = lo;
    mpfq_umul_ppmm(hi,lo,c,x[13]);
    lo += carry;
    carry = (lo<carry) + hi;
    z[13] = lo;
    lo = c*x[14] + carry;
    assert(lo >= carry);
    z[14] = lo;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_mul05) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_mod)
/* x has 29 words. z and p have 14.5 words, and the high word of p is non-zero.
 * Put x mod p in z. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mod */
/* Triggered by: 14_5_mgy_decode */
static inline
void mpfq_fixmp_14_5_mod(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t q[14+1], r[15];
    assert (p[15-1] != 0);
    mpn_tdiv_qr(q, r, 0, x, 29, p, 15);
    mpfq_copy(z, r, 15);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_mod) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_rshift)
/* a has 14.5 words. Shift it in place by cnt bits to the right.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 14_5_invmod */
static inline
void mpfq_fixmp_14_5_rshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 0; i < 15-1; ++i) {
        a[i] >>= cnt;
        a[i] |= (a[i+1] << dnt);
    }
    a[15-1] >>= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_long_rshift)
/* a has 14.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 14_5_invmod */
static inline
void mpfq_fixmp_14_5_long_rshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (int i = 0; i < 15 - off - 1; ++i) {
            a[i] = (a[i+off]>>cnt) | (a[i+off+1]<<dnt);
        }
        a[15-off-1] = a[15-1]>>cnt;
    } else {
        mpfq_copyi(a, a + off, 15 - off);
    }
    mpfq_zero(a + 15 - off, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_long_rshift) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_long_lshift)
/* a has 14.5 words. Shift it in place by off words plus cnt bits to the left.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
/* Triggered by: 14_5_invmod */
static inline
void mpfq_fixmp_14_5_long_lshift(mp_limb_t * a, int off MAYBE_UNUSED, int cnt)
{
    int i;
    if (cnt) {
        int dnt = GMP_NUMB_BITS - cnt;
        for (i = 15-1; i>off; --i) {
            a[i] = (a[i-off]<<cnt) | (a[i-off-1]>>dnt);
        }
        a[off] = a[0]<<cnt;
    } else {
        mpfq_copyd(a + off, a, 15 - off);
    }
    mpfq_zero(a, off);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_long_lshift) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_invmod)
/* x, z, and p have 14.5 words. Put inverse of x mod p in z.
 * Return non-zero if an inverse could be found.
 * If no inverse could be found, return 0 and set z to zero.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_invmod */
static inline
int mpfq_fixmp_14_5_invmod(mp_limb_t * res, const mp_limb_t * x, const mp_limb_t * p)
{
      mp_limb_t u[15], v[15], a[15], b[15], fix[15];
      int i, t, lsh;
    
      mpfq_zero(u, 15);
      mpfq_zero(v, 15);
      mpfq_copy(a, x, 15);
      mpfq_copy(b, p, 15);
      u[0] = 1UL;
      
      if (mpfq_fixmp_14_5_cmp(a, v) == 0 || mpfq_fixmp_14_5_cmp(a, b) == 0) {
        mpfq_zero(res, 15);
        return 0;
      }
    
      mpfq_fixmp_14_5_add(fix, b, u);
      mpfq_fixmp_14_5_rshift(fix, 1);
    
      assert (mpfq_fixmp_14_5_cmp(a,b) < 0);
    
      t = 0;
      
      for(i = 0 ; !a[i] ; i++) ;
      assert (i < 15);
      lsh = mpfq_ctzl(a[i]);
      mpfq_fixmp_14_5_long_rshift(a, i, lsh);
      t += lsh + i*GMP_NUMB_BITS;
      mpfq_fixmp_14_5_long_lshift(v, i, lsh);
    
      do {
        do {
          mpfq_fixmp_14_5_sub(b, b, a);
          mpfq_fixmp_14_5_add(v, v, u);
          for(i = 0 ; !b[i] ; i++) ;
          assert (i < 15);
          lsh = mpfq_ctzl(b[i]);
          mpfq_fixmp_14_5_long_rshift(b, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_14_5_long_lshift(u, i, lsh);
        } while (mpfq_fixmp_14_5_cmp(a,b) < 0);
        if (mpfq_fixmp_14_5_cmp(a, b) == 0)
          break;
        do {
          mpfq_fixmp_14_5_sub(a, a, b);
          mpfq_fixmp_14_5_add(u, u, v);
          for(i = 0 ; !a[i] ; i++) ;
          assert (i < 15);
          lsh = mpfq_ctzl(a[i]);
          mpfq_fixmp_14_5_long_rshift(a, i, lsh);
          t += lsh + i*GMP_NUMB_BITS;
          mpfq_fixmp_14_5_long_lshift(v, i, lsh);
        } while (mpfq_fixmp_14_5_cmp(b,a)<0);
      } while (mpfq_fixmp_14_5_cmp(a,b) != 0);
      {
        if (mpfq_fixmp_14_5_cmp_ui(a, 1) != 0) {
          mpfq_copy(res, a, 15);
          return 0;
        }
      }
      while (t>0) {
        mp_limb_t sig = u[0] & 1UL;
        mpfq_fixmp_14_5_rshift(u, 1);
        if (sig)
          mpfq_fixmp_14_5_add(u, u, fix);
        --t;
      }
      mpfq_copy(res, u, 15);
      return 1;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_invmod) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_redc)
/* x has 29 words, z and p have 14.5 words.
 * only one word is read from invp.
 * Assuming R=W^15 is the redc modulus, we expect that x verifies:
 *   x < R*p,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc */
static inline
void mpfq_fixmp_14_5_redc(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy;
    for(int i = 0; i < 15; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_14_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    {
        mp_limb_t ret[15] = { x[15], x[16], x[17], x[18], x[19], x[20], x[21], x[22], x[23], x[24], x[25], x[26], x[27], x[28], 0 };
        cy = mpfq_fixmp_14_5_add(x, x, ret);
    }
    /* At this point, we have (x' denotes x + cy*W^n here)
    * x' <= R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < p + p
    */
    if (cy || mpfq_fixmp_14_5_cmp(x, p) >= 0) {
        mpfq_fixmp_14_5_sub(z, x, p);
    } else {
        mpfq_copy(z, x, 15);
    }
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_redc) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_redc_ur)
/* x has 30 words, z and p have 14.5 words.
 * only one word is read from invp.
 * Assuming R=W^15 is the redc modulus, we expect that x verifies:
 *  x < W*W^14.5*p = W^0.5*R*p or the hw case, W*R*p otherwise,
 * so that we have eventually z < p, z congruent to x/R mod p.
 * The contents of the area pointed by x are clobbered by this call.
 * Note also that x may alias z.
 */
/* *Mpfq::fixmp::longlong::code_for__fixmp_redc_ur */
static inline
void mpfq_fixmp_14_5_redc_ur(mp_limb_t * z, mp_limb_t * x, const mp_limb_t * mip, const mp_limb_t * p)
{
    mp_limb_t cy, q[1];
    for(int i = 0; i < 15; ++i) {
        mp_limb_t t = x[i]*mip[0];
        cy = mpfq_fixmp_14_5_addmul1_shortz(x+i, p, t);
        assert (x[i] == 0);
        x[i] = cy;
    }
    cy = mpfq_fixmp_14_5_add(x + 15, x, x + 15);
    /* At this point, we have (x' denotes x + cy*W^(n+1) here)
    * x' <= W^0.5*R*p-1 + (W-1)*p*(1+W+...+W^{n-1}) and x = mod R.
    * x'/R < (W^0.5+1)*p
    */
    if (cy) {
        /* x'/R-p < W^0.5*p, which fits in n words. */
        mpfq_fixmp_14_5_sub(x + 15, x + 15, p);
    }
    mpn_tdiv_qr(q, z, 0, x + 15, 15, p, 15);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_redc_ur) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_mgy_encode)
/* x, z, and p have 14.5 words.
 * Assuming R=W^15 is the redc modulus, we compute z=R*x mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_encode */
static inline
void mpfq_fixmp_14_5_mgy_encode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * p)
{
    mp_limb_t t[30] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14] };
    mp_limb_t qq[15+1];
    mpn_tdiv_qr(qq, z, 0, t, 30, p, 15);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_mgy_encode) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_mgy_decode)
/* x, z, invR, and p have 14.5 words.
 * Assuming R=W^15 is the redc modulus, we compute z=x/R mod p. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_mgy_decode */
static inline
void mpfq_fixmp_14_5_mgy_decode(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * invR, const mp_limb_t * p)
{
    mp_limb_t t[30];
    mpfq_fixmp_14_5_mul(t, x, invR);
    mpfq_fixmp_14_5_mod(z, t, p);
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_mgy_decode) */

#if !defined(HAVE_native_mpfq_fixmp_14_5_lshift)
/* a has 14.5 words. Shift it in place by cnt bits to the left.
 * The shift count cnt must not exceed the word size.
 * Note that no carry is returned for the bits shifted out. */
/* *Mpfq::fixmp::longlong::code_for__fixmp_lshift */
static inline
void mpfq_fixmp_14_5_lshift(mp_limb_t * a, int cnt)
{
    if (!cnt) return;
    int i;
    int dnt = GMP_NUMB_BITS - cnt;
    for (i = 15-1; i>0; --i) {
        a[i] <<= cnt;
        a[i] |= (a[i-1] >> dnt);
    }
    a[0] <<= cnt;
}
#endif /* !defined(HAVE_native_mpfq_fixmp_14_5_lshift) */


#endif  /* MPFQ_FIXMP_LONGLONG_H_ */

/* vim:set ft=cpp: */
