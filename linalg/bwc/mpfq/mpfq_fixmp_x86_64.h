#ifndef MPFQ_FIXMP_X86_64_H_
#define MPFQ_FIXMP_X86_64_H_

/* MPFQ generated file -- do not edit */

#include <gmp.h>
#include <limits.h>
#ifdef	MPFQ_LAST_GENERATED_TAG
#undef	MPFQ_LAST_GENERATED_TAG
#endif
#define MPFQ_LAST_GENERATED_TAG      fixmp

/* Active handler: Mpfq::fixmp::x86_64 */
/* Options used:{ features={ gcc_inline_assembly=1, }, tag=fixmp, w=64, } */


#ifdef  __cplusplus
extern "C" {
#endif
#define HAVE_native_mpfq_fixmp_1_add
static inline
mp_limb_t mpfq_fixmp_1_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_1_sub
static inline
mp_limb_t mpfq_fixmp_1_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_1_add_nc
static inline
void mpfq_fixmp_1_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_1_sub_nc
static inline
void mpfq_fixmp_1_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_1_addmul1
static inline
mp_limb_t mpfq_fixmp_1_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_1_addmul1_nc
static inline
void mpfq_fixmp_1_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_1_addmul1_shortz
static inline
mp_limb_t mpfq_fixmp_1_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_1_mul
static inline
void mpfq_fixmp_1_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_1_sqr
static inline
void mpfq_fixmp_1_sqr(mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_1_mul1
static inline
void mpfq_fixmp_1_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_1_mulredc
static inline
void mpfq_fixmp_1_mulredc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_2_add
static inline
mp_limb_t mpfq_fixmp_2_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_2_sub
static inline
mp_limb_t mpfq_fixmp_2_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_2_add_nc
static inline
void mpfq_fixmp_2_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_2_sub_nc
static inline
void mpfq_fixmp_2_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_2_addmul1
static inline
mp_limb_t mpfq_fixmp_2_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_2_addmul1_nc
static inline
void mpfq_fixmp_2_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_2_addmul1_shortz
static inline
mp_limb_t mpfq_fixmp_2_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_2_mul
static inline
void mpfq_fixmp_2_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_2_sqr
static inline
void mpfq_fixmp_2_sqr(mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_2_mul1
static inline
void mpfq_fixmp_2_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_3_add
static inline
mp_limb_t mpfq_fixmp_3_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_3_sub
static inline
mp_limb_t mpfq_fixmp_3_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_3_add_nc
static inline
void mpfq_fixmp_3_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_3_sub_nc
static inline
void mpfq_fixmp_3_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_3_addmul1
static inline
mp_limb_t mpfq_fixmp_3_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_3_addmul1_nc
static inline
void mpfq_fixmp_3_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_3_addmul1_shortz
static inline
mp_limb_t mpfq_fixmp_3_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_3_mul
static inline
void mpfq_fixmp_3_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_3_sqr
static inline
void mpfq_fixmp_3_sqr(mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_3_mul1
static inline
void mpfq_fixmp_3_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_4_add
static inline
mp_limb_t mpfq_fixmp_4_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_4_sub
static inline
mp_limb_t mpfq_fixmp_4_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_4_add_nc
static inline
void mpfq_fixmp_4_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_4_sub_nc
static inline
void mpfq_fixmp_4_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_4_addmul1
static inline
mp_limb_t mpfq_fixmp_4_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_4_addmul1_nc
static inline
void mpfq_fixmp_4_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_4_addmul1_shortz
static inline
mp_limb_t mpfq_fixmp_4_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_4_mul
static inline
void mpfq_fixmp_4_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_4_sqr
static inline
void mpfq_fixmp_4_sqr(mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_4_mul1
static inline
void mpfq_fixmp_4_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_5_add
static inline
mp_limb_t mpfq_fixmp_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_5_sub
static inline
mp_limb_t mpfq_fixmp_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_5_add_nc
static inline
void mpfq_fixmp_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_5_sub_nc
static inline
void mpfq_fixmp_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_5_addmul1
static inline
mp_limb_t mpfq_fixmp_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_5_addmul1_nc
static inline
void mpfq_fixmp_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_5_addmul1_shortz
static inline
mp_limb_t mpfq_fixmp_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_5_mul
static inline
void mpfq_fixmp_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_5_sqr
static inline
void mpfq_fixmp_5_sqr(mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_5_mul1
static inline
void mpfq_fixmp_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_6_add
static inline
mp_limb_t mpfq_fixmp_6_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_6_sub
static inline
mp_limb_t mpfq_fixmp_6_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_6_add_nc
static inline
void mpfq_fixmp_6_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_6_sub_nc
static inline
void mpfq_fixmp_6_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_6_addmul1
static inline
mp_limb_t mpfq_fixmp_6_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_6_addmul1_nc
static inline
void mpfq_fixmp_6_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_6_addmul1_shortz
static inline
mp_limb_t mpfq_fixmp_6_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_6_mul
static inline
void mpfq_fixmp_6_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_6_sqr
static inline
void mpfq_fixmp_6_sqr(mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_6_mul1
static inline
void mpfq_fixmp_6_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_7_add
static inline
mp_limb_t mpfq_fixmp_7_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_7_sub
static inline
mp_limb_t mpfq_fixmp_7_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_7_add_nc
static inline
void mpfq_fixmp_7_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_7_sub_nc
static inline
void mpfq_fixmp_7_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_7_addmul1
static inline
mp_limb_t mpfq_fixmp_7_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_7_addmul1_nc
static inline
void mpfq_fixmp_7_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_7_addmul1_shortz
static inline
mp_limb_t mpfq_fixmp_7_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_7_mul
static inline
void mpfq_fixmp_7_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_7_sqr
static inline
void mpfq_fixmp_7_sqr(mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_7_mul1
static inline
void mpfq_fixmp_7_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_8_add
static inline
mp_limb_t mpfq_fixmp_8_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_8_sub
static inline
mp_limb_t mpfq_fixmp_8_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_8_add_nc
static inline
void mpfq_fixmp_8_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_8_sub_nc
static inline
void mpfq_fixmp_8_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_8_addmul1
static inline
mp_limb_t mpfq_fixmp_8_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_8_addmul1_nc
static inline
void mpfq_fixmp_8_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_8_addmul1_shortz
static inline
mp_limb_t mpfq_fixmp_8_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_8_mul
static inline
void mpfq_fixmp_8_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_8_sqr
static inline
void mpfq_fixmp_8_sqr(mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_8_mul1
static inline
void mpfq_fixmp_8_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_9_add
static inline
mp_limb_t mpfq_fixmp_9_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_9_sub
static inline
mp_limb_t mpfq_fixmp_9_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_9_add_nc
static inline
void mpfq_fixmp_9_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_9_sub_nc
static inline
void mpfq_fixmp_9_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_9_addmul1
static inline
mp_limb_t mpfq_fixmp_9_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_9_addmul1_nc
static inline
void mpfq_fixmp_9_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_9_addmul1_shortz
static inline
mp_limb_t mpfq_fixmp_9_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_9_mul
static inline
void mpfq_fixmp_9_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_9_sqr
static inline
void mpfq_fixmp_9_sqr(mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_9_mul1
static inline
void mpfq_fixmp_9_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_0_5_add
static inline
mp_limb_t mpfq_fixmp_0_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_0_5_sub
static inline
mp_limb_t mpfq_fixmp_0_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_0_5_add_nc
static inline
void mpfq_fixmp_0_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_0_5_sub_nc
static inline
void mpfq_fixmp_0_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_0_5_addmul1
static inline
mp_limb_t mpfq_fixmp_0_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_0_5_addmul1_nc
static inline
void mpfq_fixmp_0_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_0_5_addmul1_shortz
static inline
mp_limb_t mpfq_fixmp_0_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_0_5_mul
static inline
void mpfq_fixmp_0_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_0_5_sqr
static inline
void mpfq_fixmp_0_5_sqr(mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_0_5_mul1
static inline
void mpfq_fixmp_0_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_0_5_addmul05
static inline
mp_limb_t mpfq_fixmp_0_5_addmul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_0_5_addmul05_nc
static inline
void mpfq_fixmp_0_5_addmul05_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_0_5_mul05
static inline
void mpfq_fixmp_0_5_mul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_1_5_add
static inline
mp_limb_t mpfq_fixmp_1_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_1_5_sub
static inline
mp_limb_t mpfq_fixmp_1_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_1_5_add_nc
static inline
void mpfq_fixmp_1_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_1_5_sub_nc
static inline
void mpfq_fixmp_1_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_1_5_addmul1
static inline
mp_limb_t mpfq_fixmp_1_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_1_5_addmul1_nc
static inline
void mpfq_fixmp_1_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_1_5_addmul1_shortz
static inline
mp_limb_t mpfq_fixmp_1_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_1_5_mul
static inline
void mpfq_fixmp_1_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_1_5_sqr
static inline
void mpfq_fixmp_1_5_sqr(mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_1_5_mul1
static inline
void mpfq_fixmp_1_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_1_5_addmul05
static inline
mp_limb_t mpfq_fixmp_1_5_addmul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_1_5_addmul05_nc
static inline
void mpfq_fixmp_1_5_addmul05_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_1_5_mul05
static inline
void mpfq_fixmp_1_5_mul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_2_5_add
static inline
mp_limb_t mpfq_fixmp_2_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_2_5_sub
static inline
mp_limb_t mpfq_fixmp_2_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_2_5_add_nc
static inline
void mpfq_fixmp_2_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_2_5_sub_nc
static inline
void mpfq_fixmp_2_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_2_5_addmul1
static inline
mp_limb_t mpfq_fixmp_2_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_2_5_addmul1_nc
static inline
void mpfq_fixmp_2_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_2_5_addmul1_shortz
static inline
mp_limb_t mpfq_fixmp_2_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_2_5_mul
static inline
void mpfq_fixmp_2_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_2_5_sqr
static inline
void mpfq_fixmp_2_5_sqr(mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_2_5_mul1
static inline
void mpfq_fixmp_2_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_2_5_addmul05
static inline
mp_limb_t mpfq_fixmp_2_5_addmul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_2_5_addmul05_nc
static inline
void mpfq_fixmp_2_5_addmul05_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_2_5_mul05
static inline
void mpfq_fixmp_2_5_mul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_3_5_add
static inline
mp_limb_t mpfq_fixmp_3_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_3_5_sub
static inline
mp_limb_t mpfq_fixmp_3_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_3_5_add_nc
static inline
void mpfq_fixmp_3_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_3_5_sub_nc
static inline
void mpfq_fixmp_3_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_3_5_addmul1
static inline
mp_limb_t mpfq_fixmp_3_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_3_5_addmul1_nc
static inline
void mpfq_fixmp_3_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_3_5_addmul1_shortz
static inline
mp_limb_t mpfq_fixmp_3_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_3_5_mul
static inline
void mpfq_fixmp_3_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_3_5_sqr
static inline
void mpfq_fixmp_3_5_sqr(mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_3_5_mul1
static inline
void mpfq_fixmp_3_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_3_5_addmul05
static inline
mp_limb_t mpfq_fixmp_3_5_addmul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_3_5_addmul05_nc
static inline
void mpfq_fixmp_3_5_addmul05_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_3_5_mul05
static inline
void mpfq_fixmp_3_5_mul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_4_5_add
static inline
mp_limb_t mpfq_fixmp_4_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_4_5_sub
static inline
mp_limb_t mpfq_fixmp_4_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_4_5_add_nc
static inline
void mpfq_fixmp_4_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_4_5_sub_nc
static inline
void mpfq_fixmp_4_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_4_5_addmul1
static inline
mp_limb_t mpfq_fixmp_4_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_4_5_addmul1_nc
static inline
void mpfq_fixmp_4_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_4_5_addmul1_shortz
static inline
mp_limb_t mpfq_fixmp_4_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_4_5_mul
static inline
void mpfq_fixmp_4_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_4_5_sqr
static inline
void mpfq_fixmp_4_5_sqr(mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_4_5_mul1
static inline
void mpfq_fixmp_4_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_4_5_addmul05
static inline
mp_limb_t mpfq_fixmp_4_5_addmul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_4_5_addmul05_nc
static inline
void mpfq_fixmp_4_5_addmul05_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_4_5_mul05
static inline
void mpfq_fixmp_4_5_mul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_5_5_add
static inline
mp_limb_t mpfq_fixmp_5_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_5_5_sub
static inline
mp_limb_t mpfq_fixmp_5_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_5_5_add_nc
static inline
void mpfq_fixmp_5_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_5_5_sub_nc
static inline
void mpfq_fixmp_5_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_5_5_addmul1
static inline
mp_limb_t mpfq_fixmp_5_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_5_5_addmul1_nc
static inline
void mpfq_fixmp_5_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_5_5_addmul1_shortz
static inline
mp_limb_t mpfq_fixmp_5_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_5_5_mul
static inline
void mpfq_fixmp_5_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_5_5_sqr
static inline
void mpfq_fixmp_5_5_sqr(mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_5_5_mul1
static inline
void mpfq_fixmp_5_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_5_5_addmul05
static inline
mp_limb_t mpfq_fixmp_5_5_addmul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_5_5_addmul05_nc
static inline
void mpfq_fixmp_5_5_addmul05_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_5_5_mul05
static inline
void mpfq_fixmp_5_5_mul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_6_5_add
static inline
mp_limb_t mpfq_fixmp_6_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_6_5_sub
static inline
mp_limb_t mpfq_fixmp_6_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_6_5_add_nc
static inline
void mpfq_fixmp_6_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_6_5_sub_nc
static inline
void mpfq_fixmp_6_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_6_5_addmul1
static inline
mp_limb_t mpfq_fixmp_6_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_6_5_addmul1_nc
static inline
void mpfq_fixmp_6_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_6_5_addmul1_shortz
static inline
mp_limb_t mpfq_fixmp_6_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_6_5_mul
static inline
void mpfq_fixmp_6_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_6_5_sqr
static inline
void mpfq_fixmp_6_5_sqr(mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_6_5_mul1
static inline
void mpfq_fixmp_6_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_6_5_addmul05
static inline
mp_limb_t mpfq_fixmp_6_5_addmul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_6_5_addmul05_nc
static inline
void mpfq_fixmp_6_5_addmul05_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_6_5_mul05
static inline
void mpfq_fixmp_6_5_mul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_7_5_add
static inline
mp_limb_t mpfq_fixmp_7_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_7_5_sub
static inline
mp_limb_t mpfq_fixmp_7_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_7_5_add_nc
static inline
void mpfq_fixmp_7_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_7_5_sub_nc
static inline
void mpfq_fixmp_7_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_7_5_addmul1
static inline
mp_limb_t mpfq_fixmp_7_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_7_5_addmul1_nc
static inline
void mpfq_fixmp_7_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_7_5_addmul1_shortz
static inline
mp_limb_t mpfq_fixmp_7_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_7_5_mul
static inline
void mpfq_fixmp_7_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_7_5_sqr
static inline
void mpfq_fixmp_7_5_sqr(mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_7_5_mul1
static inline
void mpfq_fixmp_7_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_7_5_addmul05
static inline
mp_limb_t mpfq_fixmp_7_5_addmul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_7_5_addmul05_nc
static inline
void mpfq_fixmp_7_5_addmul05_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_7_5_mul05
static inline
void mpfq_fixmp_7_5_mul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_8_5_add
static inline
mp_limb_t mpfq_fixmp_8_5_add(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_8_5_sub
static inline
mp_limb_t mpfq_fixmp_8_5_sub(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_8_5_add_nc
static inline
void mpfq_fixmp_8_5_add_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_8_5_sub_nc
static inline
void mpfq_fixmp_8_5_sub_nc(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_8_5_addmul1
static inline
mp_limb_t mpfq_fixmp_8_5_addmul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_8_5_addmul1_nc
static inline
void mpfq_fixmp_8_5_addmul1_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_8_5_addmul1_shortz
static inline
mp_limb_t mpfq_fixmp_8_5_addmul1_shortz(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_8_5_mul
static inline
void mpfq_fixmp_8_5_mul(mp_limb_t *, const mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_8_5_sqr
static inline
void mpfq_fixmp_8_5_sqr(mp_limb_t *, const mp_limb_t *);
#define HAVE_native_mpfq_fixmp_8_5_mul1
static inline
void mpfq_fixmp_8_5_mul1(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_8_5_addmul05
static inline
mp_limb_t mpfq_fixmp_8_5_addmul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_8_5_addmul05_nc
static inline
void mpfq_fixmp_8_5_addmul05_nc(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#define HAVE_native_mpfq_fixmp_8_5_mul05
static inline
void mpfq_fixmp_8_5_mul05(mp_limb_t *, const mp_limb_t *, mp_limb_t);
#ifdef  __cplusplus
}
#endif

/* Implementations for inlines */
/* x, y, and z have 1 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add */
static inline
mp_limb_t mpfq_fixmp_1_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "addq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : [z0]"=m"(z[0]), "=&a"(carry)
    : [x0]"m"(x[0]), [y0]"m"(y[0])
    : );
    return carry;
}

/* x, y, and z have 1 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub */
static inline
mp_limb_t mpfq_fixmp_1_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "subq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : [z0]"=m"(z[0]), "=&a"(carry)
    : [x0]"m"(x[0]), [y0]"m"(y[0])
    : );
    return carry;
}

/* x, y, and z have 1 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_1_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    *z = *x + *y;
}

/* x, y, and z have 1 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_1_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    *z = *x - *y;
}

/* x has 1 words, z has 3.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_1_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "adcq    $0, %%rdx\n"
        "xorq    %%rax, %%rax\n"
        "addq    %%rdx, %[z1]\n"
        "adcq    $0, %%rax\n"
    : [z0]"+m"(z[0]), [z1]"+m"(z[1]), "=a"(carry)
    : [x0]"m"(x[0]), [mult]"m"(c)
    : "%rdx");
    return carry;
}

/* x has 1 words, z has 3.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_nc */
static inline
void mpfq_fixmp_1_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rdx, %[z1]\n"
    : [z0]"+m"(z[0]), [z1]"+m"(z[1])
    : [x0]"m"(x[0]), [mult]"m"(c)
    : "%rax", "%rdx");
}

/* x has 1 words, z has 2.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_shortz */
static inline
mp_limb_t mpfq_fixmp_1_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "adcq    $0, %%rdx\n"
    : [z0]"+m"(z[0]), [z1]"+m"(z[1]), "=d"(carry)
    : [x0]"m"(x[0]), [mult]"m"(c)
    : "%rax");
    return carry;
}

/* x and y have 1 words, z has 4. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul */
static inline
void mpfq_fixmp_1_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "mulq %[y0]\n"
              : "=a" (z[0]), "=d" (z[1])
              : "0" (x[0]), [y0] "rm1" (y[0])
              : "cc");
}

/* x has 1 words, z has 4. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_1_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    __asm__ __volatile__(
        "mulq %%rax\n"
    : "=a" (z[0]), "=d" (z[1])
    : "0" (x[0])
    : "cc");
}

/* x has 1 words, z has 3. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_1_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "mulq    %[mult]\n"
    : [z0]"=a"(z[0]), [z1]"=d"(z[1])
    : [x0]"0"(x[0]), [mult]"m"(c)
    : );
}

/* x, y, z and p have 1 words.
 * only one word is read from invp.
 * We expect that x and y are both < p.
 * Assuming R=W^2 is the redc modulus, we return x*y/R mod p in z */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mulredc */
static inline
void mpfq_fixmp_1_mulredc(mp_limb_t * pz, const mp_limb_t * x, const mp_limb_t * y, const mp_limb_t * p, const mp_limb_t * invp)
{
    mp_limb_t z;
    __asm__ __volatile__(
    "movq    %[x0], %%rax\n"
    "mulq    %[y0]\n"
    "movq    %%rdx, %[z]\n"
    "imul    %[invp0], %%rax\n"
    "mulq    %[p0]\n"
    "addq    $0xFFFFFFFFFFFFFFFF, %%rax\n" // set carry if ax is not 0
    "adcq    $0, %[z]\n"                   // this should not produce any carry
    "movq    %[z], %%rax\n"                 // ax = z
    "subq    %[p0], %%rax\n"                // ax -= p
    "addq    %%rdx, %[z]\n"                 // z += dx
    "addq    %%rdx, %%rax\n"                // ax += dx  (ax = z-p+dx)
    "cmovc   %%rax, %[z]\n"                 // z c_= ax
    : [z]"=&r"(z)
    : [x0]"rm"(x[0]), [y0]"rm"(y[0]), [p0]"rm"(p[0]), [invp0]"rm"(invp[0])
    : "%rax", "%rdx");
    *pz = z;
}

/* x, y, and z have 2 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add */
static inline
mp_limb_t mpfq_fixmp_2_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "addq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    %[y1], %%rax\n"
        "movq    %%rax, %[z1]\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), "=&a"(carry)
    : [x0]"m"(x[0]), [x1]"m"(x[1]), [y0]"m"(y[0]), [y1]"m"(y[1])
    : );
    return carry;
}

/* x, y, and z have 2 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub */
static inline
mp_limb_t mpfq_fixmp_2_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "subq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "sbbq    %[y1], %%rax\n"
        "movq    %%rax, %[z1]\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), "=&a"(carry)
    : [x0]"m"(x[0]), [x1]"m"(x[1]), [y0]"m"(y[0]), [y1]"m"(y[1])
    : );
    return carry;
}

/* x, y, and z have 2 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_2_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "addq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    %[y1], %%rax\n"
        "movq    %%rax, %[z1]\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1])
    : [x0]"m"(x[0]), [x1]"m"(x[1]), [y0]"m"(y[0]), [y1]"m"(y[1])
    : "%rax");
}

/* x, y, and z have 2 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_2_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "subq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "sbbq    %[y1], %%rax\n"
        "movq    %%rax, %[z1]\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1])
    : [x0]"m"(x[0]), [x1]"m"(x[1]), [y0]"m"(y[0]), [y1]"m"(y[1])
    : "%rax");
}

/* x has 2 words, z has 4.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_2_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, %[z1]\n"
        "adcq    $0, %%rdx\n"
        "xorq    %%rcx, %%rcx\n"
        "addq    %%rdx, %[z2]\n"
        "adcq    $0, %%rcx\n"
        "movq    %%rcx, %[carry]\n"
    : [carry]"=g"(carry), [z0]"+m"(z[0]), [z1]"+m"(z[1]), [z2]"+m"(z[2])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1])
    : "%rax", "%rcx", "%rdx");
    return carry;
}

/* x has 2 words, z has 4.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_nc */
static inline
void mpfq_fixmp_2_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, %[z1]\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rdx, %[z2]\n"
    : [z0]"+m"(z[0]), [z1]"+m"(z[1]), [z2]"+m"(z[2])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1])
    : "%rax", "%rcx", "%rdx");
}

/* x has 2 words, z has 3.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_shortz */
static inline
mp_limb_t mpfq_fixmp_2_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, %[z1]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %[carry]\n"
    : [carry]"=g"(carry), [z0]"+m"(z[0]), [z1]"+m"(z[1]), [z2]"+m"(z[2])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1])
    : "%rax", "%rcx", "%rdx");
    return carry;
}

/* x and y have 2 words, z has 6. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul */
static inline
void mpfq_fixmp_2_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "### x*y[0]\n"
        "movq    %2, %%r8\n"
        "movq    %0, %%rdi\n"
        "movq    0(%%r8), %%r9\n"
        "movq    %1, %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, 16(%%rdi)\n"
        "movq    $0, 24(%%rdi)\n"
        "### x*y[1]\n"
        "movq    8(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 8(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 24(%%rdi)\n"
      : "+m" (z)
      : "m" (x), "m" (y)
      : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", "%r9", "memory");
}

/* x has 2 words, z has 6. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_2_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    __asm__ __volatile__(
        "movq    %1, %%r8\n"
        "movq    %0, %%rdi\n"
        "movq    (%%r8), %%rax\n"
        "mulq    %%rax\n"
        "movq    %%rax, (%%rdi)\n"
        "movq    8(%%r8), %%rax\n"
        "movq    %%rdx, %%r9\n"
        "mulq    %%rax\n"
        "movq    %%rax, %%r10\n"
        "movq    (%%r8), %%rax\n"
        "movq    %%rdx, %%r11\n"
        "mulq    8(%%r8)\n"
        "addq    %%rax, %%r9\n"
        "adcq    %%rdx, %%r10\n"
        "adcq    $0, %%r11\n"
        "addq    %%rax, %%r9\n"
        "movq    %%r9, 8(%%rdi)\n"
        "adcq    %%rdx, %%r10\n"
        "movq    %%r10, 16(%%rdi)\n"
        "adcq    $0, %%r11\n"
        "movq    %%r11, 24(%%rdi)\n"
    : "+m" (z)
    : "m" (x)
    : "%rax", "%rdx", "%rdi", "%r8", "%r9", "%r10", "%r11", "memory");
}

/* x has 2 words, z has 4. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_2_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, %[z1]\n"
        "movq    %%rdx, %[z2]\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), [z2]"=m"(z[2])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1])
    : "%rax", "%rcx", "%rdx");
}

/* x, y, and z have 3 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add */
static inline
mp_limb_t mpfq_fixmp_3_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "addq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    %[y1], %%rax\n"
        "movq    %%rax, %[z1]\n"
        "movq    %[x2], %%rax\n"
        "adcq    %[y2], %%rax\n"
        "movq    %%rax, %[z2]\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), [z2]"=m"(z[2]), "=&a"(carry)
    : [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [y0]"m"(y[0]), [y1]"m"(y[1]), [y2]"m"(y[2])
    : );
    return carry;
}

/* x, y, and z have 3 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub */
static inline
mp_limb_t mpfq_fixmp_3_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "subq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "sbbq    %[y1], %%rax\n"
        "movq    %%rax, %[z1]\n"
        "movq    %[x2], %%rax\n"
        "sbbq    %[y2], %%rax\n"
        "movq    %%rax, %[z2]\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), [z2]"=m"(z[2]), "=&a"(carry)
    : [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [y0]"m"(y[0]), [y1]"m"(y[1]), [y2]"m"(y[2])
    : );
    return carry;
}

/* x, y, and z have 3 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_3_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "addq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    %[y1], %%rax\n"
        "movq    %%rax, %[z1]\n"
        "movq    %[x2], %%rax\n"
        "adcq    %[y2], %%rax\n"
        "movq    %%rax, %[z2]\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), [z2]"=m"(z[2])
    : [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [y0]"m"(y[0]), [y1]"m"(y[1]), [y2]"m"(y[2])
    : "%rax");
}

/* x, y, and z have 3 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_3_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "subq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "sbbq    %[y1], %%rax\n"
        "movq    %%rax, %[z1]\n"
        "movq    %[x2], %%rax\n"
        "sbbq    %[y2], %%rax\n"
        "movq    %%rax, %[z2]\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), [z2]"=m"(z[2])
    : [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [y0]"m"(y[0]), [y1]"m"(y[1]), [y2]"m"(y[2])
    : "%rax");
}

/* x has 3 words, z has 5.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_3_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x2], %%rax\n"
        "addq    %%rcx, %[z1]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, %[z2]\n"
        "adcq    $0, %%rdx\n"
        "xorq    %%rcx, %%rcx\n"
        "addq    %%rdx, %[z3]\n"
        "adcq    $0, %%rcx\n"
        "movq    %%rcx, %[carry]\n"
    : [carry]"=g"(carry), [z0]"+m"(z[0]), [z1]"+m"(z[1]), [z2]"+m"(z[2]), [z3]"+m"(z[3])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2])
    : "%rax", "%rcx", "%rdx");
    return carry;
}

/* x has 3 words, z has 5.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_nc */
static inline
void mpfq_fixmp_3_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x2], %%rax\n"
        "addq    %%rcx, %[z1]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, %[z2]\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rdx, %[z3]\n"
    : [z0]"+m"(z[0]), [z1]"+m"(z[1]), [z2]"+m"(z[2]), [z3]"+m"(z[3])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2])
    : "%rax", "%rcx", "%rdx");
}

/* x has 3 words, z has 4.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_shortz */
static inline
mp_limb_t mpfq_fixmp_3_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x2], %%rax\n"
        "addq    %%rcx, %[z1]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, %[z2]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %[carry]\n"
    : [carry]"=g"(carry), [z0]"+m"(z[0]), [z1]"+m"(z[1]), [z2]"+m"(z[2]), [z3]"+m"(z[3])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2])
    : "%rax", "%rcx", "%rdx");
    return carry;
}

/* x and y have 3 words, z has 8. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul */
static inline
void mpfq_fixmp_3_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "### x*y[0]\n"
        "movq    %2, %%r8\n"
        "movq    %0, %%rdi\n"
        "movq    0(%%r8), %%r9\n"
        "movq    %1, %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, 24(%%rdi)\n"
        "movq    $0, 32(%%rdi)\n"
        "movq    $0, 40(%%rdi)\n"
        "### x*y[1]\n"
        "movq    8(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 8(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 32(%%rdi)\n"
        "### x*y[2]\n"
        "movq    16(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 16(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 40(%%rdi)\n"
      : "+m" (z)
      : "m" (x), "m" (y)
      : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", "%r9", "memory");
}

/* x has 3 words, z has 8. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_3_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    __asm__ __volatile__(
        "movq    %1, %%rsi\n"
        "movq    %0, %%rdi\n"
        "### diagonal elements\n"
        "movq    (%%rsi), %%rax\n"
        "mulq    %%rax\n"
        "movq    %%rax, (%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, 8(%%rdi)\n"
        "mulq    %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rdx, 24(%%rdi)\n"
        "mulq    %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    %%rdx, 40(%%rdi)\n"
        "### precompute triangle\n"
        "### x[0]*x[1,2]\n"
        "movq    (%%rsi), %%rcx\n"
        "movq    8(%%rsi), %%rax\n"
        "mulq    %%rcx\n"
        "movq    %%rax, %%r8\n"
        "movq    %%rdx, %%r9\n"
        "movq    16(%%rsi), %%rax\n"
        "mulq    %%rcx\n"
        "addq    %%rax, %%r9\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%r10\n"
        "### x[1]*x[2]\n"
        "movq    8(%%rsi), %%rax\n"
        "mulq    16(%%rsi)\n"
        "addq    %%rax, %%r10\n"
        "adcq    $0, %%rdx\n"
        "### Shift triangle\n"
        "addq    %%r8, %%r8\n"
        "adcq    %%r9, %%r9\n"
        "adcq    %%r10, %%r10\n"
        "adcq    %%rdx, %%rdx\n"
        "adcq    $0, 40(%%rdi)\n"
        "### add shifted triangle to diagonal\n"
        "addq    %%r8, 8(%%rdi)\n"
        "adcq    %%r9, 16(%%rdi)\n"
        "adcq    %%r10, 24(%%rdi)\n"
        "adcq    %%rdx, 32(%%rdi)\n"
        "adcq    $0, 40(%%rdi)\n"
    : "+m" (z)
    : "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", "%r9", "%r10", "memory");
}

/* x has 3 words, z has 5. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_3_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x2], %%rax\n"
        "movq    %%rcx, %[z1]\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, %[z2]\n"
        "movq    %%rdx, %[z3]\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), [z2]"=m"(z[2]), [z3]"=m"(z[3])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2])
    : "%rax", "%rcx", "%rdx");
}

/* x, y, and z have 4 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add */
static inline
mp_limb_t mpfq_fixmp_4_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "addq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    %[y1], %%rax\n"
        "movq    %%rax, %[z1]\n"
        "movq    %[x2], %%rax\n"
        "adcq    %[y2], %%rax\n"
        "movq    %%rax, %[z2]\n"
        "movq    %[x3], %%rax\n"
        "adcq    %[y3], %%rax\n"
        "movq    %%rax, %[z3]\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), [z2]"=m"(z[2]), [z3]"=m"(z[3]), "=&a"(carry)
    : [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [x3]"m"(x[3]), [y0]"m"(y[0]), [y1]"m"(y[1]), [y2]"m"(y[2]), [y3]"m"(y[3])
    : );
    return carry;
}

/* x, y, and z have 4 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub */
static inline
mp_limb_t mpfq_fixmp_4_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "subq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "sbbq    %[y1], %%rax\n"
        "movq    %%rax, %[z1]\n"
        "movq    %[x2], %%rax\n"
        "sbbq    %[y2], %%rax\n"
        "movq    %%rax, %[z2]\n"
        "movq    %[x3], %%rax\n"
        "sbbq    %[y3], %%rax\n"
        "movq    %%rax, %[z3]\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), [z2]"=m"(z[2]), [z3]"=m"(z[3]), "=&a"(carry)
    : [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [x3]"m"(x[3]), [y0]"m"(y[0]), [y1]"m"(y[1]), [y2]"m"(y[2]), [y3]"m"(y[3])
    : );
    return carry;
}

/* x, y, and z have 4 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_4_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "addq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    %[y1], %%rax\n"
        "movq    %%rax, %[z1]\n"
        "movq    %[x2], %%rax\n"
        "adcq    %[y2], %%rax\n"
        "movq    %%rax, %[z2]\n"
        "movq    %[x3], %%rax\n"
        "adcq    %[y3], %%rax\n"
        "movq    %%rax, %[z3]\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), [z2]"=m"(z[2]), [z3]"=m"(z[3])
    : [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [x3]"m"(x[3]), [y0]"m"(y[0]), [y1]"m"(y[1]), [y2]"m"(y[2]), [y3]"m"(y[3])
    : "%rax");
}

/* x, y, and z have 4 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_4_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "subq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "sbbq    %[y1], %%rax\n"
        "movq    %%rax, %[z1]\n"
        "movq    %[x2], %%rax\n"
        "sbbq    %[y2], %%rax\n"
        "movq    %%rax, %[z2]\n"
        "movq    %[x3], %%rax\n"
        "sbbq    %[y3], %%rax\n"
        "movq    %%rax, %[z3]\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), [z2]"=m"(z[2]), [z3]"=m"(z[3])
    : [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [x3]"m"(x[3]), [y0]"m"(y[0]), [y1]"m"(y[1]), [y2]"m"(y[2]), [y3]"m"(y[3])
    : "%rax");
}

/* x has 4 words, z has 6.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_4_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x2], %%rax\n"
        "addq    %%rcx, %[z1]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x3], %%rax\n"
        "addq    %%rcx, %[z2]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, %[z3]\n"
        "adcq    $0, %%rdx\n"
        "xorq    %%rcx, %%rcx\n"
        "addq    %%rdx, %[z4]\n"
        "adcq    $0, %%rcx\n"
        "movq    %%rcx, %[carry]\n"
    : [carry]"=g"(carry), [z0]"+m"(z[0]), [z1]"+m"(z[1]), [z2]"+m"(z[2]), [z3]"+m"(z[3]), [z4]"+m"(z[4])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [x3]"m"(x[3])
    : "%rax", "%rcx", "%rdx");
    return carry;
}

/* x has 4 words, z has 6.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_nc */
static inline
void mpfq_fixmp_4_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x2], %%rax\n"
        "addq    %%rcx, %[z1]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x3], %%rax\n"
        "addq    %%rcx, %[z2]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, %[z3]\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rdx, %[z4]\n"
    : [z0]"+m"(z[0]), [z1]"+m"(z[1]), [z2]"+m"(z[2]), [z3]"+m"(z[3]), [z4]"+m"(z[4])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [x3]"m"(x[3])
    : "%rax", "%rcx", "%rdx");
}

/* x has 4 words, z has 5.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_shortz */
static inline
mp_limb_t mpfq_fixmp_4_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x2], %%rax\n"
        "addq    %%rcx, %[z1]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x3], %%rax\n"
        "addq    %%rcx, %[z2]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, %[z3]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %[carry]\n"
    : [carry]"=g"(carry), [z0]"+m"(z[0]), [z1]"+m"(z[1]), [z2]"+m"(z[2]), [z3]"+m"(z[3]), [z4]"+m"(z[4])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [x3]"m"(x[3])
    : "%rax", "%rcx", "%rdx");
    return carry;
}

/* x and y have 4 words, z has 10. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul */
static inline
void mpfq_fixmp_4_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "### x*y[0]\n"
        "movq    %2, %%r8\n"
        "movq    %0, %%rdi\n"
        "movq    0(%%r8), %%r9\n"
        "movq    %1, %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, 32(%%rdi)\n"
        "movq    $0, 40(%%rdi)\n"
        "movq    $0, 48(%%rdi)\n"
        "movq    $0, 56(%%rdi)\n"
        "### x*y[1]\n"
        "movq    8(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 8(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 40(%%rdi)\n"
        "### x*y[2]\n"
        "movq    16(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 16(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 48(%%rdi)\n"
        "### x*y[3]\n"
        "movq    24(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 24(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 56(%%rdi)\n"
      : "+m" (z)
      : "m" (x), "m" (y)
      : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", "%r9", "memory");
}

/* x has 4 words, z has 10. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_4_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    __asm__ __volatile__(
        "movq	%1, %%rsi\n"
        "movq	%0, %%rdi\n"
        "### diagonal elements\n"
        "movq    (%%rsi), %%rax\n"
        "mulq	%%rax\n"
        "movq    %%rax, (%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, 8(%%rdi)\n"
        "mulq	%%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq	16(%%rsi), %%rax\n"
        "movq    %%rdx, 24(%%rdi)\n"
        "mulq    %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq	24(%%rsi), %%rax\n"
        "movq    %%rdx, 40(%%rdi)\n"
        "mulq    %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
        "movq    %%rdx, 56(%%rdi)\n"
        "### precompute triangle\n"
        "### x[0]*x[1:3]\n"
        "movq	(%%rsi), %%rcx\n"
        "movq	8(%%rsi), %%rax\n"
        "mulq	%%rcx\n"
        "movq	%%rax, %%r8\n"
        "movq	%%rdx, %%r9\n"
        "movq    16(%%rsi), %%rax\n"
        "mulq	%%rcx\n"
        "addq	%%rax, %%r9\n"
        "adcq	$0, %%rdx\n"
        "movq	%%rdx, %%r10\n"
        "movq    24(%%rsi), %%rax\n"
        "mulq	%%rcx\n"
        "addq	%%rax, %%r10\n"
        "adcq	$0, %%rdx\n"
        "movq	%%rdx, %%r11\n"
        "### x[1]*x[2:3]\n"
        "movq	8(%%rsi), %%rcx\n"
        "movq	16(%%rsi), %%rax\n"
        "xorq	%%r12, %%r12\n"
        "mulq	%%rcx\n"
        "addq	%%rax, %%r10\n"
        "adcq	%%rdx, %%r11\n"
        "adcq	$0, %%r12\n"
        "movq	24(%%rsi), %%rax\n"
        "mulq	%%rcx\n"
        "addq    %%rax, %%r11\n"
        "adcq	$0, %%rdx\n"
        "addq    %%rdx, %%r12\n"
        "### x[2]*x[3]\n"
        "movq	16(%%rsi), %%rax\n"
        "mulq	24(%%rsi)\n"
        "addq	%%rax, %%r12\n"
        "adcq	$0, %%rdx\n"
        "### Shift triangle\n"
        "addq	%%r8, %%r8\n"
        "adcq	%%r9, %%r9\n"
        "adcq	%%r10, %%r10\n"
        "adcq	%%r11, %%r11\n"
        "adcq	%%r12, %%r12\n"
        "adcq	%%rdx, %%rdx\n"
        "adcq	$0, 56(%%rdi)\n"
        "### add shifted triangle to diagonal\n"
        "addq	%%r8, 8(%%rdi)\n"
        "adcq	%%r9, 16(%%rdi)\n"
        "adcq	%%r10, 24(%%rdi)\n"
        "adcq	%%r11, 32(%%rdi)\n"
        "adcq	%%r12, 40(%%rdi)\n"
        "adcq	%%rdx, 48(%%rdi)\n"
        "adcq	$0, 56(%%rdi)\n"
    : "+m" (z)
    : "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", "%r9", "%r10", "%r11", "%r12", "memory");
}

/* x has 4 words, z has 6. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_4_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x2], %%rax\n"
        "movq    %%rcx, %[z1]\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x3], %%rax\n"
        "movq    %%rcx, %[z2]\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, %[z3]\n"
        "movq    %%rdx, %[z4]\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), [z2]"=m"(z[2]), [z3]"=m"(z[3]), [z4]"=m"(z[4])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [x3]"m"(x[3])
    : "%rax", "%rcx", "%rdx");
}

/* x, y, and z have 5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add */
static inline
mp_limb_t mpfq_fixmp_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "addq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "adcq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "adcq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "adcq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : "=&a"(carry)
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x, y, and z have 5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub */
static inline
mp_limb_t mpfq_fixmp_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "subq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "sbbq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "sbbq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "sbbq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "sbbq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : "=&a"(carry)
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x, y, and z have 5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "addq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "adcq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "adcq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "adcq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
    : 
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory", "%rax");
}

/* x, y, and z have 5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "subq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "sbbq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "sbbq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "sbbq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "sbbq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
    : 
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory", "%rax");
}

/* x has 5 words, z has 7.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "xorq    %%rcx, %%rcx\n"
        "addq    %%rdx, 40(%%rdi)\n"
        "adcq    $0, %%rcx\n"
        "movq    %%rcx, %[carry]\n"
    : [carry]"=g"(carry), [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x has 5 words, z has 7.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_nc */
static inline
void mpfq_fixmp_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rdx, 40(%%rdi)\n"
    : [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x has 5 words, z has 6.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_shortz */
static inline
mp_limb_t mpfq_fixmp_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %[carry]\n"
    : [carry]"=g"(carry), [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x and y have 5 words, z has 12. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul */
static inline
void mpfq_fixmp_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "### x*y[0]\n"
        "movq    %2, %%r8\n"
        "movq    %0, %%rdi\n"
        "movq    0(%%r8), %%r9\n"
        "movq    %1, %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 32(%%rdi)\n"
        "movq    %%rdx, 40(%%rdi)\n"
        "movq    $0, 48(%%rdi)\n"
        "movq    $0, 56(%%rdi)\n"
        "movq    $0, 64(%%rdi)\n"
        "movq    $0, 72(%%rdi)\n"
        "### x*y[1]\n"
        "movq    8(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 8(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 48(%%rdi)\n"
        "### x*y[2]\n"
        "movq    16(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 16(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 56(%%rdi)\n"
        "### x*y[3]\n"
        "movq    24(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 24(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 64(%%rdi)\n"
        "### x*y[4]\n"
        "movq    32(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 32(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 72(%%rdi)\n"
      : "+m" (z)
      : "m" (x), "m" (y)
      : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", "%r9", "memory");
}

/* x has 5 words, z has 12. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mpfq_fixmp_5_mul(z, x, x);
}

/* x has 5 words, z has 7. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 32(%%rdi)\n"
        "movq    %%rdx, 40(%%rdi)\n"
    : 
    : [mult] "r" (c), [z] "m" (z), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x, y, and z have 6 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add */
static inline
mp_limb_t mpfq_fixmp_6_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "addq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "adcq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "adcq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "adcq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "adcq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : "=&a"(carry)
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x, y, and z have 6 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub */
static inline
mp_limb_t mpfq_fixmp_6_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "subq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "sbbq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "sbbq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "sbbq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "sbbq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "sbbq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : "=&a"(carry)
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x, y, and z have 6 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_6_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "addq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "adcq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "adcq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "adcq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "adcq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
    : 
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory", "%rax");
}

/* x, y, and z have 6 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_6_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "subq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "sbbq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "sbbq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "sbbq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "sbbq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "sbbq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
    : 
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory", "%rax");
}

/* x has 6 words, z has 8.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_6_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "xorq    %%rcx, %%rcx\n"
        "addq    %%rdx, 48(%%rdi)\n"
        "adcq    $0, %%rcx\n"
        "movq    %%rcx, %[carry]\n"
    : [carry]"=g"(carry), [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x has 6 words, z has 8.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_nc */
static inline
void mpfq_fixmp_6_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rdx, 48(%%rdi)\n"
    : [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x has 6 words, z has 7.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_shortz */
static inline
mp_limb_t mpfq_fixmp_6_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %[carry]\n"
    : [carry]"=g"(carry), [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x and y have 6 words, z has 14. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul */
static inline
void mpfq_fixmp_6_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "### x*y[0]\n"
        "movq    %2, %%r8\n"
        "movq    %0, %%rdi\n"
        "movq    0(%%r8), %%r9\n"
        "movq    %1, %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "movq    %%rcx, 32(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 40(%%rdi)\n"
        "movq    %%rdx, 48(%%rdi)\n"
        "movq    $0, 56(%%rdi)\n"
        "movq    $0, 64(%%rdi)\n"
        "movq    $0, 72(%%rdi)\n"
        "movq    $0, 80(%%rdi)\n"
        "movq    $0, 88(%%rdi)\n"
        "### x*y[1]\n"
        "movq    8(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 8(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 56(%%rdi)\n"
        "### x*y[2]\n"
        "movq    16(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 16(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 64(%%rdi)\n"
        "### x*y[3]\n"
        "movq    24(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 24(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 72(%%rdi)\n"
        "### x*y[4]\n"
        "movq    32(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 32(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 80(%%rdi)\n"
        "### x*y[5]\n"
        "movq    40(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 40(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 88(%%rdi)\n"
      : "+m" (z)
      : "m" (x), "m" (y)
      : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", "%r9", "memory");
}

/* x has 6 words, z has 14. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_6_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mpfq_fixmp_6_mul(z, x, x);
}

/* x has 6 words, z has 8. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_6_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "movq    %%rcx, 32(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 40(%%rdi)\n"
        "movq    %%rdx, 48(%%rdi)\n"
    : 
    : [mult] "r" (c), [z] "m" (z), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x, y, and z have 7 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add */
static inline
mp_limb_t mpfq_fixmp_7_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "addq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "adcq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "adcq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "adcq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "adcq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq    48(%%rsi), %%rax\n"
        "adcq    48(%%rdx), %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : "=&a"(carry)
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x, y, and z have 7 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub */
static inline
mp_limb_t mpfq_fixmp_7_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "subq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "sbbq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "sbbq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "sbbq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "sbbq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "sbbq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq    48(%%rsi), %%rax\n"
        "sbbq    48(%%rdx), %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : "=&a"(carry)
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x, y, and z have 7 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_7_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "addq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "adcq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "adcq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "adcq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "adcq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq    48(%%rsi), %%rax\n"
        "adcq    48(%%rdx), %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
    : 
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory", "%rax");
}

/* x, y, and z have 7 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_7_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "subq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "sbbq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "sbbq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "sbbq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "sbbq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "sbbq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq    48(%%rsi), %%rax\n"
        "sbbq    48(%%rdx), %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
    : 
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory", "%rax");
}

/* x has 7 words, z has 9.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_7_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "xorq    %%rcx, %%rcx\n"
        "addq    %%rdx, 56(%%rdi)\n"
        "adcq    $0, %%rcx\n"
        "movq    %%rcx, %[carry]\n"
    : [carry]"=g"(carry), [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x has 7 words, z has 9.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_nc */
static inline
void mpfq_fixmp_7_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rdx, 56(%%rdi)\n"
    : [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x has 7 words, z has 8.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_shortz */
static inline
mp_limb_t mpfq_fixmp_7_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %[carry]\n"
    : [carry]"=g"(carry), [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x and y have 7 words, z has 16. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul */
static inline
void mpfq_fixmp_7_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "### x*y[0]\n"
        "movq    %2, %%r8\n"
        "movq    %0, %%rdi\n"
        "movq    0(%%r8), %%r9\n"
        "movq    %1, %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "movq    %%rcx, 32(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "movq    %%rcx, 40(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 48(%%rdi)\n"
        "movq    %%rdx, 56(%%rdi)\n"
        "movq    $0, 64(%%rdi)\n"
        "movq    $0, 72(%%rdi)\n"
        "movq    $0, 80(%%rdi)\n"
        "movq    $0, 88(%%rdi)\n"
        "movq    $0, 96(%%rdi)\n"
        "movq    $0, 104(%%rdi)\n"
        "### x*y[1]\n"
        "movq    8(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 8(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 64(%%rdi)\n"
        "### x*y[2]\n"
        "movq    16(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 16(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 72(%%rdi)\n"
        "### x*y[3]\n"
        "movq    24(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 24(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 80(%%rdi)\n"
        "### x*y[4]\n"
        "movq    32(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 32(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 88(%%rdi)\n"
        "### x*y[5]\n"
        "movq    40(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 40(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 88(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 96(%%rdi)\n"
        "### x*y[6]\n"
        "movq    48(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 48(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 88(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 96(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 104(%%rdi)\n"
      : "+m" (z)
      : "m" (x), "m" (y)
      : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", "%r9", "memory");
}

/* x has 7 words, z has 16. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_7_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mpfq_fixmp_7_mul(z, x, x);
}

/* x has 7 words, z has 9. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_7_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "movq    %%rcx, 32(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "movq    %%rcx, 40(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 48(%%rdi)\n"
        "movq    %%rdx, 56(%%rdi)\n"
    : 
    : [mult] "r" (c), [z] "m" (z), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x, y, and z have 8 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add */
static inline
mp_limb_t mpfq_fixmp_8_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "addq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "adcq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "adcq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "adcq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "adcq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq    48(%%rsi), %%rax\n"
        "adcq    48(%%rdx), %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
        "movq    56(%%rsi), %%rax\n"
        "adcq    56(%%rdx), %%rax\n"
        "movq    %%rax, 56(%%rdi)\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : "=&a"(carry)
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x, y, and z have 8 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub */
static inline
mp_limb_t mpfq_fixmp_8_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "subq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "sbbq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "sbbq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "sbbq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "sbbq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "sbbq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq    48(%%rsi), %%rax\n"
        "sbbq    48(%%rdx), %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
        "movq    56(%%rsi), %%rax\n"
        "sbbq    56(%%rdx), %%rax\n"
        "movq    %%rax, 56(%%rdi)\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : "=&a"(carry)
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x, y, and z have 8 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_8_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "addq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "adcq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "adcq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "adcq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "adcq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq    48(%%rsi), %%rax\n"
        "adcq    48(%%rdx), %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
        "movq    56(%%rsi), %%rax\n"
        "adcq    56(%%rdx), %%rax\n"
        "movq    %%rax, 56(%%rdi)\n"
    : 
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory", "%rax");
}

/* x, y, and z have 8 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_8_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "subq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "sbbq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "sbbq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "sbbq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "sbbq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "sbbq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq    48(%%rsi), %%rax\n"
        "sbbq    48(%%rdx), %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
        "movq    56(%%rsi), %%rax\n"
        "sbbq    56(%%rdx), %%rax\n"
        "movq    %%rax, 56(%%rdi)\n"
    : 
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory", "%rax");
}

/* x has 8 words, z has 10.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_8_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "xorq    %%rcx, %%rcx\n"
        "addq    %%rdx, 64(%%rdi)\n"
        "adcq    $0, %%rcx\n"
        "movq    %%rcx, %[carry]\n"
    : [carry]"=g"(carry), [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x has 8 words, z has 10.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_nc */
static inline
void mpfq_fixmp_8_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rdx, 64(%%rdi)\n"
    : [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x has 8 words, z has 9.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_shortz */
static inline
mp_limb_t mpfq_fixmp_8_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %[carry]\n"
    : [carry]"=g"(carry), [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x and y have 8 words, z has 18. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul */
static inline
void mpfq_fixmp_8_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "### x*y[0]\n"
        "movq    %2, %%r8\n"
        "movq    %0, %%rdi\n"
        "movq    0(%%r8), %%r9\n"
        "movq    %1, %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "movq    %%rcx, 32(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "movq    %%rcx, 40(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "movq    %%rcx, 48(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 56(%%rdi)\n"
        "movq    %%rdx, 64(%%rdi)\n"
        "movq    $0, 72(%%rdi)\n"
        "movq    $0, 80(%%rdi)\n"
        "movq    $0, 88(%%rdi)\n"
        "movq    $0, 96(%%rdi)\n"
        "movq    $0, 104(%%rdi)\n"
        "movq    $0, 112(%%rdi)\n"
        "movq    $0, 120(%%rdi)\n"
        "### x*y[1]\n"
        "movq    8(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 8(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 72(%%rdi)\n"
        "### x*y[2]\n"
        "movq    16(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 16(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 80(%%rdi)\n"
        "### x*y[3]\n"
        "movq    24(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 24(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 88(%%rdi)\n"
        "### x*y[4]\n"
        "movq    32(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 32(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 88(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 96(%%rdi)\n"
        "### x*y[5]\n"
        "movq    40(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 40(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 88(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 96(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 104(%%rdi)\n"
        "### x*y[6]\n"
        "movq    48(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 48(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 88(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 96(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 104(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 112(%%rdi)\n"
        "### x*y[7]\n"
        "movq    56(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 56(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 88(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 96(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 104(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 112(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 120(%%rdi)\n"
      : "+m" (z)
      : "m" (x), "m" (y)
      : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", "%r9", "memory");
}

/* x has 8 words, z has 18. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_8_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mpfq_fixmp_8_mul(z, x, x);
}

/* x has 8 words, z has 10. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_8_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "movq    %%rcx, 32(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "movq    %%rcx, 40(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "movq    %%rcx, 48(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 56(%%rdi)\n"
        "movq    %%rdx, 64(%%rdi)\n"
    : 
    : [mult] "r" (c), [z] "m" (z), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x, y, and z have 9 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add */
static inline
mp_limb_t mpfq_fixmp_9_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "addq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "adcq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "adcq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "adcq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "adcq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq    48(%%rsi), %%rax\n"
        "adcq    48(%%rdx), %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
        "movq    56(%%rsi), %%rax\n"
        "adcq    56(%%rdx), %%rax\n"
        "movq    %%rax, 56(%%rdi)\n"
        "movq    64(%%rsi), %%rax\n"
        "adcq    64(%%rdx), %%rax\n"
        "movq    %%rax, 64(%%rdi)\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : "=&a"(carry)
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x, y, and z have 9 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub */
static inline
mp_limb_t mpfq_fixmp_9_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "subq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "sbbq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "sbbq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "sbbq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "sbbq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "sbbq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq    48(%%rsi), %%rax\n"
        "sbbq    48(%%rdx), %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
        "movq    56(%%rsi), %%rax\n"
        "sbbq    56(%%rdx), %%rax\n"
        "movq    %%rax, 56(%%rdi)\n"
        "movq    64(%%rsi), %%rax\n"
        "sbbq    64(%%rdx), %%rax\n"
        "movq    %%rax, 64(%%rdi)\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : "=&a"(carry)
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x, y, and z have 9 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_9_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "addq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "adcq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "adcq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "adcq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "adcq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq    48(%%rsi), %%rax\n"
        "adcq    48(%%rdx), %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
        "movq    56(%%rsi), %%rax\n"
        "adcq    56(%%rdx), %%rax\n"
        "movq    %%rax, 56(%%rdi)\n"
        "movq    64(%%rsi), %%rax\n"
        "adcq    64(%%rdx), %%rax\n"
        "movq    %%rax, 64(%%rdi)\n"
    : 
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory", "%rax");
}

/* x, y, and z have 9 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_9_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "subq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "sbbq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "sbbq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "sbbq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "sbbq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "sbbq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq    48(%%rsi), %%rax\n"
        "sbbq    48(%%rdx), %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
        "movq    56(%%rsi), %%rax\n"
        "sbbq    56(%%rdx), %%rax\n"
        "movq    %%rax, 56(%%rdi)\n"
        "movq    64(%%rsi), %%rax\n"
        "sbbq    64(%%rdx), %%rax\n"
        "movq    %%rax, 64(%%rdi)\n"
    : 
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory", "%rax");
}

/* x has 9 words, z has 11.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_9_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "xorq    %%rcx, %%rcx\n"
        "addq    %%rdx, 72(%%rdi)\n"
        "adcq    $0, %%rcx\n"
        "movq    %%rcx, %[carry]\n"
    : [carry]"=g"(carry), [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x has 9 words, z has 11.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_nc */
static inline
void mpfq_fixmp_9_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rdx, 72(%%rdi)\n"
    : [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x has 9 words, z has 10.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_shortz */
static inline
mp_limb_t mpfq_fixmp_9_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %[carry]\n"
    : [carry]"=g"(carry), [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x and y have 9 words, z has 20. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul */
static inline
void mpfq_fixmp_9_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "### x*y[0]\n"
        "movq    %2, %%r8\n"
        "movq    %0, %%rdi\n"
        "movq    0(%%r8), %%r9\n"
        "movq    %1, %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "movq    %%rcx, 32(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "movq    %%rcx, 40(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "movq    %%rcx, 48(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "movq    %%rcx, 56(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 64(%%rdi)\n"
        "movq    %%rdx, 72(%%rdi)\n"
        "movq    $0, 80(%%rdi)\n"
        "movq    $0, 88(%%rdi)\n"
        "movq    $0, 96(%%rdi)\n"
        "movq    $0, 104(%%rdi)\n"
        "movq    $0, 112(%%rdi)\n"
        "movq    $0, 120(%%rdi)\n"
        "movq    $0, 128(%%rdi)\n"
        "movq    $0, 136(%%rdi)\n"
        "### x*y[1]\n"
        "movq    8(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 8(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 80(%%rdi)\n"
        "### x*y[2]\n"
        "movq    16(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 16(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 88(%%rdi)\n"
        "### x*y[3]\n"
        "movq    24(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 24(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 88(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 96(%%rdi)\n"
        "### x*y[4]\n"
        "movq    32(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 32(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "addq    %%rcx, 88(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 96(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 104(%%rdi)\n"
        "### x*y[5]\n"
        "movq    40(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 40(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 88(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "addq    %%rcx, 96(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 104(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 112(%%rdi)\n"
        "### x*y[6]\n"
        "movq    48(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 48(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 88(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 96(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "addq    %%rcx, 104(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 112(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 120(%%rdi)\n"
        "### x*y[7]\n"
        "movq    56(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 56(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 88(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 96(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 104(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "addq    %%rcx, 112(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 120(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 128(%%rdi)\n"
        "### x*y[8]\n"
        "movq    64(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 64(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 88(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 96(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 104(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 112(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "addq    %%rcx, 120(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 128(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 136(%%rdi)\n"
      : "+m" (z)
      : "m" (x), "m" (y)
      : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", "%r9", "memory");
}

/* x has 9 words, z has 20. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_9_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mpfq_fixmp_9_mul(z, x, x);
}

/* x has 9 words, z has 11. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_9_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "movq    %%rcx, 32(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "movq    %%rcx, 40(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "movq    %%rcx, 48(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "movq    %%rcx, 56(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 64(%%rdi)\n"
        "movq    %%rdx, 72(%%rdi)\n"
    : 
    : [mult] "r" (c), [z] "m" (z), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x, y, and z have 0.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add */
static inline
mp_limb_t mpfq_fixmp_0_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "addq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : [z0]"=m"(z[0]), "=&a"(carry)
    : [x0]"m"(x[0]), [y0]"m"(y[0])
    : );
    return carry;
}

/* x, y, and z have 0.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub */
static inline
mp_limb_t mpfq_fixmp_0_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "subq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : [z0]"=m"(z[0]), "=&a"(carry)
    : [x0]"m"(x[0]), [y0]"m"(y[0])
    : );
    return carry;
}

/* x, y, and z have 0.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_0_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    *z = *x + *y;
}

/* x, y, and z have 0.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_0_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    *z = *x - *y;
}

/* x has 0.5 words, z has 2.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_0_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "adcq    $0, %%rdx\n"
        "xorq    %%rax, %%rax\n"
        "addq    %%rdx, %[z1]\n"
        "adcq    $0, %%rax\n"
    : [z0]"+m"(z[0]), [z1]"+m"(z[1]), "=a"(carry)
    : [x0]"m"(x[0]), [mult]"m"(c)
    : "%rdx");
    return carry;
}

/* x has 0.5 words, z has 2.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_nc */
static inline
void mpfq_fixmp_0_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rdx, %[z1]\n"
    : [z0]"+m"(z[0]), [z1]"+m"(z[1])
    : [x0]"m"(x[0]), [mult]"m"(c)
    : "%rax", "%rdx");
}

/* x has 0.5 words, z has 1.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_shortz */
static inline
mp_limb_t mpfq_fixmp_0_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "adcq    $0, %%rdx\n"
    : [z0]"+m"(z[0]), [z1]"+m"(z[1]), "=d"(carry)
    : [x0]"m"(x[0]), [mult]"m"(c)
    : "%rax");
    return carry;
}

/* x and y have 0.5 words, z has 1. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul */
static inline
void mpfq_fixmp_0_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "imulq %[y0], %%rax\n"
              : "=a" (z[0])
              : "0" (x[0]), [y0] "rm" (y[0])
              : "cc");
}

/* x has 0.5 words, z has 1. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_0_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    __asm__ __volatile__(
        "imulq %%rax, %%rax\n"
    : "=a" (z[0])
    : "0" (x[0])
    : "cc");
}

/* x has 0.5 words, z has 2. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_0_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "mulq    %[mult]\n"
    : [z0]"=a"(z[0]), [z1]"=d"(z[1])
    : [x0]"0"(x[0]), [mult]"m"(c)
    : );
}

/* x has 0.5 words, z has 1. c is 0.5 word.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul05 */
static inline
mp_limb_t mpfq_fixmp_0_5_addmul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "imulq    %[mult], %%rax\n"
        "xorq    %%rdx, %%rdx\n"
        "addq    %%rax, %[z0]\n"
        "adcq    $0, %%rdx\n"
    : [z0]"+m"(z[0]), "=d"(carry)
    : [x0]"m"(x[0]), [mult]"m"(c)
    : "%rax");
    return carry;
}

/* x has 0.5 words, z has 1. c is 0.5 word.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul05_nc */
static inline
void mpfq_fixmp_0_5_addmul05_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "imulq    %[mult], %%rax\n"
        "addq    %%rax, %[z0]\n"
    : [z0]"+m"(z[0])
    : [x0]"m"(x[0]), [mult]"m"(c)
    : "%rax");
}

/* x has 0.5 words, z has 1. c is 0.5 word.
 * Put (x*c) in z. No carry. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul05 */
static inline
void mpfq_fixmp_0_5_mul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "imulq    %[mult], %%rax\n"
        "movq    %%rax, %[z0]\n"
    : [z0]"=a"(z[0])
    : [x0]"0"(x[0]), [mult]"m"(c)
    : );
}

/* x, y, and z have 1.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add */
static inline
mp_limb_t mpfq_fixmp_1_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "addq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    %[y1], %%rax\n"
        "movq    %%rax, %[z1]\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), "=&a"(carry)
    : [x0]"m"(x[0]), [x1]"m"(x[1]), [y0]"m"(y[0]), [y1]"m"(y[1])
    : );
    return carry;
}

/* x, y, and z have 1.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub */
static inline
mp_limb_t mpfq_fixmp_1_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "subq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "sbbq    %[y1], %%rax\n"
        "movq    %%rax, %[z1]\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), "=&a"(carry)
    : [x0]"m"(x[0]), [x1]"m"(x[1]), [y0]"m"(y[0]), [y1]"m"(y[1])
    : );
    return carry;
}

/* x, y, and z have 1.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_1_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "addq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    %[y1], %%rax\n"
        "movq    %%rax, %[z1]\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1])
    : [x0]"m"(x[0]), [x1]"m"(x[1]), [y0]"m"(y[0]), [y1]"m"(y[1])
    : "%rax");
}

/* x, y, and z have 1.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_1_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "subq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "sbbq    %[y1], %%rax\n"
        "movq    %%rax, %[z1]\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1])
    : [x0]"m"(x[0]), [x1]"m"(x[1]), [y0]"m"(y[0]), [y1]"m"(y[1])
    : "%rax");
}

/* x has 1.5 words, z has 3.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_1_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, %[z1]\n"
        "adcq    $0, %%rdx\n"
        "xorq    %%rcx, %%rcx\n"
        "addq    %%rdx, %[z2]\n"
        "adcq    $0, %%rcx\n"
        "movq    %%rcx, %[carry]\n"
    : [carry]"=g"(carry), [z0]"+m"(z[0]), [z1]"+m"(z[1]), [z2]"+m"(z[2])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1])
    : "%rax", "%rcx", "%rdx");
    return carry;
}

/* x has 1.5 words, z has 3.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_nc */
static inline
void mpfq_fixmp_1_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, %[z1]\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rdx, %[z2]\n"
    : [z0]"+m"(z[0]), [z1]"+m"(z[1]), [z2]"+m"(z[2])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1])
    : "%rax", "%rcx", "%rdx");
}

/* x has 1.5 words, z has 2.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_shortz */
static inline
mp_limb_t mpfq_fixmp_1_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, %[z1]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %[carry]\n"
    : [carry]"=g"(carry), [z0]"+m"(z[0]), [z1]"+m"(z[1]), [z2]"+m"(z[2])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1])
    : "%rax", "%rcx", "%rdx");
    return carry;
}

/* x and y have 1.5 words, z has 3. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul */
static inline
void mpfq_fixmp_1_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "### x*y[0]\n"
        "movq    %2, %%r8\n"
        "movq    %0, %%rdi\n"
        "movq    0(%%r8), %%r9\n"
        "movq    %1, %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, 16(%%rdi)\n"
        "### x*y[1]\n"
        "movq    8(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 8(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 16(%%rdi)\n"
      : "+m" (z)
      : "m" (x), "m" (y)
      : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", "%r9", "memory");
}

/* x has 1.5 words, z has 3. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_1_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    __asm__ __volatile__(
        "movq    %1, %%r8\n"
        "movq    %0, %%rdi\n"
        "movq    (%%r8), %%rax\n"
        "mulq    %%rax\n"
        "movq    %%rax, (%%rdi)\n"
        "movq    8(%%r8), %%rax\n"
        "movq    %%rdx, %%r9\n"
        "mulq    %%rax\n"
        "movq    %%rax, %%r10\n"
        "movq    (%%r8), %%rax\n"
        "mulq    8(%%r8)\n"
        "addq    %%rax, %%r9\n"
        "adcq    %%rdx, %%r10\n"
        "addq    %%rax, %%r9\n"
        "movq    %%r9, 8(%%rdi)\n"
        "adcq    %%rdx, %%r10\n"
        "movq    %%r10, 16(%%rdi)\n"
    : "+m" (z)
    : "m" (x)
    : "%rax", "%rdx", "%rdi", "%r8", "%r9", "%r10", "%r11", "memory");
}

/* x has 1.5 words, z has 3. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_1_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, %[z1]\n"
        "movq    %%rdx, %[z2]\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), [z2]"=m"(z[2])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1])
    : "%rax", "%rcx", "%rdx");
}

/* x has 1.5 words, z has 2. c is 0.5 word.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul05 */
static inline
mp_limb_t mpfq_fixmp_1_5_addmul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "imulq    %[mult], %%rax\n"
        "xorq    %%rdx, %%rdx\n"
        "addq    %%rcx, %%rax\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rax, %[z1]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %[carry]\n"
    : [carry]"=g"(carry), [z0]"+m"(z[0]), [z1]"+m"(z[1])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1])
    : "%rax", "%rcx", "%rdx");
    return carry;
}

/* x has 1.5 words, z has 2. c is 0.5 word.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul05_nc */
static inline
void mpfq_fixmp_1_5_addmul05_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "imulq    %[mult], %%rax\n"
        "addq    %%rcx, %%rax\n"
        "addq    %%rax, %[z1]\n"
    : [z0]"+m"(z[0]), [z1]"+m"(z[1])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1])
    : "%rax", "%rcx", "%rdx");
}

/* x has 1.5 words, z has 2. c is 0.5 word.
 * Put (x*c) in z. No carry. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul05 */
static inline
void mpfq_fixmp_1_5_mul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "imulq    %[mult], %%rax\n"
        "addq    %%rcx, %%rax\n"
        "movq    %%rax, %[z1]\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1])
    : "%rax", "%rcx", "%rdx");
}

/* x, y, and z have 2.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add */
static inline
mp_limb_t mpfq_fixmp_2_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "addq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    %[y1], %%rax\n"
        "movq    %%rax, %[z1]\n"
        "movq    %[x2], %%rax\n"
        "adcq    %[y2], %%rax\n"
        "movq    %%rax, %[z2]\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), [z2]"=m"(z[2]), "=&a"(carry)
    : [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [y0]"m"(y[0]), [y1]"m"(y[1]), [y2]"m"(y[2])
    : );
    return carry;
}

/* x, y, and z have 2.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub */
static inline
mp_limb_t mpfq_fixmp_2_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "subq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "sbbq    %[y1], %%rax\n"
        "movq    %%rax, %[z1]\n"
        "movq    %[x2], %%rax\n"
        "sbbq    %[y2], %%rax\n"
        "movq    %%rax, %[z2]\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), [z2]"=m"(z[2]), "=&a"(carry)
    : [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [y0]"m"(y[0]), [y1]"m"(y[1]), [y2]"m"(y[2])
    : );
    return carry;
}

/* x, y, and z have 2.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_2_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "addq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    %[y1], %%rax\n"
        "movq    %%rax, %[z1]\n"
        "movq    %[x2], %%rax\n"
        "adcq    %[y2], %%rax\n"
        "movq    %%rax, %[z2]\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), [z2]"=m"(z[2])
    : [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [y0]"m"(y[0]), [y1]"m"(y[1]), [y2]"m"(y[2])
    : "%rax");
}

/* x, y, and z have 2.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_2_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "subq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "sbbq    %[y1], %%rax\n"
        "movq    %%rax, %[z1]\n"
        "movq    %[x2], %%rax\n"
        "sbbq    %[y2], %%rax\n"
        "movq    %%rax, %[z2]\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), [z2]"=m"(z[2])
    : [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [y0]"m"(y[0]), [y1]"m"(y[1]), [y2]"m"(y[2])
    : "%rax");
}

/* x has 2.5 words, z has 4.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_2_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x2], %%rax\n"
        "addq    %%rcx, %[z1]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, %[z2]\n"
        "adcq    $0, %%rdx\n"
        "xorq    %%rcx, %%rcx\n"
        "addq    %%rdx, %[z3]\n"
        "adcq    $0, %%rcx\n"
        "movq    %%rcx, %[carry]\n"
    : [carry]"=g"(carry), [z0]"+m"(z[0]), [z1]"+m"(z[1]), [z2]"+m"(z[2]), [z3]"+m"(z[3])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2])
    : "%rax", "%rcx", "%rdx");
    return carry;
}

/* x has 2.5 words, z has 4.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_nc */
static inline
void mpfq_fixmp_2_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x2], %%rax\n"
        "addq    %%rcx, %[z1]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, %[z2]\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rdx, %[z3]\n"
    : [z0]"+m"(z[0]), [z1]"+m"(z[1]), [z2]"+m"(z[2]), [z3]"+m"(z[3])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2])
    : "%rax", "%rcx", "%rdx");
}

/* x has 2.5 words, z has 3.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_shortz */
static inline
mp_limb_t mpfq_fixmp_2_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x2], %%rax\n"
        "addq    %%rcx, %[z1]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, %[z2]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %[carry]\n"
    : [carry]"=g"(carry), [z0]"+m"(z[0]), [z1]"+m"(z[1]), [z2]"+m"(z[2]), [z3]"+m"(z[3])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2])
    : "%rax", "%rcx", "%rdx");
    return carry;
}

/* x and y have 2.5 words, z has 5. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul */
static inline
void mpfq_fixmp_2_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "### x*y[0]\n"
        "movq    %2, %%r8\n"
        "movq    %0, %%rdi\n"
        "movq    0(%%r8), %%r9\n"
        "movq    %1, %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, 24(%%rdi)\n"
        "movq    $0, 32(%%rdi)\n"
        "### x*y[1]\n"
        "movq    8(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 8(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 32(%%rdi)\n"
        "### x*y[2]\n"
        "movq    16(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 16(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 32(%%rdi)\n"
      : "+m" (z)
      : "m" (x), "m" (y)
      : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", "%r9", "memory");
}

/* x has 2.5 words, z has 5. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_2_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    __asm__ __volatile__(
        "movq    %1, %%rsi\n"
        "movq    %0, %%rdi\n"
        "### diagonal elements\n"
        "movq    (%%rsi), %%rax\n"
        "mulq    %%rax\n"
        "movq    %%rax, (%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, 8(%%rdi)\n"
        "mulq    %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rdx, 24(%%rdi)\n"
        "mulq    %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "### precompute triangle\n"
        "### x[0]*x[1,2]\n"
        "movq    (%%rsi), %%rcx\n"
        "movq    8(%%rsi), %%rax\n"
        "mulq    %%rcx\n"
        "movq    %%rax, %%r8\n"
        "movq    %%rdx, %%r9\n"
        "movq    16(%%rsi), %%rax\n"
        "mulq    %%rcx\n"
        "addq    %%rax, %%r9\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%r10\n"
        "### x[1]*x[2]\n"
        "movq    8(%%rsi), %%rax\n"
        "mulq    16(%%rsi)\n"
        "addq    %%rax, %%r10\n"
        "adcq    $0, %%rdx\n"
        "### Shift triangle\n"
        "addq    %%r8, %%r8\n"
        "adcq    %%r9, %%r9\n"
        "adcq    %%r10, %%r10\n"
        "adcq    %%rdx, %%rdx\n"
        "### add shifted triangle to diagonal\n"
        "addq    %%r8, 8(%%rdi)\n"
        "adcq    %%r9, 16(%%rdi)\n"
        "adcq    %%r10, 24(%%rdi)\n"
        "adcq    %%rdx, 32(%%rdi)\n"
    : "+m" (z)
    : "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", "%r9", "%r10", "memory");
}

/* x has 2.5 words, z has 4. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_2_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x2], %%rax\n"
        "movq    %%rcx, %[z1]\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, %[z2]\n"
        "movq    %%rdx, %[z3]\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), [z2]"=m"(z[2]), [z3]"=m"(z[3])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2])
    : "%rax", "%rcx", "%rdx");
}

/* x has 2.5 words, z has 3. c is 0.5 word.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul05 */
static inline
mp_limb_t mpfq_fixmp_2_5_addmul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x2], %%rax\n"
        "addq    %%rcx, %[z1]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "imulq    %[mult], %%rax\n"
        "xorq    %%rdx, %%rdx\n"
        "addq    %%rcx, %%rax\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rax, %[z2]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %[carry]\n"
    : [carry]"=g"(carry), [z0]"+m"(z[0]), [z1]"+m"(z[1]), [z2]"+m"(z[2])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2])
    : "%rax", "%rcx", "%rdx");
    return carry;
}

/* x has 2.5 words, z has 3. c is 0.5 word.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul05_nc */
static inline
void mpfq_fixmp_2_5_addmul05_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x2], %%rax\n"
        "addq    %%rcx, %[z1]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "imulq    %[mult], %%rax\n"
        "addq    %%rcx, %%rax\n"
        "addq    %%rax, %[z2]\n"
    : [z0]"+m"(z[0]), [z1]"+m"(z[1]), [z2]"+m"(z[2])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2])
    : "%rax", "%rcx", "%rdx");
}

/* x has 2.5 words, z has 3. c is 0.5 word.
 * Put (x*c) in z. No carry. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul05 */
static inline
void mpfq_fixmp_2_5_mul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x2], %%rax\n"
        "movq    %%rcx, %[z1]\n"
        "movq    %%rdx, %%rcx\n"
        "imulq    %[mult], %%rax\n"
        "addq    %%rcx, %%rax\n"
        "movq    %%rax, %[z2]\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), [z2]"=m"(z[2])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2])
    : "%rax", "%rcx", "%rdx");
}

/* x, y, and z have 3.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add */
static inline
mp_limb_t mpfq_fixmp_3_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "addq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    %[y1], %%rax\n"
        "movq    %%rax, %[z1]\n"
        "movq    %[x2], %%rax\n"
        "adcq    %[y2], %%rax\n"
        "movq    %%rax, %[z2]\n"
        "movq    %[x3], %%rax\n"
        "adcq    %[y3], %%rax\n"
        "movq    %%rax, %[z3]\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), [z2]"=m"(z[2]), [z3]"=m"(z[3]), "=&a"(carry)
    : [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [x3]"m"(x[3]), [y0]"m"(y[0]), [y1]"m"(y[1]), [y2]"m"(y[2]), [y3]"m"(y[3])
    : );
    return carry;
}

/* x, y, and z have 3.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub */
static inline
mp_limb_t mpfq_fixmp_3_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "subq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "sbbq    %[y1], %%rax\n"
        "movq    %%rax, %[z1]\n"
        "movq    %[x2], %%rax\n"
        "sbbq    %[y2], %%rax\n"
        "movq    %%rax, %[z2]\n"
        "movq    %[x3], %%rax\n"
        "sbbq    %[y3], %%rax\n"
        "movq    %%rax, %[z3]\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), [z2]"=m"(z[2]), [z3]"=m"(z[3]), "=&a"(carry)
    : [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [x3]"m"(x[3]), [y0]"m"(y[0]), [y1]"m"(y[1]), [y2]"m"(y[2]), [y3]"m"(y[3])
    : );
    return carry;
}

/* x, y, and z have 3.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_3_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "addq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    %[y1], %%rax\n"
        "movq    %%rax, %[z1]\n"
        "movq    %[x2], %%rax\n"
        "adcq    %[y2], %%rax\n"
        "movq    %%rax, %[z2]\n"
        "movq    %[x3], %%rax\n"
        "adcq    %[y3], %%rax\n"
        "movq    %%rax, %[z3]\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), [z2]"=m"(z[2]), [z3]"=m"(z[3])
    : [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [x3]"m"(x[3]), [y0]"m"(y[0]), [y1]"m"(y[1]), [y2]"m"(y[2]), [y3]"m"(y[3])
    : "%rax");
}

/* x, y, and z have 3.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_3_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "subq    %[y0], %%rax\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "sbbq    %[y1], %%rax\n"
        "movq    %%rax, %[z1]\n"
        "movq    %[x2], %%rax\n"
        "sbbq    %[y2], %%rax\n"
        "movq    %%rax, %[z2]\n"
        "movq    %[x3], %%rax\n"
        "sbbq    %[y3], %%rax\n"
        "movq    %%rax, %[z3]\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), [z2]"=m"(z[2]), [z3]"=m"(z[3])
    : [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [x3]"m"(x[3]), [y0]"m"(y[0]), [y1]"m"(y[1]), [y2]"m"(y[2]), [y3]"m"(y[3])
    : "%rax");
}

/* x has 3.5 words, z has 5.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_3_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x2], %%rax\n"
        "addq    %%rcx, %[z1]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x3], %%rax\n"
        "addq    %%rcx, %[z2]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, %[z3]\n"
        "adcq    $0, %%rdx\n"
        "xorq    %%rcx, %%rcx\n"
        "addq    %%rdx, %[z4]\n"
        "adcq    $0, %%rcx\n"
        "movq    %%rcx, %[carry]\n"
    : [carry]"=g"(carry), [z0]"+m"(z[0]), [z1]"+m"(z[1]), [z2]"+m"(z[2]), [z3]"+m"(z[3]), [z4]"+m"(z[4])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [x3]"m"(x[3])
    : "%rax", "%rcx", "%rdx");
    return carry;
}

/* x has 3.5 words, z has 5.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_nc */
static inline
void mpfq_fixmp_3_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x2], %%rax\n"
        "addq    %%rcx, %[z1]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x3], %%rax\n"
        "addq    %%rcx, %[z2]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, %[z3]\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rdx, %[z4]\n"
    : [z0]"+m"(z[0]), [z1]"+m"(z[1]), [z2]"+m"(z[2]), [z3]"+m"(z[3]), [z4]"+m"(z[4])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [x3]"m"(x[3])
    : "%rax", "%rcx", "%rdx");
}

/* x has 3.5 words, z has 4.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_shortz */
static inline
mp_limb_t mpfq_fixmp_3_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x2], %%rax\n"
        "addq    %%rcx, %[z1]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x3], %%rax\n"
        "addq    %%rcx, %[z2]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, %[z3]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %[carry]\n"
    : [carry]"=g"(carry), [z0]"+m"(z[0]), [z1]"+m"(z[1]), [z2]"+m"(z[2]), [z3]"+m"(z[3]), [z4]"+m"(z[4])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [x3]"m"(x[3])
    : "%rax", "%rcx", "%rdx");
    return carry;
}

/* x and y have 3.5 words, z has 7. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul */
static inline
void mpfq_fixmp_3_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "### x*y[0]\n"
        "movq    %2, %%r8\n"
        "movq    %0, %%rdi\n"
        "movq    0(%%r8), %%r9\n"
        "movq    %1, %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, 32(%%rdi)\n"
        "movq    $0, 40(%%rdi)\n"
        "movq    $0, 48(%%rdi)\n"
        "### x*y[1]\n"
        "movq    8(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 8(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 40(%%rdi)\n"
        "### x*y[2]\n"
        "movq    16(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 16(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 48(%%rdi)\n"
        "### x*y[3]\n"
        "movq    24(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 24(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 48(%%rdi)\n"
      : "+m" (z)
      : "m" (x), "m" (y)
      : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", "%r9", "memory");
}

/* x has 3.5 words, z has 7. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_3_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    __asm__ __volatile__(
        "movq    %1, %%rsi\n"
        "movq    %0, %%rdi\n"
        "### diagonal elements\n"
        "movq    (%%rsi), %%rax\n"
        "mulq    %%rax\n"
        "movq    %%rax, (%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, 8(%%rdi)\n"
        "mulq    %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rdx, 24(%%rdi)\n"
        "mulq    %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rdx, 40(%%rdi)\n"
        "mulq    %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
        "### precompute triangle\n"
        "### x[0]*x[1:3]\n"
        "movq    (%%rsi), %%rcx\n"
        "movq    8(%%rsi), %%rax\n"
        "mulq    %%rcx\n"
        "movq    %%rax, %%r8\n"
        "movq    %%rdx, %%r9\n"
        "movq    16(%%rsi), %%rax\n"
        "mulq    %%rcx\n"
        "addq    %%rax, %%r9\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%r10\n"
        "movq    24(%%rsi), %%rax\n"
        "mulq    %%rcx\n"
        "addq    %%rax, %%r10\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%r11\n"
        "### x[1]*x[2:3]\n"
        "movq    8(%%rsi), %%rcx\n"
        "movq    16(%%rsi), %%rax\n"
        "xorq    %%r12, %%r12\n"
        "mulq    %%rcx\n"
        "addq    %%rax, %%r10\n"
        "adcq    %%rdx, %%r11\n"
        "adcq    $0, %%r12\n"
        "movq    24(%%rsi), %%rax\n"
        "mulq    %%rcx\n"
        "addq    %%rax, %%r11\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rdx, %%r12\n"
        "### x[2]*x[3]\n"
        "movq    16(%%rsi), %%rax\n"
        "mulq    24(%%rsi)\n"
        "addq    %%rax, %%r12\n"
        "adcq    $0, %%rdx\n"
        "### Shift triangle\n"
        "addq    %%r8, %%r8\n"
        "adcq    %%r9, %%r9\n"
        "adcq    %%r10, %%r10\n"
        "adcq    %%r11, %%r11\n"
        "adcq    %%r12, %%r12\n"
        "adcq    %%rdx, %%rdx\n"
        "### add shifted triangle to diagonal\n"
        "addq    %%r8, 8(%%rdi)\n"
        "adcq    %%r9, 16(%%rdi)\n"
        "adcq    %%r10, 24(%%rdi)\n"
        "adcq    %%r11, 32(%%rdi)\n"
        "adcq    %%r12, 40(%%rdi)\n"
        "adcq    %%rdx, 48(%%rdi)\n"
    : "+m" (z)
    : "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", "%r9", "%r10", "%r11", "%r12", "memory");
}

/* x has 3.5 words, z has 5. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_3_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x2], %%rax\n"
        "movq    %%rcx, %[z1]\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x3], %%rax\n"
        "movq    %%rcx, %[z2]\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, %[z3]\n"
        "movq    %%rdx, %[z4]\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), [z2]"=m"(z[2]), [z3]"=m"(z[3]), [z4]"=m"(z[4])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [x3]"m"(x[3])
    : "%rax", "%rcx", "%rdx");
}

/* x has 3.5 words, z has 4. c is 0.5 word.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul05 */
static inline
mp_limb_t mpfq_fixmp_3_5_addmul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x2], %%rax\n"
        "addq    %%rcx, %[z1]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x3], %%rax\n"
        "addq    %%rcx, %[z2]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "imulq    %[mult], %%rax\n"
        "xorq    %%rdx, %%rdx\n"
        "addq    %%rcx, %%rax\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rax, %[z3]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %[carry]\n"
    : [carry]"=g"(carry), [z0]"+m"(z[0]), [z1]"+m"(z[1]), [z2]"+m"(z[2]), [z3]"+m"(z[3])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [x3]"m"(x[3])
    : "%rax", "%rcx", "%rdx");
    return carry;
}

/* x has 3.5 words, z has 4. c is 0.5 word.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul05_nc */
static inline
void mpfq_fixmp_3_5_addmul05_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x2], %%rax\n"
        "addq    %%rcx, %[z1]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x3], %%rax\n"
        "addq    %%rcx, %[z2]\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "imulq    %[mult], %%rax\n"
        "addq    %%rcx, %%rax\n"
        "addq    %%rax, %[z3]\n"
    : [z0]"+m"(z[0]), [z1]"+m"(z[1]), [z2]"+m"(z[2]), [z3]"+m"(z[3])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [x3]"m"(x[3])
    : "%rax", "%rcx", "%rdx");
}

/* x has 3.5 words, z has 4. c is 0.5 word.
 * Put (x*c) in z. No carry. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul05 */
static inline
void mpfq_fixmp_3_5_mul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[x0], %%rax\n"
        "mulq    %[mult]\n"
        "movq    %%rax, %[z0]\n"
        "movq    %[x1], %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x2], %%rax\n"
        "movq    %%rcx, %[z1]\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %[x3], %%rax\n"
        "movq    %%rcx, %[z2]\n"
        "movq    %%rdx, %%rcx\n"
        "imulq    %[mult], %%rax\n"
        "addq    %%rcx, %%rax\n"
        "movq    %%rax, %[z3]\n"
    : [z0]"=m"(z[0]), [z1]"=m"(z[1]), [z2]"=m"(z[2]), [z3]"=m"(z[3])
    : [mult] "r" (c), [x0]"m"(x[0]), [x1]"m"(x[1]), [x2]"m"(x[2]), [x3]"m"(x[3])
    : "%rax", "%rcx", "%rdx");
}

/* x, y, and z have 4.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add */
static inline
mp_limb_t mpfq_fixmp_4_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "addq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "adcq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "adcq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "adcq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : "=&a"(carry)
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x, y, and z have 4.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub */
static inline
mp_limb_t mpfq_fixmp_4_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "subq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "sbbq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "sbbq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "sbbq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "sbbq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : "=&a"(carry)
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x, y, and z have 4.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_4_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "addq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "adcq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "adcq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "adcq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
    : 
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory", "%rax");
}

/* x, y, and z have 4.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_4_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "subq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "sbbq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "sbbq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "sbbq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "sbbq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
    : 
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory", "%rax");
}

/* x has 4.5 words, z has 6.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_4_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "xorq    %%rcx, %%rcx\n"
        "addq    %%rdx, 40(%%rdi)\n"
        "adcq    $0, %%rcx\n"
        "movq    %%rcx, %[carry]\n"
    : [carry]"=g"(carry), [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x has 4.5 words, z has 6.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_nc */
static inline
void mpfq_fixmp_4_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rdx, 40(%%rdi)\n"
    : [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x has 4.5 words, z has 5.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_shortz */
static inline
mp_limb_t mpfq_fixmp_4_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %[carry]\n"
    : [carry]"=g"(carry), [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x and y have 4.5 words, z has 9. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul */
static inline
void mpfq_fixmp_4_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "### x*y[0]\n"
        "movq    %2, %%r8\n"
        "movq    %0, %%rdi\n"
        "movq    0(%%r8), %%r9\n"
        "movq    %1, %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 32(%%rdi)\n"
        "movq    %%rdx, 40(%%rdi)\n"
        "movq    $0, 48(%%rdi)\n"
        "movq    $0, 56(%%rdi)\n"
        "movq    $0, 64(%%rdi)\n"
        "### x*y[1]\n"
        "movq    8(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 8(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 48(%%rdi)\n"
        "### x*y[2]\n"
        "movq    16(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 16(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 56(%%rdi)\n"
        "### x*y[3]\n"
        "movq    24(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 24(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 64(%%rdi)\n"
        "### x*y[4]\n"
        "movq    32(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 32(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 64(%%rdi)\n"
      : "+m" (z)
      : "m" (x), "m" (y)
      : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", "%r9", "memory");
}

/* x has 4.5 words, z has 9. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_4_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mpfq_fixmp_4_5_mul(z, x, x);
}

/* x has 4.5 words, z has 6. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_4_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 32(%%rdi)\n"
        "movq    %%rdx, 40(%%rdi)\n"
    : 
    : [mult] "r" (c), [z] "m" (z), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x has 4.5 words, z has 5. c is 0.5 word.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul05 */
static inline
mp_limb_t mpfq_fixmp_4_5_addmul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "imulq    %[mult], %%rax\n"
        "xorq    %%rdx, %%rdx\n"
        "addq    %%rcx, %%rax\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rax, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %[carry]\n"
    : [carry]"=g"(carry), [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x has 4.5 words, z has 5. c is 0.5 word.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul05_nc */
static inline
void mpfq_fixmp_4_5_addmul05_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "imulq    %[mult], %%rax\n"
        "addq    %%rcx, %%rax\n"
        "addq    %%rax, 32(%%rdi)\n"
    : [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x has 4.5 words, z has 5. c is 0.5 word.
 * Put (x*c) in z. No carry. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul05 */
static inline
void mpfq_fixmp_4_5_mul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "imulq    %[mult], %%rax\n"
        "addq    %%rcx, %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
    : 
    : [mult] "r" (c), [z] "m" (z), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x, y, and z have 5.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add */
static inline
mp_limb_t mpfq_fixmp_5_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "addq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "adcq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "adcq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "adcq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "adcq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : "=&a"(carry)
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x, y, and z have 5.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub */
static inline
mp_limb_t mpfq_fixmp_5_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "subq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "sbbq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "sbbq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "sbbq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "sbbq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "sbbq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : "=&a"(carry)
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x, y, and z have 5.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_5_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "addq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "adcq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "adcq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "adcq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "adcq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
    : 
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory", "%rax");
}

/* x, y, and z have 5.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_5_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "subq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "sbbq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "sbbq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "sbbq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "sbbq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "sbbq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
    : 
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory", "%rax");
}

/* x has 5.5 words, z has 7.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_5_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "xorq    %%rcx, %%rcx\n"
        "addq    %%rdx, 48(%%rdi)\n"
        "adcq    $0, %%rcx\n"
        "movq    %%rcx, %[carry]\n"
    : [carry]"=g"(carry), [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x has 5.5 words, z has 7.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_nc */
static inline
void mpfq_fixmp_5_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rdx, 48(%%rdi)\n"
    : [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x has 5.5 words, z has 6.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_shortz */
static inline
mp_limb_t mpfq_fixmp_5_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %[carry]\n"
    : [carry]"=g"(carry), [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x and y have 5.5 words, z has 11. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul */
static inline
void mpfq_fixmp_5_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "### x*y[0]\n"
        "movq    %2, %%r8\n"
        "movq    %0, %%rdi\n"
        "movq    0(%%r8), %%r9\n"
        "movq    %1, %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "movq    %%rcx, 32(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 40(%%rdi)\n"
        "movq    %%rdx, 48(%%rdi)\n"
        "movq    $0, 56(%%rdi)\n"
        "movq    $0, 64(%%rdi)\n"
        "movq    $0, 72(%%rdi)\n"
        "movq    $0, 80(%%rdi)\n"
        "### x*y[1]\n"
        "movq    8(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 8(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 56(%%rdi)\n"
        "### x*y[2]\n"
        "movq    16(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 16(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 64(%%rdi)\n"
        "### x*y[3]\n"
        "movq    24(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 24(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 72(%%rdi)\n"
        "### x*y[4]\n"
        "movq    32(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 32(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 80(%%rdi)\n"
        "### x*y[5]\n"
        "movq    40(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 40(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 80(%%rdi)\n"
      : "+m" (z)
      : "m" (x), "m" (y)
      : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", "%r9", "memory");
}

/* x has 5.5 words, z has 11. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_5_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mpfq_fixmp_5_5_mul(z, x, x);
}

/* x has 5.5 words, z has 7. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_5_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "movq    %%rcx, 32(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 40(%%rdi)\n"
        "movq    %%rdx, 48(%%rdi)\n"
    : 
    : [mult] "r" (c), [z] "m" (z), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x has 5.5 words, z has 6. c is 0.5 word.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul05 */
static inline
mp_limb_t mpfq_fixmp_5_5_addmul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "imulq    %[mult], %%rax\n"
        "xorq    %%rdx, %%rdx\n"
        "addq    %%rcx, %%rax\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rax, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %[carry]\n"
    : [carry]"=g"(carry), [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x has 5.5 words, z has 6. c is 0.5 word.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul05_nc */
static inline
void mpfq_fixmp_5_5_addmul05_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "imulq    %[mult], %%rax\n"
        "addq    %%rcx, %%rax\n"
        "addq    %%rax, 40(%%rdi)\n"
    : [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x has 5.5 words, z has 6. c is 0.5 word.
 * Put (x*c) in z. No carry. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul05 */
static inline
void mpfq_fixmp_5_5_mul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "movq    %%rcx, 32(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "imulq    %[mult], %%rax\n"
        "addq    %%rcx, %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
    : 
    : [mult] "r" (c), [z] "m" (z), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x, y, and z have 6.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add */
static inline
mp_limb_t mpfq_fixmp_6_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "addq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "adcq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "adcq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "adcq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "adcq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq    48(%%rsi), %%rax\n"
        "adcq    48(%%rdx), %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : "=&a"(carry)
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x, y, and z have 6.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub */
static inline
mp_limb_t mpfq_fixmp_6_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "subq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "sbbq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "sbbq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "sbbq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "sbbq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "sbbq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq    48(%%rsi), %%rax\n"
        "sbbq    48(%%rdx), %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : "=&a"(carry)
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x, y, and z have 6.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_6_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "addq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "adcq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "adcq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "adcq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "adcq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq    48(%%rsi), %%rax\n"
        "adcq    48(%%rdx), %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
    : 
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory", "%rax");
}

/* x, y, and z have 6.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_6_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "subq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "sbbq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "sbbq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "sbbq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "sbbq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "sbbq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq    48(%%rsi), %%rax\n"
        "sbbq    48(%%rdx), %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
    : 
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory", "%rax");
}

/* x has 6.5 words, z has 8.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_6_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "xorq    %%rcx, %%rcx\n"
        "addq    %%rdx, 56(%%rdi)\n"
        "adcq    $0, %%rcx\n"
        "movq    %%rcx, %[carry]\n"
    : [carry]"=g"(carry), [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x has 6.5 words, z has 8.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_nc */
static inline
void mpfq_fixmp_6_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rdx, 56(%%rdi)\n"
    : [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x has 6.5 words, z has 7.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_shortz */
static inline
mp_limb_t mpfq_fixmp_6_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %[carry]\n"
    : [carry]"=g"(carry), [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x and y have 6.5 words, z has 13. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul */
static inline
void mpfq_fixmp_6_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "### x*y[0]\n"
        "movq    %2, %%r8\n"
        "movq    %0, %%rdi\n"
        "movq    0(%%r8), %%r9\n"
        "movq    %1, %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "movq    %%rcx, 32(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "movq    %%rcx, 40(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 48(%%rdi)\n"
        "movq    %%rdx, 56(%%rdi)\n"
        "movq    $0, 64(%%rdi)\n"
        "movq    $0, 72(%%rdi)\n"
        "movq    $0, 80(%%rdi)\n"
        "movq    $0, 88(%%rdi)\n"
        "movq    $0, 96(%%rdi)\n"
        "### x*y[1]\n"
        "movq    8(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 8(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 64(%%rdi)\n"
        "### x*y[2]\n"
        "movq    16(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 16(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 72(%%rdi)\n"
        "### x*y[3]\n"
        "movq    24(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 24(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 80(%%rdi)\n"
        "### x*y[4]\n"
        "movq    32(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 32(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 88(%%rdi)\n"
        "### x*y[5]\n"
        "movq    40(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 40(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 88(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 96(%%rdi)\n"
        "### x*y[6]\n"
        "movq    48(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 48(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 88(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 96(%%rdi)\n"
      : "+m" (z)
      : "m" (x), "m" (y)
      : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", "%r9", "memory");
}

/* x has 6.5 words, z has 13. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_6_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mpfq_fixmp_6_5_mul(z, x, x);
}

/* x has 6.5 words, z has 8. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_6_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "movq    %%rcx, 32(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "movq    %%rcx, 40(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 48(%%rdi)\n"
        "movq    %%rdx, 56(%%rdi)\n"
    : 
    : [mult] "r" (c), [z] "m" (z), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x has 6.5 words, z has 7. c is 0.5 word.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul05 */
static inline
mp_limb_t mpfq_fixmp_6_5_addmul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "imulq    %[mult], %%rax\n"
        "xorq    %%rdx, %%rdx\n"
        "addq    %%rcx, %%rax\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rax, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %[carry]\n"
    : [carry]"=g"(carry), [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x has 6.5 words, z has 7. c is 0.5 word.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul05_nc */
static inline
void mpfq_fixmp_6_5_addmul05_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "imulq    %[mult], %%rax\n"
        "addq    %%rcx, %%rax\n"
        "addq    %%rax, 48(%%rdi)\n"
    : [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x has 6.5 words, z has 7. c is 0.5 word.
 * Put (x*c) in z. No carry. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul05 */
static inline
void mpfq_fixmp_6_5_mul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "movq    %%rcx, 32(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "movq    %%rcx, 40(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "imulq    %[mult], %%rax\n"
        "addq    %%rcx, %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
    : 
    : [mult] "r" (c), [z] "m" (z), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x, y, and z have 7.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add */
static inline
mp_limb_t mpfq_fixmp_7_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "addq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "adcq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "adcq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "adcq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "adcq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq    48(%%rsi), %%rax\n"
        "adcq    48(%%rdx), %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
        "movq    56(%%rsi), %%rax\n"
        "adcq    56(%%rdx), %%rax\n"
        "movq    %%rax, 56(%%rdi)\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : "=&a"(carry)
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x, y, and z have 7.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub */
static inline
mp_limb_t mpfq_fixmp_7_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "subq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "sbbq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "sbbq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "sbbq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "sbbq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "sbbq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq    48(%%rsi), %%rax\n"
        "sbbq    48(%%rdx), %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
        "movq    56(%%rsi), %%rax\n"
        "sbbq    56(%%rdx), %%rax\n"
        "movq    %%rax, 56(%%rdi)\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : "=&a"(carry)
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x, y, and z have 7.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_7_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "addq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "adcq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "adcq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "adcq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "adcq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq    48(%%rsi), %%rax\n"
        "adcq    48(%%rdx), %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
        "movq    56(%%rsi), %%rax\n"
        "adcq    56(%%rdx), %%rax\n"
        "movq    %%rax, 56(%%rdi)\n"
    : 
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory", "%rax");
}

/* x, y, and z have 7.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_7_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "subq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "sbbq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "sbbq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "sbbq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "sbbq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "sbbq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq    48(%%rsi), %%rax\n"
        "sbbq    48(%%rdx), %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
        "movq    56(%%rsi), %%rax\n"
        "sbbq    56(%%rdx), %%rax\n"
        "movq    %%rax, 56(%%rdi)\n"
    : 
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory", "%rax");
}

/* x has 7.5 words, z has 9.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_7_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "xorq    %%rcx, %%rcx\n"
        "addq    %%rdx, 64(%%rdi)\n"
        "adcq    $0, %%rcx\n"
        "movq    %%rcx, %[carry]\n"
    : [carry]"=g"(carry), [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x has 7.5 words, z has 9.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_nc */
static inline
void mpfq_fixmp_7_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rdx, 64(%%rdi)\n"
    : [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x has 7.5 words, z has 8.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_shortz */
static inline
mp_limb_t mpfq_fixmp_7_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %[carry]\n"
    : [carry]"=g"(carry), [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x and y have 7.5 words, z has 15. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul */
static inline
void mpfq_fixmp_7_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "### x*y[0]\n"
        "movq    %2, %%r8\n"
        "movq    %0, %%rdi\n"
        "movq    0(%%r8), %%r9\n"
        "movq    %1, %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "movq    %%rcx, 32(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "movq    %%rcx, 40(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "movq    %%rcx, 48(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 56(%%rdi)\n"
        "movq    %%rdx, 64(%%rdi)\n"
        "movq    $0, 72(%%rdi)\n"
        "movq    $0, 80(%%rdi)\n"
        "movq    $0, 88(%%rdi)\n"
        "movq    $0, 96(%%rdi)\n"
        "movq    $0, 104(%%rdi)\n"
        "movq    $0, 112(%%rdi)\n"
        "### x*y[1]\n"
        "movq    8(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 8(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 72(%%rdi)\n"
        "### x*y[2]\n"
        "movq    16(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 16(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 80(%%rdi)\n"
        "### x*y[3]\n"
        "movq    24(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 24(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 88(%%rdi)\n"
        "### x*y[4]\n"
        "movq    32(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 32(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 88(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 96(%%rdi)\n"
        "### x*y[5]\n"
        "movq    40(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 40(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 88(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 96(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 104(%%rdi)\n"
        "### x*y[6]\n"
        "movq    48(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 48(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 88(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 96(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 104(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 112(%%rdi)\n"
        "### x*y[7]\n"
        "movq    56(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 56(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 88(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 96(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 104(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 112(%%rdi)\n"
      : "+m" (z)
      : "m" (x), "m" (y)
      : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", "%r9", "memory");
}

/* x has 7.5 words, z has 15. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_7_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mpfq_fixmp_7_5_mul(z, x, x);
}

/* x has 7.5 words, z has 9. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_7_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "movq    %%rcx, 32(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "movq    %%rcx, 40(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "movq    %%rcx, 48(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 56(%%rdi)\n"
        "movq    %%rdx, 64(%%rdi)\n"
    : 
    : [mult] "r" (c), [z] "m" (z), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x has 7.5 words, z has 8. c is 0.5 word.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul05 */
static inline
mp_limb_t mpfq_fixmp_7_5_addmul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "imulq    %[mult], %%rax\n"
        "xorq    %%rdx, %%rdx\n"
        "addq    %%rcx, %%rax\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rax, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %[carry]\n"
    : [carry]"=g"(carry), [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x has 7.5 words, z has 8. c is 0.5 word.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul05_nc */
static inline
void mpfq_fixmp_7_5_addmul05_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "imulq    %[mult], %%rax\n"
        "addq    %%rcx, %%rax\n"
        "addq    %%rax, 56(%%rdi)\n"
    : [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x has 7.5 words, z has 8. c is 0.5 word.
 * Put (x*c) in z. No carry. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul05 */
static inline
void mpfq_fixmp_7_5_mul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "movq    %%rcx, 32(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "movq    %%rcx, 40(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "movq    %%rcx, 48(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "imulq    %[mult], %%rax\n"
        "addq    %%rcx, %%rax\n"
        "movq    %%rax, 56(%%rdi)\n"
    : 
    : [mult] "r" (c), [z] "m" (z), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x, y, and z have 8.5 words. Result in z. Return carry bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add */
static inline
mp_limb_t mpfq_fixmp_8_5_add(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "addq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "adcq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "adcq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "adcq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "adcq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq    48(%%rsi), %%rax\n"
        "adcq    48(%%rdx), %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
        "movq    56(%%rsi), %%rax\n"
        "adcq    56(%%rdx), %%rax\n"
        "movq    %%rax, 56(%%rdi)\n"
        "movq    64(%%rsi), %%rax\n"
        "adcq    64(%%rdx), %%rax\n"
        "movq    %%rax, 64(%%rdi)\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : "=&a"(carry)
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x, y, and z have 8.5 words. Result in z. Return borrow bit */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub */
static inline
mp_limb_t mpfq_fixmp_8_5_sub(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "subq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "sbbq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "sbbq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "sbbq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "sbbq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "sbbq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq    48(%%rsi), %%rax\n"
        "sbbq    48(%%rdx), %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
        "movq    56(%%rsi), %%rax\n"
        "sbbq    56(%%rdx), %%rax\n"
        "movq    %%rax, 56(%%rdi)\n"
        "movq    64(%%rsi), %%rax\n"
        "sbbq    64(%%rdx), %%rax\n"
        "movq    %%rax, 64(%%rdi)\n"
        "movq $0, %%rax\n"
        "adcq $0, %%rax\n"
    : "=&a"(carry)
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x, y, and z have 8.5 words. Result in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_add_nc */
static inline
void mpfq_fixmp_8_5_add_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "addq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "adcq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "adcq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "adcq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "adcq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq    48(%%rsi), %%rax\n"
        "adcq    48(%%rdx), %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
        "movq    56(%%rsi), %%rax\n"
        "adcq    56(%%rdx), %%rax\n"
        "movq    %%rax, 56(%%rdi)\n"
        "movq    64(%%rsi), %%rax\n"
        "adcq    64(%%rdx), %%rax\n"
        "movq    %%rax, 64(%%rdi)\n"
    : 
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory", "%rax");
}

/* x, y, and z have 8.5 words. Result in z. Borrow bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sub_nc */
static inline
void mpfq_fixmp_8_5_sub_nc(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    %[y], %%rdx\n"
        "movq    0(%%rsi), %%rax\n"
        "subq    0(%%rdx), %%rax\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "sbbq    8(%%rdx), %%rax\n"
        "movq    %%rax, 8(%%rdi)\n"
        "movq    16(%%rsi), %%rax\n"
        "sbbq    16(%%rdx), %%rax\n"
        "movq    %%rax, 16(%%rdi)\n"
        "movq    24(%%rsi), %%rax\n"
        "sbbq    24(%%rdx), %%rax\n"
        "movq    %%rax, 24(%%rdi)\n"
        "movq    32(%%rsi), %%rax\n"
        "sbbq    32(%%rdx), %%rax\n"
        "movq    %%rax, 32(%%rdi)\n"
        "movq    40(%%rsi), %%rax\n"
        "sbbq    40(%%rdx), %%rax\n"
        "movq    %%rax, 40(%%rdi)\n"
        "movq    48(%%rsi), %%rax\n"
        "sbbq    48(%%rdx), %%rax\n"
        "movq    %%rax, 48(%%rdi)\n"
        "movq    56(%%rsi), %%rax\n"
        "sbbq    56(%%rdx), %%rax\n"
        "movq    %%rax, 56(%%rdi)\n"
        "movq    64(%%rsi), %%rax\n"
        "sbbq    64(%%rdx), %%rax\n"
        "movq    %%rax, 64(%%rdi)\n"
    : 
    : [z]"m"(z), [x]"m"(x), [y]"m"(y)
    : "%rdx", "%rsi", "%rdi", "memory", "%rax");
}

/* x has 8.5 words, z has 10.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1 */
static inline
mp_limb_t mpfq_fixmp_8_5_addmul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "xorq    %%rcx, %%rcx\n"
        "addq    %%rdx, 72(%%rdi)\n"
        "adcq    $0, %%rcx\n"
        "movq    %%rcx, %[carry]\n"
    : [carry]"=g"(carry), [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x has 8.5 words, z has 10.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_nc */
static inline
void mpfq_fixmp_8_5_addmul1_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rdx, 72(%%rdi)\n"
    : [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x has 8.5 words, z has 9.
 * Put (z+x*c) in z. Return carry word. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul1_shortz */
static inline
mp_limb_t mpfq_fixmp_8_5_addmul1_shortz(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %[carry]\n"
    : [carry]"=g"(carry), [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x and y have 8.5 words, z has 17. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul */
static inline
void mpfq_fixmp_8_5_mul(mp_limb_t * z, const mp_limb_t * x, const mp_limb_t * y)
{
    __asm__ __volatile__(
        "### x*y[0]\n"
        "movq    %2, %%r8\n"
        "movq    %0, %%rdi\n"
        "movq    0(%%r8), %%r9\n"
        "movq    %1, %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "movq    %%rcx, 32(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "movq    %%rcx, 40(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "movq    %%rcx, 48(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "movq    %%rcx, 56(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 64(%%rdi)\n"
        "movq    %%rdx, 72(%%rdi)\n"
        "movq    $0, 80(%%rdi)\n"
        "movq    $0, 88(%%rdi)\n"
        "movq    $0, 96(%%rdi)\n"
        "movq    $0, 104(%%rdi)\n"
        "movq    $0, 112(%%rdi)\n"
        "movq    $0, 120(%%rdi)\n"
        "movq    $0, 128(%%rdi)\n"
        "### x*y[1]\n"
        "movq    8(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 8(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 80(%%rdi)\n"
        "### x*y[2]\n"
        "movq    16(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 16(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 88(%%rdi)\n"
        "### x*y[3]\n"
        "movq    24(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 24(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 88(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 96(%%rdi)\n"
        "### x*y[4]\n"
        "movq    32(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 32(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "addq    %%rcx, 88(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 96(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 104(%%rdi)\n"
        "### x*y[5]\n"
        "movq    40(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 40(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 88(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "addq    %%rcx, 96(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 104(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 112(%%rdi)\n"
        "### x*y[6]\n"
        "movq    48(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 48(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 88(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 96(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "addq    %%rcx, 104(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 112(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 120(%%rdi)\n"
        "### x*y[7]\n"
        "movq    56(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 56(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 88(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 96(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 104(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "addq    %%rcx, 112(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 120(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, 128(%%rdi)\n"
        "### x*y[8]\n"
        "movq    64(%%r8), %%r9\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %%r9\n"
        "addq    %%rax, 64(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 72(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 80(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 88(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 96(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 104(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 112(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "addq    %%rcx, 120(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %%r9\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rcx, 128(%%rdi)\n"
      : "+m" (z)
      : "m" (x), "m" (y)
      : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "%r8", "%r9", "memory");
}

/* x has 8.5 words, z has 17. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_sqr */
static inline
void mpfq_fixmp_8_5_sqr(mp_limb_t * z, const mp_limb_t * x)
{
    mpfq_fixmp_8_5_mul(z, x, x);
}

/* x has 8.5 words, z has 10. Put x*y in z. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul1 */
static inline
void mpfq_fixmp_8_5_mul1(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "movq    %%rcx, 32(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "movq    %%rcx, 40(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "movq    %%rcx, 48(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "movq    %%rcx, 56(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rcx, 64(%%rdi)\n"
        "movq    %%rdx, 72(%%rdi)\n"
    : 
    : [mult] "r" (c), [z] "m" (z), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x has 8.5 words, z has 9. c is 0.5 word.
 * Put (z+x*c) in z. Return carry bit. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul05 */
static inline
mp_limb_t mpfq_fixmp_8_5_addmul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    mp_limb_t carry;
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "imulq    %[mult], %%rax\n"
        "xorq    %%rdx, %%rdx\n"
        "addq    %%rcx, %%rax\n"
        "adcq    $0, %%rdx\n"
        "addq    %%rax, 64(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %[carry]\n"
    : [carry]"=g"(carry), [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
    return carry;
}

/* x has 8.5 words, z has 9. c is 0.5 word.
 * Put (z+x*c) in z. Carry bit is lost. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_addmul05_nc */
static inline
void mpfq_fixmp_8_5_addmul05_nc(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "addq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "addq    %%rcx, 8(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "addq    %%rcx, 16(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "addq    %%rcx, 24(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "addq    %%rcx, 32(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "addq    %%rcx, 40(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "addq    %%rcx, 48(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "addq    %%rcx, 56(%%rdi)\n"
        "adcq    $0, %%rdx\n"
        "movq    %%rdx, %%rcx\n"
        "imulq    %[mult], %%rax\n"
        "addq    %%rcx, %%rax\n"
        "addq    %%rax, 64(%%rdi)\n"
    : [z] "+m" (z)
    : [mult] "r" (c), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}

/* x has 8.5 words, z has 9. c is 0.5 word.
 * Put (x*c) in z. No carry. */
/* *Mpfq::fixmp::x86_64::code_for__fixmp_mul05 */
static inline
void mpfq_fixmp_8_5_mul05(mp_limb_t * z, const mp_limb_t * x, mp_limb_t c)
{
    __asm__ __volatile__(
        "movq    %[z], %%rdi\n"
        "movq    %[x], %%rsi\n"
        "movq    0(%%rsi), %%rax\n"
        "mulq    %[mult]\n"
        "movq    %%rax, 0(%%rdi)\n"
        "movq    8(%%rsi), %%rax\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    16(%%rsi), %%rax\n"
        "movq    %%rcx, 8(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    24(%%rsi), %%rax\n"
        "movq    %%rcx, 16(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    32(%%rsi), %%rax\n"
        "movq    %%rcx, 24(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    40(%%rsi), %%rax\n"
        "movq    %%rcx, 32(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    48(%%rsi), %%rax\n"
        "movq    %%rcx, 40(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    56(%%rsi), %%rax\n"
        "movq    %%rcx, 48(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "mulq    %[mult]\n"
        "addq    %%rax, %%rcx\n"
        "adcq    $0, %%rdx\n"
        "movq    64(%%rsi), %%rax\n"
        "movq    %%rcx, 56(%%rdi)\n"
        "movq    %%rdx, %%rcx\n"
        "imulq    %[mult], %%rax\n"
        "addq    %%rcx, %%rax\n"
        "movq    %%rax, 64(%%rdi)\n"
    : 
    : [mult] "r" (c), [z] "m" (z), [x] "m" (x)
    : "%rax", "%rcx", "%rdx", "%rsi", "%rdi", "memory");
}


#endif  /* MPFQ_FIXMP_X86_64_H_ */

/* vim:set ft=cpp: */
