#ifndef ABASE_PZ_H_
#define ABASE_PZ_H_

/* MPFQ generated file -- do not edit */

#include "mpfq.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <ctype.h>
#include "mpfq_gfp_common.h"
#include <stddef.h>
#include <stdio.h>
#include "assert.h"
#include "select_mpi.h"
#include "abase_vbase.h"
#ifdef	MPFQ_LAST_GENERATED_TAG
#undef	MPFQ_LAST_GENERATED_TAG
#endif
#define MPFQ_LAST_GENERATED_TAG      pz

/* Active handler: simd_pz */
/* Automatically generated code  */
/* Active handler: pz */
/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::poly */
/* Active handler: Mpfq::defaults::mpi_flat */
/* Options used:{
   family=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=pz, }, ],
   fieldtype=prime,
   tag=pz,
   type=plain,
   vbase_stuff={
    choose_byfeatures=<code>,
    families=[
     [ u64k1, u64k2, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_1, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_2, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_3, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_4, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_8, }, ],
     [ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=pz, }, ],
     ],
    member_templates_restrict={
     p_1=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_1, }, ],
     p_2=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_2, }, ],
     p_3=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_3, }, ],
     p_4=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_4, }, ],
     p_8=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_8, }, ],
     pz=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=pz, }, ],
     u64k1=[ u64k1, u64k2, ],
     u64k2=[ u64k1, u64k2, ],
     },
    vc:includes=[ <stdarg.h>, ],
    },
   virtual_base={
    filebase=abase_vbase,
    global_prefix=abase_,
    name=abase_vbase,
    substitutions=[
     [ (?^:abase_pz_elt \*), void *, ],
     [ (?^:abase_pz_src_elt\b), const void *, ],
     [ (?^:abase_pz_elt\b), void *, ],
     [ (?^:abase_pz_dst_elt\b), void *, ],
     [ (?^:abase_pz_elt_ur \*), void *, ],
     [ (?^:abase_pz_src_elt_ur\b), const void *, ],
     [ (?^:abase_pz_elt_ur\b), void *, ],
     [ (?^:abase_pz_dst_elt_ur\b), void *, ],
     [ (?^:abase_pz_vec \*), void *, ],
     [ (?^:abase_pz_src_vec\b), const void *, ],
     [ (?^:abase_pz_vec\b), void *, ],
     [ (?^:abase_pz_dst_vec\b), void *, ],
     [ (?^:abase_pz_vec_ur \*), void *, ],
     [ (?^:abase_pz_src_vec_ur\b), const void *, ],
     [ (?^:abase_pz_vec_ur\b), void *, ],
     [ (?^:abase_pz_dst_vec_ur\b), void *, ],
     [ (?^:abase_pz_poly \*), void *, ],
     [ (?^:abase_pz_src_poly\b), const void *, ],
     [ (?^:abase_pz_poly\b), void *, ],
     [ (?^:abase_pz_dst_poly\b), void *, ],
     ],
    },
   vtag=pz,
   w=64,
   } */

typedef mpfq_p_field abase_pz_field;
typedef mpfq_p_dst_field abase_pz_dst_field;

typedef mp_limb_t * abase_pz_elt;
typedef mp_limb_t * abase_pz_dst_elt;
typedef const mp_limb_t * abase_pz_src_elt;

typedef mp_limb_t * abase_pz_elt_ur;
typedef mp_limb_t * abase_pz_dst_elt_ur;
typedef const mp_limb_t * abase_pz_src_elt_ur;

typedef mp_limb_t * abase_pz_vec;
typedef mp_limb_t * abase_pz_dst_vec;
typedef const mp_limb_t * abase_pz_src_vec;

typedef mp_limb_t * abase_pz_vec_ur;
typedef mp_limb_t * abase_pz_dst_vec_ur;
typedef const mp_limb_t * abase_pz_src_vec_ur;

typedef struct {
  abase_pz_vec c;
  unsigned int alloc;
  unsigned int size;
} abase_pz_poly_struct;
typedef abase_pz_poly_struct abase_pz_poly [1];
typedef abase_pz_poly_struct * abase_pz_dst_poly;
typedef abase_pz_poly_struct * abase_pz_src_poly;

#ifdef  __cplusplus
extern "C" {
#endif
/* *Mpfq::defaults::code_for_impl_name, pz */
#define abase_pz_impl_name()	"pz"
/* *pz::code_for_impl_max_characteristic_bits */
#define abase_pz_impl_max_characteristic_bits()	UINT_MAX
/* *pz::code_for_impl_max_degree */
#define abase_pz_impl_max_degree()	1

/* Functions operating on the field structure */
void abase_pz_field_characteristic(abase_pz_dst_field, mpz_t);
/* *pz::code_for_field_degree */
#define abase_pz_field_degree(k)	1
void abase_pz_field_init(abase_pz_dst_field);
void abase_pz_field_clear(abase_pz_dst_field);
void abase_pz_field_specify(abase_pz_dst_field, unsigned long, void *);
/* *pz::code_for_field_setopt */
#define abase_pz_field_setopt(k, x, y)	/**/

/* Element allocation functions */
void abase_pz_init(abase_pz_dst_field, abase_pz_elt *);
void abase_pz_clear(abase_pz_dst_field, abase_pz_elt *);

/* Elementary assignment functions */
static inline
void abase_pz_set(abase_pz_dst_field, abase_pz_dst_elt, abase_pz_src_elt);
static inline
void abase_pz_set_ui(abase_pz_dst_field, abase_pz_dst_elt, unsigned long);
static inline
void abase_pz_set_zero(abase_pz_dst_field, abase_pz_dst_elt);
static inline
unsigned long abase_pz_get_ui(abase_pz_dst_field, abase_pz_src_elt);
void abase_pz_set_mpn(abase_pz_dst_field, abase_pz_dst_elt, mp_limb_t *, size_t);
void abase_pz_set_mpz(abase_pz_dst_field, abase_pz_dst_elt, mpz_t);
static inline
void abase_pz_get_mpn(abase_pz_dst_field, mp_limb_t *, abase_pz_src_elt);
static inline
void abase_pz_get_mpz(abase_pz_dst_field, mpz_t, abase_pz_src_elt);

/* Assignment of random values */
void abase_pz_random(abase_pz_dst_field, abase_pz_dst_elt, gmp_randstate_t);
void abase_pz_random2(abase_pz_dst_field, abase_pz_dst_elt, gmp_randstate_t);

/* Arithmetic operations on elements */
void abase_pz_add(abase_pz_dst_field, abase_pz_dst_elt, abase_pz_src_elt, abase_pz_src_elt);
void abase_pz_sub(abase_pz_dst_field, abase_pz_dst_elt, abase_pz_src_elt, abase_pz_src_elt);
void abase_pz_neg(abase_pz_dst_field, abase_pz_dst_elt, abase_pz_src_elt);
void abase_pz_mul(abase_pz_dst_field, abase_pz_dst_elt, abase_pz_src_elt, abase_pz_src_elt);
void abase_pz_sqr(abase_pz_dst_field, abase_pz_dst_elt, abase_pz_src_elt);
int abase_pz_is_sqr(abase_pz_dst_field, abase_pz_src_elt);
/* missing sqrt */
void abase_pz_pow(abase_pz_dst_field, abase_pz_dst_elt, abase_pz_src_elt, unsigned long *, size_t);
/* *pz::code_for_frobenius */
#define abase_pz_frobenius(k, x, y)	abase_pz_set(k,x,y)
void abase_pz_add_ui(abase_pz_dst_field, abase_pz_dst_elt, abase_pz_src_elt, unsigned long);
void abase_pz_sub_ui(abase_pz_dst_field, abase_pz_dst_elt, abase_pz_src_elt, unsigned long);
void abase_pz_mul_ui(abase_pz_dst_field, abase_pz_dst_elt, abase_pz_src_elt, unsigned long);
int abase_pz_inv(abase_pz_dst_field, abase_pz_dst_elt, abase_pz_src_elt);

/* Operations involving unreduced elements */
void abase_pz_elt_ur_init(abase_pz_dst_field, abase_pz_elt_ur *);
void abase_pz_elt_ur_clear(abase_pz_dst_field, abase_pz_elt_ur *);
static inline
void abase_pz_elt_ur_set(abase_pz_dst_field, abase_pz_dst_elt_ur, abase_pz_src_elt_ur);
static inline
void abase_pz_elt_ur_set_elt(abase_pz_dst_field, abase_pz_dst_elt_ur, abase_pz_src_elt);
static inline
void abase_pz_elt_ur_set_zero(abase_pz_dst_field, abase_pz_dst_elt_ur);
static inline
void abase_pz_elt_ur_set_ui(abase_pz_dst_field, abase_pz_dst_elt_ur, unsigned long);
void abase_pz_elt_ur_add(abase_pz_dst_field, abase_pz_dst_elt_ur, abase_pz_src_elt_ur, abase_pz_src_elt_ur);
void abase_pz_elt_ur_neg(abase_pz_dst_field, abase_pz_dst_elt_ur, abase_pz_src_elt_ur);
void abase_pz_elt_ur_sub(abase_pz_dst_field, abase_pz_dst_elt_ur, abase_pz_src_elt_ur, abase_pz_src_elt_ur);
void abase_pz_mul_ur(abase_pz_dst_field, abase_pz_dst_elt_ur, abase_pz_src_elt, abase_pz_src_elt);
void abase_pz_sqr_ur(abase_pz_dst_field, abase_pz_dst_elt_ur, abase_pz_src_elt);
void abase_pz_reduce(abase_pz_dst_field, abase_pz_dst_elt, abase_pz_dst_elt_ur);
#define HAVE_abase_pz_normalize
void abase_pz_normalize(abase_pz_dst_field, abase_pz_dst_elt);
#define HAVE_abase_pz_addmul_si_ur
static inline
void abase_pz_addmul_si_ur(abase_pz_dst_field, abase_pz_dst_elt_ur, abase_pz_src_elt, long);

/* Comparison functions */
static inline
int abase_pz_cmp(abase_pz_dst_field, abase_pz_src_elt, abase_pz_src_elt);
static inline
int abase_pz_cmp_ui(abase_pz_dst_field, abase_pz_src_elt, unsigned long);
static inline
int abase_pz_is_zero(abase_pz_dst_field, abase_pz_src_elt);

/* Input/output functions */
void abase_pz_asprint(abase_pz_dst_field, char * *, abase_pz_src_elt);
void abase_pz_fprint(abase_pz_dst_field, FILE *, abase_pz_src_elt);
/* *Mpfq::defaults::code_for_print, pz */
#define abase_pz_print(k, x)	abase_pz_fprint(k,stdout,x)
int abase_pz_sscan(abase_pz_dst_field, abase_pz_dst_elt, const char *);
int abase_pz_fscan(abase_pz_dst_field, FILE *, abase_pz_dst_elt);
/* *Mpfq::defaults::code_for_scan, pz */
#define abase_pz_scan(k, x)	abase_pz_fscan(k,stdin,x)

/* Vector functions */
void abase_pz_vec_init(abase_pz_dst_field, abase_pz_vec *, unsigned int);
void abase_pz_vec_reinit(abase_pz_dst_field, abase_pz_vec *, unsigned int, unsigned int);
void abase_pz_vec_clear(abase_pz_dst_field, abase_pz_vec *, unsigned int);
void abase_pz_vec_set(abase_pz_dst_field, abase_pz_dst_vec, abase_pz_src_vec, unsigned int);
void abase_pz_vec_set_partial(abase_pz_dst_field, abase_pz_dst_vec, abase_pz_src_vec, unsigned int, unsigned int, unsigned int);
void abase_pz_vec_set_zero(abase_pz_dst_field, abase_pz_dst_vec, unsigned int);
void abase_pz_vec_setcoef(abase_pz_dst_field, abase_pz_dst_vec, abase_pz_src_elt, unsigned int);
void abase_pz_vec_setcoef_ui(abase_pz_dst_field, abase_pz_dst_vec, unsigned long, unsigned int);
void abase_pz_vec_getcoef(abase_pz_dst_field, abase_pz_dst_elt, abase_pz_src_vec, unsigned int);
void abase_pz_vec_add(abase_pz_dst_field, abase_pz_dst_vec, abase_pz_src_vec, abase_pz_src_vec, unsigned int);
void abase_pz_vec_neg(abase_pz_dst_field, abase_pz_dst_vec, abase_pz_src_vec, unsigned int);
void abase_pz_vec_rev(abase_pz_dst_field, abase_pz_dst_vec, abase_pz_src_vec, unsigned int);
void abase_pz_vec_sub(abase_pz_dst_field, abase_pz_dst_vec, abase_pz_src_vec, abase_pz_src_vec, unsigned int);
void abase_pz_vec_scal_mul(abase_pz_dst_field, abase_pz_dst_vec, abase_pz_src_vec, abase_pz_src_elt, unsigned int);
void abase_pz_vec_conv(abase_pz_dst_field, abase_pz_dst_vec, abase_pz_src_vec, unsigned int, abase_pz_src_vec, unsigned int);
void abase_pz_vec_random(abase_pz_dst_field, abase_pz_dst_vec, unsigned int, gmp_randstate_t);
void abase_pz_vec_random2(abase_pz_dst_field, abase_pz_dst_vec, unsigned int, gmp_randstate_t);
int abase_pz_vec_cmp(abase_pz_dst_field, abase_pz_src_vec, abase_pz_src_vec, unsigned int);
int abase_pz_vec_is_zero(abase_pz_dst_field, abase_pz_src_vec, unsigned int);
static inline
abase_pz_dst_vec abase_pz_vec_subvec(abase_pz_dst_field, abase_pz_dst_vec, int);
static inline
abase_pz_src_vec abase_pz_vec_subvec_const(abase_pz_dst_field, abase_pz_src_vec, int);
static inline
abase_pz_dst_elt abase_pz_vec_coeff_ptr(abase_pz_dst_field, abase_pz_dst_vec, int);
static inline
abase_pz_src_elt abase_pz_vec_coeff_ptr_const(abase_pz_dst_field, abase_pz_src_vec, int);
void abase_pz_vec_asprint(abase_pz_dst_field, char * *, abase_pz_src_vec, unsigned int);
void abase_pz_vec_fprint(abase_pz_dst_field, FILE *, abase_pz_src_vec, unsigned int);
void abase_pz_vec_print(abase_pz_dst_field, abase_pz_src_vec, unsigned int);
int abase_pz_vec_sscan(abase_pz_dst_field, abase_pz_vec *, unsigned int *, const char *);
int abase_pz_vec_fscan(abase_pz_dst_field, FILE *, abase_pz_vec *, unsigned int *);
/* missing vec_scan */
void abase_pz_vec_ur_init(abase_pz_dst_field, abase_pz_vec_ur *, unsigned int);
void abase_pz_vec_ur_set_zero(abase_pz_dst_field, abase_pz_dst_vec_ur, unsigned int);
void abase_pz_vec_ur_set_vec(abase_pz_dst_field, abase_pz_dst_vec_ur, abase_pz_src_vec, unsigned int);
void abase_pz_vec_ur_reinit(abase_pz_dst_field, abase_pz_vec_ur *, unsigned int, unsigned int);
void abase_pz_vec_ur_clear(abase_pz_dst_field, abase_pz_vec_ur *, unsigned int);
void abase_pz_vec_ur_set(abase_pz_dst_field, abase_pz_dst_vec_ur, abase_pz_src_vec_ur, unsigned int);
void abase_pz_vec_ur_setcoef(abase_pz_dst_field, abase_pz_dst_vec_ur, abase_pz_src_elt_ur, unsigned int);
void abase_pz_vec_ur_getcoef(abase_pz_dst_field, abase_pz_dst_elt_ur, abase_pz_src_vec_ur, unsigned int);
void abase_pz_vec_ur_add(abase_pz_dst_field, abase_pz_dst_vec_ur, abase_pz_src_vec_ur, abase_pz_src_vec_ur, unsigned int);
void abase_pz_vec_ur_sub(abase_pz_dst_field, abase_pz_dst_vec_ur, abase_pz_src_vec_ur, abase_pz_src_vec_ur, unsigned int);
void abase_pz_vec_ur_neg(abase_pz_dst_field, abase_pz_dst_vec_ur, abase_pz_src_vec_ur, unsigned int);
void abase_pz_vec_ur_rev(abase_pz_dst_field, abase_pz_dst_vec_ur, abase_pz_src_vec_ur, unsigned int);
void abase_pz_vec_scal_mul_ur(abase_pz_dst_field, abase_pz_dst_vec_ur, abase_pz_src_vec, abase_pz_src_elt, unsigned int);
void abase_pz_vec_conv_ur(abase_pz_dst_field, abase_pz_dst_vec_ur, abase_pz_src_vec, unsigned int, abase_pz_src_vec, unsigned int);
void abase_pz_vec_reduce(abase_pz_dst_field, abase_pz_dst_vec, abase_pz_dst_vec_ur, unsigned int);
static inline
abase_pz_dst_vec_ur abase_pz_vec_ur_subvec(abase_pz_dst_field, abase_pz_dst_vec_ur, int);
static inline
abase_pz_src_vec_ur abase_pz_vec_ur_subvec_const(abase_pz_dst_field, abase_pz_src_vec_ur, int);
static inline
abase_pz_dst_elt abase_pz_vec_ur_coeff_ptr(abase_pz_dst_field, abase_pz_dst_vec_ur, int);
static inline
abase_pz_src_elt abase_pz_vec_ur_coeff_ptr_const(abase_pz_dst_field, abase_pz_src_vec_ur, int);
ptrdiff_t abase_pz_vec_elt_stride(abase_pz_dst_field, int);

/* Polynomial functions */
static inline
void abase_pz_poly_init(abase_pz_dst_field, abase_pz_poly, unsigned int);
static inline
void abase_pz_poly_clear(abase_pz_dst_field, abase_pz_poly);
static inline
void abase_pz_poly_set(abase_pz_dst_field, abase_pz_dst_poly, abase_pz_src_poly);
void abase_pz_poly_setmonic(abase_pz_dst_field, abase_pz_dst_poly, abase_pz_src_poly);
static inline
void abase_pz_poly_setcoef(abase_pz_dst_field, abase_pz_dst_poly, abase_pz_src_elt, unsigned int);
static inline
void abase_pz_poly_setcoef_ui(abase_pz_dst_field, abase_pz_dst_poly, unsigned long, unsigned int);
static inline
void abase_pz_poly_getcoef(abase_pz_dst_field, abase_pz_dst_elt, abase_pz_src_poly, unsigned int);
static inline
int abase_pz_poly_deg(abase_pz_dst_field, abase_pz_src_poly);
static inline
void abase_pz_poly_add(abase_pz_dst_field, abase_pz_dst_poly, abase_pz_src_poly, abase_pz_src_poly);
static inline
void abase_pz_poly_sub(abase_pz_dst_field, abase_pz_dst_poly, abase_pz_src_poly, abase_pz_src_poly);
static inline
void abase_pz_poly_add_ui(abase_pz_dst_field, abase_pz_dst_poly, abase_pz_src_poly, unsigned long);
static inline
void abase_pz_poly_sub_ui(abase_pz_dst_field, abase_pz_dst_poly, abase_pz_src_poly, unsigned long);
static inline
void abase_pz_poly_neg(abase_pz_dst_field, abase_pz_dst_poly, abase_pz_src_poly);
static inline
void abase_pz_poly_scal_mul(abase_pz_dst_field, abase_pz_dst_poly, abase_pz_src_poly, abase_pz_src_elt);
static inline
void abase_pz_poly_mul(abase_pz_dst_field, abase_pz_dst_poly, abase_pz_src_poly, abase_pz_src_poly);
void abase_pz_poly_divmod(abase_pz_dst_field, abase_pz_dst_poly, abase_pz_dst_poly, abase_pz_src_poly, abase_pz_src_poly);
void abase_pz_poly_precomp_mod(abase_pz_dst_field, abase_pz_dst_poly, abase_pz_src_poly);
void abase_pz_poly_mod_pre(abase_pz_dst_field, abase_pz_dst_poly, abase_pz_src_poly, abase_pz_src_poly, abase_pz_src_poly);
static inline
void abase_pz_poly_gcd(abase_pz_dst_field, abase_pz_dst_poly, abase_pz_src_poly, abase_pz_src_poly);
static inline
void abase_pz_poly_xgcd(abase_pz_dst_field, abase_pz_dst_poly, abase_pz_dst_poly, abase_pz_dst_poly, abase_pz_src_poly, abase_pz_src_poly);
static inline
void abase_pz_poly_random(abase_pz_dst_field, abase_pz_dst_poly, unsigned int, gmp_randstate_t);
static inline
void abase_pz_poly_random2(abase_pz_dst_field, abase_pz_dst_poly, unsigned int, gmp_randstate_t);
static inline
int abase_pz_poly_cmp(abase_pz_dst_field, abase_pz_src_poly, abase_pz_src_poly);
static inline
void abase_pz_poly_asprint(abase_pz_dst_field, char * *, abase_pz_src_poly);
static inline
void abase_pz_poly_fprint(abase_pz_dst_field, FILE *, abase_pz_src_poly);
static inline
void abase_pz_poly_print(abase_pz_dst_field, abase_pz_src_poly);
static inline
int abase_pz_poly_sscan(abase_pz_dst_field, abase_pz_dst_poly, const char *);
static inline
int abase_pz_poly_fscan(abase_pz_dst_field, FILE *, abase_pz_dst_poly);
/* missing poly_scan */

/* Functions related to SIMD operation */
/* *simd_pz::code_for_groupsize */
#define abase_pz_groupsize(K)	1
/* *simd_pz::code_for_offset */
#define abase_pz_offset(K, n)	n /* TO BE DEPRECATED */
/* *simd_pz::code_for_stride */
#define abase_pz_stride(K)	1 /* TO BE DEPRECATED */
static inline
void abase_pz_set_ui_at(abase_pz_dst_field, abase_pz_dst_elt, int, unsigned long);
static inline
void abase_pz_set_ui_all(abase_pz_dst_field, abase_pz_dst_elt, unsigned long);
static inline
void abase_pz_elt_ur_set_ui_at(abase_pz_dst_field, abase_pz_dst_elt, int, unsigned long);
static inline
void abase_pz_elt_ur_set_ui_all(abase_pz_dst_field, abase_pz_dst_elt, unsigned long);
void abase_pz_dotprod(abase_pz_dst_field, abase_pz_dst_vec, abase_pz_src_vec, abase_pz_src_vec, unsigned int);

/* Member templates related to SIMD operation */

/* MPI interface */
void abase_pz_mpi_ops_init(abase_pz_dst_field);
MPI_Datatype abase_pz_mpi_datatype(abase_pz_dst_field);
MPI_Datatype abase_pz_mpi_datatype_ur(abase_pz_dst_field);
MPI_Op abase_pz_mpi_addition_op(abase_pz_dst_field);
MPI_Op abase_pz_mpi_addition_op_ur(abase_pz_dst_field);
void abase_pz_mpi_ops_clear(abase_pz_dst_field);

/* Object-oriented interface */
static inline
void abase_pz_oo_field_clear(abase_vbase_ptr);
void abase_pz_oo_field_init(abase_vbase_ptr);
#ifdef  __cplusplus
}
#endif

/* Implementations for inlines */
/* *pz::code_for_set */
static inline
void abase_pz_set(abase_pz_dst_field k, abase_pz_dst_elt z, abase_pz_src_elt x)
{
        memcpy(z, x, k->kl * sizeof(mp_limb_t));
}

/* *pz::code_for_set_ui */
static inline
void abase_pz_set_ui(abase_pz_dst_field k, abase_pz_dst_elt z, unsigned long x0)
{
        z[0] = x0;
        memset(z + 1, 0, (k->kl - 1) * sizeof(mp_limb_t));
}

/* *pz::code_for_set_zero */
static inline
void abase_pz_set_zero(abase_pz_dst_field k, abase_pz_dst_elt z)
{
        memset(z, 0, (k->kl) * sizeof(mp_limb_t));
}

/* *pz::code_for_get_ui */
static inline
unsigned long abase_pz_get_ui(abase_pz_dst_field k MAYBE_UNUSED, abase_pz_src_elt x)
{
        return x[0];
}

/* *pz::code_for_get_mpn */
static inline
void abase_pz_get_mpn(abase_pz_dst_field k, mp_limb_t * z, abase_pz_src_elt x)
{
        memcpy(z, x, k->kl * sizeof(mp_limb_t));
}

/* *pz::code_for_get_mpz */
static inline
void abase_pz_get_mpz(abase_pz_dst_field k, mpz_t z, abase_pz_src_elt x)
{
        int n = k->kl;
        mpz_set_ui(z, x[n - 1]);
        for (int i = n - 2; i >= 0; --i) {
        mpz_mul_2exp(z, z, 64);
        mpz_add_ui(z, z, x[i]);
        }
}

/* *pz::code_for_elt_ur_set */
static inline
void abase_pz_elt_ur_set(abase_pz_dst_field k, abase_pz_dst_elt_ur z, abase_pz_src_elt_ur x)
{
        memcpy(z, x, k->url * sizeof(mp_limb_t));
}

/* *pz::code_for_elt_ur_set_elt */
static inline
void abase_pz_elt_ur_set_elt(abase_pz_dst_field k, abase_pz_dst_elt_ur z, abase_pz_src_elt x)
{
        memcpy(z, x, k->kl * sizeof(mp_limb_t));
        memset(z + k->kl, 0, (k->url - k->kl) * sizeof(mp_limb_t));
}

/* *pz::code_for_elt_ur_set_zero */
static inline
void abase_pz_elt_ur_set_zero(abase_pz_dst_field k, abase_pz_dst_elt_ur z)
{
        memset(z, 0, (k->url) * sizeof(mp_limb_t));
}

/* *pz::code_for_elt_ur_set_ui */
static inline
void abase_pz_elt_ur_set_ui(abase_pz_dst_field k, abase_pz_dst_elt_ur z, unsigned long x0)
{
        z[0] = x0;
        memset(z + 1, 0, (k->url - 1) * sizeof(mp_limb_t));
}

/* *simd_pz::code_for_addmul_si_ur */
static inline
void abase_pz_addmul_si_ur(abase_pz_dst_field K MAYBE_UNUSED, abase_pz_dst_elt_ur w, abase_pz_src_elt u, long v)
{
        abase_pz_elt_ur s;
        abase_pz_elt vx;
        abase_pz_elt_ur_init(K, &s);
        abase_pz_init(K, &vx);
        if (v>0) {
            abase_pz_set_ui(K, vx, v);
            abase_pz_mul_ur(K, s, u, vx);
            abase_pz_elt_ur_add(K, w, w, s);
        } else {
            abase_pz_set_ui(K, vx, -v);
            abase_pz_mul_ur(K, s, u, vx);
            abase_pz_elt_ur_sub(K, w, w, s);
        }
        abase_pz_clear(K, &vx);
        abase_pz_elt_ur_clear(K, &s);
}

/* *pz::code_for_cmp */
static inline
int abase_pz_cmp(abase_pz_dst_field k, abase_pz_src_elt x, abase_pz_src_elt y)
{
        return mpn_cmp(x, y, k->kl);
}

/* *pz::code_for_cmp_ui */
static inline
int abase_pz_cmp_ui(abase_pz_dst_field k, abase_pz_src_elt x, unsigned long y0)
{
        for (int i = k->kl - 1; i > 0; --i) {
        if (x[i] != 0)
            return 1;
        }
        if (x[0] > y0)
        return 1;
        if (x[0] < y0)
        return -1;
        return 0;
}

/* *pz::code_for_is_zero */
static inline
int abase_pz_is_zero(abase_pz_dst_field k, abase_pz_src_elt x)
{
        for (unsigned int i = 0; i < k->kl; ++i)
        if (x[i] != 0)
            return 0;
        return 1;
}

/* *pz::code_for_vec_subvec */
static inline
abase_pz_dst_vec abase_pz_vec_subvec(abase_pz_dst_field K MAYBE_UNUSED, abase_pz_dst_vec v, int i)
{
    return v + i * K->kl;
}

/* *pz::code_for_vec_subvec_const */
static inline
abase_pz_src_vec abase_pz_vec_subvec_const(abase_pz_dst_field K MAYBE_UNUSED, abase_pz_src_vec v, int i)
{
    return v + i * K->kl;
}

/* *pz::code_for_vec_coeff_ptr */
static inline
abase_pz_dst_elt abase_pz_vec_coeff_ptr(abase_pz_dst_field K MAYBE_UNUSED, abase_pz_dst_vec v, int i)
{
    return v + i*K->kl;
}

/* *pz::code_for_vec_coeff_ptr_const */
static inline
abase_pz_src_elt abase_pz_vec_coeff_ptr_const(abase_pz_dst_field K MAYBE_UNUSED, abase_pz_src_vec v, int i)
{
    return v + i*K->kl;
}

/* *pz::code_for_vec_ur_subvec */
static inline
abase_pz_dst_vec_ur abase_pz_vec_ur_subvec(abase_pz_dst_field K MAYBE_UNUSED, abase_pz_dst_vec_ur v, int i)
{
    return v + i * K->url;
}

/* *pz::code_for_vec_ur_subvec_const */
static inline
abase_pz_src_vec_ur abase_pz_vec_ur_subvec_const(abase_pz_dst_field K MAYBE_UNUSED, abase_pz_src_vec_ur v, int i)
{
    return v + i * K->url;
}

/* *pz::code_for_vec_ur_coeff_ptr */
static inline
abase_pz_dst_elt abase_pz_vec_ur_coeff_ptr(abase_pz_dst_field K MAYBE_UNUSED, abase_pz_dst_vec_ur v, int i)
{
    return v + i*K->url;
}

/* *pz::code_for_vec_ur_coeff_ptr_const */
static inline
abase_pz_src_elt abase_pz_vec_ur_coeff_ptr_const(abase_pz_dst_field K MAYBE_UNUSED, abase_pz_src_vec_ur v, int i)
{
    return v + i*K->url;
}

/* *Mpfq::defaults::poly::code_for_poly_init, pz */
static inline
void abase_pz_poly_init(abase_pz_dst_field k MAYBE_UNUSED, abase_pz_poly p, unsigned int n)
{
    abase_pz_vec_init(k, &(p->c), n);
    p->alloc=n;
    p->size=0;
}

/* *Mpfq::defaults::poly::code_for_poly_clear, pz */
static inline
void abase_pz_poly_clear(abase_pz_dst_field k MAYBE_UNUSED, abase_pz_poly p)
{
    abase_pz_vec_clear(k, &(p->c), p->alloc);
}

/* *Mpfq::defaults::poly::code_for_poly_set, pz */
static inline
void abase_pz_poly_set(abase_pz_dst_field k MAYBE_UNUSED, abase_pz_dst_poly w, abase_pz_src_poly u)
{
    if (w->alloc < u->size) {
        abase_pz_vec_reinit(k, &(w->c), w->alloc, u->size);
        w->alloc = u->size;
    }
    abase_pz_vec_set(k, w->c, u->c, u->size);
    w->size = u->size;
}

/* *Mpfq::defaults::poly::code_for_poly_setcoef, pz */
static inline
void abase_pz_poly_setcoef(abase_pz_dst_field k MAYBE_UNUSED, abase_pz_dst_poly w, abase_pz_src_elt x, unsigned int i)
{
    unsigned long j;
    if (w->alloc < (i+1)) {
        abase_pz_vec_reinit(k, &(w->c), w->alloc, i+1);
        w->alloc = i+1;
    }
    abase_pz_vec_setcoef(k, w->c, x, i);
    if (w->size < (i+1)) {
        for (j = w->size; j < i; ++j) {
            abase_pz_vec_setcoef_ui(k, w->c, 0, j);
        }  
        w->size = i+1;
    }
}

/* *Mpfq::defaults::poly::code_for_poly_setcoef_ui, pz */
static inline
void abase_pz_poly_setcoef_ui(abase_pz_dst_field k MAYBE_UNUSED, abase_pz_dst_poly w, unsigned long x, unsigned int i)
{
    unsigned long j;
    if (w->alloc < (i+1)) {
        abase_pz_vec_reinit(k, &(w->c), w->alloc, i+1);
        w->alloc = i+1;
    }
    abase_pz_vec_setcoef_ui(k, w->c, x, i);
    if (w->size < (i+1)) {
        for (j = w->size; j < i; ++j) {
            abase_pz_vec_setcoef_ui(k, w->c, 0, j);
        }  
        w->size = i+1;
    }
}

/* *Mpfq::defaults::poly::code_for_poly_getcoef, pz */
static inline
void abase_pz_poly_getcoef(abase_pz_dst_field k MAYBE_UNUSED, abase_pz_dst_elt x, abase_pz_src_poly w, unsigned int i)
{
    if (w->size < (i+1)) {
       abase_pz_set_ui(k,x,0);
    } else {
       abase_pz_vec_getcoef(k, x, w->c, i);
    }
}

/* *Mpfq::defaults::poly::code_for_poly_deg, pz */
static inline
int abase_pz_poly_deg(abase_pz_dst_field K MAYBE_UNUSED, abase_pz_src_poly w)
{
    if (w->size == 0)
        return -1;
    int deg = w->size-1;
    abase_pz_elt temp;
    abase_pz_init(K, &temp);
    abase_pz_vec_getcoef(K, temp, w->c, deg);
    int comp=abase_pz_cmp_ui(K, temp, 0);
    while ((deg >= 0) && (comp == 0)){
        deg--;
        if (deg>=0) {
           abase_pz_vec_getcoef(K, temp, w->c, deg);
           comp=abase_pz_cmp_ui(K, temp, 0);
        }
    }
    abase_pz_clear(K, &temp);
    return deg;
}

/* *Mpfq::defaults::poly::code_for_poly_add, pz */
static inline
void abase_pz_poly_add(abase_pz_dst_field k MAYBE_UNUSED, abase_pz_dst_poly w, abase_pz_src_poly u, abase_pz_src_poly v)
{
    unsigned int minsize MAYBE_UNUSED = MIN(u->size, v->size);
    unsigned int maxsize MAYBE_UNUSED = MAX(u->size, v->size);
    if (w->alloc < maxsize) {
        abase_pz_vec_reinit(k, &(w->c), w->alloc, maxsize);
        w->alloc = maxsize;
    }
    if (u->size <= v->size) {
        abase_pz_vec_add(k, w->c, u->c, v->c, u->size);
        abase_pz_vec_set_partial(k, (w->c), (v->c), u->size, u->size, v->size-u->size);
    } else {
        abase_pz_vec_add(k, w->c, u->c, v->c, v->size);
        abase_pz_vec_set_partial(k, (w->c), (u->c), v->size, v->size, u->size-v->size);
    }
    w->size=maxsize;
    unsigned int wdeg = abase_pz_poly_deg(k, w);
    w->size=wdeg+1;
}

/* *Mpfq::defaults::poly::code_for_poly_sub, pz */
static inline
void abase_pz_poly_sub(abase_pz_dst_field k MAYBE_UNUSED, abase_pz_dst_poly w, abase_pz_src_poly u, abase_pz_src_poly v)
{
    unsigned int minsize MAYBE_UNUSED = MIN(u->size, v->size);
    unsigned int maxsize MAYBE_UNUSED = MAX(u->size, v->size);
    if (w->alloc < maxsize) {
        abase_pz_vec_reinit(k, &(w->c), w->alloc, maxsize);
        w->alloc = maxsize;
    }
    if (u->size <= v->size) {
        abase_pz_vec_sub(k, w->c, u->c, v->c, u->size);
        unsigned int i;
        abase_pz_elt temp;
        abase_pz_init(k, &temp);
        for (i = u->size; i< v->size; ++i) {
            abase_pz_poly_getcoef(k, temp, v, i);
            abase_pz_neg(k, temp, temp);
            abase_pz_poly_setcoef(k, w, temp, i);
        }
    } else {
        abase_pz_vec_sub(k, w->c, u->c, v->c, v->size);
        abase_pz_vec_set_partial(k, (w->c), (u->c), v->size, v->size, u->size-v->size);
    }
    w->size=maxsize;
    //abase_pz_poly_neg(k, w, v);
    //abase_pz_poly_add(k, w, u, w);
    unsigned int wdeg = abase_pz_poly_deg(k, w);
    w->size=wdeg+1;
}

/* *Mpfq::defaults::poly::code_for_poly_add_ui, pz */
static inline
void abase_pz_poly_add_ui(abase_pz_dst_field k MAYBE_UNUSED, abase_pz_dst_poly w, abase_pz_src_poly u, unsigned long x)
{
    if (u->size == 0) {
        if (x == 0) {
            w->size = 0;
            return;
        }
        if (w->alloc == 0) {
            abase_pz_vec_reinit(k, &(w->c), w->alloc, 1);
            w->alloc = 1;
        }
        w->size = 1;
        abase_pz_vec_setcoef_ui(k, w->c, x, 0);
        return;
    }
    if (w->alloc < u->size) {
        abase_pz_vec_reinit(k, &(w->c), w->alloc, u->size);
        w->alloc = u->size;
    }
    abase_pz_elt temp;
    abase_pz_init(k, &temp);
    abase_pz_poly_getcoef(k, temp, u, 0);
    abase_pz_add_ui(k, temp, temp, x);
    abase_pz_vec_setcoef(k, w->c, temp, 0);
    abase_pz_vec_set_partial(k, w->c, u->c, 1, 1, u->size-1);
    w->size=u->size;
    abase_pz_clear(k, &temp);
}

/* *Mpfq::defaults::poly::code_for_poly_sub_ui, pz */
static inline
void abase_pz_poly_sub_ui(abase_pz_dst_field k MAYBE_UNUSED, abase_pz_dst_poly w, abase_pz_src_poly u, unsigned long x)
{
    if (u->size == 0) {
        if (x == 0) {
            w->size = 0;
            return;
        }
        if (w->alloc == 0) {
            abase_pz_vec_reinit(k, &(w->c), w->alloc, 1);
            w->alloc = 1;
        }
        w->size = 1;
        abase_pz_elt temp;
        abase_pz_init(k, &temp);
        abase_pz_set_ui(k, temp, x);
        abase_pz_neg(k, temp, temp);
        abase_pz_vec_setcoef(k, w->c, temp, 0);
        return;
    }
    if (w->alloc < u->size) {
        abase_pz_vec_reinit(k, &(w->c), w->alloc, u->size);
        w->alloc = u->size;
    }
    abase_pz_elt temp;
    abase_pz_init(k, &temp);
    abase_pz_poly_getcoef(k, temp, u, 0);
    abase_pz_sub_ui(k, temp, temp, x);
    abase_pz_vec_setcoef(k, w->c, temp, 0);
    abase_pz_vec_set_partial(k, w->c, u->c, 1, 1, u->size-1);
    w->size=u->size;
    abase_pz_clear(k, &temp);
}

/* *Mpfq::defaults::poly::code_for_poly_neg, pz */
static inline
void abase_pz_poly_neg(abase_pz_dst_field k MAYBE_UNUSED, abase_pz_dst_poly w, abase_pz_src_poly u)
{
    if (w->alloc < u->size) {
        abase_pz_vec_reinit(k, &(w->c), w->alloc, u->size);
        w->alloc = u->size;
    }
    abase_pz_vec_neg(k, w->c, u->c, u->size);
    w->size = u->size;
}

/* *Mpfq::defaults::poly::code_for_poly_scal_mul, pz */
static inline
void abase_pz_poly_scal_mul(abase_pz_dst_field k MAYBE_UNUSED, abase_pz_dst_poly w, abase_pz_src_poly u, abase_pz_src_elt x)
{
    if (abase_pz_cmp_ui(k, x, 0) == 0) {
        w->size = 0;
        return;
    }
    unsigned int n = u->size;
    if (w->alloc < n) {
        abase_pz_vec_reinit(k, &(w->c), w->alloc, n);
        w->alloc = n;
    }
    abase_pz_vec_scal_mul(k, w->c, u->c, x, n);
    w->size=n;
}

/* *Mpfq::defaults::poly::code_for_poly_mul, pz */
static inline
void abase_pz_poly_mul(abase_pz_dst_field K MAYBE_UNUSED, abase_pz_dst_poly w, abase_pz_src_poly u, abase_pz_src_poly v)
{
    unsigned int usize = abase_pz_poly_deg(K, u)+1;
    unsigned int vsize = abase_pz_poly_deg(K, v)+1;
    if ((usize == 0) || (vsize == 0)) {
        w->size = 0;
        return;
    }
    unsigned int wsize = usize + vsize - 1;
    if (w->alloc < wsize) {
        abase_pz_vec_reinit(K, &(w->c), w->alloc, wsize);
        w->alloc = wsize;
    }
    abase_pz_vec_conv(K, w->c, u->c, usize, v->c, vsize);
    w->size=wsize;
}

/* *Mpfq::defaults::polygcd::code_for_poly_gcd, Mpfq::defaults::poly, pz */
static inline
void abase_pz_poly_gcd(abase_pz_dst_field k MAYBE_UNUSED, abase_pz_dst_poly g, abase_pz_src_poly a0, abase_pz_src_poly b0)
{
    abase_pz_poly a,b,q,r;
    int da0=abase_pz_poly_deg(k,a0), db0=abase_pz_poly_deg(k,b0);
    if (db0==-1)
     abase_pz_poly_set(k,g,a0);
    else {
     abase_pz_poly_init(k,a,da0+1);
     abase_pz_poly_init(k,b,db0+1);
     abase_pz_poly_init(k,q,1);
     abase_pz_poly_init(k,r,db0);
     abase_pz_poly_set(k,a,a0);
     abase_pz_poly_set(k,b,b0);
     while (abase_pz_poly_deg(k,b)>=0) {
      abase_pz_poly_divmod(k,q,r,a,b);
      abase_pz_poly_set(k,a,b);
      abase_pz_poly_set(k,b,r); 
     }
     abase_pz_poly_setmonic(k,g,a);
    abase_pz_poly_clear(k,a);
    abase_pz_poly_clear(k,b);
    abase_pz_poly_clear(k,q);
    abase_pz_poly_clear(k,r);
    }
}

/* *Mpfq::defaults::polygcd::code_for_poly_xgcd, Mpfq::defaults::poly, pz */
static inline
void abase_pz_poly_xgcd(abase_pz_dst_field k MAYBE_UNUSED, abase_pz_dst_poly g, abase_pz_dst_poly u0, abase_pz_dst_poly v0, abase_pz_src_poly a0, abase_pz_src_poly b0)
{
    abase_pz_poly a,b,u,v,w,x,q,r;
    abase_pz_elt c;
    abase_pz_init(k,&c);
    int da0=abase_pz_poly_deg(k,a0), db0=abase_pz_poly_deg(k,b0), dega;
    if (db0==-1) {
     if (da0==-1) {
      abase_pz_poly_set(k,u0,a0);
      abase_pz_poly_set(k,v0,b0);
      abase_pz_poly_set(k,g,a0);
     } else {
      abase_pz_poly_getcoef(k,c,a0,da0);
      abase_pz_inv(k,c,c);
      abase_pz_poly_scal_mul(k,g,a0,c);
      abase_pz_poly_set(k,v0,b0);
      abase_pz_poly_set(k,u0,b0);
      abase_pz_poly_setcoef(k,u0,c,0);
     }
    }
    else {
     abase_pz_poly_init(k,a,da0+1);
     abase_pz_poly_init(k,b,db0+1);
     abase_pz_poly_init(k,q,1);
     abase_pz_poly_init(k,r,db0);
     abase_pz_poly_set(k,a,a0);
     abase_pz_poly_set(k,b,b0);
     abase_pz_poly_init(k,u,1);
     abase_pz_poly_init(k,v,1);
     abase_pz_poly_init(k,w,1);
     abase_pz_poly_init(k,x,1);
     abase_pz_poly_setcoef_ui(k,u,1,0);
     abase_pz_poly_setcoef_ui(k,x,1,0);
     /* u*a_initial + v*b_initial = a */
     /* w*a_initial + x*b_initial = b */
     while (abase_pz_poly_deg(k,b)>=0) {
      abase_pz_poly_divmod(k,q,r,a,b);
      abase_pz_poly_set(k,a,b);  /* a,b <- b,a-qb=r */
      abase_pz_poly_set(k,b,r);
      abase_pz_poly_mul(k,r,q,w);
      abase_pz_poly_sub(k,r,u,r);
      abase_pz_poly_set(k,u,w);   /* u,w <- w,u-qw */
      abase_pz_poly_set(k,w,r);
      abase_pz_poly_mul(k,r,q,x); /* v,x <- x,v-qx */
      abase_pz_poly_sub(k,r,v,r);
      abase_pz_poly_set(k,v,x);
      abase_pz_poly_set(k,x,r);
     }
     dega=abase_pz_poly_deg(k,a);
     abase_pz_poly_getcoef(k,c,a,dega);
     abase_pz_inv(k,c,c);
     abase_pz_poly_scal_mul(k,g,a,c);
     abase_pz_poly_scal_mul(k,u0,u,c);
     abase_pz_poly_scal_mul(k,v0,v,c);
     abase_pz_poly_clear(k,a);
     abase_pz_poly_clear(k,b);
     abase_pz_poly_clear(k,u);
     abase_pz_poly_clear(k,v);
     abase_pz_poly_clear(k,w);
     abase_pz_poly_clear(k,x);
     abase_pz_poly_clear(k,q);
     abase_pz_poly_clear(k,r);
    }
    abase_pz_clear(k,&c);
}

/* *Mpfq::defaults::poly::code_for_poly_random, pz */
static inline
void abase_pz_poly_random(abase_pz_dst_field k MAYBE_UNUSED, abase_pz_dst_poly w, unsigned int n, gmp_randstate_t state)
{
    n++;
    if (w->alloc < n) {
        abase_pz_vec_reinit(k, &(w->c), w->alloc, n);
        w->alloc = n;
    }
    abase_pz_vec_random(k, w->c, n,state);
    w->size=n;
    int wdeg = abase_pz_poly_deg(k, w);
    w->size=wdeg+1;
}

/* *Mpfq::defaults::poly::code_for_poly_random2, pz */
static inline
void abase_pz_poly_random2(abase_pz_dst_field k MAYBE_UNUSED, abase_pz_dst_poly w, unsigned int n, gmp_randstate_t state)
{
    n++;
    if (w->alloc < n) {
        abase_pz_vec_reinit(k, &(w->c), w->alloc, n);
        w->alloc = n;
    }
    abase_pz_vec_random2(k, w->c, n,state);
    w->size=n;
    int wdeg = abase_pz_poly_deg(k, w);
    w->size=wdeg+1;
}

/* *Mpfq::defaults::poly::code_for_poly_cmp, pz */
static inline
int abase_pz_poly_cmp(abase_pz_dst_field k MAYBE_UNUSED, abase_pz_src_poly u, abase_pz_src_poly v)
{
    if (u->size != v->size)
        return (int)(u->size) - (int)(v->size);
    else
        return abase_pz_vec_cmp(k, u->c, v->c, u->size);
}

/* *Mpfq::defaults::poly::code_for_poly_asprint, pz */
static inline
void abase_pz_poly_asprint(abase_pz_dst_field k MAYBE_UNUSED, char * * pstr, abase_pz_src_poly w)
{
    abase_pz_vec_asprint(k, pstr, w->c, w->size);
}

/* *Mpfq::defaults::poly::code_for_poly_fprint, pz */
static inline
void abase_pz_poly_fprint(abase_pz_dst_field k MAYBE_UNUSED, FILE * file, abase_pz_src_poly w)
{
    abase_pz_vec_fprint(k, file, w->c, w->size);
}

/* *Mpfq::defaults::poly::code_for_poly_print, pz */
static inline
void abase_pz_poly_print(abase_pz_dst_field k MAYBE_UNUSED, abase_pz_src_poly w)
{
    abase_pz_vec_print(k, w->c, w->size);
}

/* *Mpfq::defaults::poly::code_for_poly_sscan, pz */
static inline
int abase_pz_poly_sscan(abase_pz_dst_field k MAYBE_UNUSED, abase_pz_dst_poly w, const char * str)
{
    int ret;
    ret = abase_pz_vec_sscan(k, &(w->c), &(w->alloc), str);
    w->size = w->alloc;
    return ret;
}

/* *Mpfq::defaults::poly::code_for_poly_fscan, pz */
static inline
int abase_pz_poly_fscan(abase_pz_dst_field k MAYBE_UNUSED, FILE * file, abase_pz_dst_poly w)
{
    int ret;
    ret = abase_pz_vec_fscan(k, file, &(w->c), &(w->alloc));
    w->size = w->alloc;
    return ret;
}

/* *simd_pz::code_for_set_ui_at */
static inline
void abase_pz_set_ui_at(abase_pz_dst_field K MAYBE_UNUSED, abase_pz_dst_elt p, int k MAYBE_UNUSED, unsigned long v)
{
    abase_pz_set_ui(K,p,v);
}

/* *simd_pz::code_for_set_ui_all */
static inline
void abase_pz_set_ui_all(abase_pz_dst_field K MAYBE_UNUSED, abase_pz_dst_elt p, unsigned long v)
{
    abase_pz_set_ui(K,p,v);
}

/* *simd_pz::code_for_elt_ur_set_ui_at */
static inline
void abase_pz_elt_ur_set_ui_at(abase_pz_dst_field K MAYBE_UNUSED, abase_pz_dst_elt p, int k MAYBE_UNUSED, unsigned long v)
{
    abase_pz_set_ui(K,p,v);
}

/* *simd_pz::code_for_elt_ur_set_ui_all */
static inline
void abase_pz_elt_ur_set_ui_all(abase_pz_dst_field K MAYBE_UNUSED, abase_pz_dst_elt p, unsigned long v)
{
    abase_pz_set_ui(K,p,v);
}

static inline
void abase_pz_oo_field_clear(abase_vbase_ptr f)
{
    abase_pz_field_clear((abase_pz_dst_field)(f->obj));
    free(f->obj);
    f->obj = NULL;
}


#endif  /* ABASE_PZ_H_ */

/* vim:set ft=cpp: */
