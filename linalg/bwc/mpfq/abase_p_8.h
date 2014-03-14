#ifndef ABASE_P_8_H_
#define ABASE_P_8_H_

/* MPFQ generated file -- do not edit */

#include "mpfq.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include "assert.h"
#include <limits.h>
#include "fixmp.h"
#include "mpfq_gfp_common.h"
#include "select_mpi.h"
#include "abase_vbase.h"
#ifdef	MPFQ_LAST_GENERATED_TAG
#undef	MPFQ_LAST_GENERATED_TAG
#endif
#define MPFQ_LAST_GENERATED_TAG      p_8

/* Active handler: simd_gfp */
/* Automatically generated code  */
/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::vec */
/* Active handler: Mpfq::defaults::poly */
/* Active handler: Mpfq::gfp::field */
/* Active handler: Mpfq::gfp::elt */
/* Active handler: Mpfq::defaults::mpi_flat */
/* Options used:{
   family=[ { cpp_ifdef=COMPILE_MPFQ_PRIME_FIELDS, tag=p_8, }, ],
   fieldtype=prime,
   n=8,
   nn=17,
   opthw=,
   tag=p_8,
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
     [ (?^:abase_p_8_elt \*), void *, ],
     [ (?^:abase_p_8_src_elt\b), const void *, ],
     [ (?^:abase_p_8_elt\b), void *, ],
     [ (?^:abase_p_8_dst_elt\b), void *, ],
     [ (?^:abase_p_8_elt_ur \*), void *, ],
     [ (?^:abase_p_8_src_elt_ur\b), const void *, ],
     [ (?^:abase_p_8_elt_ur\b), void *, ],
     [ (?^:abase_p_8_dst_elt_ur\b), void *, ],
     [ (?^:abase_p_8_vec \*), void *, ],
     [ (?^:abase_p_8_src_vec\b), const void *, ],
     [ (?^:abase_p_8_vec\b), void *, ],
     [ (?^:abase_p_8_dst_vec\b), void *, ],
     [ (?^:abase_p_8_vec_ur \*), void *, ],
     [ (?^:abase_p_8_src_vec_ur\b), const void *, ],
     [ (?^:abase_p_8_vec_ur\b), void *, ],
     [ (?^:abase_p_8_dst_vec_ur\b), void *, ],
     [ (?^:abase_p_8_poly \*), void *, ],
     [ (?^:abase_p_8_src_poly\b), const void *, ],
     [ (?^:abase_p_8_poly\b), void *, ],
     [ (?^:abase_p_8_dst_poly\b), void *, ],
     ],
    },
   vtag=p_8,
   w=64,
   } */

typedef mpfq_p_field abase_p_8_field;
typedef mpfq_p_dst_field abase_p_8_dst_field;

typedef unsigned long abase_p_8_elt[8];
typedef unsigned long * abase_p_8_dst_elt;
typedef const unsigned long * abase_p_8_src_elt;

typedef unsigned long abase_p_8_elt_ur[17];
typedef unsigned long * abase_p_8_dst_elt_ur;
typedef const unsigned long * abase_p_8_src_elt_ur;

typedef abase_p_8_elt * abase_p_8_vec;
typedef abase_p_8_elt * abase_p_8_dst_vec;
typedef abase_p_8_elt * abase_p_8_src_vec;

typedef abase_p_8_elt_ur * abase_p_8_vec_ur;
typedef abase_p_8_elt_ur * abase_p_8_dst_vec_ur;
typedef abase_p_8_elt_ur * abase_p_8_src_vec_ur;

typedef struct {
  abase_p_8_vec c;
  unsigned int alloc;
  unsigned int size;
} abase_p_8_poly_struct;
typedef abase_p_8_poly_struct abase_p_8_poly [1];
typedef abase_p_8_poly_struct * abase_p_8_dst_poly;
typedef abase_p_8_poly_struct * abase_p_8_src_poly;

#ifdef  __cplusplus
extern "C" {
#endif
/* *Mpfq::defaults::code_for_impl_name, Mpfq::gfp */
#define abase_p_8_impl_name()	"p_8"
/* *Mpfq::gfp::field::code_for_impl_max_characteristic_bits, Mpfq::gfp */
#define abase_p_8_impl_max_characteristic_bits()	512
/* *Mpfq::gfp::field::code_for_impl_max_degree, Mpfq::gfp */
#define abase_p_8_impl_max_degree()	1

/* Functions operating on the field structure */
void abase_p_8_field_characteristic(abase_p_8_dst_field, mpz_t);
/* *Mpfq::gfp::field::code_for_field_degree, Mpfq::gfp */
#define abase_p_8_field_degree(K)	1
static inline
void abase_p_8_field_init(abase_p_8_dst_field);
void abase_p_8_field_clear(abase_p_8_dst_field);
void abase_p_8_field_specify(abase_p_8_dst_field, unsigned long, void *);
/* *Mpfq::gfp::field::code_for_field_setopt, Mpfq::gfp */
#define abase_p_8_field_setopt(f, x, y)	/**/

/* Element allocation functions */
static inline
void abase_p_8_init(abase_p_8_dst_field, abase_p_8_elt *);
static inline
void abase_p_8_clear(abase_p_8_dst_field, abase_p_8_elt *);

/* Elementary assignment functions */
static inline
void abase_p_8_set(abase_p_8_dst_field, abase_p_8_dst_elt, abase_p_8_src_elt);
static inline
void abase_p_8_set_ui(abase_p_8_dst_field, abase_p_8_dst_elt, unsigned long);
static inline
void abase_p_8_set_zero(abase_p_8_dst_field, abase_p_8_dst_elt);
static inline
unsigned long abase_p_8_get_ui(abase_p_8_dst_field, abase_p_8_src_elt);
static inline
void abase_p_8_set_mpn(abase_p_8_dst_field, abase_p_8_dst_elt, mp_limb_t *, size_t);
static inline
void abase_p_8_set_mpz(abase_p_8_dst_field, abase_p_8_dst_elt, mpz_t);
static inline
void abase_p_8_get_mpn(abase_p_8_dst_field, mp_limb_t *, abase_p_8_src_elt);
static inline
void abase_p_8_get_mpz(abase_p_8_dst_field, mpz_t, abase_p_8_src_elt);

/* Assignment of random values */
static inline
void abase_p_8_random(abase_p_8_dst_field, abase_p_8_dst_elt, gmp_randstate_t);
static inline
void abase_p_8_random2(abase_p_8_dst_field, abase_p_8_dst_elt, gmp_randstate_t);

/* Arithmetic operations on elements */
static inline
void abase_p_8_add(abase_p_8_dst_field, abase_p_8_dst_elt, abase_p_8_src_elt, abase_p_8_src_elt);
static inline
void abase_p_8_sub(abase_p_8_dst_field, abase_p_8_dst_elt, abase_p_8_src_elt, abase_p_8_src_elt);
static inline
void abase_p_8_neg(abase_p_8_dst_field, abase_p_8_dst_elt, abase_p_8_src_elt);
static inline
void abase_p_8_mul(abase_p_8_dst_field, abase_p_8_dst_elt, abase_p_8_src_elt, abase_p_8_src_elt);
static inline
void abase_p_8_sqr(abase_p_8_dst_field, abase_p_8_dst_elt, abase_p_8_src_elt);
static inline
int abase_p_8_is_sqr(abase_p_8_dst_field, abase_p_8_src_elt);
int abase_p_8_sqrt(abase_p_8_dst_field, abase_p_8_dst_elt, abase_p_8_src_elt);
static inline
void abase_p_8_pow(abase_p_8_dst_field, abase_p_8_dst_elt, abase_p_8_src_elt, unsigned long *, size_t);
/* *Mpfq::gfp::elt::code_for_frobenius, Mpfq::gfp */
#define abase_p_8_frobenius(k, x, y)	abase_p_8_set(k, x, y)
static inline
void abase_p_8_add_ui(abase_p_8_dst_field, abase_p_8_dst_elt, abase_p_8_src_elt, unsigned long);
static inline
void abase_p_8_sub_ui(abase_p_8_dst_field, abase_p_8_dst_elt, abase_p_8_src_elt, unsigned long);
static inline
void abase_p_8_mul_ui(abase_p_8_dst_field, abase_p_8_dst_elt, abase_p_8_src_elt, unsigned long);
static inline
int abase_p_8_inv(abase_p_8_dst_field, abase_p_8_dst_elt, abase_p_8_src_elt);
#define HAVE_abase_p_8_hadamard
static inline
void abase_p_8_hadamard(abase_p_8_dst_field, abase_p_8_dst_elt, abase_p_8_dst_elt, abase_p_8_dst_elt, abase_p_8_dst_elt);

/* Operations involving unreduced elements */
static inline
void abase_p_8_elt_ur_init(abase_p_8_dst_field, abase_p_8_elt_ur *);
static inline
void abase_p_8_elt_ur_clear(abase_p_8_dst_field, abase_p_8_elt_ur *);
static inline
void abase_p_8_elt_ur_set(abase_p_8_dst_field, abase_p_8_dst_elt_ur, abase_p_8_src_elt_ur);
static inline
void abase_p_8_elt_ur_set_elt(abase_p_8_dst_field, abase_p_8_dst_elt_ur, abase_p_8_src_elt);
static inline
void abase_p_8_elt_ur_set_zero(abase_p_8_dst_field, abase_p_8_dst_elt_ur);
static inline
void abase_p_8_elt_ur_set_ui(abase_p_8_dst_field, abase_p_8_dst_elt_ur, unsigned long);
static inline
void abase_p_8_elt_ur_add(abase_p_8_dst_field, abase_p_8_dst_elt_ur, abase_p_8_src_elt_ur, abase_p_8_src_elt_ur);
static inline
void abase_p_8_elt_ur_neg(abase_p_8_dst_field, abase_p_8_dst_elt_ur, abase_p_8_src_elt_ur);
static inline
void abase_p_8_elt_ur_sub(abase_p_8_dst_field, abase_p_8_dst_elt_ur, abase_p_8_src_elt_ur, abase_p_8_src_elt_ur);
static inline
void abase_p_8_mul_ur(abase_p_8_dst_field, abase_p_8_dst_elt_ur, abase_p_8_src_elt, abase_p_8_src_elt);
static inline
void abase_p_8_sqr_ur(abase_p_8_dst_field, abase_p_8_dst_elt_ur, abase_p_8_src_elt);
static inline
void abase_p_8_reduce(abase_p_8_dst_field, abase_p_8_dst_elt, abase_p_8_dst_elt_ur);
#define HAVE_abase_p_8_normalize
static inline
void abase_p_8_normalize(abase_p_8_dst_field, abase_p_8_dst_elt);
#define HAVE_abase_p_8_addmul_si_ur
static inline
void abase_p_8_addmul_si_ur(abase_p_8_dst_field, abase_p_8_dst_elt_ur, abase_p_8_src_elt, long);

/* Comparison functions */
static inline
int abase_p_8_cmp(abase_p_8_dst_field, abase_p_8_src_elt, abase_p_8_src_elt);
static inline
int abase_p_8_cmp_ui(abase_p_8_dst_field, abase_p_8_src_elt, unsigned long);
static inline
int abase_p_8_is_zero(abase_p_8_dst_field, abase_p_8_src_elt);

/* Input/output functions */
void abase_p_8_asprint(abase_p_8_dst_field, char * *, abase_p_8_src_elt);
void abase_p_8_fprint(abase_p_8_dst_field, FILE *, abase_p_8_src_elt);
/* *Mpfq::defaults::code_for_print, Mpfq::gfp */
#define abase_p_8_print(k, x)	abase_p_8_fprint(k,stdout,x)
int abase_p_8_sscan(abase_p_8_dst_field, abase_p_8_dst_elt, const char *);
int abase_p_8_fscan(abase_p_8_dst_field, FILE *, abase_p_8_dst_elt);
/* *Mpfq::defaults::code_for_scan, Mpfq::gfp */
#define abase_p_8_scan(k, x)	abase_p_8_fscan(k,stdin,x)

/* Vector functions */
void abase_p_8_vec_init(abase_p_8_dst_field, abase_p_8_vec *, unsigned int);
void abase_p_8_vec_reinit(abase_p_8_dst_field, abase_p_8_vec *, unsigned int, unsigned int);
void abase_p_8_vec_clear(abase_p_8_dst_field, abase_p_8_vec *, unsigned int);
static inline
void abase_p_8_vec_set(abase_p_8_dst_field, abase_p_8_dst_vec, abase_p_8_src_vec, unsigned int);
static inline
void abase_p_8_vec_set_partial(abase_p_8_dst_field, abase_p_8_dst_vec, abase_p_8_src_vec, unsigned int, unsigned int, unsigned int);
static inline
void abase_p_8_vec_set_zero(abase_p_8_dst_field, abase_p_8_dst_vec, unsigned int);
static inline
void abase_p_8_vec_setcoef(abase_p_8_dst_field, abase_p_8_dst_vec, abase_p_8_src_elt, unsigned int);
static inline
void abase_p_8_vec_setcoef_ui(abase_p_8_dst_field, abase_p_8_dst_vec, unsigned long, unsigned int);
static inline
void abase_p_8_vec_getcoef(abase_p_8_dst_field, abase_p_8_dst_elt, abase_p_8_src_vec, unsigned int);
static inline
void abase_p_8_vec_add(abase_p_8_dst_field, abase_p_8_dst_vec, abase_p_8_src_vec, abase_p_8_src_vec, unsigned int);
static inline
void abase_p_8_vec_neg(abase_p_8_dst_field, abase_p_8_dst_vec, abase_p_8_src_vec, unsigned int);
static inline
void abase_p_8_vec_rev(abase_p_8_dst_field, abase_p_8_dst_vec, abase_p_8_src_vec, unsigned int);
static inline
void abase_p_8_vec_sub(abase_p_8_dst_field, abase_p_8_dst_vec, abase_p_8_src_vec, abase_p_8_src_vec, unsigned int);
static inline
void abase_p_8_vec_scal_mul(abase_p_8_dst_field, abase_p_8_dst_vec, abase_p_8_src_vec, abase_p_8_src_elt, unsigned int);
static inline
void abase_p_8_vec_conv(abase_p_8_dst_field, abase_p_8_dst_vec, abase_p_8_src_vec, unsigned int, abase_p_8_src_vec, unsigned int);
static inline
void abase_p_8_vec_random(abase_p_8_dst_field, abase_p_8_dst_vec, unsigned int, gmp_randstate_t);
static inline
void abase_p_8_vec_random2(abase_p_8_dst_field, abase_p_8_dst_vec, unsigned int, gmp_randstate_t);
static inline
int abase_p_8_vec_cmp(abase_p_8_dst_field, abase_p_8_src_vec, abase_p_8_src_vec, unsigned int);
static inline
int abase_p_8_vec_is_zero(abase_p_8_dst_field, abase_p_8_src_vec, unsigned int);
static inline
abase_p_8_dst_vec abase_p_8_vec_subvec(abase_p_8_dst_field, abase_p_8_dst_vec, int);
static inline
abase_p_8_src_vec abase_p_8_vec_subvec_const(abase_p_8_dst_field, abase_p_8_src_vec, int);
static inline
abase_p_8_dst_elt abase_p_8_vec_coeff_ptr(abase_p_8_dst_field, abase_p_8_dst_vec, int);
static inline
abase_p_8_src_elt abase_p_8_vec_coeff_ptr_const(abase_p_8_dst_field, abase_p_8_src_vec, int);
void abase_p_8_vec_asprint(abase_p_8_dst_field, char * *, abase_p_8_src_vec, unsigned int);
void abase_p_8_vec_fprint(abase_p_8_dst_field, FILE *, abase_p_8_src_vec, unsigned int);
void abase_p_8_vec_print(abase_p_8_dst_field, abase_p_8_src_vec, unsigned int);
int abase_p_8_vec_sscan(abase_p_8_dst_field, abase_p_8_vec *, unsigned int *, const char *);
int abase_p_8_vec_fscan(abase_p_8_dst_field, FILE *, abase_p_8_vec *, unsigned int *);
/* *Mpfq::defaults::vec::io::code_for_vec_scan, Mpfq::defaults::vec, Mpfq::gfp */
#define abase_p_8_vec_scan(K, w, n)	abase_p_8_vec_fscan(K,stdout,w,n)
void abase_p_8_vec_ur_init(abase_p_8_dst_field, abase_p_8_vec_ur *, unsigned int);
static inline
void abase_p_8_vec_ur_set_zero(abase_p_8_dst_field, abase_p_8_dst_vec_ur, unsigned int);
static inline
void abase_p_8_vec_ur_set_vec(abase_p_8_dst_field, abase_p_8_dst_vec_ur, abase_p_8_src_vec, unsigned int);
void abase_p_8_vec_ur_reinit(abase_p_8_dst_field, abase_p_8_vec_ur *, unsigned int, unsigned int);
void abase_p_8_vec_ur_clear(abase_p_8_dst_field, abase_p_8_vec_ur *, unsigned int);
static inline
void abase_p_8_vec_ur_set(abase_p_8_dst_field, abase_p_8_dst_vec_ur, abase_p_8_src_vec_ur, unsigned int);
static inline
void abase_p_8_vec_ur_setcoef(abase_p_8_dst_field, abase_p_8_dst_vec_ur, abase_p_8_src_elt_ur, unsigned int);
static inline
void abase_p_8_vec_ur_getcoef(abase_p_8_dst_field, abase_p_8_dst_elt_ur, abase_p_8_src_vec_ur, unsigned int);
static inline
void abase_p_8_vec_ur_add(abase_p_8_dst_field, abase_p_8_dst_vec_ur, abase_p_8_src_vec_ur, abase_p_8_src_vec_ur, unsigned int);
static inline
void abase_p_8_vec_ur_sub(abase_p_8_dst_field, abase_p_8_dst_vec_ur, abase_p_8_src_vec_ur, abase_p_8_src_vec_ur, unsigned int);
static inline
void abase_p_8_vec_ur_neg(abase_p_8_dst_field, abase_p_8_dst_vec_ur, abase_p_8_src_vec_ur, unsigned int);
static inline
void abase_p_8_vec_ur_rev(abase_p_8_dst_field, abase_p_8_dst_vec_ur, abase_p_8_src_vec_ur, unsigned int);
static inline
void abase_p_8_vec_scal_mul_ur(abase_p_8_dst_field, abase_p_8_dst_vec_ur, abase_p_8_src_vec, abase_p_8_src_elt, unsigned int);
static inline
void abase_p_8_vec_conv_ur_n(abase_p_8_dst_field, abase_p_8_dst_vec_ur, abase_p_8_src_vec, abase_p_8_src_vec, unsigned int);
void abase_p_8_vec_conv_ur_ks(abase_p_8_dst_field, abase_p_8_dst_vec_ur, abase_p_8_src_vec, unsigned int, abase_p_8_src_vec, unsigned int);
static inline
void abase_p_8_vec_conv_ur(abase_p_8_dst_field, abase_p_8_dst_vec_ur, abase_p_8_src_vec, unsigned int, abase_p_8_src_vec, unsigned int);
static inline
void abase_p_8_vec_reduce(abase_p_8_dst_field, abase_p_8_dst_vec, abase_p_8_dst_vec_ur, unsigned int);
static inline
abase_p_8_dst_vec_ur abase_p_8_vec_ur_subvec(abase_p_8_dst_field, abase_p_8_dst_vec_ur, int);
static inline
abase_p_8_src_vec_ur abase_p_8_vec_ur_subvec_const(abase_p_8_dst_field, abase_p_8_src_vec_ur, int);
static inline
abase_p_8_dst_elt abase_p_8_vec_ur_coeff_ptr(abase_p_8_dst_field, abase_p_8_dst_vec_ur, int);
static inline
abase_p_8_src_elt abase_p_8_vec_ur_coeff_ptr_const(abase_p_8_dst_field, abase_p_8_src_vec_ur, int);
/* *Mpfq::defaults::flatdata::code_for_vec_elt_stride, Mpfq::gfp::elt, Mpfq::gfp */
#define abase_p_8_vec_elt_stride(K, n)	((n)*sizeof(abase_p_8_elt))
/* *Mpfq::defaults::flatdata::code_for_vec_ur_elt_stride, Mpfq::gfp::elt, Mpfq::gfp */
#define abase_p_8_vec_ur_elt_stride(K, n)	((n)*sizeof(abase_p_8_elt_ur))

/* Polynomial functions */
static inline
void abase_p_8_poly_init(abase_p_8_dst_field, abase_p_8_poly, unsigned int);
static inline
void abase_p_8_poly_clear(abase_p_8_dst_field, abase_p_8_poly);
static inline
void abase_p_8_poly_set(abase_p_8_dst_field, abase_p_8_dst_poly, abase_p_8_src_poly);
void abase_p_8_poly_setmonic(abase_p_8_dst_field, abase_p_8_dst_poly, abase_p_8_src_poly);
static inline
void abase_p_8_poly_setcoef(abase_p_8_dst_field, abase_p_8_dst_poly, abase_p_8_src_elt, unsigned int);
static inline
void abase_p_8_poly_setcoef_ui(abase_p_8_dst_field, abase_p_8_dst_poly, unsigned long, unsigned int);
static inline
void abase_p_8_poly_getcoef(abase_p_8_dst_field, abase_p_8_dst_elt, abase_p_8_src_poly, unsigned int);
static inline
int abase_p_8_poly_deg(abase_p_8_dst_field, abase_p_8_src_poly);
static inline
void abase_p_8_poly_add(abase_p_8_dst_field, abase_p_8_dst_poly, abase_p_8_src_poly, abase_p_8_src_poly);
static inline
void abase_p_8_poly_sub(abase_p_8_dst_field, abase_p_8_dst_poly, abase_p_8_src_poly, abase_p_8_src_poly);
static inline
void abase_p_8_poly_add_ui(abase_p_8_dst_field, abase_p_8_dst_poly, abase_p_8_src_poly, unsigned long);
static inline
void abase_p_8_poly_sub_ui(abase_p_8_dst_field, abase_p_8_dst_poly, abase_p_8_src_poly, unsigned long);
static inline
void abase_p_8_poly_neg(abase_p_8_dst_field, abase_p_8_dst_poly, abase_p_8_src_poly);
static inline
void abase_p_8_poly_scal_mul(abase_p_8_dst_field, abase_p_8_dst_poly, abase_p_8_src_poly, abase_p_8_src_elt);
static inline
void abase_p_8_poly_mul(abase_p_8_dst_field, abase_p_8_dst_poly, abase_p_8_src_poly, abase_p_8_src_poly);
void abase_p_8_poly_divmod(abase_p_8_dst_field, abase_p_8_dst_poly, abase_p_8_dst_poly, abase_p_8_src_poly, abase_p_8_src_poly);
void abase_p_8_poly_precomp_mod(abase_p_8_dst_field, abase_p_8_dst_poly, abase_p_8_src_poly);
void abase_p_8_poly_mod_pre(abase_p_8_dst_field, abase_p_8_dst_poly, abase_p_8_src_poly, abase_p_8_src_poly, abase_p_8_src_poly);
static inline
void abase_p_8_poly_gcd(abase_p_8_dst_field, abase_p_8_dst_poly, abase_p_8_src_poly, abase_p_8_src_poly);
static inline
void abase_p_8_poly_xgcd(abase_p_8_dst_field, abase_p_8_dst_poly, abase_p_8_dst_poly, abase_p_8_dst_poly, abase_p_8_src_poly, abase_p_8_src_poly);
static inline
void abase_p_8_poly_random(abase_p_8_dst_field, abase_p_8_dst_poly, unsigned int, gmp_randstate_t);
static inline
void abase_p_8_poly_random2(abase_p_8_dst_field, abase_p_8_dst_poly, unsigned int, gmp_randstate_t);
static inline
int abase_p_8_poly_cmp(abase_p_8_dst_field, abase_p_8_src_poly, abase_p_8_src_poly);
static inline
void abase_p_8_poly_asprint(abase_p_8_dst_field, char * *, abase_p_8_src_poly);
static inline
void abase_p_8_poly_fprint(abase_p_8_dst_field, FILE *, abase_p_8_src_poly);
static inline
void abase_p_8_poly_print(abase_p_8_dst_field, abase_p_8_src_poly);
static inline
int abase_p_8_poly_sscan(abase_p_8_dst_field, abase_p_8_dst_poly, const char *);
static inline
int abase_p_8_poly_fscan(abase_p_8_dst_field, FILE *, abase_p_8_dst_poly);
static inline
int abase_p_8_poly_scan(abase_p_8_dst_field, abase_p_8_dst_poly);

/* Functions related to SIMD operation */
/* *simd_gfp::code_for_groupsize */
#define abase_p_8_groupsize(K)	1
/* *simd_gfp::code_for_offset */
#define abase_p_8_offset(K, n)	n /* TO BE DEPRECATED */
/* *simd_gfp::code_for_stride */
#define abase_p_8_stride(K)	1 /* TO BE DEPRECATED */
static inline
void abase_p_8_set_ui_at(abase_p_8_dst_field, abase_p_8_dst_elt, int, unsigned long);
static inline
void abase_p_8_set_ui_all(abase_p_8_dst_field, abase_p_8_dst_elt, unsigned long);
static inline
void abase_p_8_elt_ur_set_ui_at(abase_p_8_dst_field, abase_p_8_dst_elt, int, unsigned long);
static inline
void abase_p_8_elt_ur_set_ui_all(abase_p_8_dst_field, abase_p_8_dst_elt, unsigned long);
void abase_p_8_dotprod(abase_p_8_dst_field, abase_p_8_dst_vec, abase_p_8_src_vec, abase_p_8_src_vec, unsigned int);

/* Member templates related to SIMD operation */

/* MPI interface */
void abase_p_8_mpi_ops_init(abase_p_8_dst_field);
MPI_Datatype abase_p_8_mpi_datatype(abase_p_8_dst_field);
MPI_Datatype abase_p_8_mpi_datatype_ur(abase_p_8_dst_field);
MPI_Op abase_p_8_mpi_addition_op(abase_p_8_dst_field);
MPI_Op abase_p_8_mpi_addition_op_ur(abase_p_8_dst_field);
void abase_p_8_mpi_ops_clear(abase_p_8_dst_field);

/* Object-oriented interface */
static inline
void abase_p_8_oo_field_clear(abase_vbase_ptr);
void abase_p_8_oo_field_init(abase_vbase_ptr);
#ifdef  __cplusplus
}
#endif

/* Implementations for inlines */
/* *Mpfq::gfp::field::code_for_field_init, Mpfq::gfp */
static inline
void abase_p_8_field_init(abase_p_8_dst_field k)
{
    k->p = NULL;
    k->bigmul_p = NULL;
    k->io_base = 10;
    mpz_init(k->factor);
    k->ts_info.e=0;
}

/* *Mpfq::gfp::elt::code_for_init, Mpfq::gfp */
static inline
void abase_p_8_init(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_elt * x MAYBE_UNUSED)
{
    assert(k);
    assert(k->p);
    assert(*x);
}

/* *Mpfq::gfp::elt::code_for_clear, Mpfq::gfp */
static inline
void abase_p_8_clear(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_elt * x MAYBE_UNUSED)
{
    assert(k);
    assert(*x);
}

/* *Mpfq::defaults::flatdata::code_for_set, Mpfq::gfp::elt, Mpfq::gfp */
static inline
void abase_p_8_set(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_elt r, abase_p_8_src_elt s)
{
    if (r != s) memcpy(r,s,sizeof(abase_p_8_elt));
}

/* *Mpfq::gfp::elt::code_for_set_ui, Mpfq::gfp */
static inline
void abase_p_8_set_ui(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_dst_elt r, unsigned long x)
{
    int i; 
    assert (r);
    r[0] = x;
    for (i = 1; i < 8; ++i)
        r[i] = 0;
}

/* *Mpfq::defaults::flatdata::code_for_set_zero, Mpfq::gfp::elt, Mpfq::gfp */
static inline
void abase_p_8_set_zero(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_elt r)
{
    abase_p_8_vec_set_zero(K,(abase_p_8_dst_vec)r,1);
}

/* *Mpfq::gfp::elt::code_for_get_ui, Mpfq::gfp */
static inline
unsigned long abase_p_8_get_ui(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_src_elt x)
{
    return x[0];
}

/* *Mpfq::gfp::elt::code_for_set_mpn, Mpfq::gfp */
static inline
void abase_p_8_set_mpn(abase_p_8_dst_field k, abase_p_8_dst_elt r, mp_limb_t * x, size_t n)
{
    int i;
    if (n < 8) {
        for (i = 0; i < (int)n; ++i)
            r[i] = x[i];
        for (i = n; i < 8; ++i)
            r[i] = 0;
    } else {
        mp_limb_t tmp[n-8+1];
        mpn_tdiv_qr(tmp, r, 0, x, n, k->p, 8);
    }
}

/* *Mpfq::gfp::elt::code_for_set_mpz, Mpfq::gfp */
static inline
void abase_p_8_set_mpz(abase_p_8_dst_field k, abase_p_8_dst_elt r, mpz_t z)
{
    if (z->_mp_size < 0) {
        abase_p_8_set_mpn(k, r, z->_mp_d, -z->_mp_size);
        abase_p_8_neg(k, r, r);
    } else {
        abase_p_8_set_mpn(k, r, z->_mp_d, z->_mp_size);
    }
}

/* *Mpfq::gfp::elt::code_for_get_mpn, Mpfq::gfp */
static inline
void abase_p_8_get_mpn(abase_p_8_dst_field k MAYBE_UNUSED, mp_limb_t * r, abase_p_8_src_elt x)
{
    int i; 
    assert (r);
    assert (x);
    for (i = 0; i < 8; ++i)
        r[i] = x[i];
}

/* *Mpfq::gfp::elt::code_for_get_mpz, Mpfq::gfp */
static inline
void abase_p_8_get_mpz(abase_p_8_dst_field k MAYBE_UNUSED, mpz_t z, abase_p_8_src_elt y)
{
    int i; 
    mpz_realloc2(z, 8*64);
    for (i = 0; i < 8; ++i)
        z->_mp_d[i] = y[i];
    i = 8;
    while (i>=1 && z->_mp_d[i-1] == 0)
        i--;
    z->_mp_size = i;
}

/* *Mpfq::gfp::elt::code_for_random, Mpfq::gfp */
static inline
void abase_p_8_random(abase_p_8_dst_field k, abase_p_8_dst_elt x, gmp_randstate_t state)
{
      mpz_t z;
      mpz_init(z);
      mpz_urandomb(z, state, 8 * GMP_LIMB_BITS);
      memcpy(x, z->_mp_d, 8 * sizeof(mp_limb_t));  /* UGLY */
      mpz_clear(z);
    abase_p_8_normalize(k, x);
}

/* *Mpfq::gfp::elt::code_for_random2, Mpfq::gfp */
static inline
void abase_p_8_random2(abase_p_8_dst_field k, abase_p_8_dst_elt x, gmp_randstate_t state)
{
      mpz_t z;
      mpz_init(z);
      mpz_rrandomb(z, state, 8 * GMP_LIMB_BITS);
      memcpy(x, z->_mp_d, 8 * sizeof(mp_limb_t));  /* UGLY */
      mpz_clear(z);
    abase_p_8_normalize(k, x);
}

/* *Mpfq::gfp::elt::code_for_add, Mpfq::gfp */
static inline
void abase_p_8_add(abase_p_8_dst_field k, abase_p_8_dst_elt z, abase_p_8_src_elt x, abase_p_8_src_elt y)
{
    mp_limb_t cy;
    cy = add_8(z, x, y);
    if (cy || (cmp_8(z, k->p) >= 0))
        sub_8(z, z, k->p);
}

/* *Mpfq::gfp::elt::code_for_sub, Mpfq::gfp */
static inline
void abase_p_8_sub(abase_p_8_dst_field k, abase_p_8_dst_elt z, abase_p_8_src_elt x, abase_p_8_src_elt y)
{
    mp_limb_t cy;
    cy = sub_8(z, x, y);
    if (cy) // negative result
        add_8(z, z, k->p);
}

/* *Mpfq::gfp::elt::code_for_neg, Mpfq::gfp */
static inline
void abase_p_8_neg(abase_p_8_dst_field k, abase_p_8_dst_elt z, abase_p_8_src_elt x)
{
    if (cmp_ui_8(x, 0))
        sub_8(z, k->p, x);
    else {
        int i;
        for (i = 0; i < 8; ++i)
            z[i] = 0;
        }
}

/* *Mpfq::gfp::elt::code_for_mul, Mpfq::gfp */
static inline
void abase_p_8_mul(abase_p_8_dst_field k, abase_p_8_dst_elt z, abase_p_8_src_elt x, abase_p_8_src_elt y)
{
    mp_limb_t tmp[2*8];
    mul_8(tmp, x, y);
    mod_8(z, tmp, k->p);
}

/* *Mpfq::gfp::elt::code_for_sqr, Mpfq::gfp */
static inline
void abase_p_8_sqr(abase_p_8_dst_field k, abase_p_8_dst_elt z, abase_p_8_src_elt x)
{
    mp_limb_t tmp[2*8];
    sqr_8(tmp, x);
    mod_8(z, tmp, k->p);
}

/* *Mpfq::gfp::elt::code_for_is_sqr, Mpfq::gfp */
static inline
int abase_p_8_is_sqr(abase_p_8_dst_field k, abase_p_8_src_elt x)
{
    mp_limb_t pp[8];
    abase_p_8_elt y;
    sub_ui_nc_8(pp, k->p, 1);
    rshift_8(pp, 1);
    abase_p_8_init(k, &y);
    abase_p_8_pow(k, y, x, pp, 8);
    int res = cmp_ui_8(y, 1);
    abase_p_8_clear(k, &y);
    if (res == 0)
        return 1;
    else 
        return 0;
}

/* *Mpfq::defaults::pow::code_for_pow, Mpfq::gfp::elt, Mpfq::gfp */
static inline
void abase_p_8_pow(abase_p_8_dst_field k, abase_p_8_dst_elt res, abase_p_8_src_elt r, unsigned long * x, size_t n)
{
    abase_p_8_elt u, a;
    long i, j, lead;     /* it is a signed type */
    unsigned long mask;
    
    assert (n>0);
    
    /* get the correct (i,j) position of the most significant bit in x */
    for(i = n-1; i>=0 && x[i]==0; i--)
        ;
    if (i < 0) {
        abase_p_8_set_ui(k, res, 0);
        return;
    }
    j = 64 - 1;
    mask = (1UL<<j);
    for( ; (x[i]&mask)==0 ;j--, mask>>=1)
        ;
    lead = i*64+j;      /* Ensured. */
    
    abase_p_8_init(k, &u);
    abase_p_8_init(k, &a);
    abase_p_8_set(k, a, r);
    for( ; lead > 0; lead--) {
        if (j-- == 0) {
            i--;
            j = 64-1;
            mask = (1UL<<j);
        } else {
            mask >>= 1;
        }
        if (x[i]&mask) {
            abase_p_8_sqr(k, u, a);
            abase_p_8_mul(k, a, u, r);
        } else {
            abase_p_8_sqr(k, a,a);
        }
    }
    abase_p_8_set(k, res, a);
    abase_p_8_clear(k, &u);
    abase_p_8_clear(k, &a);
}

/* *Mpfq::gfp::elt::code_for_add_ui, Mpfq::gfp */
static inline
void abase_p_8_add_ui(abase_p_8_dst_field k, abase_p_8_dst_elt z, abase_p_8_src_elt x, unsigned long y)
{
    mp_limb_t cy;
    cy = add_ui_8(z, x, y);
    if (cy || (cmp_8(z, k->p) >= 0))
        sub_8(z, z, k->p);
}

/* *Mpfq::gfp::elt::code_for_sub_ui, Mpfq::gfp */
static inline
void abase_p_8_sub_ui(abase_p_8_dst_field k, abase_p_8_dst_elt z, abase_p_8_src_elt x, unsigned long y)
{
    mp_limb_t cy;
    cy = sub_ui_8(z, x, y);
    if (cy) // negative result
        add_8(z, z, k->p);
}

/* *Mpfq::gfp::elt::code_for_mul_ui, Mpfq::gfp */
static inline
void abase_p_8_mul_ui(abase_p_8_dst_field k, abase_p_8_dst_elt z, abase_p_8_src_elt x, unsigned long y)
{
    mp_limb_t tmp[8+1], q[2];
    mul1_8(tmp,x,y);
    mpn_tdiv_qr(q, z, 0, tmp, 8+1, k->p, 8);
}

/* *Mpfq::gfp::elt::code_for_inv, Mpfq::gfp */
static inline
int abase_p_8_inv(abase_p_8_dst_field k, abase_p_8_dst_elt z, abase_p_8_src_elt x)
{
    int ret=invmod_8(z, x, k->p);
    if (!ret)
        abase_p_8_get_mpz(k, k->factor, z);
    return ret;
}

/* *Mpfq::gfp::elt::code_for_hadamard, Mpfq::gfp */
static inline
void abase_p_8_hadamard(abase_p_8_dst_field k, abase_p_8_dst_elt x, abase_p_8_dst_elt y, abase_p_8_dst_elt z, abase_p_8_dst_elt t)
{
    abase_p_8_elt tmp;
    abase_p_8_init(k, &tmp);
    abase_p_8_add(k, tmp, x, y);
    abase_p_8_sub(k, y, x, y);
    abase_p_8_set(k, x, tmp);
    abase_p_8_add(k, tmp, z, t);
    abase_p_8_sub(k, t, z, t);
    abase_p_8_set(k, z, tmp);
    abase_p_8_sub(k, tmp, x, z);
    abase_p_8_add(k, x, x, z);
    abase_p_8_add(k, z, y, t);
    abase_p_8_sub(k, t, y, t);
    abase_p_8_set(k, y, tmp);
    abase_p_8_clear(k, &tmp); 
}

/* *Mpfq::gfp::elt::code_for_elt_ur_init, Mpfq::gfp */
static inline
void abase_p_8_elt_ur_init(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_elt_ur * x MAYBE_UNUSED)
{
    assert(k);
    assert(k->p);
    assert(*x);
}

/* *Mpfq::gfp::elt::code_for_elt_ur_clear, Mpfq::gfp */
static inline
void abase_p_8_elt_ur_clear(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_elt_ur * x MAYBE_UNUSED)
{
    assert(k);
    assert(*x);
}

/* *Mpfq::gfp::elt::code_for_elt_ur_set, Mpfq::gfp */
static inline
void abase_p_8_elt_ur_set(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_dst_elt_ur z, abase_p_8_src_elt_ur x)
{
    int i;
    for (i = 0; i < 17; ++i) 
        z[i] = x[i];
}

/* *Mpfq::defaults::flatdata::code_for_elt_ur_set_elt, Mpfq::gfp::elt, Mpfq::gfp */
static inline
void abase_p_8_elt_ur_set_elt(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_elt_ur r, abase_p_8_src_elt s)
{
    memset(r, 0, sizeof(abase_p_8_elt_ur)); memcpy(r,s,sizeof(abase_p_8_elt));
}

/* *Mpfq::defaults::flatdata::code_for_elt_ur_set_zero, Mpfq::gfp::elt, Mpfq::gfp */
static inline
void abase_p_8_elt_ur_set_zero(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_elt_ur r)
{
    memset(r, 0, sizeof(abase_p_8_elt_ur));
}

/* *Mpfq::gfp::elt::code_for_elt_ur_set_ui, Mpfq::gfp */
static inline
void abase_p_8_elt_ur_set_ui(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_dst_elt_ur r, unsigned long x)
{
    int i; 
    assert (r); 
    r[0] = x;
    for (i = 1; i < 17; ++i)
        r[i] = 0;
}

/* *Mpfq::gfp::elt::code_for_elt_ur_add, Mpfq::gfp */
static inline
void abase_p_8_elt_ur_add(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_dst_elt_ur z, abase_p_8_src_elt_ur x, abase_p_8_src_elt_ur y)
{
    mpn_add_n(z, x, y, 17);
}

/* *Mpfq::gfp::elt::code_for_elt_ur_neg, Mpfq::gfp */
static inline
void abase_p_8_elt_ur_neg(abase_p_8_dst_field k, abase_p_8_dst_elt_ur z, abase_p_8_src_elt_ur x)
{
    abase_p_8_elt_ur tmp;
    abase_p_8_elt_ur_init(k, &tmp);
    int i;
    for (i = 0; i < 17; ++i) 
        tmp[i] = 0;
    mpn_sub_n(z, tmp, x, 17);
    abase_p_8_elt_ur_clear(k, &tmp);
}

/* *Mpfq::gfp::elt::code_for_elt_ur_sub, Mpfq::gfp */
static inline
void abase_p_8_elt_ur_sub(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_dst_elt_ur z, abase_p_8_src_elt_ur x, abase_p_8_src_elt_ur y)
{
    mpn_sub_n(z, x, y, 17);
}

/* *Mpfq::gfp::elt::code_for_mul_ur, Mpfq::gfp */
static inline
void abase_p_8_mul_ur(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_dst_elt_ur z, abase_p_8_src_elt x, abase_p_8_src_elt y)
{
    mul_8(z, x, y);
    int i;
    for (i = 16; i < 17; ++i) {
        z[i] = 0;
    }
}

/* *Mpfq::gfp::elt::code_for_sqr_ur, Mpfq::gfp */
static inline
void abase_p_8_sqr_ur(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_dst_elt_ur z, abase_p_8_src_elt x)
{
    sqr_8(z, x);
    int i;
    for (i = 16; i < 17; ++i) {
        z[i] = 0;
    }
}

/* *Mpfq::gfp::elt::code_for_reduce, Mpfq::gfp */
static inline
void abase_p_8_reduce(abase_p_8_dst_field k, abase_p_8_dst_elt z, abase_p_8_dst_elt_ur x)
{
    mp_limb_t q[17+1];
    if (x[17-1]>>(64-1)) {
        // negative number, add bigmul_p to make it positive before reduction
        mpn_add_n(x, x, k->bigmul_p, 17);
    }
    mpn_tdiv_qr(q, z, 0, x, 17, k->p, 8);
}

/* *Mpfq::gfp::elt::code_for_normalize, Mpfq::gfp */
static inline
void abase_p_8_normalize(abase_p_8_dst_field k, abase_p_8_dst_elt x)
{
    if (cmp_8(x,k->p)>=0) {
      mp_limb_t q[8+1];
      abase_p_8_elt r;
      mpn_tdiv_qr(q, r, 0, x, 8, k->p, 8);
      abase_p_8_set(k, x, r);
    }
}

/* *simd_gfp::code_for_addmul_si_ur */
static inline
void abase_p_8_addmul_si_ur(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_elt_ur w, abase_p_8_src_elt u, long v)
{
        abase_p_8_elt_ur s;
        abase_p_8_elt vx;
        abase_p_8_elt_ur_init(K, &s);
        abase_p_8_init(K, &vx);
        if (v>0) {
            abase_p_8_set_ui(K, vx, v);
            abase_p_8_mul_ur(K, s, u, vx);
            abase_p_8_elt_ur_add(K, w, w, s);
        } else {
            abase_p_8_set_ui(K, vx, -v);
            abase_p_8_mul_ur(K, s, u, vx);
            abase_p_8_elt_ur_sub(K, w, w, s);
        }
        abase_p_8_clear(K, &vx);
        abase_p_8_elt_ur_clear(K, &s);
}

/* *Mpfq::gfp::elt::code_for_cmp, Mpfq::gfp */
static inline
int abase_p_8_cmp(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_src_elt x, abase_p_8_src_elt y)
{
    return cmp_8(x,y);
}

/* *Mpfq::gfp::elt::code_for_cmp_ui, Mpfq::gfp */
static inline
int abase_p_8_cmp_ui(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_src_elt x, unsigned long y)
{
    return cmp_ui_8(x,y);
}

/* *Mpfq::defaults::flatdata::code_for_is_zero, Mpfq::gfp::elt, Mpfq::gfp */
static inline
int abase_p_8_is_zero(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_src_elt r)
{
        unsigned int i;
        for(i = 0 ; i < sizeof(abase_p_8_elt)/sizeof(r[0]) ; i++) {
            if (r[i]) return 0;
        }
        return 1;
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_set, Mpfq::defaults::flatdata, Mpfq::gfp::elt, Mpfq::gfp */
static inline
void abase_p_8_vec_set(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec r, abase_p_8_src_vec s, unsigned int n)
{
    if (r != s) memmove(r, s, n*sizeof(abase_p_8_elt));
}

/* *Mpfq::defaults::vec::getset::code_for_vec_set_partial, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void abase_p_8_vec_set_partial(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_dst_vec w, abase_p_8_src_vec u, unsigned int bw, unsigned int bu, unsigned int l)
{
    unsigned int i;
    for(i = 0; i < l; ++i)
        abase_p_8_set(k, w[bw+i], u[bu+i]);
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_set_zero, Mpfq::defaults::flatdata, Mpfq::gfp::elt, Mpfq::gfp */
static inline
void abase_p_8_vec_set_zero(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec r, unsigned int n)
{
    memset(r, 0, n*sizeof(abase_p_8_elt));
}

/* *Mpfq::defaults::vec::getset::code_for_vec_setcoef, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void abase_p_8_vec_setcoef(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec w, abase_p_8_src_elt x, unsigned int i)
{
    abase_p_8_set(K, w[i], x);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_setcoef_ui, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void abase_p_8_vec_setcoef_ui(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec w, unsigned long x, unsigned int i)
{
    abase_p_8_set_ui(K, w[i], x);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_getcoef, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void abase_p_8_vec_getcoef(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_elt x, abase_p_8_src_vec w, unsigned int i)
{
    abase_p_8_set(K, x, w[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_add, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void abase_p_8_vec_add(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec w, abase_p_8_src_vec u, abase_p_8_src_vec v, unsigned int n)
{
        unsigned int i;
    for(i = 0; i < n; i+=1)
        abase_p_8_add(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_neg, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void abase_p_8_vec_neg(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec w, abase_p_8_src_vec u, unsigned int n)
{
        unsigned int i;
    for(i = 0; i < n; ++i)
        abase_p_8_neg(K, w[i], u[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_rev, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void abase_p_8_vec_rev(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec w, abase_p_8_src_vec u, unsigned int n)
{
    unsigned int nn = n >> 1;
    abase_p_8_elt tmp[1];
    abase_p_8_init(K, tmp);
    unsigned int i;
    for(i = 0; i < nn; ++i) {
        abase_p_8_set(K, tmp[0], u[i]);
        abase_p_8_set(K, w[i], u[n-1-i]);
        abase_p_8_set(K, w[n-1-i], tmp[0]);
    }
    if (n & 1)
        abase_p_8_set(K, w[nn], u[nn]);
    abase_p_8_clear(K, tmp);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_sub, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void abase_p_8_vec_sub(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec w, abase_p_8_src_vec u, abase_p_8_src_vec v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; ++i)
        abase_p_8_sub(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::mul::code_for_vec_scal_mul, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void abase_p_8_vec_scal_mul(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec w, abase_p_8_src_vec u, abase_p_8_src_elt x, unsigned int n)
{
        unsigned int i;
    for(i = 0; i < n; i+=1)
        abase_p_8_mul(K, w[i], u[i], x);
}

/* *Mpfq::defaults::vec::conv::code_for_vec_conv, Mpfq::gfp */
static inline
void abase_p_8_vec_conv(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec w, abase_p_8_src_vec u, unsigned int n, abase_p_8_src_vec v, unsigned int m)
{
    abase_p_8_vec_ur tmp;
    abase_p_8_vec_ur_init(K, &tmp, m+n-1);
    abase_p_8_vec_conv_ur(K, tmp, u, n, v, m);
    abase_p_8_vec_reduce(K, w, tmp, m+n-1);
    abase_p_8_vec_ur_clear(K, &tmp, m+n-1);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_random, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void abase_p_8_vec_random(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec w, unsigned int n, gmp_randstate_t state)
{
    unsigned int i;
    for(i = 0; i < n; ++i)
        abase_p_8_random(K, w[i], state);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_random2, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void abase_p_8_vec_random2(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec w, unsigned int n, gmp_randstate_t state)
{
    unsigned int i;
    for(i = 0; i < n; ++i)
        abase_p_8_random2(K, w[i],state);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_cmp, Mpfq::defaults::vec, Mpfq::gfp */
static inline
int abase_p_8_vec_cmp(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_src_vec u, abase_p_8_src_vec v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; ++i) {
        int ret = abase_p_8_cmp(K, u[i], v[i]);
        if (ret != 0)
            return ret;
    }
    return 0;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_is_zero, Mpfq::defaults::vec, Mpfq::gfp */
static inline
int abase_p_8_vec_is_zero(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_src_vec r, unsigned int n)
{
    unsigned int i;
    for(i = 0 ; i < n ; i+=1) {
        if (!abase_p_8_is_zero(K,r[i])) return 0;
    }
    return 1;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_subvec, Mpfq::defaults::vec, Mpfq::gfp */
static inline
abase_p_8_dst_vec abase_p_8_vec_subvec(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec v, int i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_subvec_const, Mpfq::defaults::vec, Mpfq::gfp */
static inline
abase_p_8_src_vec abase_p_8_vec_subvec_const(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_src_vec v, int i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_coeff_ptr, Mpfq::defaults::vec, Mpfq::gfp */
static inline
abase_p_8_dst_elt abase_p_8_vec_coeff_ptr(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec v, int i)
{
    return v[i];
}

/* *Mpfq::defaults::vec::getset::code_for_vec_coeff_ptr_const, Mpfq::defaults::vec, Mpfq::gfp */
static inline
abase_p_8_src_elt abase_p_8_vec_coeff_ptr_const(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_src_vec v, int i)
{
    return v[i];
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_ur_set_zero, Mpfq::defaults::flatdata, Mpfq::gfp::elt, Mpfq::gfp */
static inline
void abase_p_8_vec_ur_set_zero(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec_ur r, unsigned int n)
{
    memset(r, 0, n*sizeof(abase_p_8_elt_ur));
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_set_vec, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void abase_p_8_vec_ur_set_vec(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec_ur w, abase_p_8_src_vec u, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        abase_p_8_elt_ur_set_elt(K, w[i], u[i]);
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_ur_set, Mpfq::defaults::flatdata, Mpfq::gfp::elt, Mpfq::gfp */
static inline
void abase_p_8_vec_ur_set(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec_ur r, abase_p_8_src_vec_ur s, unsigned int n)
{
    if (r != s) memmove(r, s, n*sizeof(abase_p_8_elt_ur));
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_setcoef, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void abase_p_8_vec_ur_setcoef(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec_ur w, abase_p_8_src_elt_ur x, unsigned int i)
{
    abase_p_8_elt_ur_set(K, w[i], x);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_getcoef, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void abase_p_8_vec_ur_getcoef(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_elt_ur x, abase_p_8_src_vec_ur w, unsigned int i)
{
    abase_p_8_elt_ur_set(K, x, w[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_add, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void abase_p_8_vec_ur_add(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec_ur w, abase_p_8_src_vec_ur u, abase_p_8_src_vec_ur v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        abase_p_8_elt_ur_add(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_sub, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void abase_p_8_vec_ur_sub(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec_ur w, abase_p_8_src_vec_ur u, abase_p_8_src_vec_ur v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        abase_p_8_elt_ur_sub(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_neg, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void abase_p_8_vec_ur_neg(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec_ur w, abase_p_8_src_vec_ur u, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; ++i)
        abase_p_8_elt_ur_neg(K, w[i], u[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_rev, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void abase_p_8_vec_ur_rev(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec_ur w, abase_p_8_src_vec_ur u, unsigned int n)
{
    unsigned int nn = n >> 1;
    abase_p_8_elt_ur tmp[1];
    abase_p_8_elt_ur_init(K, tmp);
    unsigned int i;
    for(i = 0; i < nn; ++i) {
        abase_p_8_elt_ur_set(K, tmp[0], u[i]);
        abase_p_8_elt_ur_set(K, w[i], u[n-1-i]);
        abase_p_8_elt_ur_set(K, w[n-1-i], tmp[0]);
    }
    if (n & 1)
        abase_p_8_elt_ur_set(K, w[nn], u[nn]);
    abase_p_8_elt_ur_clear(K, tmp);
}

/* *Mpfq::defaults::vec::mul::code_for_vec_scal_mul_ur, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void abase_p_8_vec_scal_mul_ur(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec_ur w, abase_p_8_src_vec u, abase_p_8_src_elt x, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        abase_p_8_mul_ur(K, w[i], u[i], x);
}

/* *Mpfq::defaults::vec::conv::code_for_vec_conv_ur, Mpfq::gfp */
static inline
void abase_p_8_vec_conv_ur_n(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec_ur w, abase_p_8_src_vec u, abase_p_8_src_vec v, unsigned int n)
{
    if (n == 0)
        return;
    if (n == 1) {
        abase_p_8_mul_ur(K, w[0], u[0], v[0]);
        return;
    }
    if (n == 2) {  // Kara 2
        abase_p_8_elt t1, t2;
        abase_p_8_init(K, &t1);
        abase_p_8_init(K, &t2);
        abase_p_8_mul_ur(K, w[0], u[0], v[0]);
        abase_p_8_mul_ur(K, w[2], u[1], v[1]);
        abase_p_8_add(K, t1, u[0], u[1]);
        abase_p_8_add(K, t2, v[0], v[1]);
        abase_p_8_mul_ur(K, w[1], t1, t2);
        abase_p_8_elt_ur_sub(K, w[1], w[1], w[0]);
        abase_p_8_elt_ur_sub(K, w[1], w[1], w[2]);
        abase_p_8_clear(K, &t1);
        abase_p_8_clear(K, &t2);
        return;
    }
    if (n == 3) {  // do it in 6
        abase_p_8_elt t1, t2;
        abase_p_8_elt_ur s;
        abase_p_8_init(K, &t1);
        abase_p_8_init(K, &t2);
        abase_p_8_elt_ur_init(K, &s);
        // a0*b0*(1 - X)
        abase_p_8_mul_ur(K, w[0], u[0], v[0]);
        abase_p_8_elt_ur_neg(K, w[1], w[0]);
        // a1*b1*(-X + 2*X^2 - X^3)
        abase_p_8_mul_ur(K, w[2], u[1], v[1]);
        abase_p_8_elt_ur_neg(K, w[3], w[2]);
        abase_p_8_elt_ur_add(K, w[2], w[2], w[2]);
        abase_p_8_elt_ur_add(K, w[1], w[1], w[3]);
        // a2*b2*(-X^3+X^4)
        abase_p_8_mul_ur(K, w[4], u[2], v[2]);
        abase_p_8_elt_ur_sub(K, w[3], w[3], w[4]);
        // (a0+a1)*(b0+b1)*(X - X^2)
        abase_p_8_add(K, t1, u[0], u[1]);
        abase_p_8_add(K, t2, v[0], v[1]);
        abase_p_8_mul_ur(K, s, t1, t2);
        abase_p_8_elt_ur_add(K, w[1], w[1], s);
        abase_p_8_elt_ur_sub(K, w[2], w[2], s);
        // (a1+a2)*(b1+b2)*(X^3 - X^2)
        abase_p_8_add(K, t1, u[1], u[2]);
        abase_p_8_add(K, t2, v[1], v[2]);
        abase_p_8_mul_ur(K, s, t1, t2);
        abase_p_8_elt_ur_add(K, w[3], w[3], s);
        abase_p_8_elt_ur_sub(K, w[2], w[2], s);
        // (a0+a1+a2)*(b0+b1+b2)* X^2
        abase_p_8_add(K, t1, u[0], t1);
        abase_p_8_add(K, t2, v[0], t2);
        abase_p_8_mul_ur(K, s, t1, t2);
        abase_p_8_elt_ur_add(K, w[2], w[2], s);
        return;
    }
    unsigned int n0, n1;
    n0 = n / 2;
    n1 = n - n0;
    abase_p_8_vec_conv_ur_n(K, w, u, v, n0);
    abase_p_8_vec_conv_ur_n(K, w + 2*n0, u + n0, v + n0, n1);
    abase_p_8_elt_ur_set_ui(K, w[2*n0-1], 0);
    
    abase_p_8_vec tmpu, tmpv;
    abase_p_8_vec_ur tmpw;
    abase_p_8_vec_init(K, &tmpu, n1);
    abase_p_8_vec_init(K, &tmpv, n1);
    abase_p_8_vec_ur_init(K, &tmpw, 2*n1-1);
    
    abase_p_8_vec_set(K, tmpu, u, n0);
    if (n1 != n0) 
        abase_p_8_set_ui(K, tmpu[n0], 0);
    abase_p_8_vec_add(K, tmpu, tmpu, u+n0, n1);
    abase_p_8_vec_set(K, tmpv, v, n0);
    if (n1 != n0) 
        abase_p_8_set_ui(K, tmpv[n0], 0);
    abase_p_8_vec_add(K, tmpv, tmpv, v+n0, n1);
    abase_p_8_vec_conv_ur_n(K, tmpw, tmpu, tmpv, n1);
    abase_p_8_vec_ur_sub(K, tmpw, tmpw, w, 2*n0-1);
    abase_p_8_vec_ur_sub(K, tmpw, tmpw, w + 2*n0, 2*n1-1);
    abase_p_8_vec_ur_add(K, w + n0, w + n0, tmpw, 2*n1-1);
    
    abase_p_8_vec_clear(K, &tmpu, n1);
    abase_p_8_vec_clear(K, &tmpv, n1);
    abase_p_8_vec_ur_clear(K, &tmpw, 2*n1-1);
    return;
}

/* *Mpfq::defaults::vec::conv::code_for_vec_conv_ur, Mpfq::gfp */
static inline
void abase_p_8_vec_conv_ur(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec_ur w, abase_p_8_src_vec u, unsigned int n, abase_p_8_src_vec v, unsigned int m)
{
    if ((n > 1) && (m > 1) && (n+m > 15)) {
        abase_p_8_vec_conv_ur_ks(K, w, u, n, v, m);
        return;
    }
    if (n == m) {
        abase_p_8_vec_conv_ur_n(K, w, u, v, n);
        return;
    }
    unsigned int i, j MAYBE_UNUSED, k;
    abase_p_8_elt_ur acc, z;
    abase_p_8_elt_ur_init(K, &acc);
    abase_p_8_elt_ur_init(K, &z);
    // swap pointers to have n <= m
    abase_p_8_src_vec uu, vv;
    if (n <= m) {
        uu = u; vv = v;
    } else {
        uu = v; vv = u;
        unsigned int tmp = n;
        n = m; m = tmp;
    }
    for(k = 0; k < n; ++k) {
        abase_p_8_mul_ur(K, acc, uu[0], vv[k]);
        for(i = 1; i <= k; ++i) {
            abase_p_8_mul_ur(K, z, uu[i], vv[k-i]);
            abase_p_8_elt_ur_add(K, acc, acc, z);
        }
        abase_p_8_elt_ur_set(K, w[k], acc);
    }
    for(k = n; k < m; ++k) {
        abase_p_8_mul_ur(K, acc, uu[0], vv[k]);
        for(i = 1; i < n; ++i) {
            abase_p_8_mul_ur(K, z, uu[i], vv[k-i]);
            abase_p_8_elt_ur_add(K, acc, acc, z);
        }
        abase_p_8_elt_ur_set(K, w[k], acc);
    }
    for(k = m; k < n+m-1; ++k) {
        abase_p_8_mul_ur(K, acc, uu[k-m+1], vv[m-1]);
        for(i = k-m+2; i < n; ++i) {
            abase_p_8_mul_ur(K, z, uu[i], vv[k-i]);
            abase_p_8_elt_ur_add(K, acc, acc, z);
        }
        abase_p_8_elt_ur_set(K, w[k], acc);
    }
    abase_p_8_elt_ur_clear(K, &acc);
    abase_p_8_elt_ur_clear(K, &z);
}

/* *Mpfq::defaults::vec::mul::code_for_vec_reduce, Mpfq::defaults::vec, Mpfq::gfp */
static inline
void abase_p_8_vec_reduce(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec w, abase_p_8_dst_vec_ur u, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        abase_p_8_reduce(K, w[i], u[i]);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_subvec, Mpfq::defaults::vec, Mpfq::gfp */
static inline
abase_p_8_dst_vec_ur abase_p_8_vec_ur_subvec(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec_ur v, int i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_subvec_const, Mpfq::defaults::vec, Mpfq::gfp */
static inline
abase_p_8_src_vec_ur abase_p_8_vec_ur_subvec_const(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_src_vec_ur v, int i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_coeff_ptr, Mpfq::defaults::vec, Mpfq::gfp */
static inline
abase_p_8_dst_elt abase_p_8_vec_ur_coeff_ptr(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_vec_ur v, int i)
{
    return v[i];
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_coeff_ptr_const, Mpfq::defaults::vec, Mpfq::gfp */
static inline
abase_p_8_src_elt abase_p_8_vec_ur_coeff_ptr_const(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_src_vec_ur v, int i)
{
    return v[i];
}

/* *Mpfq::defaults::poly::code_for_poly_init, Mpfq::gfp */
static inline
void abase_p_8_poly_init(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_poly p, unsigned int n)
{
    abase_p_8_vec_init(k, &(p->c), n);
    p->alloc=n;
    p->size=0;
}

/* *Mpfq::defaults::poly::code_for_poly_clear, Mpfq::gfp */
static inline
void abase_p_8_poly_clear(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_poly p)
{
    abase_p_8_vec_clear(k, &(p->c), p->alloc);
}

/* *Mpfq::defaults::poly::code_for_poly_set, Mpfq::gfp */
static inline
void abase_p_8_poly_set(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_dst_poly w, abase_p_8_src_poly u)
{
    if (w->alloc < u->size) {
        abase_p_8_vec_reinit(k, &(w->c), w->alloc, u->size);
        w->alloc = u->size;
    }
    abase_p_8_vec_set(k, w->c, u->c, u->size);
    w->size = u->size;
}

/* *Mpfq::defaults::poly::code_for_poly_setcoef, Mpfq::gfp */
static inline
void abase_p_8_poly_setcoef(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_dst_poly w, abase_p_8_src_elt x, unsigned int i)
{
    unsigned long j;
    if (w->alloc < (i+1)) {
        abase_p_8_vec_reinit(k, &(w->c), w->alloc, i+1);
        w->alloc = i+1;
    }
    abase_p_8_vec_setcoef(k, w->c, x, i);
    if (w->size < (i+1)) {
        for (j = w->size; j < i; ++j) {
            abase_p_8_vec_setcoef_ui(k, w->c, 0, j);
        }  
        w->size = i+1;
    }
}

/* *Mpfq::defaults::poly::code_for_poly_setcoef_ui, Mpfq::gfp */
static inline
void abase_p_8_poly_setcoef_ui(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_dst_poly w, unsigned long x, unsigned int i)
{
    unsigned long j;
    if (w->alloc < (i+1)) {
        abase_p_8_vec_reinit(k, &(w->c), w->alloc, i+1);
        w->alloc = i+1;
    }
    abase_p_8_vec_setcoef_ui(k, w->c, x, i);
    if (w->size < (i+1)) {
        for (j = w->size; j < i; ++j) {
            abase_p_8_vec_setcoef_ui(k, w->c, 0, j);
        }  
        w->size = i+1;
    }
}

/* *Mpfq::defaults::poly::code_for_poly_getcoef, Mpfq::gfp */
static inline
void abase_p_8_poly_getcoef(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_dst_elt x, abase_p_8_src_poly w, unsigned int i)
{
    if (w->size < (i+1)) {
       abase_p_8_set_ui(k,x,0);
    } else {
       abase_p_8_vec_getcoef(k, x, w->c, i);
    }
}

/* *Mpfq::defaults::poly::code_for_poly_deg, Mpfq::gfp */
static inline
int abase_p_8_poly_deg(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_src_poly w)
{
    if (w->size == 0)
        return -1;
    int deg = w->size-1;
    abase_p_8_elt temp;
    abase_p_8_init(K, &temp);
    abase_p_8_vec_getcoef(K, temp, w->c, deg);
    int comp=abase_p_8_cmp_ui(K, temp, 0);
    while ((deg >= 0) && (comp == 0)){
        deg--;
        if (deg>=0) {
           abase_p_8_vec_getcoef(K, temp, w->c, deg);
           comp=abase_p_8_cmp_ui(K, temp, 0);
        }
    }
    abase_p_8_clear(K, &temp);
    return deg;
}

/* *Mpfq::defaults::poly::code_for_poly_add, Mpfq::gfp */
static inline
void abase_p_8_poly_add(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_dst_poly w, abase_p_8_src_poly u, abase_p_8_src_poly v)
{
    unsigned int minsize MAYBE_UNUSED = MIN(u->size, v->size);
    unsigned int maxsize MAYBE_UNUSED = MAX(u->size, v->size);
    if (w->alloc < maxsize) {
        abase_p_8_vec_reinit(k, &(w->c), w->alloc, maxsize);
        w->alloc = maxsize;
    }
    if (u->size <= v->size) {
        abase_p_8_vec_add(k, w->c, u->c, v->c, u->size);
        abase_p_8_vec_set_partial(k, (w->c), (v->c), u->size, u->size, v->size-u->size);
    } else {
        abase_p_8_vec_add(k, w->c, u->c, v->c, v->size);
        abase_p_8_vec_set_partial(k, (w->c), (u->c), v->size, v->size, u->size-v->size);
    }
    w->size=maxsize;
    unsigned int wdeg = abase_p_8_poly_deg(k, w);
    w->size=wdeg+1;
}

/* *Mpfq::defaults::poly::code_for_poly_sub, Mpfq::gfp */
static inline
void abase_p_8_poly_sub(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_dst_poly w, abase_p_8_src_poly u, abase_p_8_src_poly v)
{
    unsigned int minsize MAYBE_UNUSED = MIN(u->size, v->size);
    unsigned int maxsize MAYBE_UNUSED = MAX(u->size, v->size);
    if (w->alloc < maxsize) {
        abase_p_8_vec_reinit(k, &(w->c), w->alloc, maxsize);
        w->alloc = maxsize;
    }
    if (u->size <= v->size) {
        abase_p_8_vec_sub(k, w->c, u->c, v->c, u->size);
        unsigned int i;
        abase_p_8_elt temp;
        abase_p_8_init(k, &temp);
        for (i = u->size; i< v->size; ++i) {
            abase_p_8_poly_getcoef(k, temp, v, i);
            abase_p_8_neg(k, temp, temp);
            abase_p_8_poly_setcoef(k, w, temp, i);
        }
    } else {
        abase_p_8_vec_sub(k, w->c, u->c, v->c, v->size);
        abase_p_8_vec_set_partial(k, (w->c), (u->c), v->size, v->size, u->size-v->size);
    }
    w->size=maxsize;
    //abase_p_8_poly_neg(k, w, v);
    //abase_p_8_poly_add(k, w, u, w);
    unsigned int wdeg = abase_p_8_poly_deg(k, w);
    w->size=wdeg+1;
}

/* *Mpfq::defaults::poly::code_for_poly_add_ui, Mpfq::gfp */
static inline
void abase_p_8_poly_add_ui(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_dst_poly w, abase_p_8_src_poly u, unsigned long x)
{
    if (u->size == 0) {
        if (x == 0) {
            w->size = 0;
            return;
        }
        if (w->alloc == 0) {
            abase_p_8_vec_reinit(k, &(w->c), w->alloc, 1);
            w->alloc = 1;
        }
        w->size = 1;
        abase_p_8_vec_setcoef_ui(k, w->c, x, 0);
        return;
    }
    if (w->alloc < u->size) {
        abase_p_8_vec_reinit(k, &(w->c), w->alloc, u->size);
        w->alloc = u->size;
    }
    abase_p_8_elt temp;
    abase_p_8_init(k, &temp);
    abase_p_8_poly_getcoef(k, temp, u, 0);
    abase_p_8_add_ui(k, temp, temp, x);
    abase_p_8_vec_setcoef(k, w->c, temp, 0);
    abase_p_8_vec_set_partial(k, w->c, u->c, 1, 1, u->size-1);
    w->size=u->size;
    abase_p_8_clear(k, &temp);
}

/* *Mpfq::defaults::poly::code_for_poly_sub_ui, Mpfq::gfp */
static inline
void abase_p_8_poly_sub_ui(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_dst_poly w, abase_p_8_src_poly u, unsigned long x)
{
    if (u->size == 0) {
        if (x == 0) {
            w->size = 0;
            return;
        }
        if (w->alloc == 0) {
            abase_p_8_vec_reinit(k, &(w->c), w->alloc, 1);
            w->alloc = 1;
        }
        w->size = 1;
        abase_p_8_elt temp;
        abase_p_8_init(k, &temp);
        abase_p_8_set_ui(k, temp, x);
        abase_p_8_neg(k, temp, temp);
        abase_p_8_vec_setcoef(k, w->c, temp, 0);
        return;
    }
    if (w->alloc < u->size) {
        abase_p_8_vec_reinit(k, &(w->c), w->alloc, u->size);
        w->alloc = u->size;
    }
    abase_p_8_elt temp;
    abase_p_8_init(k, &temp);
    abase_p_8_poly_getcoef(k, temp, u, 0);
    abase_p_8_sub_ui(k, temp, temp, x);
    abase_p_8_vec_setcoef(k, w->c, temp, 0);
    abase_p_8_vec_set_partial(k, w->c, u->c, 1, 1, u->size-1);
    w->size=u->size;
    abase_p_8_clear(k, &temp);
}

/* *Mpfq::defaults::poly::code_for_poly_neg, Mpfq::gfp */
static inline
void abase_p_8_poly_neg(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_dst_poly w, abase_p_8_src_poly u)
{
    if (w->alloc < u->size) {
        abase_p_8_vec_reinit(k, &(w->c), w->alloc, u->size);
        w->alloc = u->size;
    }
    abase_p_8_vec_neg(k, w->c, u->c, u->size);
    w->size = u->size;
}

/* *Mpfq::defaults::poly::code_for_poly_scal_mul, Mpfq::gfp */
static inline
void abase_p_8_poly_scal_mul(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_dst_poly w, abase_p_8_src_poly u, abase_p_8_src_elt x)
{
    if (abase_p_8_cmp_ui(k, x, 0) == 0) {
        w->size = 0;
        return;
    }
    unsigned int n = u->size;
    if (w->alloc < n) {
        abase_p_8_vec_reinit(k, &(w->c), w->alloc, n);
        w->alloc = n;
    }
    abase_p_8_vec_scal_mul(k, w->c, u->c, x, n);
    w->size=n;
}

/* *Mpfq::defaults::poly::code_for_poly_mul, Mpfq::gfp */
static inline
void abase_p_8_poly_mul(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_poly w, abase_p_8_src_poly u, abase_p_8_src_poly v)
{
    unsigned int usize = abase_p_8_poly_deg(K, u)+1;
    unsigned int vsize = abase_p_8_poly_deg(K, v)+1;
    if ((usize == 0) || (vsize == 0)) {
        w->size = 0;
        return;
    }
    unsigned int wsize = usize + vsize - 1;
    if (w->alloc < wsize) {
        abase_p_8_vec_reinit(K, &(w->c), w->alloc, wsize);
        w->alloc = wsize;
    }
    abase_p_8_vec_conv(K, w->c, u->c, usize, v->c, vsize);
    w->size=wsize;
}

/* *Mpfq::defaults::polygcd::code_for_poly_gcd, Mpfq::defaults::poly, Mpfq::gfp */
static inline
void abase_p_8_poly_gcd(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_dst_poly g, abase_p_8_src_poly a0, abase_p_8_src_poly b0)
{
    abase_p_8_poly a,b,q,r;
    int da0=abase_p_8_poly_deg(k,a0), db0=abase_p_8_poly_deg(k,b0);
    if (db0==-1)
     abase_p_8_poly_set(k,g,a0);
    else {
     abase_p_8_poly_init(k,a,da0+1);
     abase_p_8_poly_init(k,b,db0+1);
     abase_p_8_poly_init(k,q,1);
     abase_p_8_poly_init(k,r,db0);
     abase_p_8_poly_set(k,a,a0);
     abase_p_8_poly_set(k,b,b0);
     while (abase_p_8_poly_deg(k,b)>=0) {
      abase_p_8_poly_divmod(k,q,r,a,b);
      abase_p_8_poly_set(k,a,b);
      abase_p_8_poly_set(k,b,r); 
     }
     abase_p_8_poly_setmonic(k,g,a);
    abase_p_8_poly_clear(k,a);
    abase_p_8_poly_clear(k,b);
    abase_p_8_poly_clear(k,q);
    abase_p_8_poly_clear(k,r);
    }
}

/* *Mpfq::defaults::polygcd::code_for_poly_xgcd, Mpfq::defaults::poly, Mpfq::gfp */
static inline
void abase_p_8_poly_xgcd(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_dst_poly g, abase_p_8_dst_poly u0, abase_p_8_dst_poly v0, abase_p_8_src_poly a0, abase_p_8_src_poly b0)
{
    abase_p_8_poly a,b,u,v,w,x,q,r;
    abase_p_8_elt c;
    abase_p_8_init(k,&c);
    int da0=abase_p_8_poly_deg(k,a0), db0=abase_p_8_poly_deg(k,b0), dega;
    if (db0==-1) {
     if (da0==-1) {
      abase_p_8_poly_set(k,u0,a0);
      abase_p_8_poly_set(k,v0,b0);
      abase_p_8_poly_set(k,g,a0);
     } else {
      abase_p_8_poly_getcoef(k,c,a0,da0);
      abase_p_8_inv(k,c,c);
      abase_p_8_poly_scal_mul(k,g,a0,c);
      abase_p_8_poly_set(k,v0,b0);
      abase_p_8_poly_set(k,u0,b0);
      abase_p_8_poly_setcoef(k,u0,c,0);
     }
    }
    else {
     abase_p_8_poly_init(k,a,da0+1);
     abase_p_8_poly_init(k,b,db0+1);
     abase_p_8_poly_init(k,q,1);
     abase_p_8_poly_init(k,r,db0);
     abase_p_8_poly_set(k,a,a0);
     abase_p_8_poly_set(k,b,b0);
     abase_p_8_poly_init(k,u,1);
     abase_p_8_poly_init(k,v,1);
     abase_p_8_poly_init(k,w,1);
     abase_p_8_poly_init(k,x,1);
     abase_p_8_poly_setcoef_ui(k,u,1,0);
     abase_p_8_poly_setcoef_ui(k,x,1,0);
     /* u*a_initial + v*b_initial = a */
     /* w*a_initial + x*b_initial = b */
     while (abase_p_8_poly_deg(k,b)>=0) {
      abase_p_8_poly_divmod(k,q,r,a,b);
      abase_p_8_poly_set(k,a,b);  /* a,b <- b,a-qb=r */
      abase_p_8_poly_set(k,b,r);
      abase_p_8_poly_mul(k,r,q,w);
      abase_p_8_poly_sub(k,r,u,r);
      abase_p_8_poly_set(k,u,w);   /* u,w <- w,u-qw */
      abase_p_8_poly_set(k,w,r);
      abase_p_8_poly_mul(k,r,q,x); /* v,x <- x,v-qx */
      abase_p_8_poly_sub(k,r,v,r);
      abase_p_8_poly_set(k,v,x);
      abase_p_8_poly_set(k,x,r);
     }
     dega=abase_p_8_poly_deg(k,a);
     abase_p_8_poly_getcoef(k,c,a,dega);
     abase_p_8_inv(k,c,c);
     abase_p_8_poly_scal_mul(k,g,a,c);
     abase_p_8_poly_scal_mul(k,u0,u,c);
     abase_p_8_poly_scal_mul(k,v0,v,c);
     abase_p_8_poly_clear(k,a);
     abase_p_8_poly_clear(k,b);
     abase_p_8_poly_clear(k,u);
     abase_p_8_poly_clear(k,v);
     abase_p_8_poly_clear(k,w);
     abase_p_8_poly_clear(k,x);
     abase_p_8_poly_clear(k,q);
     abase_p_8_poly_clear(k,r);
    }
    abase_p_8_clear(k,&c);
}

/* *Mpfq::defaults::poly::code_for_poly_random, Mpfq::gfp */
static inline
void abase_p_8_poly_random(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_dst_poly w, unsigned int n, gmp_randstate_t state)
{
    n++;
    if (w->alloc < n) {
        abase_p_8_vec_reinit(k, &(w->c), w->alloc, n);
        w->alloc = n;
    }
    abase_p_8_vec_random(k, w->c, n,state);
    w->size=n;
    int wdeg = abase_p_8_poly_deg(k, w);
    w->size=wdeg+1;
}

/* *Mpfq::defaults::poly::code_for_poly_random2, Mpfq::gfp */
static inline
void abase_p_8_poly_random2(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_dst_poly w, unsigned int n, gmp_randstate_t state)
{
    n++;
    if (w->alloc < n) {
        abase_p_8_vec_reinit(k, &(w->c), w->alloc, n);
        w->alloc = n;
    }
    abase_p_8_vec_random2(k, w->c, n,state);
    w->size=n;
    int wdeg = abase_p_8_poly_deg(k, w);
    w->size=wdeg+1;
}

/* *Mpfq::defaults::poly::code_for_poly_cmp, Mpfq::gfp */
static inline
int abase_p_8_poly_cmp(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_src_poly u, abase_p_8_src_poly v)
{
    if (u->size != v->size)
        return (int)(u->size) - (int)(v->size);
    else
        return abase_p_8_vec_cmp(k, u->c, v->c, u->size);
}

/* *Mpfq::defaults::poly::code_for_poly_asprint, Mpfq::gfp */
static inline
void abase_p_8_poly_asprint(abase_p_8_dst_field k MAYBE_UNUSED, char * * pstr, abase_p_8_src_poly w)
{
    abase_p_8_vec_asprint(k, pstr, w->c, w->size);
}

/* *Mpfq::defaults::poly::code_for_poly_fprint, Mpfq::gfp */
static inline
void abase_p_8_poly_fprint(abase_p_8_dst_field k MAYBE_UNUSED, FILE * file, abase_p_8_src_poly w)
{
    abase_p_8_vec_fprint(k, file, w->c, w->size);
}

/* *Mpfq::defaults::poly::code_for_poly_print, Mpfq::gfp */
static inline
void abase_p_8_poly_print(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_src_poly w)
{
    abase_p_8_vec_print(k, w->c, w->size);
}

/* *Mpfq::defaults::poly::code_for_poly_sscan, Mpfq::gfp */
static inline
int abase_p_8_poly_sscan(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_dst_poly w, const char * str)
{
    int ret;
    ret = abase_p_8_vec_sscan(k, &(w->c), &(w->alloc), str);
    w->size = w->alloc;
    return ret;
}

/* *Mpfq::defaults::poly::code_for_poly_fscan, Mpfq::gfp */
static inline
int abase_p_8_poly_fscan(abase_p_8_dst_field k MAYBE_UNUSED, FILE * file, abase_p_8_dst_poly w)
{
    int ret;
    ret = abase_p_8_vec_fscan(k, file, &(w->c), &(w->alloc));
    w->size = w->alloc;
    return ret;
}

/* *Mpfq::defaults::poly::code_for_poly_scan, Mpfq::gfp */
static inline
int abase_p_8_poly_scan(abase_p_8_dst_field k MAYBE_UNUSED, abase_p_8_dst_poly w)
{
    int ret;
    ret = abase_p_8_vec_scan(k, &(w->c), &(w->alloc));
    w->size = w->alloc;
    return ret;
}

/* *simd_gfp::code_for_set_ui_at */
static inline
void abase_p_8_set_ui_at(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_elt p, int k MAYBE_UNUSED, unsigned long v)
{
    abase_p_8_set_ui(K,p,v);
}

/* *simd_gfp::code_for_set_ui_all */
static inline
void abase_p_8_set_ui_all(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_elt p, unsigned long v)
{
    abase_p_8_set_ui(K,p,v);
}

/* *simd_gfp::code_for_elt_ur_set_ui_at */
static inline
void abase_p_8_elt_ur_set_ui_at(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_elt p, int k MAYBE_UNUSED, unsigned long v)
{
    abase_p_8_set_ui(K,p,v);
}

/* *simd_gfp::code_for_elt_ur_set_ui_all */
static inline
void abase_p_8_elt_ur_set_ui_all(abase_p_8_dst_field K MAYBE_UNUSED, abase_p_8_dst_elt p, unsigned long v)
{
    abase_p_8_set_ui(K,p,v);
}

static inline
void abase_p_8_oo_field_clear(abase_vbase_ptr f)
{
    abase_p_8_field_clear((abase_p_8_dst_field)(f->obj));
    free(f->obj);
    f->obj = NULL;
}


#endif  /* ABASE_P_8_H_ */

/* vim:set ft=cpp: */
