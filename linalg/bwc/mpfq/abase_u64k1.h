#ifndef ABASE_U64K1_H_
#define ABASE_U64K1_H_

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
#include "select_mpi.h"
#include "abase_vbase.h"
#ifdef	MPFQ_LAST_GENERATED_TAG
#undef	MPFQ_LAST_GENERATED_TAG
#endif
#define MPFQ_LAST_GENERATED_TAG      u64k1

/* Active handler: simd_u64k */
/* Automatically generated code  */
/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::vec */
/* Active handler: Mpfq::defaults::mpi_flat */
/* Active handler: simd_dotprod */
/* Active handler: io */
/* Active handler: trivialities */
/* Active handler: simd_char2 */
/* Options used:{
   family=[ u64k1, u64k2, ],
   k=1,
   tag=u64k1,
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
     [ (?^:abase_u64k1_elt \*), void *, ],
     [ (?^:abase_u64k1_src_elt\b), const void *, ],
     [ (?^:abase_u64k1_elt\b), void *, ],
     [ (?^:abase_u64k1_dst_elt\b), void *, ],
     [ (?^:abase_u64k1_elt_ur \*), void *, ],
     [ (?^:abase_u64k1_src_elt_ur\b), const void *, ],
     [ (?^:abase_u64k1_elt_ur\b), void *, ],
     [ (?^:abase_u64k1_dst_elt_ur\b), void *, ],
     [ (?^:abase_u64k1_vec \*), void *, ],
     [ (?^:abase_u64k1_src_vec\b), const void *, ],
     [ (?^:abase_u64k1_vec\b), void *, ],
     [ (?^:abase_u64k1_dst_vec\b), void *, ],
     [ (?^:abase_u64k1_vec_ur \*), void *, ],
     [ (?^:abase_u64k1_src_vec_ur\b), const void *, ],
     [ (?^:abase_u64k1_vec_ur\b), void *, ],
     [ (?^:abase_u64k1_dst_vec_ur\b), void *, ],
     [ (?^:abase_u64k1_poly \*), void *, ],
     [ (?^:abase_u64k1_src_poly\b), const void *, ],
     [ (?^:abase_u64k1_poly\b), void *, ],
     [ (?^:abase_u64k1_dst_poly\b), void *, ],
     ],
    },
   w=64,
   } */

typedef void * abase_u64k1_field[1];
typedef void * abase_u64k1_dst_field;

typedef uint64_t abase_u64k1_elt[1];
typedef uint64_t * abase_u64k1_dst_elt;
typedef const uint64_t * abase_u64k1_src_elt;

typedef uint64_t abase_u64k1_elt_ur[1];
typedef uint64_t * abase_u64k1_dst_elt_ur;
typedef const uint64_t * abase_u64k1_src_elt_ur;

typedef abase_u64k1_elt * abase_u64k1_vec;
typedef abase_u64k1_elt * abase_u64k1_dst_vec;
typedef abase_u64k1_elt * abase_u64k1_src_vec;

typedef abase_u64k1_elt_ur * abase_u64k1_vec_ur;
typedef abase_u64k1_elt_ur * abase_u64k1_dst_vec_ur;
typedef abase_u64k1_elt_ur * abase_u64k1_src_vec_ur;

typedef struct {
  abase_u64k1_vec c;
  unsigned int alloc;
  unsigned int size;
} abase_u64k1_poly_struct;
typedef abase_u64k1_poly_struct abase_u64k1_poly [1];
typedef abase_u64k1_poly_struct * abase_u64k1_dst_poly;
typedef abase_u64k1_poly_struct * abase_u64k1_src_poly;

#ifdef  __cplusplus
extern "C" {
#endif
/* *Mpfq::defaults::code_for_impl_name */
#define abase_u64k1_impl_name()	"u64k1"
/* missing impl_max_characteristic_bits */
/* missing impl_max_degree */

/* Functions operating on the field structure */
/* *simd_char2::code_for_field_characteristic */
#define abase_u64k1_field_characteristic(K, z)	mpz_set_ui(z,2)
/* missing field_characteristic_bits */
/* *simd_u64k::code_for_field_degree */
#define abase_u64k1_field_degree(f)	1
static inline
void abase_u64k1_field_init(abase_u64k1_dst_field);
/* *simd_u64k::code_for_field_clear */
#define abase_u64k1_field_clear(K)	/**/
void abase_u64k1_field_specify(abase_u64k1_dst_field, unsigned long, void *);
/* *simd_u64k::code_for_field_setopt */
#define abase_u64k1_field_setopt(f, x, y)	/**/

/* Element allocation functions */
/* *Mpfq::defaults::flatdata::code_for_init, simd_flat */
#define abase_u64k1_init(f, px)	/**/
/* *Mpfq::defaults::flatdata::code_for_clear, simd_flat */
#define abase_u64k1_clear(f, px)	/**/

/* Elementary assignment functions */
static inline
void abase_u64k1_set(abase_u64k1_dst_field, abase_u64k1_dst_elt, abase_u64k1_src_elt);
/* missing set_ui */
static inline
void abase_u64k1_set_zero(abase_u64k1_dst_field, abase_u64k1_dst_elt);
/* missing get_ui */
/* missing set_mpn */
/* missing set_mpz */
/* missing get_mpn */
/* missing get_mpz */

/* Assignment of random values */
static inline
void abase_u64k1_random(abase_u64k1_dst_field, abase_u64k1_dst_elt, gmp_randstate_t);
/* missing random2 */

/* Arithmetic operations on elements */
static inline
void abase_u64k1_add(abase_u64k1_dst_field, abase_u64k1_dst_elt, abase_u64k1_src_elt, abase_u64k1_src_elt);
/* *simd_char2::code_for_sub */
#define abase_u64k1_sub(K, r, s1, s2)	abase_u64k1_add(K,r,s1,s2)
/* *simd_char2::code_for_neg */
#define abase_u64k1_neg(K, r, s)	abase_u64k1_set(K,r,s)
/* missing mul */
/* missing sqr */
/* missing is_sqr */
/* missing sqrt */
/* missing pow */
/* missing powz */
/* missing frobenius */
/* missing add_ui */
/* missing sub_ui */
/* missing mul_ui */
/* missing inv */

/* Operations involving unreduced elements */
/* *Mpfq::defaults::flatdata::code_for_elt_ur_init, simd_flat */
#define abase_u64k1_elt_ur_init(f, px)	/**/
/* *Mpfq::defaults::flatdata::code_for_elt_ur_clear, simd_flat */
#define abase_u64k1_elt_ur_clear(f, px)	/**/
static inline
void abase_u64k1_elt_ur_set(abase_u64k1_dst_field, abase_u64k1_dst_elt_ur, abase_u64k1_src_elt_ur);
static inline
void abase_u64k1_elt_ur_set_elt(abase_u64k1_dst_field, abase_u64k1_dst_elt_ur, abase_u64k1_src_elt);
static inline
void abase_u64k1_elt_ur_set_zero(abase_u64k1_dst_field, abase_u64k1_dst_elt_ur);
/* missing elt_ur_set_ui */
static inline
void abase_u64k1_elt_ur_add(abase_u64k1_dst_field, abase_u64k1_dst_elt_ur, abase_u64k1_src_elt_ur, abase_u64k1_src_elt_ur);
/* *simd_char2::code_for_elt_ur_neg */
#define abase_u64k1_elt_ur_neg(K, r, s)	abase_u64k1_elt_ur_set(K,r,s)
/* *simd_char2::code_for_elt_ur_sub */
#define abase_u64k1_elt_ur_sub(K, r, s1, s2)	abase_u64k1_elt_ur_add(K,r,s1,s2)
/* missing mul_ur */
/* missing sqr_ur */
/* missing reduce */

/* Comparison functions */
static inline
int abase_u64k1_cmp(abase_u64k1_dst_field, abase_u64k1_src_elt, abase_u64k1_src_elt);
/* missing cmp_ui */
static inline
int abase_u64k1_is_zero(abase_u64k1_dst_field, abase_u64k1_src_elt);

/* Input/output functions */
void abase_u64k1_asprint(abase_u64k1_dst_field, char * *, abase_u64k1_src_elt);
void abase_u64k1_fprint(abase_u64k1_dst_field, FILE *, abase_u64k1_src_elt);
/* *io::code_for_print */
#define abase_u64k1_print(k, x)	abase_u64k1_fprint(k,stdout,x)
int abase_u64k1_sscan(abase_u64k1_dst_field, abase_u64k1_dst_elt, const char *);
int abase_u64k1_fscan(abase_u64k1_dst_field, FILE *, abase_u64k1_dst_elt);
/* *Mpfq::defaults::code_for_scan */
#define abase_u64k1_scan(k, x)	abase_u64k1_fscan(k,stdin,x)

/* Vector functions */
void abase_u64k1_vec_init(abase_u64k1_dst_field, abase_u64k1_vec *, unsigned int);
void abase_u64k1_vec_reinit(abase_u64k1_dst_field, abase_u64k1_vec *, unsigned int, unsigned int);
void abase_u64k1_vec_clear(abase_u64k1_dst_field, abase_u64k1_vec *, unsigned int);
static inline
void abase_u64k1_vec_set(abase_u64k1_dst_field, abase_u64k1_dst_vec, abase_u64k1_src_vec, unsigned int);
static inline
void abase_u64k1_vec_set_zero(abase_u64k1_dst_field, abase_u64k1_dst_vec, unsigned int);
static inline
void abase_u64k1_vec_setcoeff(abase_u64k1_dst_field, abase_u64k1_dst_vec, abase_u64k1_src_elt, unsigned int);
/* missing vec_setcoeff_ui */
static inline
void abase_u64k1_vec_getcoeff(abase_u64k1_dst_field, abase_u64k1_dst_elt, abase_u64k1_src_vec, unsigned int);
static inline
void abase_u64k1_vec_add(abase_u64k1_dst_field, abase_u64k1_dst_vec, abase_u64k1_src_vec, abase_u64k1_src_vec, unsigned int);
static inline
void abase_u64k1_vec_neg(abase_u64k1_dst_field, abase_u64k1_dst_vec, abase_u64k1_src_vec, unsigned int);
static inline
void abase_u64k1_vec_rev(abase_u64k1_dst_field, abase_u64k1_dst_vec, abase_u64k1_src_vec, unsigned int);
static inline
void abase_u64k1_vec_sub(abase_u64k1_dst_field, abase_u64k1_dst_vec, abase_u64k1_src_vec, abase_u64k1_src_vec, unsigned int);
/* missing vec_scal_mul */
/* missing vec_conv */
static inline
void abase_u64k1_vec_random(abase_u64k1_dst_field, abase_u64k1_dst_vec, unsigned int, gmp_randstate_t);
/* missing vec_random2 */
static inline
int abase_u64k1_vec_cmp(abase_u64k1_dst_field, abase_u64k1_src_vec, abase_u64k1_src_vec, unsigned int);
static inline
int abase_u64k1_vec_is_zero(abase_u64k1_dst_field, abase_u64k1_src_vec, unsigned int);
static inline
abase_u64k1_dst_vec abase_u64k1_vec_subvec(abase_u64k1_dst_field, abase_u64k1_dst_vec, int);
static inline
abase_u64k1_src_vec abase_u64k1_vec_subvec_const(abase_u64k1_dst_field, abase_u64k1_src_vec, int);
static inline
abase_u64k1_dst_elt abase_u64k1_vec_coeff_ptr(abase_u64k1_dst_field, abase_u64k1_dst_vec, int);
static inline
abase_u64k1_src_elt abase_u64k1_vec_coeff_ptr_const(abase_u64k1_dst_field, abase_u64k1_src_vec, int);
void abase_u64k1_vec_asprint(abase_u64k1_dst_field, char * *, abase_u64k1_src_vec, unsigned int);
void abase_u64k1_vec_fprint(abase_u64k1_dst_field, FILE *, abase_u64k1_src_vec, unsigned int);
void abase_u64k1_vec_print(abase_u64k1_dst_field, abase_u64k1_src_vec, unsigned int);
int abase_u64k1_vec_sscan(abase_u64k1_dst_field, abase_u64k1_vec *, unsigned int *, const char *);
int abase_u64k1_vec_fscan(abase_u64k1_dst_field, FILE *, abase_u64k1_vec *, unsigned int *);
/* *Mpfq::defaults::vec::io::code_for_vec_scan, Mpfq::defaults::vec */
#define abase_u64k1_vec_scan(K, w, n)	abase_u64k1_vec_fscan(K,stdout,w,n)
void abase_u64k1_vec_ur_init(abase_u64k1_dst_field, abase_u64k1_vec_ur *, unsigned int);
static inline
void abase_u64k1_vec_ur_set_zero(abase_u64k1_dst_field, abase_u64k1_dst_vec_ur, unsigned int);
static inline
void abase_u64k1_vec_ur_set_vec(abase_u64k1_dst_field, abase_u64k1_dst_vec_ur, abase_u64k1_src_vec, unsigned int);
void abase_u64k1_vec_ur_reinit(abase_u64k1_dst_field, abase_u64k1_vec_ur *, unsigned int, unsigned int);
void abase_u64k1_vec_ur_clear(abase_u64k1_dst_field, abase_u64k1_vec_ur *, unsigned int);
static inline
void abase_u64k1_vec_ur_set(abase_u64k1_dst_field, abase_u64k1_dst_vec_ur, abase_u64k1_src_vec_ur, unsigned int);
static inline
void abase_u64k1_vec_ur_setcoeff(abase_u64k1_dst_field, abase_u64k1_dst_vec_ur, abase_u64k1_src_elt_ur, unsigned int);
static inline
void abase_u64k1_vec_ur_getcoeff(abase_u64k1_dst_field, abase_u64k1_dst_elt_ur, abase_u64k1_src_vec_ur, unsigned int);
static inline
void abase_u64k1_vec_ur_add(abase_u64k1_dst_field, abase_u64k1_dst_vec_ur, abase_u64k1_src_vec_ur, abase_u64k1_src_vec_ur, unsigned int);
static inline
void abase_u64k1_vec_ur_sub(abase_u64k1_dst_field, abase_u64k1_dst_vec_ur, abase_u64k1_src_vec_ur, abase_u64k1_src_vec_ur, unsigned int);
static inline
void abase_u64k1_vec_ur_neg(abase_u64k1_dst_field, abase_u64k1_dst_vec_ur, abase_u64k1_src_vec_ur, unsigned int);
static inline
void abase_u64k1_vec_ur_rev(abase_u64k1_dst_field, abase_u64k1_dst_vec_ur, abase_u64k1_src_vec_ur, unsigned int);
/* missing vec_scal_mul_ur */
/* missing vec_conv_ur */
/* missing vec_reduce */
static inline
abase_u64k1_dst_vec_ur abase_u64k1_vec_ur_subvec(abase_u64k1_dst_field, abase_u64k1_dst_vec_ur, int);
static inline
abase_u64k1_src_vec_ur abase_u64k1_vec_ur_subvec_const(abase_u64k1_dst_field, abase_u64k1_src_vec_ur, int);
static inline
abase_u64k1_dst_elt abase_u64k1_vec_ur_coeff_ptr(abase_u64k1_dst_field, abase_u64k1_dst_vec_ur, int);
static inline
abase_u64k1_src_elt abase_u64k1_vec_ur_coeff_ptr_const(abase_u64k1_dst_field, abase_u64k1_src_vec_ur, int);
/* *Mpfq::defaults::flatdata::code_for_vec_elt_stride, simd_flat */
#define abase_u64k1_vec_elt_stride(K, n)	((n)*sizeof(abase_u64k1_elt))
/* *Mpfq::defaults::flatdata::code_for_vec_ur_elt_stride, simd_flat */
#define abase_u64k1_vec_ur_elt_stride(K, n)	((n)*sizeof(abase_u64k1_elt_ur))

/* Polynomial functions */
/* missing poly_init */
/* missing poly_clear */
/* missing poly_set */
/* missing poly_setmonic */
/* missing poly_setcoeff */
/* missing poly_setcoeff_ui */
/* missing poly_getcoeff */
/* missing poly_deg */
/* missing poly_add */
/* missing poly_sub */
/* missing poly_add_ui */
/* missing poly_sub_ui */
/* missing poly_neg */
/* missing poly_scal_mul */
/* missing poly_mul */
/* missing poly_divmod */
/* missing poly_precomp_mod */
/* missing poly_mod_pre */
/* missing poly_gcd */
/* missing poly_xgcd */
/* missing poly_random */
/* missing poly_random2 */
/* missing poly_cmp */
/* missing poly_asprint */
/* missing poly_fprint */
/* missing poly_print */
/* missing poly_sscan */
/* missing poly_fscan */
/* missing poly_scan */

/* Functions related to SIMD operation */
/* *simd_u64k::code_for_groupsize */
#define abase_u64k1_groupsize(K)	64
/* *trivialities::code_for_offset */
#define abase_u64k1_offset(K, n)	n /* TO BE DEPRECATED */
/* *trivialities::code_for_stride */
#define abase_u64k1_stride(K)	1 /* TO BE DEPRECATED */
static inline
void abase_u64k1_set_ui_at(abase_u64k1_dst_field, abase_u64k1_dst_elt, int, unsigned long);
static inline
void abase_u64k1_set_ui_all(abase_u64k1_dst_field, abase_u64k1_dst_elt, unsigned long);
static inline
void abase_u64k1_elt_ur_set_ui_at(abase_u64k1_dst_field, abase_u64k1_dst_elt, int, unsigned long);
static inline
void abase_u64k1_elt_ur_set_ui_all(abase_u64k1_dst_field, abase_u64k1_dst_elt, unsigned long);
void abase_u64k1_dotprod(abase_u64k1_dst_field, abase_u64k1_dst_vec, abase_u64k1_src_vec, abase_u64k1_src_vec, unsigned int);

/* Member templates related to SIMD operation */

/* MPI interface */
void abase_u64k1_mpi_ops_init(abase_u64k1_dst_field);
MPI_Datatype abase_u64k1_mpi_datatype(abase_u64k1_dst_field);
MPI_Datatype abase_u64k1_mpi_datatype_ur(abase_u64k1_dst_field);
MPI_Op abase_u64k1_mpi_addition_op(abase_u64k1_dst_field);
MPI_Op abase_u64k1_mpi_addition_op_ur(abase_u64k1_dst_field);
void abase_u64k1_mpi_ops_clear(abase_u64k1_dst_field);

/* Object-oriented interface */
void abase_u64k1_oo_field_init(abase_vbase_ptr);
static inline
void abase_u64k1_oo_field_clear(abase_vbase_ptr);
#ifdef  __cplusplus
}
#endif

/* Implementations for inlines */
/* *simd_u64k::code_for_field_init */
static inline
void abase_u64k1_field_init(abase_u64k1_dst_field f MAYBE_UNUSED)
{
}

/* *Mpfq::defaults::flatdata::code_for_set, simd_flat */
static inline
void abase_u64k1_set(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_elt r, abase_u64k1_src_elt s)
{
    if (r != s) memcpy(r,s,sizeof(abase_u64k1_elt));
}

/* *simd_flat::code_for_set_zero */
static inline
void abase_u64k1_set_zero(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_elt r)
{
    memset(r, 0, sizeof(abase_u64k1_elt));
}

/* *simd_flat::code_for_random */
static inline
void abase_u64k1_random(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_elt r, gmp_randstate_t state)
{
        for(unsigned int i = 0 ; i < sizeof(abase_u64k1_elt) ; i++) {
            ((unsigned char*)r)[i] = gmp_urandomb_ui(state, 8);
        }
}

/* *simd_char2::code_for_add */
static inline
void abase_u64k1_add(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_elt r, abase_u64k1_src_elt s1, abase_u64k1_src_elt s2)
{
        for(unsigned int i = 0 ; i < sizeof(abase_u64k1_elt)/sizeof(*r) ; i++) {
            r[i] = s1[i] ^ s2[i];
        }
}

/* *Mpfq::defaults::flatdata::code_for_elt_ur_set, simd_flat */
static inline
void abase_u64k1_elt_ur_set(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_elt_ur r, abase_u64k1_src_elt_ur s)
{
    if (r != s) memcpy(r,s,sizeof(abase_u64k1_elt_ur));
}

/* *Mpfq::defaults::flatdata::code_for_elt_ur_set_elt, simd_flat */
static inline
void abase_u64k1_elt_ur_set_elt(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_elt_ur r, abase_u64k1_src_elt s)
{
    memset(r, 0, sizeof(abase_u64k1_elt_ur)); memcpy(r,s,sizeof(abase_u64k1_elt));
}

/* *Mpfq::defaults::flatdata::code_for_elt_ur_set_zero, simd_flat */
static inline
void abase_u64k1_elt_ur_set_zero(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_elt_ur r)
{
    memset(r, 0, sizeof(abase_u64k1_elt_ur));
}

/* *simd_char2::code_for_elt_ur_add */
static inline
void abase_u64k1_elt_ur_add(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_elt_ur r, abase_u64k1_src_elt_ur s1, abase_u64k1_src_elt_ur s2)
{
        for(unsigned int i = 0 ; i < sizeof(abase_u64k1_elt)/sizeof(*r) ; i++) {
            r[i] = s1[i] ^ s2[i];
        }
}

/* *Mpfq::defaults::flatdata::code_for_cmp, simd_flat */
static inline
int abase_u64k1_cmp(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_src_elt r, abase_u64k1_src_elt s)
{
    return memcmp(r,s,sizeof(abase_u64k1_elt));
}

/* *Mpfq::defaults::flatdata::code_for_is_zero, simd_flat */
static inline
int abase_u64k1_is_zero(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_src_elt r)
{
        unsigned int i;
        for(i = 0 ; i < sizeof(abase_u64k1_elt)/sizeof(r[0]) ; i++) {
            if (r[i]) return 0;
        }
        return 1;
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_set, Mpfq::defaults::flatdata, simd_flat */
static inline
void abase_u64k1_vec_set(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_vec r, abase_u64k1_src_vec s, unsigned int n)
{
    if (r != s) memmove(r, s, n*sizeof(abase_u64k1_elt));
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_set_zero, Mpfq::defaults::flatdata, simd_flat */
static inline
void abase_u64k1_vec_set_zero(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_vec r, unsigned int n)
{
    memset(r, 0, n*sizeof(abase_u64k1_elt));
}

/* *Mpfq::defaults::vec::getset::code_for_vec_setcoeff, Mpfq::defaults::vec */
static inline
void abase_u64k1_vec_setcoeff(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_vec w, abase_u64k1_src_elt x, unsigned int i)
{
    abase_u64k1_set(K, w[i], x);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_getcoeff, Mpfq::defaults::vec */
static inline
void abase_u64k1_vec_getcoeff(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_elt x, abase_u64k1_src_vec w, unsigned int i)
{
    abase_u64k1_set(K, x, w[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_add, Mpfq::defaults::vec */
static inline
void abase_u64k1_vec_add(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_vec w, abase_u64k1_src_vec u, abase_u64k1_src_vec v, unsigned int n)
{
        unsigned int i;
    for(i = 0; i < n; i+=1)
        abase_u64k1_add(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_neg, Mpfq::defaults::vec */
static inline
void abase_u64k1_vec_neg(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_vec w, abase_u64k1_src_vec u, unsigned int n)
{
        unsigned int i;
    for(i = 0; i < n; ++i)
        abase_u64k1_neg(K, w[i], u[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_rev, Mpfq::defaults::vec */
static inline
void abase_u64k1_vec_rev(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_vec w, abase_u64k1_src_vec u, unsigned int n)
{
    unsigned int nn = n >> 1;
    abase_u64k1_elt tmp[1];
    abase_u64k1_init(K, tmp);
    unsigned int i;
    for(i = 0; i < nn; ++i) {
        abase_u64k1_set(K, tmp[0], u[i]);
        abase_u64k1_set(K, w[i], u[n-1-i]);
        abase_u64k1_set(K, w[n-1-i], tmp[0]);
    }
    if (n & 1)
        abase_u64k1_set(K, w[nn], u[nn]);
    abase_u64k1_clear(K, tmp);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_sub, Mpfq::defaults::vec */
static inline
void abase_u64k1_vec_sub(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_vec w, abase_u64k1_src_vec u, abase_u64k1_src_vec v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; ++i)
        abase_u64k1_sub(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_random, Mpfq::defaults::vec */
static inline
void abase_u64k1_vec_random(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_vec w, unsigned int n, gmp_randstate_t state)
{
    unsigned int i;
    for(i = 0; i < n; ++i)
        abase_u64k1_random(K, w[i], state);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_cmp, Mpfq::defaults::vec */
static inline
int abase_u64k1_vec_cmp(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_src_vec u, abase_u64k1_src_vec v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; ++i) {
        int ret = abase_u64k1_cmp(K, u[i], v[i]);
        if (ret != 0)
            return ret;
    }
    return 0;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_is_zero, Mpfq::defaults::vec */
static inline
int abase_u64k1_vec_is_zero(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_src_vec r, unsigned int n)
{
    unsigned int i;
    for(i = 0 ; i < n ; i+=1) {
        if (!abase_u64k1_is_zero(K,r[i])) return 0;
    }
    return 1;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_subvec, Mpfq::defaults::vec */
static inline
abase_u64k1_dst_vec abase_u64k1_vec_subvec(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_vec v, int i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_subvec_const, Mpfq::defaults::vec */
static inline
abase_u64k1_src_vec abase_u64k1_vec_subvec_const(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_src_vec v, int i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_coeff_ptr, Mpfq::defaults::vec */
static inline
abase_u64k1_dst_elt abase_u64k1_vec_coeff_ptr(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_vec v, int i)
{
    return v[i];
}

/* *Mpfq::defaults::vec::getset::code_for_vec_coeff_ptr_const, Mpfq::defaults::vec */
static inline
abase_u64k1_src_elt abase_u64k1_vec_coeff_ptr_const(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_src_vec v, int i)
{
    return v[i];
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_ur_set_zero, Mpfq::defaults::flatdata, simd_flat */
static inline
void abase_u64k1_vec_ur_set_zero(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_vec_ur r, unsigned int n)
{
    memset(r, 0, n*sizeof(abase_u64k1_elt_ur));
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_set_vec, Mpfq::defaults::vec */
static inline
void abase_u64k1_vec_ur_set_vec(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_vec_ur w, abase_u64k1_src_vec u, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        abase_u64k1_elt_ur_set_elt(K, w[i], u[i]);
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_ur_set, Mpfq::defaults::flatdata, simd_flat */
static inline
void abase_u64k1_vec_ur_set(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_vec_ur r, abase_u64k1_src_vec_ur s, unsigned int n)
{
    if (r != s) memmove(r, s, n*sizeof(abase_u64k1_elt_ur));
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_setcoeff, Mpfq::defaults::vec */
static inline
void abase_u64k1_vec_ur_setcoeff(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_vec_ur w, abase_u64k1_src_elt_ur x, unsigned int i)
{
    abase_u64k1_elt_ur_set(K, w[i], x);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_getcoeff, Mpfq::defaults::vec */
static inline
void abase_u64k1_vec_ur_getcoeff(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_elt_ur x, abase_u64k1_src_vec_ur w, unsigned int i)
{
    abase_u64k1_elt_ur_set(K, x, w[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_add, Mpfq::defaults::vec */
static inline
void abase_u64k1_vec_ur_add(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_vec_ur w, abase_u64k1_src_vec_ur u, abase_u64k1_src_vec_ur v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        abase_u64k1_elt_ur_add(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_sub, Mpfq::defaults::vec */
static inline
void abase_u64k1_vec_ur_sub(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_vec_ur w, abase_u64k1_src_vec_ur u, abase_u64k1_src_vec_ur v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        abase_u64k1_elt_ur_sub(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_neg, Mpfq::defaults::vec */
static inline
void abase_u64k1_vec_ur_neg(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_vec_ur w, abase_u64k1_src_vec_ur u, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; ++i)
        abase_u64k1_elt_ur_neg(K, w[i], u[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_rev, Mpfq::defaults::vec */
static inline
void abase_u64k1_vec_ur_rev(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_vec_ur w, abase_u64k1_src_vec_ur u, unsigned int n)
{
    unsigned int nn = n >> 1;
    abase_u64k1_elt_ur tmp[1];
    abase_u64k1_elt_ur_init(K, tmp);
    unsigned int i;
    for(i = 0; i < nn; ++i) {
        abase_u64k1_elt_ur_set(K, tmp[0], u[i]);
        abase_u64k1_elt_ur_set(K, w[i], u[n-1-i]);
        abase_u64k1_elt_ur_set(K, w[n-1-i], tmp[0]);
    }
    if (n & 1)
        abase_u64k1_elt_ur_set(K, w[nn], u[nn]);
    abase_u64k1_elt_ur_clear(K, tmp);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_subvec, Mpfq::defaults::vec */
static inline
abase_u64k1_dst_vec_ur abase_u64k1_vec_ur_subvec(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_vec_ur v, int i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_subvec_const, Mpfq::defaults::vec */
static inline
abase_u64k1_src_vec_ur abase_u64k1_vec_ur_subvec_const(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_src_vec_ur v, int i)
{
    return v+i;
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_coeff_ptr, Mpfq::defaults::vec */
static inline
abase_u64k1_dst_elt abase_u64k1_vec_ur_coeff_ptr(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_vec_ur v, int i)
{
    return v[i];
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_coeff_ptr_const, Mpfq::defaults::vec */
static inline
abase_u64k1_src_elt abase_u64k1_vec_ur_coeff_ptr_const(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_src_vec_ur v, int i)
{
    return v[i];
}

/* *simd_flat::code_for_set_ui_at */
static inline
void abase_u64k1_set_ui_at(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_elt p, int k, unsigned long v)
{
        assert(k < abase_u64k1_groupsize(K));
        uint64_t * xp = (uint64_t *) p;
        uint64_t mask = ((uint64_t)1) << (k%64);
        xp[k/64] = (xp[k/64] & ~mask) | ((((uint64_t)v) << (k%64))&mask);
}

/* *simd_flat::code_for_set_ui_all */
static inline
void abase_u64k1_set_ui_all(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_elt r, unsigned long v)
{
        for(unsigned int i = 0 ; i < sizeof(abase_u64k1_elt)/sizeof(*r) ; i++) {
            r[i] = ~v;
        }
}

/* *simd_flat::code_for_elt_ur_set_ui_at */
static inline
void abase_u64k1_elt_ur_set_ui_at(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_elt p, int k, unsigned long v)
{
        assert(k < abase_u64k1_groupsize(K));
        uint64_t * xp = (uint64_t *) p;
        uint64_t mask = ((uint64_t)1) << (k%64);
        xp[k/64] = (xp[k/64] & ~mask) | ((((uint64_t)v) << (k%64))&mask);
}

/* *simd_flat::code_for_elt_ur_set_ui_all */
static inline
void abase_u64k1_elt_ur_set_ui_all(abase_u64k1_dst_field K MAYBE_UNUSED, abase_u64k1_dst_elt r, unsigned long v)
{
        for(unsigned int i = 0 ; i < sizeof(abase_u64k1_elt)/sizeof(*r) ; i++) {
            r[i] = ~v;
        }
}

static inline
void abase_u64k1_oo_field_clear(abase_vbase_ptr f)
{
    abase_u64k1_field_clear((abase_u64k1_dst_field)(f->obj));
    free(f->obj);
    f->obj = NULL;
}


#endif  /* ABASE_U64K1_H_ */

/* vim:set ft=cpp: */
