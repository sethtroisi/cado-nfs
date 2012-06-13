#ifndef ABASE_P16_H_
#define ABASE_P16_H_

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
#include "random_generation.h"
#include "abase_vbase.h"
#ifdef	MPFQ_LAST_GENERATED_TAG
#undef	MPFQ_LAST_GENERATED_TAG
#endif
#define MPFQ_LAST_GENERATED_TAG      p16

/* Active handler: simd_p16 */
/* Automatically generated code  */
/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::vec */
/* Active handler: Mpfq::defaults::mpi_flat */
/* Active handler: io */
/* Active handler: trivialities */
/* Active handler: simd_flat */
/* Options used: vtag=p16 tag=p16 vbase_stuff={
                 'vc:includes' => [
                                    '<stdarg.h>'
                                  ],
                 'member_templates_restrict' => {
                                                  'u64k2' => [
                                                               'u64k1',
                                                               'u64k2'
                                                             ],
                                                  'p16' => [
                                                             'p16'
                                                           ],
                                                  'u64k1' => $vbase_stuff->{'member_templates_restrict'}{'u64k2'}
                                                },
                 'families' => [
                                 $vbase_stuff->{'member_templates_restrict'}{'u64k2'},
                                 $vbase_stuff->{'member_templates_restrict'}{'p16'}
                               ],
                 'choose_byfeatures' => sub { "DUMMY" }
               };
 family=[p16] virtual_base={
                  'filebase' => 'abase_vbase',
                  'substitutions' => [
                                       [
                                         qr/(?^:abase_p16_elt \*)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p16_src_elt\b)/,
                                         'const void *'
                                       ],
                                       [
                                         qr/(?^:abase_p16_elt\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p16_dst_elt\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p16_elt_ur \*)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p16_src_elt_ur\b)/,
                                         'const void *'
                                       ],
                                       [
                                         qr/(?^:abase_p16_elt_ur\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p16_dst_elt_ur\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p16_vec \*)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p16_src_vec\b)/,
                                         'const void *'
                                       ],
                                       [
                                         qr/(?^:abase_p16_vec\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p16_dst_vec\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p16_vec_ur \*)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p16_src_vec_ur\b)/,
                                         'const void *'
                                       ],
                                       [
                                         qr/(?^:abase_p16_vec_ur\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p16_dst_vec_ur\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p16_poly \*)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p16_src_poly\b)/,
                                         'const void *'
                                       ],
                                       [
                                         qr/(?^:abase_p16_poly\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?^:abase_p16_dst_poly\b)/,
                                         'void *'
                                       ]
                                     ],
                  'name' => 'abase_vbase',
                  'global_prefix' => 'abase_'
                };
 */

typedef int32_t abase_p16_field[1];
typedef int32_t * abase_p16_dst_field;

typedef int32_t abase_p16_elt[1];
typedef int32_t * abase_p16_dst_elt;
typedef const int32_t * abase_p16_src_elt;

typedef int32_t abase_p16_elt_ur[1];
typedef int32_t * abase_p16_dst_elt_ur;
typedef const int32_t * abase_p16_src_elt_ur;

typedef abase_p16_elt * abase_p16_vec;
typedef abase_p16_elt * abase_p16_dst_vec;
typedef abase_p16_elt * abase_p16_src_vec;

typedef abase_p16_elt_ur * abase_p16_vec_ur;
typedef abase_p16_elt_ur * abase_p16_dst_vec_ur;
typedef abase_p16_elt_ur * abase_p16_src_vec_ur;

#ifdef  __cplusplus
extern "C" {
#endif

/* Functions operating on the field structure */
static inline
void abase_p16_field_characteristic(abase_p16_dst_field, mpz_t);
/* *simd_p16::code_for_field_degree */
#define abase_p16_field_degree(K)	1
/* *simd_p16::code_for_field_init */
#define abase_p16_field_init(f)	/**/
/* *simd_p16::code_for_field_clear */
#define abase_p16_field_clear(f)	/**/
void abase_p16_field_specify(abase_p16_dst_field, unsigned long, void *);
/* *simd_p16::code_for_field_setopt */
#define abase_p16_field_setopt(f, x, y)	/**/

/* Element allocation functions */
/* *Mpfq::defaults::flatdata::code_for_init, simd_flat */
#define abase_p16_init(f, px)	/**/
/* *Mpfq::defaults::flatdata::code_for_clear, simd_flat */
#define abase_p16_clear(f, px)	/**/

/* Elementary assignment functions */
static inline
void abase_p16_set(abase_p16_dst_field, abase_p16_dst_elt, abase_p16_src_elt);
/* missing set_ui */
static inline
void abase_p16_set_zero(abase_p16_dst_field, abase_p16_dst_elt);
/* missing get_ui */
/* missing set_mpn */
/* missing set_mpz */
/* missing get_mpn */
/* missing get_mpz */

/* Assignment of random values */
static inline
void abase_p16_random(abase_p16_dst_field, abase_p16_dst_elt);

/* Arithmetic operations on elements */
/* *simd_p16::code_for_add */
#define abase_p16_add(K, r, s1, s2)	*r = *s1 + *s2
/* *simd_p16::code_for_sub */
#define abase_p16_sub(K, r, s1, s2)	*r = *s1 - *s2
/* *simd_p16::code_for_neg */
#define abase_p16_neg(K, r, s)	*r=-*s
/* missing mul */
/* missing sqr */
/* missing is_sqr */
/* missing sqrt */
/* missing pow */
/* missing add_ui */
/* missing sub_ui */
/* missing mul_ui */

/* Operations involving unreduced elements */
/* *Mpfq::defaults::flatdata::code_for_elt_ur_init, simd_flat */
#define abase_p16_elt_ur_init(f, px)	/**/
/* *Mpfq::defaults::flatdata::code_for_elt_ur_clear, simd_flat */
#define abase_p16_elt_ur_clear(f, px)	/**/
static inline
void abase_p16_elt_ur_set(abase_p16_dst_field, abase_p16_dst_elt_ur, abase_p16_src_elt_ur);
static inline
void abase_p16_elt_ur_set_zero(abase_p16_dst_field, abase_p16_dst_elt_ur);
/* missing elt_ur_set_ui */
/* *simd_p16::code_for_elt_ur_add */
#define abase_p16_elt_ur_add(K, r, s1, s2)	*r = *s1 + *s2
/* *simd_p16::code_for_elt_ur_neg */
#define abase_p16_elt_ur_neg(K, r, s)	*r=-*s
/* *simd_p16::code_for_elt_ur_sub */
#define abase_p16_elt_ur_sub(K, r, s1, s2)	*r = *s1 - *s2
/* missing mul_ur */
/* missing sqr_ur */
static inline
void abase_p16_reduce(abase_p16_dst_field, abase_p16_dst_elt, abase_p16_dst_elt_ur);
#define HAVE_abase_p16_addmul1
/* *simd_p16::code_for_addmul1 */
#define abase_p16_addmul1(K, r, s1, v)	*r += *s1 * v

/* Comparison functions */
static inline
int abase_p16_cmp(abase_p16_dst_field, abase_p16_src_elt, abase_p16_src_elt);
/* missing cmp_ui */
static inline
int abase_p16_is_zero(abase_p16_dst_field, abase_p16_src_elt);

/* Input/output functions */
void abase_p16_asprint(abase_p16_dst_field, char * *, abase_p16_src_elt);
void abase_p16_fprint(abase_p16_dst_field, FILE *, abase_p16_src_elt);
/* *io::code_for_print */
#define abase_p16_print(k, x)	abase_p16_fprint(k,stdout,x)
int abase_p16_sscan(abase_p16_dst_field, abase_p16_dst_elt, const char *);
int abase_p16_fscan(abase_p16_dst_field, FILE *, abase_p16_dst_elt);
/* *simd_p16::code_for_scan */
#define abase_p16_scan(k, x)	abase_p16_fscan(k,stdin,x)

/* Vector functions */
void abase_p16_vec_init(abase_p16_dst_field, abase_p16_vec *, unsigned int);
void abase_p16_vec_reinit(abase_p16_dst_field, abase_p16_vec *, unsigned int, unsigned int);
void abase_p16_vec_clear(abase_p16_dst_field, abase_p16_vec *, unsigned int);
static inline
void abase_p16_vec_set(abase_p16_dst_field, abase_p16_dst_vec, abase_p16_src_vec, unsigned int);
static inline
void abase_p16_vec_set_zero(abase_p16_dst_field, abase_p16_dst_vec, unsigned int);
static inline
void abase_p16_vec_setcoef(abase_p16_dst_field, abase_p16_dst_vec, abase_p16_src_elt, unsigned int);
/* missing vec_setcoef_ui */
static inline
void abase_p16_vec_getcoef(abase_p16_dst_field, abase_p16_dst_elt, abase_p16_src_vec, unsigned int);
static inline
void abase_p16_vec_add(abase_p16_dst_field, abase_p16_dst_vec, abase_p16_src_vec, abase_p16_src_vec, unsigned int);
static inline
void abase_p16_vec_neg(abase_p16_dst_field, abase_p16_dst_vec, abase_p16_src_vec, unsigned int);
static inline
void abase_p16_vec_rev(abase_p16_dst_field, abase_p16_dst_vec, abase_p16_src_vec, unsigned int);
static inline
void abase_p16_vec_sub(abase_p16_dst_field, abase_p16_dst_vec, abase_p16_src_vec, abase_p16_src_vec, unsigned int);
/* missing vec_scal_mul */
/* missing vec_conv */
static inline
void abase_p16_vec_random(abase_p16_dst_field, abase_p16_dst_vec, unsigned int);
static inline
int abase_p16_vec_cmp(abase_p16_dst_field, abase_p16_src_vec, abase_p16_src_vec, unsigned int);
static inline
int abase_p16_vec_is_zero(abase_p16_dst_field, abase_p16_src_vec, unsigned int);
void abase_p16_vec_asprint(abase_p16_dst_field, char * *, abase_p16_src_vec, unsigned int);
void abase_p16_vec_fprint(abase_p16_dst_field, FILE *, abase_p16_src_vec, unsigned int);
void abase_p16_vec_print(abase_p16_dst_field, abase_p16_src_vec, unsigned int);
int abase_p16_vec_sscan(abase_p16_dst_field, abase_p16_vec *, unsigned int *, const char *);
int abase_p16_vec_fscan(abase_p16_dst_field, FILE *, abase_p16_vec *, unsigned int *);
/* *Mpfq::defaults::vec::io::code_for_vec_scan, Mpfq::defaults::vec, Mpfq::defaults */
#define abase_p16_vec_scan(K, w, n)	abase_p16_vec_fscan(K,stdout,w,n)
void abase_p16_vec_ur_init(abase_p16_dst_field, abase_p16_vec_ur *, unsigned int);
void abase_p16_vec_ur_reinit(abase_p16_dst_field, abase_p16_vec_ur *, unsigned int, unsigned int);
void abase_p16_vec_ur_clear(abase_p16_dst_field, abase_p16_vec_ur *, unsigned int);
static inline
void abase_p16_vec_ur_set(abase_p16_dst_field, abase_p16_dst_vec_ur, abase_p16_src_vec_ur, unsigned int);
static inline
void abase_p16_vec_ur_setcoef(abase_p16_dst_field, abase_p16_dst_vec_ur, abase_p16_src_elt_ur, unsigned int);
static inline
void abase_p16_vec_ur_getcoef(abase_p16_dst_field, abase_p16_dst_elt_ur, abase_p16_src_vec_ur, unsigned int);
static inline
void abase_p16_vec_ur_add(abase_p16_dst_field, abase_p16_dst_vec_ur, abase_p16_src_vec_ur, abase_p16_src_vec_ur, unsigned int);
static inline
void abase_p16_vec_ur_sub(abase_p16_dst_field, abase_p16_dst_vec_ur, abase_p16_src_vec_ur, abase_p16_src_vec_ur, unsigned int);
/* missing vec_scal_mul_ur */
/* missing vec_conv_ur */
static inline
void abase_p16_vec_reduce(abase_p16_dst_field, abase_p16_dst_vec, abase_p16_dst_vec_ur, unsigned int);
/* *Mpfq::defaults::flatdata::code_for_vec_elt_stride, simd_flat */
#define abase_p16_vec_elt_stride(K, n)	((n)*sizeof(abase_p16_elt))

/* Functions related to SIMD operation */
/* *simd_p16::code_for_groupsize */
#define abase_p16_groupsize(K)	1
/* *simd_p16::code_for_offset */
#define abase_p16_offset(K, n)	n /* TO BE DEPRECATED */
/* *simd_p16::code_for_stride */
#define abase_p16_stride(K)	1 /* TO BE DEPRECATED */
static inline
void abase_p16_set_ui_at(abase_p16_dst_field, abase_p16_dst_elt, int, unsigned long);
static inline
void abase_p16_set_ui_all(abase_p16_dst_field, abase_p16_dst_elt, unsigned long);
static inline
void abase_p16_elt_ur_set_ui_at(abase_p16_dst_field, abase_p16_dst_elt, int, unsigned long);
static inline
void abase_p16_elt_ur_set_ui_all(abase_p16_dst_field, abase_p16_dst_elt, unsigned long);
void abase_p16_dotprod(abase_p16_dst_field, abase_p16_dst_vec, abase_p16_src_vec, abase_p16_src_vec, unsigned int);

/* Member templates related to SIMD operation */

/* MPI interface */
void abase_p16_mpi_ops_init(abase_p16_dst_field);
MPI_Datatype abase_p16_mpi_datatype(abase_p16_dst_field);
MPI_Datatype abase_p16_mpi_datatype_ur(abase_p16_dst_field);
MPI_Op abase_p16_mpi_addition_op(abase_p16_dst_field);
MPI_Op abase_p16_mpi_addition_op_ur(abase_p16_dst_field);
void abase_p16_mpi_ops_clear(abase_p16_dst_field);

/* Object-oriented interface */
#define abase_p16_oo_impl_name(v)	"p16"
static inline
void abase_p16_oo_field_clear(abase_vbase_ptr);
void abase_p16_oo_field_init(abase_vbase_ptr);
#ifdef  __cplusplus
}
#endif

/* Implementations for inlines */
/* *simd_p16::code_for_field_characteristic */
static inline
void abase_p16_field_characteristic(abase_p16_dst_field K, mpz_t z)
{
    mpz_set_si(z,*K);
}

/* *Mpfq::defaults::flatdata::code_for_set, simd_flat */
static inline
void abase_p16_set(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_dst_elt r, abase_p16_src_elt s)
{
    if (r != s) memcpy(r,s,sizeof(abase_p16_elt));
}

/* *simd_flat::code_for_set_zero */
static inline
void abase_p16_set_zero(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_dst_elt r)
{
    memset(r, 0, sizeof(abase_p16_elt));
}

/* *simd_p16::code_for_random */
static inline
void abase_p16_random(abase_p16_dst_field K, abase_p16_dst_elt r)
{
    *r = rand() % *K;
}

/* *Mpfq::defaults::flatdata::code_for_elt_ur_set, simd_flat */
static inline
void abase_p16_elt_ur_set(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_dst_elt_ur r, abase_p16_src_elt_ur s)
{
    if (r != s) memcpy(r,s,sizeof(abase_p16_elt_ur));
}

/* *Mpfq::defaults::flatdata::code_for_elt_ur_set_zero, simd_flat */
static inline
void abase_p16_elt_ur_set_zero(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_dst_elt_ur r)
{
    memset(r, 0, sizeof(abase_p16_elt_ur));
}

/* *simd_p16::code_for_reduce */
static inline
void abase_p16_reduce(abase_p16_dst_field K, abase_p16_dst_elt x, abase_p16_dst_elt_ur y)
{
    *x = *y % *K; *x += *K & -(*x < 0);
}

/* *Mpfq::defaults::flatdata::code_for_cmp, simd_flat */
static inline
int abase_p16_cmp(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_src_elt r, abase_p16_src_elt s)
{
    return memcmp(r,s,sizeof(abase_p16_elt));
}

/* *Mpfq::defaults::flatdata::code_for_is_zero, simd_flat */
static inline
int abase_p16_is_zero(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_src_elt r)
{
        unsigned int i;
        for(i = 0 ; i < sizeof(abase_p16_elt)/sizeof(r[0]) ; i++) {
            if (r[i]) return 0;
        }
        return 1;
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_set, Mpfq::defaults::flatdata, simd_flat */
static inline
void abase_p16_vec_set(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_dst_vec r, abase_p16_src_vec s, unsigned int n)
{
    if (r != s) memcpy(r, s, n*sizeof(abase_p16_elt));
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_set_zero, Mpfq::defaults::flatdata, simd_flat */
static inline
void abase_p16_vec_set_zero(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_dst_vec r, unsigned int n)
{
    memset(r, 0, n*sizeof(abase_p16_elt));
}

/* *Mpfq::defaults::vec::getset::code_for_vec_setcoef, Mpfq::defaults::vec, Mpfq::defaults */
static inline
void abase_p16_vec_setcoef(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_dst_vec w, abase_p16_src_elt x, unsigned int i)
{
    abase_p16_set(K, w[i], x);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_getcoef, Mpfq::defaults::vec, Mpfq::defaults */
static inline
void abase_p16_vec_getcoef(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_dst_elt x, abase_p16_src_vec w, unsigned int i)
{
    abase_p16_set(K, x, w[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_add, Mpfq::defaults::vec, Mpfq::defaults */
static inline
void abase_p16_vec_add(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_dst_vec w, abase_p16_src_vec u, abase_p16_src_vec v, unsigned int n)
{
        unsigned int i;
    for(i = 0; i < n; i+=1)
        abase_p16_add(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_neg, Mpfq::defaults::vec, Mpfq::defaults */
static inline
void abase_p16_vec_neg(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_dst_vec w, abase_p16_src_vec u, unsigned int n)
{
        unsigned int i;
    for(i = 0; i < n; ++i)
        abase_p16_neg(K, w[i], u[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_rev, Mpfq::defaults::vec, Mpfq::defaults */
static inline
void abase_p16_vec_rev(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_dst_vec w, abase_p16_src_vec u, unsigned int n)
{
    unsigned int nn = n >> 1;
    abase_p16_elt tmp[1];
    abase_p16_init(K, tmp);
    unsigned int i;
    for(i = 0; i < nn; ++i) {
        abase_p16_set(K, tmp[0], u[i]);
        abase_p16_set(K, w[i], u[n-1-i]);
        abase_p16_set(K, w[n-1-i], tmp[0]);
    }
    if (n & 1)
        abase_p16_set(K, w[nn], u[nn]);
    abase_p16_clear(K, tmp);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_sub, Mpfq::defaults::vec, Mpfq::defaults */
static inline
void abase_p16_vec_sub(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_dst_vec w, abase_p16_src_vec u, abase_p16_src_vec v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; ++i)
        abase_p16_sub(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_random, Mpfq::defaults::vec, Mpfq::defaults */
static inline
void abase_p16_vec_random(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_dst_vec w, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; ++i)
        abase_p16_random(K, w[i]);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_cmp, Mpfq::defaults::vec, Mpfq::defaults */
static inline
int abase_p16_vec_cmp(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_src_vec u, abase_p16_src_vec v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; ++i) {
        int ret = abase_p16_cmp(K, u[i], v[i]);
        if (ret != 0)
            return ret;
    }
    return 0;
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_is_zero, Mpfq::defaults::flatdata, simd_flat */
static inline
int abase_p16_vec_is_zero(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_src_vec r, unsigned int n)
{
    unsigned int i;
    for(i = 0 ; i < n ; i+=1) {
        if (!abase_p16_is_zero(K,r[i])) return 0;
    }
    return 1;
}

/* *Mpfq::defaults::vec::flatdata::code_for_vec_ur_set, Mpfq::defaults::flatdata, simd_flat */
static inline
void abase_p16_vec_ur_set(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_dst_vec_ur r, abase_p16_src_vec_ur s, unsigned int n)
{
    if (r != s) memcpy(r, s, n*sizeof(abase_p16_elt_ur));
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_setcoef, Mpfq::defaults::vec, Mpfq::defaults */
static inline
void abase_p16_vec_ur_setcoef(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_dst_vec_ur w, abase_p16_src_elt_ur x, unsigned int i)
{
    abase_p16_elt_ur_set(K, w[i], x);
}

/* *Mpfq::defaults::vec::getset::code_for_vec_ur_getcoef, Mpfq::defaults::vec, Mpfq::defaults */
static inline
void abase_p16_vec_ur_getcoef(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_dst_elt_ur x, abase_p16_src_vec_ur w, unsigned int i)
{
    abase_p16_elt_ur_set(K, x, w[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_add, Mpfq::defaults::vec, Mpfq::defaults */
static inline
void abase_p16_vec_ur_add(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_dst_vec_ur w, abase_p16_src_vec_ur u, abase_p16_src_vec_ur v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        abase_p16_elt_ur_add(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::addsub::code_for_vec_ur_sub, Mpfq::defaults::vec, Mpfq::defaults */
static inline
void abase_p16_vec_ur_sub(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_dst_vec_ur w, abase_p16_src_vec_ur u, abase_p16_src_vec_ur v, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        abase_p16_elt_ur_sub(K, w[i], u[i], v[i]);
}

/* *Mpfq::defaults::vec::mul::code_for_vec_reduce, Mpfq::defaults::vec, Mpfq::defaults */
static inline
void abase_p16_vec_reduce(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_dst_vec w, abase_p16_dst_vec_ur u, unsigned int n)
{
    unsigned int i;
    for(i = 0; i < n; i+=1)
        abase_p16_reduce(K, w[i], u[i]);
}

/* *simd_p16::code_for_set_ui_at */
static inline
void abase_p16_set_ui_at(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_dst_elt p, int k MAYBE_UNUSED, unsigned long v)
{
        assert(k < abase_p16_groupsize(K));
        *p = v;
}

/* *simd_p16::code_for_set_ui_all */
static inline
void abase_p16_set_ui_all(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_dst_elt p, unsigned long v)
{
    *p=v;
}

/* *simd_p16::code_for_elt_ur_set_ui_at */
static inline
void abase_p16_elt_ur_set_ui_at(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_dst_elt p, int k MAYBE_UNUSED, unsigned long v)
{
        assert(k < abase_p16_groupsize(K));
        *p = v;
}

/* *simd_p16::code_for_elt_ur_set_ui_all */
static inline
void abase_p16_elt_ur_set_ui_all(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_dst_elt p, unsigned long v)
{
    *p=v;
}

static inline
void abase_p16_oo_field_clear(abase_vbase_ptr f)
{
    abase_p16_field_clear((abase_p16_dst_field)(f->obj));
    free(f->obj);
    f->obj = NULL;
}


#endif  /* ABASE_P16_H_ */

/* vim:set ft=cpp: */
