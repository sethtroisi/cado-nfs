/* MPFQ generated file -- do not edit */

#define _POSIX_C_SOURCE 200112L
#include "abase_p16_t.h"

/* Active handler: Mpfq::defaults::vec */
/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::engine::defaults::mpi_flat */
/* Active handler: simd_flat */
/* Active handler: io */
/* Active handler: simd_p16 */
/* Automatically generated code  */
/* Options used: vtag=p16 tag=p16 choose_by_groupsize=<code> prefix=abase_ virtual_base={
                  'filebase' => 'abase_vbase',
                  'substitutions' => [
                                       [
                                         qr/(?-xism:abase_p16_elt \*)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?-xism:abase_p16_src_elt\b)/,
                                         'const void *'
                                       ],
                                       [
                                         qr/(?-xism:abase_p16_elt\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?-xism:abase_p16_dst_elt\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?-xism:abase_p16_elt_ur \*)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?-xism:abase_p16_src_elt_ur\b)/,
                                         'const void *'
                                       ],
                                       [
                                         qr/(?-xism:abase_p16_elt_ur\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?-xism:abase_p16_dst_elt_ur\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?-xism:abase_p16_vec \*)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?-xism:abase_p16_src_vec\b)/,
                                         'const void *'
                                       ],
                                       [
                                         qr/(?-xism:abase_p16_vec\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?-xism:abase_p16_dst_vec\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?-xism:abase_p16_vec_ur \*)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?-xism:abase_p16_src_vec_ur\b)/,
                                         'const void *'
                                       ],
                                       [
                                         qr/(?-xism:abase_p16_vec_ur\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?-xism:abase_p16_dst_vec_ur\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?-xism:abase_p16_poly \*)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?-xism:abase_p16_src_poly\b)/,
                                         'const void *'
                                       ],
                                       [
                                         qr/(?-xism:abase_p16_poly\b)/,
                                         'void *'
                                       ],
                                       [
                                         qr/(?-xism:abase_p16_dst_poly\b)/,
                                         'void *'
                                       ]
                                     ],
                  'name' => 'abase_vbase'
                };
 family=[p16] */


/* Functions operating on the field structure */

/* Element allocation functions */

/* Elementary assignment functions */

/* Assignment of random values */

/* Arithmetic operations on elements */

/* Operations involving unreduced elements */

/* Comparison functions */

/* Input/output functions */

/* Vector functions */

/* Functions related to SIMD operation */

/* Member templates related to SIMD operation */

/* MPI interface */

/* Object-oriented interface */
/* *simd_p16::code_for_member_template_dotprod */
void abase_p16_p16_dotprod(abase_p16_dst_field K0 MAYBE_UNUSED, abase_p16_dst_field K1 MAYBE_UNUSED, abase_p16_dst_vec xw, abase_p16_src_vec xu1, abase_p16_src_vec xu0, unsigned int n)
{
        int32_t s = 0;
        for(unsigned int i = 0 ; i < n ; i++) {
            s+=xu0[i][0] * xu1[i][0];
        }
        xw[0][0] =s % *K0;
}

/* *simd_p16::code_for_member_template_addmul_tiny */
void abase_p16_p16_addmul_tiny(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_dst_field L MAYBE_UNUSED, abase_p16_dst_vec w, abase_p16_src_vec u, abase_p16_dst_vec v, unsigned int n)
{
        for(unsigned int i = 0 ; i < n ; i++) {
            w[i][0] += u[i][0] * v[0][0];
        }
}

/* *simd_p16::code_for_member_template_transpose */
void abase_p16_p16_transpose(abase_p16_dst_field K MAYBE_UNUSED, abase_p16_dst_field L MAYBE_UNUSED, abase_p16_dst_vec w, abase_p16_src_vec u)
{
    w[0][0]=u[0][0];
}


/* vim:set ft=cpp: */
