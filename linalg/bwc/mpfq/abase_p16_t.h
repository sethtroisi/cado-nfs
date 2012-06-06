#ifndef ABASE_P16_T_H_
#define ABASE_P16_T_H_

/* MPFQ generated file -- do not edit */

#include "abase_p16.h"
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
void abase_p16_p16_dotprod(abase_p16_dst_field, abase_p16_dst_field, abase_p16_dst_vec, abase_p16_src_vec, abase_p16_src_vec, unsigned int);
void abase_p16_p16_addmul_tiny(abase_p16_dst_field, abase_p16_dst_field, abase_p16_dst_vec, abase_p16_src_vec, abase_p16_dst_vec, unsigned int);
void abase_p16_p16_transpose(abase_p16_dst_field, abase_p16_dst_field, abase_p16_dst_vec, abase_p16_src_vec);

#endif  /* ABASE_P16_T_H_ */

/* vim:set ft=cpp: */
