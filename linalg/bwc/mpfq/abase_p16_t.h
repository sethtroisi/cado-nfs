#ifndef ABASE_P16_T_H_
#define ABASE_P16_T_H_

/* MPFQ generated file -- do not edit */

#include "abase_p16.h"
/* Active handler: simd_p16 */
/* Automatically generated code  */
/* Active handler: p16 */
/* Active handler: Mpfq::defaults */
/* Active handler: Mpfq::defaults::vec */
/* Active handler: Mpfq::defaults::mpi_flat */
/* Active handler: io */
/* Active handler: trivialities */
/* Active handler: simd_flat */
/* Options used: w=64 vtag=p16 tag=p16 vbase_stuff={
                 'vc:includes' => [
                                    '<stdarg.h>'
                                  ],
                 'member_templates_restrict' => {
                                                  'p_4' => [
                                                             'p_4'
                                                           ],
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
                                 $vbase_stuff->{'member_templates_restrict'}{'p16'},
                                 $vbase_stuff->{'member_templates_restrict'}{'p_4'}
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
